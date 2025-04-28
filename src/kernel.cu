#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include "const.h"
#include "kernel.h"

__global__ void computeCellIdSoA(int particleCount, float *posX, float *posY, float *posZ, int *cellId, float cellSize,
                                 float xLen, float yLen, float zLen, int gridDimX, int gridDimY, int gridDimZ) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particleCount) return;

  float px = posX[i];
  float py = posY[i];
  float pz = posZ[i];

  float halfX = xLen / 2.0f;
  float halfY = yLen / 2.0f;
  float halfZ = zLen / 2.0f;
  int ix = (int)floorf((px + halfX) / cellSize);
  int iy = (int)floorf((py + halfY) / cellSize);
  int iz = (int)floorf((pz + halfZ) / cellSize);

  ix = min(max(ix, 0), gridDimX - 1);
  iy = min(max(iy, 0), gridDimY - 1);
  iz = min(max(iz, 0), gridDimZ - 1);

  cellId[i] = ix + iy * gridDimX + iz * gridDimX * gridDimY;
}

__global__ void initCells(int *cellStart, int *cellEnd, int totalCells) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= totalCells) return;
  cellStart[idx] = -1;
  cellEnd[idx] = -1;
}

__global__ void findCellStartEndSoA(int particleCount, const int *particleIndices, const int *cellId, int *cellStart,
                                    int *cellEnd, int totalCells) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= particleCount) return;
  if (particleCount == 0) return;

  int originalIndex = particleIndices[idx];
  int currentCellId = cellId[originalIndex];

  if (idx == 0) {
    cellStart[currentCellId] = 0;
  } else {
    int prevOriginalIndex = particleIndices[idx - 1];
    int prevCellId = cellId[prevOriginalIndex];
    if (currentCellId != prevCellId) {
      cellEnd[prevCellId] = idx;
      cellStart[currentCellId] = idx;
    }
  }

  if (idx == particleCount - 1) {
    // Last particle marks the end of its cell
    cellEnd[currentCellId] = particleCount;
  }
}

__global__ void computeDensityPressureSoA(int particleCount, const int *particleIndices, float *posX, float *posY,
                                          float *posZ, float *density, float *pressure, float mass,
                                          const int *cellStart, const int *cellEnd, float cellSize, int gridDimX,
                                          int gridDimY, int gridDimZ, float xLen, float yLen, float zLen, float POLY6,
                                          float WEIGHT_AT_0) {
  int i_sorted = blockDim.x * blockIdx.x + threadIdx.x;
  if (i_sorted >= particleCount) return;

  int i = particleIndices[i_sorted];

  float px_i = posX[i];
  float py_i = posY[i];
  float pz_i = posZ[i];
  float currentDensity = 0.0f;

  float h2 = KERNEL_RADIUS * KERNEL_RADIUS;
  float C = mass * POLY6;

  float halfX = xLen / 2.0f;
  float halfY = yLen / 2.0f;
  float halfZ = zLen / 2.0f;
  int ix = min(max((int)floorf((px_i + halfX) / cellSize), 0), gridDimX - 1);
  int iy = min(max((int)floorf((py_i + halfY) / cellSize), 0), gridDimY - 1);
  int iz = min(max((int)floorf((pz_i + halfZ) / cellSize), 0), gridDimZ - 1);

  // Loop over neighboring cells (3×3×3)
  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {
        int nx = ix + dx;
        int ny = iy + dy;
        int nz = iz + dz;
        if (nx < 0 || nx >= gridDimX || ny < 0 || ny >= gridDimY || nz < 0 || nz >= gridDimZ) continue;

        int neighborCellIdx = nx + ny * gridDimX + nz * gridDimX * gridDimY;
        int start = cellStart[neighborCellIdx];
        int end = cellEnd[neighborCellIdx];

        if (start != -1) {
          for (int j_sorted = start; j_sorted < end; j_sorted++) {
            int j = particleIndices[j_sorted];
            if (j == i) continue;

            float dx_ij = px_i - posX[j];
            float dy_ij = py_i - posY[j];
            float dz_ij = pz_i - posZ[j];
            float r2 = dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij;
            float zVal = h2 - r2;

            if (zVal > 0.0f) {
              float rho_contrib = C * zVal * zVal * zVal;
              currentDensity += rho_contrib;
            }
          }
        }
      }
    }
  }

  currentDensity += mass * WEIGHT_AT_0;  // Self-density contribution
  density[i] = currentDensity;
  pressure[i] = (powf(currentDensity / REST_DENSITY, 7.0f) - 1.0f) * STIFFNESS;
}

__global__ void computeAccelSoA(int particleCount, const int *particleIndices, float *posX, float *posY, float *posZ,
                                float *velX, float *velY, float *velZ, float *avgVelX, float *avgVelY, float *avgVelZ,
                                float *accX, float *accY, float *accZ, const float *density, const float *pressure,
                                float mass, const int *cellStart, const int *cellEnd, float cellSize, int gridDimX,
                                int gridDimY, int gridDimZ, float xLen, float yLen, float zLen,
                                float VISCOSITY_LAPLACIAN) {
  int i_sorted = blockDim.x * blockIdx.x + threadIdx.x;
  if (i_sorted >= particleCount) return;

  int i = particleIndices[i_sorted];

  float px_i = posX[i];
  float py_i = posY[i];
  float pz_i = posZ[i];
  float density_i = density[i];
  float pressure_i = pressure[i];
  float avgVelX_i = avgVelX[i];
  float avgVelY_i = avgVelY[i];
  float avgVelZ_i = avgVelZ[i];

  float forceX = 0.0f;
  float forceY = 0.0f;
  float forceZ = 0.0f;

  float h = KERNEL_RADIUS;
  float h2 = h * h;

  float halfX = xLen / 2.0f;
  float halfY = yLen / 2.0f;
  float halfZ = zLen / 2.0f;
  int ix = min(max((int)floorf((px_i + halfX) / cellSize), 0), gridDimX - 1);
  int iy = min(max((int)floorf((py_i + halfY) / cellSize), 0), gridDimY - 1);
  int iz = min(max((int)floorf((pz_i + halfZ) / cellSize), 0), gridDimZ - 1);

  // Loop over neighboring cells (3×3×3)
  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {
        int nx = ix + dx;
        int ny = iy + dy;
        int nz = iz + dz;
        if (nx < 0 || nx >= gridDimX || ny < 0 || ny >= gridDimY || nz < 0 || nz >= gridDimZ) continue;

        int neighborCellIdx = nx + ny * gridDimX + nz * gridDimX * gridDimY;
        int start = cellStart[neighborCellIdx];
        int end = cellEnd[neighborCellIdx];

        if (start != -1) {
          for (int j_sorted = start; j_sorted < end; j_sorted++) {
            int j = particleIndices[j_sorted];
            if (i == j) continue;

            float dx_ij = px_i - posX[j];
            float dy_ij = py_i - posY[j];
            float dz_ij = pz_i - posZ[j];
            float r2 = dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij;

            if (r2 < h2 && r2 > 1e-12f) {
              float r = sqrtf(r2);
              float h_minus_r = h - r;

              float density_j = density[j];
              float pressure_j = pressure[j];

              float pressure_term = mass * (pressure_i + pressure_j) / (2.0f * density_j);
              float spiky_grad_factor = -VISCOSITY_LAPLACIAN * h_minus_r * h_minus_r / r;
              forceX += pressure_term * spiky_grad_factor * dx_ij;
              forceY += pressure_term * spiky_grad_factor * dy_ij;
              forceZ += pressure_term * spiky_grad_factor * dz_ij;

              float viscosity_term = VISCOSITY * mass / density_j * VISCOSITY_LAPLACIAN * h_minus_r;
              forceX += viscosity_term * (avgVelX[j] - avgVelX_i);
              forceY += viscosity_term * (avgVelY[j] - avgVelY_i);
              forceZ += viscosity_term * (avgVelZ[j] - avgVelZ_i);
            }
          }
        }
      }
    }
  }

  if (density_i > 1e-12f) {
    accX[i] = forceX / density_i;
    accY[i] = forceY / density_i;
    accZ[i] = forceZ / density_i;
  } else {
    accX[i] = 0.0f;
    accY[i] = 0.0f;
    accZ[i] = 0.0f;
  }
}

__device__ void reflectInSinkSoA(float &px, float &py, float &pz, float &vx, float &vy, float &vz, float xLen,
                                 float yLen, float zLen) {
  float halfX = xLen / 2.0f;
  float halfY = yLen / 2.0f;
  float halfZ = zLen / 2.0f;
  float damping = REFLECT_DAMP;

  // X boundary
  if (px > halfX) {
    px = halfX - (px - halfX);
    vx = -vx * damping;
    vy *= damping;
    vz *= damping;
  } else if (px < -halfX) {
    px = -halfX + (-halfX - px);
    vx = -vx * damping;
    vy *= damping;
    vz *= damping;
  }

  // Y boundary
  if (py > halfY) {
    py = halfY - (py - halfY);
    vy = -vy * damping;
    vx *= damping;
    vz *= damping;
  } else if (py < -halfY) {
    py = -halfY + (-halfY - py);
    vy = -vy * damping;
    vx *= damping;
    vz *= damping;
  }

  // Z boundary
  if (pz > halfZ) {
    pz = halfZ - (pz - halfZ);
    vz = -vz * damping;
    vx *= damping;
    vy *= damping;
  } else if (pz < -halfZ) {
    pz = -halfZ + (-halfZ - pz);
    vz = -vz * damping;
    vx *= damping;
    vy *= damping;
  }
}

__device__ void reflectInTroughSoA(float &px, float &py, float &pz, float &vx, float &vy, float &vz, float zLen,
                                   float slope, float intercept, Vec3 normal) {
  float halfZ = zLen / 2.0f;
  float damping = REFLECT_DAMP;

  float planeY = slope * px + intercept;
  if (py < planeY) {
    py = planeY + (planeY - py);
    float dotVN = vx * normal.x + vy * normal.y + vz * normal.z;
    vx = (vx - 2.0f * dotVN * normal.x) * damping;
    vy = (vy - 2.0f * dotVN * normal.y) * damping;
    vz = (vz - 2.0f * dotVN * normal.z) * damping;
    py = planeY + 0.001f;
  }

  if (pz > halfZ) {
    pz = halfZ - (pz - halfZ);
    vz = -vz * damping;
    vx *= damping;
    vy *= damping;
  } else if (pz < -halfZ) {
    pz = -halfZ + (-halfZ - pz);
    vz = -vz * damping;
    vx *= damping;
    vy *= damping;
  }
}

__global__ void integrationSoA(int particleCount, const int *particleIndices, float *posX, float *posY, float *posZ,
                               float *velX, float *velY, float *velZ, float *avgVelX, float *avgVelY, float *avgVelZ,
                               float *accX, float *accY, float *accZ, char *inSink, float sinkXLen, float sinkYLen,
                               float sinkZLen, float troughZLen, float slope, float intercept, Vec3 normal) {
  int i_sorted = blockDim.x * blockIdx.x + threadIdx.x;
  if (i_sorted >= particleCount) return;

  int i = particleIndices[i_sorted];

  float vx_i = velX[i] + accX[i] * DELTA_T;
  float vy_i = velY[i] + (accY[i] + GRAVITY) * DELTA_T;
  float vz_i = velZ[i] + accZ[i] * DELTA_T;

  // Update position
  float px_i = posX[i] + vx_i * DELTA_T;
  float py_i = posY[i] + vy_i * DELTA_T;
  float pz_i = posZ[i] + vz_i * DELTA_T;

  char inSink_i = inSink[i];
  if (inSink_i == 0) {
    bool nowInSink = (px_i > -sinkXLen / 2.0f && px_i < sinkXLen / 2.0f) &&
                     (py_i > -sinkYLen / 2.0f && py_i < sinkYLen / 2.0f) &&
                     (pz_i > -sinkZLen / 2.0f && pz_i < sinkZLen / 2.0f);
    if (nowInSink) {
      inSink_i = 1;
    }
  }

  if (inSink_i == 1) {
    reflectInSinkSoA(px_i, py_i, pz_i, vx_i, vy_i, vz_i, sinkXLen, sinkYLen, sinkZLen);
  } else {
    reflectInTroughSoA(px_i, py_i, pz_i, vx_i, vy_i, vz_i, troughZLen, slope, intercept, normal);
  }

  avgVelX[i] = (avgVelX[i] + vx_i) / 2.0f;
  avgVelY[i] = (avgVelY[i] + vy_i) / 2.0f;
  avgVelZ[i] = (avgVelZ[i] + vz_i) / 2.0f;

  posX[i] = px_i;
  posY[i] = py_i;
  posZ[i] = pz_i;
  velX[i] = vx_i;
  velY[i] = vy_i;
  velZ[i] = vz_i;
  inSink[i] = inSink_i;
}

__global__ void coordTransformSoA(int particleCount, const float *posX, const float *posY, const float *posZ,
                                  const float *transformMat, Vec2 *screenPosOnGPU) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particleCount) {
    return;
  }

  float worldPos[4];
  worldPos[0] = posX[i];
  worldPos[1] = posY[i];
  worldPos[2] = posZ[i];
  worldPos[3] = 1.0f;

  float ndcPos[4];
  for (int row = 0; row < 4; row++) {
    ndcPos[row] = 0.0f;
    for (int col = 0; col < 4; col++) {
      ndcPos[row] += transformMat[row * 4 + col] * worldPos[col];
    }
  }

  float invW = (ndcPos[3] == 0.0f) ? 1.0f : 1.0f / ndcPos[3];
  float ndcX = ndcPos[0] * invW;
  float ndcY = ndcPos[1] * invW;

  if (ndcPos[3] <= 0.0f || ndcX < -1.0f || ndcX > 1.0f || ndcY < -1.0f || ndcY > 1.0f) {
    screenPosOnGPU[i].x = -1.0f;
    screenPosOnGPU[i].y = -1.0f;
  } else {
    ndcX = fmaxf(-1.0f, fminf(1.0f, ndcX));
    ndcY = fmaxf(-1.0f, fminf(1.0f, ndcY));

    float screenX = (ndcX + 1.0f) * 0.5f * SCREEN_WIDTH;
    float screenY = (1.0f - ndcY) * 0.5f * SCREEN_HEIGHT;

    screenPosOnGPU[i].x = screenX;
    screenPosOnGPU[i].y = screenY;
  }
}