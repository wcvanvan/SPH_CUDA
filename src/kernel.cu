#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include "const.h"
#include "kernel.h"

// Assign a cell id to each particle based on its position.
__global__ void computeCellId(Particle *particles, int particleCount, float cellSize, float xLen, float yLen,
                              float zLen, int gridDimX, int gridDimY, int gridDimZ) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particleCount) return;
  Particle &particle = particles[i];

  float halfX = xLen / 2.0f;
  float halfY = yLen / 2.0f;
  float halfZ = zLen / 2.0f;
  int ix = (int)floor((particle.position.x + halfX) / cellSize);
  int iy = (int)floor((particle.position.y + halfY) / cellSize);
  int iz = (int)floor((particle.position.z + halfZ) / cellSize);

  ix = min(max(ix, 0), gridDimX - 1);
  iy = min(max(iy, 0), gridDimY - 1);
  iz = min(max(iz, 0), gridDimZ - 1);

  particle.cellId = ix + iy * gridDimX + iz * gridDimX * gridDimY;
}

__global__ void initCells(int *cellStart, int *cellEnd, int totalCells) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= totalCells) return;
  cellStart[idx] = -1;
  cellEnd[idx] = -1;
}

// cellStart[i] = the first particle in cell i; cellEnd[i] = the first particle in cell i+1
__global__ void findCellStartEnd(Particle *particles, int particleCount, int *cellStart, int *cellEnd, int totalCells) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= particleCount) return;
  if (particleCount == 0) return;
  if (idx == 0) {
    cellStart[particles[0].cellId] = 0;
  } else {
    int cid = particles[idx].cellId;
    int prevCid = particles[idx - 1].cellId;
    if (cid != prevCid) {
      // as particles are sorted, no data race would happen
      cellEnd[prevCid] = idx;
      cellStart[cid] = idx;
    }
  }
  if (idx == particleCount - 1) {
    cellEnd[particles[particleCount - 1].cellId] = particleCount;
  }
}

// Update: each thread loops only over particles in its own and neighboring grid cells.
__global__ void computeDensityPressureSorted(Particle *particles, int particleCount, float mass, int *cellStart,
                                             int *cellEnd, float cellSize, int gridDimX, int gridDimY, int gridDimZ,
                                             float xLen, float yLen, float zLen, float POLY6, float WEIGHT_AT_0) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particleCount) return;
  Particle &particle = particles[i];

  particle.density = 0.0f;
  float h2 = KERNEL_RADIUS * KERNEL_RADIUS;
  float C = mass * POLY6;

  float halfX = xLen / 2.0f;
  float halfY = yLen / 2.0f;
  float halfZ = zLen / 2.0f;
  int ix = min(max((int)floor((particle.position.x + halfX) / cellSize), 0), gridDimX - 1);
  int iy = min(max((int)floor((particle.position.y + halfY) / cellSize), 0), gridDimY - 1);
  int iz = min(max((int)floor((particle.position.z + halfZ) / cellSize), 0), gridDimZ - 1);

  // Loop over neighboring cells (3×3×3)
  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {
        int nx = ix + dx;
        int ny = iy + dy;
        int nz = iz + dz;
        if (nx < 0 || nx >= gridDimX || ny < 0 || ny >= gridDimY || nz < 0 || nz >= gridDimZ) continue;
        int neighborCell = nx + ny * gridDimX + nz * gridDimX * gridDimY;
        int start = cellStart[neighborCell];
        int end = cellEnd[neighborCell];
        if (start == -1) continue;
        for (int j = start; j < end; j++) {
          if (j == i) continue;
          float dx = particle.position.x - particles[j].position.x;
          float dy = particle.position.y - particles[j].position.y;
          float dz = particle.position.z - particles[j].position.z;
          float r2 = dx * dx + dy * dy + dz * dz;
          float zVal = h2 - r2;
          if (zVal <= 0 || r2 < 1e-12) continue;
          float rho = C * zVal * zVal * zVal;
          particle.density += rho;
        }
      }
    }
  }
  particle.density += mass * WEIGHT_AT_0;  // contributing to the density of itself
  particle.pressure = (pow(particle.density / REST_DENSITY, 7) - 1.0f) * STIFFNESS;
}

// Update: each thread loops only over particles in its own and neighboring grid cells.
__global__ void computeAccelSorted(Particle *particles, int particleCount, float mass, int *cellStart, int *cellEnd,
                                   float cellSize, int gridDimX, int gridDimY, int gridDimZ, float xLen, float yLen,
                                   float zLen, float VISCOSITY_LAPLACIAN) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particleCount) return;
  Particle &particle = particles[i];

  particle.acceleration.x = 0.0f;
  particle.acceleration.y = 0.0f;
  particle.acceleration.z = 0.0f;

  float h2 = KERNEL_RADIUS * KERNEL_RADIUS;

  float halfX = xLen / 2.0f;
  float halfY = yLen / 2.0f;
  float halfZ = zLen / 2.0f;
  int ix = min(max((int)floor((particle.position.x + halfX) / cellSize), 0), gridDimX - 1);
  int iy = min(max((int)floor((particle.position.y + halfY) / cellSize), 0), gridDimY - 1);
  int iz = min(max((int)floor((particle.position.z + halfZ) / cellSize), 0), gridDimZ - 1);

  // Loop over neighboring cells (3×3×3)
  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {
        int nx = ix + dx;
        int ny = iy + dy;
        int nz = iz + dz;
        if (nx < 0 || nx >= gridDimX || ny < 0 || ny >= gridDimY || nz < 0 || nz >= gridDimZ) continue;
        int neighborCell = nx + ny * gridDimX + nz * gridDimX * gridDimY;
        int start = cellStart[neighborCell];
        int end = cellEnd[neighborCell];
        if (start == -1) continue;
        for (int j = start; j < end; j++) {
          if (j == i) continue;
          float dx = particles[i].position.x - particles[j].position.x;
          float dy = particles[i].position.y - particles[j].position.y;
          float dz = particles[i].position.z - particles[j].position.z;
          float r2 = dx * dx + dy * dy + dz * dz;

          if (r2 >= h2 || r2 <= 1e-12) continue;
          float r = sqrtf(r2);

          // pressure force push particles away
          float V = mass / particles[j].density / 2.0f;
          float Kr = KERNEL_RADIUS - r;
          float Kp = (-VISCOSITY_LAPLACIAN) * Kr * Kr;
          float pressureForce = V * (particle.pressure + particles[j].pressure) * Kp;
          particle.acceleration.x -= dx * pressureForce / r;
          particle.acceleration.y -= dy * pressureForce / r;
          particle.acceleration.z -= dz * pressureForce / r;

          // viscosity force pulls particles closer
          float Kv = VISCOSITY_LAPLACIAN * (KERNEL_RADIUS - r);
          float viscosityForce = V * VISCOSITY * Kv;
          float dvx = particles[j].averageVelocity.x - particle.averageVelocity.x;
          float dvy = particles[j].averageVelocity.y - particle.averageVelocity.y;
          float dvz = particles[j].averageVelocity.z - particle.averageVelocity.z;
          particle.acceleration.x += dvx * viscosityForce;
          particle.acceleration.y += dvy * viscosityForce;
          particle.acceleration.z += dvz * viscosityForce;
        }
      }
    }
  }
  particle.acceleration.x /= particle.density;
  particle.acceleration.y /= particle.density;
  particle.acceleration.z /= particle.density;
}

__device__ void reflectInSink(Particle &particle, float xLen, float yLen, float zLen) {
  float tbounce = 0.0f;
  if (particle.velocity.x != 0 && (particle.position.x > xLen / 2 || particle.position.x < -xLen / 2)) {
    if (particle.position.x > xLen / 2) {
      tbounce = (particle.position.x - xLen / 2) / particle.velocity.x;
      particle.position.x = xLen - particle.position.x;
    } else {
      tbounce = (particle.position.x + xLen / 2) / particle.velocity.x;
      particle.position.x = -xLen - particle.position.x;
    }
    // revert the movement for the period
    particle.position.y -= particle.velocity.y * (1 - REFLECT_DAMP) * tbounce;
    particle.position.z -= particle.velocity.z * (1 - REFLECT_DAMP) * tbounce;
    particle.velocity.x = -particle.velocity.x;
    particle.velocity.x *= REFLECT_DAMP;
    particle.velocity.y *= REFLECT_DAMP;
    particle.velocity.z *= REFLECT_DAMP;
  }
  if (particle.velocity.y != 0 && (particle.position.y > yLen / 2 || particle.position.y < -yLen / 2)) {
    // bounce back
    if (particle.position.y > yLen / 2) {
      tbounce = (particle.position.y - yLen / 2) / particle.velocity.y;
      particle.position.y = yLen - particle.position.y;
    } else {
      tbounce = (particle.position.y + yLen / 2) / particle.velocity.y;
      particle.position.y = -yLen - particle.position.y;
    }
    // revert the movement for the period
    particle.position.x -= particle.velocity.x * (1 - REFLECT_DAMP) * tbounce;
    particle.position.z -= particle.velocity.z * (1 - REFLECT_DAMP) * tbounce;
    particle.velocity.y = -particle.velocity.y;
    particle.velocity.x *= REFLECT_DAMP;
    particle.velocity.y *= REFLECT_DAMP;
    particle.velocity.z *= REFLECT_DAMP;
  }
  if (particle.velocity.z != 0 && (particle.position.z > zLen / 2 || particle.position.z < -zLen / 2)) {
    // bounce back
    if (particle.position.z > zLen / 2) {
      tbounce = (particle.position.z - zLen / 2) / particle.velocity.z;
      particle.position.z = zLen - particle.position.z;
    } else {
      tbounce = (particle.position.z + zLen / 2) / particle.velocity.z;
      particle.position.z = -zLen - particle.position.z;
    }
    // revert the movement for the period
    particle.position.x -= particle.velocity.x * (1 - REFLECT_DAMP) * tbounce;
    particle.position.y -= particle.velocity.y * (1 - REFLECT_DAMP) * tbounce;
    particle.velocity.z = -particle.velocity.z;
    particle.velocity.x *= REFLECT_DAMP;
    particle.velocity.y *= REFLECT_DAMP;
    particle.velocity.z *= REFLECT_DAMP;
  }
}

__device__ void reflectInTrough(Particle &particle, float zLen, float slope, float intercept, Vec3 normal) {
  float y = particle.position.x * slope + intercept;
  if (y > particle.position.y) {
    // hitting the bottom of the trough: v' = v - 2(v·N)N
    float dotV = particle.velocity.x * normal.x + particle.velocity.y * normal.y + particle.velocity.z * normal.z;
    float newVx = particle.velocity.x - 2 * dotV * normal.x;
    float newVy = particle.velocity.y - 2 * dotV * normal.y;
    float newVz = particle.velocity.z - 2 * dotV * normal.z;
    newVx *= REFLECT_DAMP;
    newVy *= REFLECT_DAMP;
    newVz *= REFLECT_DAMP;
    particle.velocity.x = newVx;
    particle.velocity.y = newVy;
    particle.velocity.z = newVz;
    particle.position.y = y + 0.001;  // simple method putting the particle back on trough
  }
  if (particle.velocity.z != 0 && (particle.position.z > zLen / 2 || particle.position.z < -zLen / 2)) {
    // hitting the side of the trough
    float tbounce = 0.0f;
    if (particle.position.z > zLen / 2) {
      tbounce = (particle.position.z - zLen / 2) / particle.velocity.z;
      particle.position.z = zLen - particle.position.z;
    } else {
      tbounce = (particle.position.z + zLen / 2) / particle.velocity.z;
      particle.position.z = -zLen - particle.position.z;
    }
    // revert the movement for the period
    particle.position.x -= particle.velocity.x * (1 - REFLECT_DAMP) * tbounce;
    particle.position.y -= particle.velocity.y * (1 - REFLECT_DAMP) * tbounce;
    particle.velocity.z = -particle.velocity.z;
    particle.velocity.x *= REFLECT_DAMP;
    particle.velocity.y *= REFLECT_DAMP;
    particle.velocity.z *= REFLECT_DAMP;
  }
}

__global__ void integration(Particle *particles, int particleCount, float sinkXLen, float sinkYLen, float sinkZLen,
                            float troughZLen, float slope, float intercept, Vec3 normal) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particleCount) {
    return;
  }
  Particle &particle = particles[i];
  particle.velocity.x += particle.acceleration.x * DELTA_T;
  particle.velocity.y += particle.acceleration.y * DELTA_T + GRAVITY * DELTA_T;
  particle.velocity.z += particle.acceleration.z * DELTA_T;
  particle.position.x += particle.velocity.x * DELTA_T;
  particle.position.y += particle.velocity.y * DELTA_T;
  particle.position.z += particle.velocity.z * DELTA_T;
  if (particle.inSink == false && (particle.position.x > -sinkXLen / 2.0f && particle.position.x < sinkXLen / 2.0f) &&
      (particle.position.y > -sinkYLen / 2.0f && particle.position.y < sinkYLen / 2.0f) &&
      (particle.position.z > -sinkZLen / 2.0f && particle.position.z < sinkZLen / 2.0f)) {
    particle.inSink = true;
  }
  if (particle.inSink) {
    reflectInSink(particle, sinkXLen, sinkYLen, sinkZLen);
  } else {
    reflectInTrough(particle, troughZLen, slope, intercept, normal);
  }
  particle.averageVelocity.x = (particle.averageVelocity.x + particle.velocity.x) / 2.0f;
  particle.averageVelocity.y = (particle.averageVelocity.y + particle.velocity.y) / 2.0f;
  particle.averageVelocity.z = (particle.averageVelocity.z + particle.velocity.z) / 2.0f;
}

__global__ void coordTransform(Particle *particles, int particleCount, float *transformMat, Vec2 *screenPosOnGPU) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particleCount) {
    return;
  }
  Particle &particle = particles[i];
  Vec2 &screenPos = screenPosOnGPU[i];
  float worldPos[4], result[4];
  worldPos[0] = particle.position.x;
  worldPos[1] = particle.position.y;
  worldPos[2] = particle.position.z;
  worldPos[3] = 1.0f;
  for (int i = 0; i < 4; i++) {
    result[i] = 0.0f;
    for (int j = 0; j < 4; j++) {
      result[i] += transformMat[i * 4 + j] * worldPos[j];
    }
  }
  float x = result[0] / result[3];
  float y = result[1] / result[3];

  if (x < -1.0f || x > 1.0f || y < -1.0f || y > 1.0f) {
    screenPos.x = -1.0f;
    screenPos.y = -1.0f;
  } else {
    float screenX = fmaxf(0.0f, fminf(1.0f, (x + 1.0f) * 0.5f)) * SCREEN_WIDTH;
    float screenY = fmaxf(0.0f, fminf(1.0f, (1.0f - y) * 0.5f)) * SCREEN_HEIGHT;
    screenPos.x = screenX;
    screenPos.y = screenY;
  }
}