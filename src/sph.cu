#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "const.h"
#include "sph.h"

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

// Self-defined comparator for thrust::sort
struct ParticleComparator {
  __host__ __device__ bool operator()(const Particle &a, const Particle &b) const { return a.cellId < b.cellId; }
};

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

float *allocateMatOnGPU(Mat4 &mat) {
  float data[16];
  int count = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      data[count++] = mat[i][j];
    }
  }
  float *matOnGPU;
  cudaMalloc((void **)&matOnGPU, 16 * sizeof(float));
  cudaMemcpy(matOnGPU, data, 16 * sizeof(float), cudaMemcpyHostToDevice);
  return matOnGPU;
}

// Sort particles based on their cell IDs and find the start and end indices of each cell.
void sortParticles(Particle *particles, int particleCount, int *&cellStart, int *&cellEnd, float cellSize, float xLen,
                   float yLen, float zLen, int gridDimX, int gridDimY, int gridDimZ) {
  int totalCells = gridDimX * gridDimY * gridDimZ;

  int threads = 128;
  int blocks = (particleCount + threads - 1) / threads;
  computeCellId<<<blocks, threads>>>(particles, particleCount, cellSize, xLen, yLen, zLen, gridDimX, gridDimY,
                                     gridDimZ);
  cudaDeviceSynchronize();

  thrust::device_ptr<Particle> dev_ptr(particles);
  thrust::sort(dev_ptr, dev_ptr + particleCount, ParticleComparator());

  initCells<<<(totalCells + threads - 1) / threads, threads>>>(cellStart, cellEnd, totalCells);
  cudaDeviceSynchronize();
  findCellStartEnd<<<(particleCount + (threads - 1)) / threads, threads>>>(particles, particleCount, cellStart, cellEnd,
                                                                           totalCells);
  cudaDeviceSynchronize();
}

Particle *placeParticles(int &particleCount, int &droppingparticleCount, Sink &sink, Trough &trough) {
  float h = KERNEL_RADIUS;
  float hh = h / 2.0f;
  std::cout << "hh: " << hh << std::endl;
  int particlesInSink = 0;
  for (float x = -sink.xLen / 2.0f; x <= sink.xLen / 2.0f; x += hh) {
    for (float z = -sink.zLen / 2.0f; z <= sink.zLen / 2.0f; z += hh) {
      for (float y = -sink.yLen / 2.0f; y <= sink.yLen / 4.0f; y += hh) {
        particlesInSink++;
      }
    }
  }
  int dropping = 0;
  for (float x = trough.vertices[0].x + 0.001f; x <= trough.vertices[1].x - 0.001f; x += hh) {
    float y0 = trough.slope * x + trough.intercept;
    for (float z = -trough.zLen / 2.0f + 0.001f; z <= trough.zLen / 2.0f - 0.001f; z += hh) {
      for (float y = y0; y <= y0 + trough.yLen; y += hh) {
        dropping++;
      }
    }
  }
  droppingparticleCount = dropping;
  std::cout << "dropping particle count: " << droppingparticleCount << std::endl;
  particleCount = particlesInSink + droppingparticleCount;
  std::cout << "Particle Count: " << particleCount << std::endl;
  Particle *particlesOnGPU;
  cudaMalloc(&particlesOnGPU, particleCount * sizeof(Particle));
  Particle *particlesOnCPU = new Particle[particleCount];
  int count = 0;
  // particles generated in sink
  for (float x = -sink.xLen / 2.0f; x <= sink.xLen / 2.0f; x += hh) {
    for (float z = -sink.zLen / 2.0f; z <= sink.zLen / 2.0f; z += hh) {
      for (float y = -sink.yLen / 2.0f; y <= sink.yLen / 4.0f; y += hh) {
        particlesOnCPU[count].position = {x, y, z};
        particlesOnCPU[count].density = 0.0f;
        particlesOnCPU[count].inSink = true;
        particlesOnCPU[count].velocity = {0.0f, 0.0f, 0.0f};
        particlesOnCPU[count].averageVelocity = {0.0f, 0.0f, 0.0f};
        particlesOnCPU[count].acceleration = {0.0f, 0.0f, 0.0f};
        count++;
      }
    }
  }

  // particles that will fall on trough
  float vx = 2.5f;
  float vy = vx * trough.slope;
  for (float x = trough.vertices[0].x + 0.001f; x <= trough.vertices[1].x - 0.001f; x += hh) {
    float y0 = trough.slope * x + trough.intercept;
    for (float z = -trough.zLen / 2.0f + 0.001f; z <= trough.zLen / 2.0f - 0.001f; z += hh) {
      for (float y = y0; y <= y0 + trough.yLen; y += hh) {
        particlesOnCPU[count].position = {x, y, z};
        particlesOnCPU[count].density = 0.0f;
        particlesOnCPU[count].inSink = false;
        particlesOnCPU[count].velocity = {vx, vy, 0.0f};
        particlesOnCPU[count].averageVelocity = {vx, vy, 0.0f};
        particlesOnCPU[count].acceleration = {0.0f, 0.0f, 0.0f};
        count++;
      }
    }
  }
  assert(count == particleCount);
  cudaError_t copyErr =
      cudaMemcpy(particlesOnGPU, particlesOnCPU, particleCount * sizeof(Particle), cudaMemcpyHostToDevice);
  if (copyErr != cudaSuccess) {
    std::cerr << "cudaMemcpy particles error: " << cudaGetErrorString(copyErr) << std::endl;
  }
  return particlesOnGPU;
}

int *initCellStart(int totalCells) {
  int *cellStart;
  cudaMalloc(&cellStart, totalCells * sizeof(int));
  return cellStart;
}

int *initCellEnd(int totalCells) {
  int *cellEnd;
  cudaMalloc(&cellEnd, totalCells * sizeof(int));
  return cellEnd;
}

// Normalize mass based on the density of particles in the sink
float normalizeMass(Particle *particles, int particleCount, const Sink &sink, int *cellStart, int *cellEnd,
                    float POLY6, float WEIGHT_AT_0) {
  float mass = 1.0f;
  int blockDim = 32;
  int gridDim = (particleCount + (blockDim - 1)) / blockDim;

  // Use Sink dimensions as simulation domain
  float xLen = sink.xLen;
  float yLen = sink.yLen;
  float zLen = sink.zLen;
  float cellSize = KERNEL_RADIUS;
  int gridDimX = (int)ceil(sink.xLen / cellSize);
  int gridDimY = (int)ceil(sink.yLen / cellSize);
  int gridDimZ = (int)ceil(sink.zLen / cellSize);
  sortParticles(particles, particleCount, cellStart, cellEnd, cellSize, xLen, yLen, zLen, gridDimX, gridDimY, gridDimZ);

  computeDensityPressureSorted<<<gridDim, blockDim>>>(particles, particleCount, mass, cellStart, cellEnd, cellSize,
                                                      gridDimX, gridDimY, gridDimZ, xLen, yLen, zLen, POLY6, WEIGHT_AT_0);
  cudaDeviceSynchronize();

  cudaError_t err;
  if ((err = cudaGetLastError()) != cudaSuccess)
    std::cerr << "Kernel error (computeDensityPressureSorted): " << cudaGetErrorString(err) << std::endl;

  float rho0 = REST_DENSITY;
  float rho2s = 0.0f;
  float rhos = 0.0f;
  Particle *particlesOnCPU = new Particle[particleCount];
  cudaError_t copyErr = cudaMemcpy(particlesOnCPU, particles, particleCount * sizeof(Particle), cudaMemcpyDeviceToHost);
  if (copyErr != cudaSuccess) {
    std::cerr << "cudaMemcpy Error: " << cudaGetErrorString(copyErr) << std::endl;
  }
  for (int i = 0; i < particleCount; i++) {
    rho2s += particlesOnCPU[i].density * particlesOnCPU[i].density;
    rhos += particlesOnCPU[i].density;
  }
  mass *= (rho0 * rhos / rho2s);
  std::cout << "Mass: " << mass << std::endl;

  return mass;
}

Particle *initParticles(int &particleCount, float &mass, Sink &sink, Trough &trough, int *cellStart, int *cellEnd,
                        float POLY6, float WEIGHT_AT_0) {
  int droppingparticleCount = 0;
  Particle *particlesOnGPU = placeParticles(particleCount, droppingparticleCount, sink, trough);
  mass = normalizeMass(particlesOnGPU, particleCount - droppingparticleCount, sink, cellStart, cellEnd, POLY6, WEIGHT_AT_0);
  return particlesOnGPU;
}

Vec2 *initScreenPos(int particleCount) {
  Vec2 *screenPos;
  cudaMalloc(&screenPos, particleCount * sizeof(Vec2));
  return screenPos;
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

void updateSimulation(Particle *particles, int particleCount, const Sink &sink, const Trough &trough, float mass,
                      float *transformMat, int *cellStart, int *cellEnd, Vec2 *screenPosOnGPU, Vec2 *screenPosOnCPU,
                      float POLY6, float VISCOSITY_LAPLACIAN, float WEIGHT_AT_0) {
  int blockDim = 32;
  int gridDim = (particleCount + (blockDim - 1)) / blockDim;

  float cellSize = KERNEL_RADIUS;
  float xLen = sink.xLen;
  float yLen = sink.yLen;
  float zLen = sink.zLen;
  int gridDimX = (int)ceil(xLen / cellSize);
  int gridDimY = (int)ceil(yLen / cellSize);
  int gridDimZ = (int)ceil(zLen / cellSize);

  sortParticles(particles, particleCount, cellStart, cellEnd, cellSize, xLen, yLen, zLen, gridDimX, gridDimY, gridDimZ);

  cudaError_t err;
  computeDensityPressureSorted<<<gridDim, blockDim>>>(particles, particleCount, mass, cellStart, cellEnd, cellSize,
                                                      gridDimX, gridDimY, gridDimZ, xLen, yLen, zLen, POLY6, WEIGHT_AT_0);
  cudaDeviceSynchronize();
  if ((err = cudaGetLastError()) != cudaSuccess)
    std::cerr << "Kernel error (computeDensityPressureSorted): " << cudaGetErrorString(err) << std::endl;

  computeAccelSorted<<<gridDim, blockDim>>>(particles, particleCount, mass, cellStart, cellEnd, cellSize, gridDimX,
                                            gridDimY, gridDimZ, xLen, yLen, zLen, VISCOSITY_LAPLACIAN);
  if ((err = cudaGetLastError()) != cudaSuccess)
    std::cerr << "Kernel error (computeAccelSorted): " << cudaGetErrorString(err) << std::endl;

  integration<<<gridDim, blockDim>>>(particles, particleCount, sink.xLen, sink.yLen, sink.zLen, trough.zLen,
                                     trough.slope, trough.intercept, trough.normal);
  coordTransform<<<gridDim, blockDim>>>(particles, particleCount, transformMat, screenPosOnGPU);
  cudaDeviceSynchronize();
  if ((err = cudaGetLastError()) != cudaSuccess)
    std::cerr << "Kernel error (integration or coordTransform): " << cudaGetErrorString(err) << std::endl;
  cudaError_t copyErr =
      cudaMemcpy(screenPosOnCPU, screenPosOnGPU, particleCount * sizeof(Vec2), cudaMemcpyDeviceToHost);
  if (copyErr != cudaSuccess) {
    std::cerr << "cudaMemcpy Error: " << cudaGetErrorString(copyErr) << std::endl;
  }
}
