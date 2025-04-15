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
__global__ void computeCellId(Particle *particles, int particlesCount, float cellSize, float xLen, float yLen,
                              float zLen, int gridDimX, int gridDimY, int gridDimZ) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particlesCount) return;

  float halfX = xLen / 2.0f;
  float halfY = yLen / 2.0f;
  float halfZ = zLen / 2.0f;
  int ix = (int)floor((particles[i].position.x + halfX) / cellSize);
  int iy = (int)floor((particles[i].position.y + halfY) / cellSize);
  int iz = (int)floor((particles[i].position.z + halfZ) / cellSize);

  ix = min(max(ix, 0), gridDimX - 1);
  iy = min(max(iy, 0), gridDimY - 1);
  iz = min(max(iz, 0), gridDimZ - 1);

  particles[i].cellId = ix + iy * gridDimX + iz * gridDimX * gridDimY;
}

// cellStart[i] = the first particle in cell i; cellEnd[i] = the first particle in cell i+1
__global__ void findCellStartEnd(Particle *particles, int particlesCount, int *cellStart, int *cellEnd,
                                 int totalCells) {
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    for (int i = 0; i < totalCells; i++) {
      cellStart[i] = -1;
      cellEnd[i] = -1;
    }
    if (particlesCount == 0) return;
    int currentCell = particles[0].cellId;
    cellStart[currentCell] = 0;
    for (int i = 1; i < particlesCount; i++) {
      int cid = particles[i].cellId;
      int prevCid = particles[i - 1].cellId;
      if (cid != prevCid) {
        cellEnd[prevCid] = i;
        cellStart[cid] = i;
      }
    }
    cellEnd[particles[particlesCount - 1].cellId] = particlesCount;
  }
}

// Self-defined comparator for thrust::sort
struct ParticleComparator {
  __host__ __device__ bool operator()(const Particle &a, const Particle &b) const { return a.cellId < b.cellId; }
};

// Update: each thread loops only over particles in its own and neighboring grid cells.
__global__ void computeDensitySorted(Particle *particles, int particlesCount, float mass, int *cellStart, int *cellEnd,
                                     float cellSize, int gridDimX, int gridDimY, int gridDimZ, float xLen, float yLen,
                                     float zLen) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particlesCount) return;

  particles[i].density = 0.0f;
  float h2 = KERNEL_RADIUS * KERNEL_RADIUS;
  float h4 = h2 * h2;
  float h8 = h4 * h4;
  float C = 4 * mass / M_PI / h8;
  particles[i].density += 4 * mass / M_PI / h2;

  float halfX = xLen / 2.0f;
  float halfY = yLen / 2.0f;
  float halfZ = zLen / 2.0f;
  int ix = min(max((int)floor((particles[i].position.x + halfX) / cellSize), 0), gridDimX - 1);
  int iy = min(max((int)floor((particles[i].position.y + halfY) / cellSize), 0), gridDimY - 1);
  int iz = min(max((int)floor((particles[i].position.z + halfZ) / cellSize), 0), gridDimZ - 1);

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
          float zVal = h2 - r2;
          if (zVal <= 0) continue;
          float rho = C * zVal * zVal * zVal;
          particles[i].density += rho;
        }
      }
    }
  }
}

// Update: each thread loops only over particles in its own and neighboring grid cells.
__global__ void computeAccelSorted(Particle *particles, int particlesCount, float mass, int *cellStart, int *cellEnd,
                                   float cellSize, int gridDimX, int gridDimY, int gridDimZ, float xLen, float yLen,
                                   float zLen) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particlesCount) return;

  particles[i].acceleration.x = 0.0f;
  particles[i].acceleration.y = GRAVITY;
  particles[i].acceleration.z = 0.0f;

  float h2 = KERNEL_RADIUS * KERNEL_RADIUS;
  float C0 = mass / M_PI / (h2 * h2);
  float Cp = 15 * STIFFNESS;
  float Cv = -40.0f * VISCOSITY;

  float halfX = xLen / 2.0f;
  float halfY = yLen / 2.0f;
  float halfZ = zLen / 2.0f;
  int ix = min(max((int)floor((particles[i].position.x + halfX) / cellSize), 0), gridDimX - 1);
  int iy = min(max((int)floor((particles[i].position.y + halfY) / cellSize), 0), gridDimY - 1);
  int iz = min(max((int)floor((particles[i].position.z + halfZ) / cellSize), 0), gridDimZ - 1);

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
          if (r2 >= h2) continue;
          float r = sqrtf(r2);
          float q = r / KERNEL_RADIUS;
          float u = 1 - q;
          float w0 = C0 * u / (particles[i].density * particles[j].density);
          float wp = w0 * Cp * (particles[i].density + particles[j].density - 2 * REST_DENSITY) * u / (q + 1e-6f);
          float wv = w0 * Cv;
          float dvx = particles[i].velocity.x - particles[j].velocity.x;
          float dvy = particles[i].velocity.y - particles[j].velocity.y;
          float dvz = particles[i].velocity.z - particles[j].velocity.z;
          particles[i].acceleration.x += (wp * dx + wv * dvx);
          particles[i].acceleration.y += (wp * dy + wv * dvy);
          particles[i].acceleration.z += (wp * dz + wv * dvz);
        }
      }
    }
  }
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
void sortParticles(Particle *particles, int particlesCount, int *&cellStart, int *&cellEnd, float cellSize, float xLen,
                   float yLen, float zLen, int gridDimX, int gridDimY, int gridDimZ) {
  int totalCells = gridDimX * gridDimY * gridDimZ;
  cudaMalloc(&cellStart, totalCells * sizeof(int));
  cudaMalloc(&cellEnd, totalCells * sizeof(int));

  int threads = 128;
  int blocks = (particlesCount + threads - 1) / threads;
  computeCellId<<<blocks, threads>>>(particles, particlesCount, cellSize, xLen, yLen, zLen, gridDimX, gridDimY,
                                     gridDimZ);
  cudaDeviceSynchronize();

  thrust::device_ptr<Particle> dev_ptr(particles);
  thrust::sort(dev_ptr, dev_ptr + particlesCount, ParticleComparator());

  findCellStartEnd<<<1, 1>>>(particles, particlesCount, cellStart, cellEnd, totalCells);
  cudaDeviceSynchronize();
}

// Normalize mass based on the density of particles in the sink
float normalizeMass(Particle *particles, int particlesCount, const Sink &sink) {
  float mass = 1.0f;
  int blockDim = 32;
  int gridDim = (particlesCount + (blockDim - 1)) / blockDim;

  // Use Sink dimensions as simulation domain
  float xLen = sink.xLen;
  float yLen = sink.yLen;
  float zLen = sink.zLen;

  float cellSize = KERNEL_RADIUS;
  int gridDimX = (int)ceil(xLen / cellSize);
  int gridDimY = (int)ceil(yLen / cellSize);
  int gridDimZ = (int)ceil(zLen / cellSize);

  int *cellStart, *cellEnd;
  sortParticles(particles, particlesCount, cellStart, cellEnd, cellSize, xLen, yLen, zLen, gridDimX, gridDimY,
                gridDimZ);

  computeDensitySorted<<<gridDim, blockDim>>>(particles, particlesCount, mass, cellStart, cellEnd, cellSize, gridDimX,
                                              gridDimY, gridDimZ, xLen, yLen, zLen);
  cudaDeviceSynchronize();

  cudaError_t err;
  if ((err = cudaGetLastError()) != cudaSuccess)
    std::cerr << "Kernel error (computeDensitySorted): " << cudaGetErrorString(err) << std::endl;

  float rho0 = REST_DENSITY;
  float rho2s = 0.0f;
  float rhos = 0.0f;
  for (int i = 0; i < particlesCount; i++) {
    rho2s += particles[i].density * particles[i].density;
    rhos += particles[i].density;
  }
  mass *= (rho0 * rhos / rho2s);
  std::cout << "Mass: " << mass << std::endl;

  cudaFree(cellStart);
  cudaFree(cellEnd);
  return mass;
}

Particle *placeParticles(int &particlesCount, int &droppingParticlesCount, Sink &sink, Trough &trough) {
  float h = KERNEL_RADIUS;
  float hh = h / 1.3f;  // don't change when not necessary
  std::cout << "hh: " << hh << std::endl;
  int particlesInSink = 0;
  for (float x = -sink.xLen / 2.0f; x <= sink.xLen / 2.0f; x += hh) {
    for (float z = -sink.zLen / 2.0f; z <= sink.zLen / 2.0f; z += hh) {
      for (float y = -sink.yLen / 2.0f; y <= 0; y += hh) {
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
  droppingParticlesCount = dropping;
  std::cout << "dropping particle count: " << droppingParticlesCount << std::endl;
  particlesCount = particlesInSink + droppingParticlesCount;
  std::cout << "Particle Count: " << particlesCount << std::endl;
  Particle *particles;
  cudaMallocManaged(&particles, particlesCount * sizeof(Particle));
  int count = 0;
  // particles generated in sink
  for (float x = -sink.xLen / 2.0f; x <= sink.xLen / 2.0f; x += hh) {
    for (float z = -sink.zLen / 2.0f; z <= sink.zLen / 2.0f; z += hh) {
      for (float y = -sink.yLen / 2.0f; y <= 0; y += hh) {
        particles[count].position = {x, y, z};
        particles[count].density = 0.0f;
        particles[count].inSink = true;
        particles[count].velocity = {0.0f, 0.0f, 0.0f};
        particles[count].velocityHalf = {0.0f, 0.0f, 0.0f};
        particles[count].acceleration = {0.0f, 0.0f, 0.0f};
        count++;
      }
    }
  }
  // particles that will fall on trough
  float vx = 10.0f;
  float vy = vx * trough.slope;
  for (float x = trough.vertices[0].x + 0.001f; x <= trough.vertices[1].x - 0.001f; x += hh) {
    float y0 = trough.slope * x + trough.intercept;
    for (float z = -trough.zLen / 2.0f + 0.001f; z <= trough.zLen / 2.0f - 0.001f; z += hh) {
      for (float y = y0; y <= y0 + trough.yLen; y += hh) {
        particles[count].position = {x, y, z};
        particles[count].density = 0.0f;
        particles[count].inSink = false;
        particles[count].velocity = {vx, vy, 0.0f};
        particles[count].velocityHalf = {vx, vy, 0.0f};
        particles[count].acceleration = {0.0f, 0.0f, 0.0f};
        count++;
      }
    }
  }
  assert(count == particlesCount);
  return particles;
}

Particle *initParticles(int &particlesCount, float &mass, Sink &sink, Trough &trough) {
  int droppingParticlesCount = 0;
  Particle *particles = placeParticles(particlesCount, droppingParticlesCount, sink, trough);
  mass = normalizeMass(particles, particlesCount - droppingParticlesCount, sink);
  return particles;
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
    particle.velocityHalf.x = -particle.velocityHalf.x;
    particle.velocity.x *= REFLECT_DAMP;
    particle.velocity.y *= REFLECT_DAMP;
    particle.velocity.z *= REFLECT_DAMP;
    particle.velocityHalf.x *= REFLECT_DAMP;
    particle.velocityHalf.y *= REFLECT_DAMP;
    particle.velocityHalf.z *= REFLECT_DAMP;
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
    particle.velocityHalf.y = -particle.velocityHalf.y;
    particle.velocity.x *= REFLECT_DAMP;
    particle.velocity.y *= REFLECT_DAMP;
    particle.velocity.z *= REFLECT_DAMP;
    particle.velocityHalf.x *= REFLECT_DAMP;
    particle.velocityHalf.y *= REFLECT_DAMP;
    particle.velocityHalf.z *= REFLECT_DAMP;
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
    particle.velocityHalf.z = -particle.velocityHalf.z;
    particle.velocity.x *= REFLECT_DAMP;
    particle.velocity.y *= REFLECT_DAMP;
    particle.velocity.z *= REFLECT_DAMP;
    particle.velocityHalf.x *= REFLECT_DAMP;
    particle.velocityHalf.y *= REFLECT_DAMP;
    particle.velocityHalf.z *= REFLECT_DAMP;
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
    float dotVH =
        particle.velocityHalf.x * normal.x + particle.velocityHalf.y * normal.y + particle.velocityHalf.z * normal.z;
    float newVHx = particle.velocityHalf.x - 2 * dotVH * normal.x;
    float newVHy = particle.velocityHalf.y - 2 * dotVH * normal.y;
    float newVHz = particle.velocityHalf.z - 2 * dotVH * normal.z;
    newVHx *= REFLECT_DAMP;
    newVHy *= REFLECT_DAMP;
    newVHz *= REFLECT_DAMP;
    particle.velocityHalf.x = newVHx;
    particle.velocityHalf.y = newVHy;
    particle.velocityHalf.z = newVHz;
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
    particle.velocityHalf.z = -particle.velocityHalf.z;
    particle.velocity.x *= REFLECT_DAMP;
    particle.velocity.y *= REFLECT_DAMP;
    particle.velocity.z *= REFLECT_DAMP;
    particle.velocityHalf.x *= REFLECT_DAMP;
    particle.velocityHalf.y *= REFLECT_DAMP;
    particle.velocityHalf.z *= REFLECT_DAMP;
  }
}

__global__ void leapfrogStart(Particle *particles, int particlesCount, float sinkXLen, float sinkYLen, float sinkZLen,
                              float troughZLen, float slope, float intercept, Vec3 normal) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particlesCount) {
    return;
  }
  Particle &particle = particles[i];
  particle.velocityHalf.x = particle.velocity.x + particle.acceleration.x * DELTA_T / 2;
  particle.velocityHalf.y = particle.velocity.y + particle.acceleration.y * DELTA_T / 2;
  particle.velocityHalf.z = particle.velocity.z + particle.acceleration.z * DELTA_T / 2;
  particle.velocity.x = particle.velocity.x + particle.acceleration.x * DELTA_T;
  particle.velocity.y = particle.velocity.y + particle.acceleration.y * DELTA_T;
  particle.velocity.z = particle.velocity.z + particle.acceleration.z * DELTA_T;
  particle.position.x = particle.position.x + particle.velocityHalf.x * DELTA_T;
  particle.position.y = particle.position.y + particle.velocityHalf.y * DELTA_T;
  particle.position.z = particle.position.z + particle.velocityHalf.z * DELTA_T;
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
}

__global__ void leapfrogStep(Particle *particles, int particlesCount, float sinkXLen, float sinkYLen, float sinkZLen,
                             float troughZLen, float slope, float intercept, Vec3 normal) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particlesCount) {
    return;
  }
  Particle &particle = particles[i];
  particle.velocityHalf.x = particle.velocityHalf.x + particle.acceleration.x * DELTA_T;
  particle.velocityHalf.y = particle.velocityHalf.y + particle.acceleration.y * DELTA_T;
  particle.velocityHalf.z = particle.velocityHalf.z + particle.acceleration.z * DELTA_T;
  particle.velocity.x = particle.velocityHalf.x + particle.acceleration.x * DELTA_T / 2;
  particle.velocity.y = particle.velocityHalf.y + particle.acceleration.y * DELTA_T / 2;
  particle.velocity.z = particle.velocityHalf.z + particle.acceleration.z * DELTA_T / 2;
  particle.position.x = particle.position.x + particle.velocityHalf.x * DELTA_T;
  particle.position.y = particle.position.y + particle.velocityHalf.y * DELTA_T;
  particle.position.z = particle.position.z + particle.velocityHalf.z * DELTA_T;
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
}

__global__ void coordTransform(Particle *particles, int particlesCount, float *transformMat) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= particlesCount) {
    return;
  }
  float worldPos[4], result[4];
  worldPos[0] = particles[i].position.x;
  worldPos[1] = particles[i].position.y;
  worldPos[2] = particles[i].position.z;
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
    particles[i].screenPos.x = -1.0f;
    particles[i].screenPos.y = -1.0f;
  } else {
    float screenX = fmaxf(0.0f, fminf(1.0f, (x + 1.0f) * 0.5f)) * SCREEN_WIDTH;
    float screenY = fmaxf(0.0f, fminf(1.0f, (1.0f - y) * 0.5f)) * SCREEN_HEIGHT;
    particles[i].screenPos.x = screenX;
    particles[i].screenPos.y = screenY;
  }
}

void initSimulation(Particle *particles, int particlesCount, const Sink &sink, const Trough &trough, float mass,
                    float *transformMat) {
  int blockDim = 32;
  int gridDim = (particlesCount + (blockDim - 1)) / blockDim;

  // Set up uniform grid parameters using the sink dimensions
  float cellSize = KERNEL_RADIUS;
  float xLen = sink.xLen;
  float yLen = sink.yLen;
  float zLen = sink.zLen;
  int gridDimX = (int)ceil(xLen / cellSize);
  int gridDimY = (int)ceil(yLen / cellSize);
  int gridDimZ = (int)ceil(zLen / cellSize);

  int *cellStart, *cellEnd;
  sortParticles(particles, particlesCount, cellStart, cellEnd, cellSize, xLen, yLen, zLen, gridDimX, gridDimY,
                gridDimZ);

  cudaError_t err;
  computeDensitySorted<<<gridDim, blockDim>>>(particles, particlesCount, mass, cellStart, cellEnd, cellSize, gridDimX,
                                              gridDimY, gridDimZ, xLen, yLen, zLen);
  cudaDeviceSynchronize();
  if ((err = cudaGetLastError()) != cudaSuccess)
    std::cerr << "Kernel error (computeDensitySorted): " << cudaGetErrorString(err) << std::endl;

  computeAccelSorted<<<gridDim, blockDim>>>(particles, particlesCount, mass, cellStart, cellEnd, cellSize, gridDimX,
                                            gridDimY, gridDimZ, xLen, yLen, zLen);
  if ((err = cudaGetLastError()) != cudaSuccess)
    std::cerr << "Kernel error (computeAccelSorted): " << cudaGetErrorString(err) << std::endl;

  leapfrogStart<<<gridDim, blockDim>>>(particles, particlesCount, sink.xLen, sink.yLen, sink.zLen, trough.zLen,
                                       trough.slope, trough.intercept, trough.normal);
  coordTransform<<<gridDim, blockDim>>>(particles, particlesCount, transformMat);
  cudaDeviceSynchronize();
  if ((err = cudaGetLastError()) != cudaSuccess)
    std::cerr << "Kernel error (leapfrogStart or coordTransform): " << cudaGetErrorString(err) << std::endl;

  cudaFree(cellStart);
  cudaFree(cellEnd);
}

void updateSimulation(Particle *particles, int particlesCount, const Sink &sink, const Trough &trough, float mass,
                      float *transformMat) {
  int blockDim = 32;
  int gridDim = (particlesCount + (blockDim - 1)) / blockDim;

  float cellSize = KERNEL_RADIUS;
  float xLen = sink.xLen;
  float yLen = sink.yLen;
  float zLen = sink.zLen;
  int gridDimX = (int)ceil(xLen / cellSize);
  int gridDimY = (int)ceil(yLen / cellSize);
  int gridDimZ = (int)ceil(zLen / cellSize);

  int *cellStart, *cellEnd;
  sortParticles(particles, particlesCount, cellStart, cellEnd, cellSize, xLen, yLen, zLen, gridDimX, gridDimY,
                gridDimZ);

  cudaError_t err;
  computeDensitySorted<<<gridDim, blockDim>>>(particles, particlesCount, mass, cellStart, cellEnd, cellSize, gridDimX,
                                              gridDimY, gridDimZ, xLen, yLen, zLen);
  cudaDeviceSynchronize();
  if ((err = cudaGetLastError()) != cudaSuccess)
    std::cerr << "Kernel error (computeDensitySorted): " << cudaGetErrorString(err) << std::endl;

  computeAccelSorted<<<gridDim, blockDim>>>(particles, particlesCount, mass, cellStart, cellEnd, cellSize, gridDimX,
                                            gridDimY, gridDimZ, xLen, yLen, zLen);
  if ((err = cudaGetLastError()) != cudaSuccess)
    std::cerr << "Kernel error (computeAccelSorted): " << cudaGetErrorString(err) << std::endl;

  leapfrogStep<<<gridDim, blockDim>>>(particles, particlesCount, sink.xLen, sink.yLen, sink.zLen, trough.zLen,
                                      trough.slope, trough.intercept, trough.normal);
  coordTransform<<<gridDim, blockDim>>>(particles, particlesCount, transformMat);
  cudaDeviceSynchronize();
  if ((err = cudaGetLastError()) != cudaSuccess)
    std::cerr << "Kernel error (leapfrogStep or coordTransform): " << cudaGetErrorString(err) << std::endl;

  cudaFree(cellStart);
  cudaFree(cellEnd);
}
