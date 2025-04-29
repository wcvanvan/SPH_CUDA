#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include "const.h"
#include "kernel.h"
#include "sph.h"

// Helper to check CUDA errors
#define CUDA_CHECK(err)                                                                                        \
  do {                                                                                                         \
    cudaError_t err_ = (err);                                                                                  \
    if (err_ != cudaSuccess) {                                                                                 \
      std::cerr << "CUDA error in " << __FILE__ << " at line " << __LINE__ << ": " << cudaGetErrorString(err_) \
                << std::endl;                                                                                  \
      exit(EXIT_FAILURE);                                                                                      \
    }                                                                                                          \
    err_ = cudaGetLastError();                                                                                 \
    if (err_ != cudaSuccess) {                                                                                 \
      std::cerr << "CUDA kernel launch error in " << __FILE__ << " at line " << __LINE__ << ": "               \
                << cudaGetErrorString(err_) << std::endl;                                                      \
      exit(EXIT_FAILURE);                                                                                      \
    }                                                                                                          \
  } while (0)

float *allocateMatOnGPU(Mat4 &mat) {
  float data[16];
  int count = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      data[count++] = mat[i][j];
    }
  }
  float *matOnGPU;
  CUDA_CHECK(cudaMalloc((void **)&matOnGPU, 16 * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(matOnGPU, data, 16 * sizeof(float), cudaMemcpyHostToDevice));
  return matOnGPU;
}

void cleanupMatOnGPU(float *matOnGPU) {
  if (matOnGPU) CUDA_CHECK(cudaFree(matOnGPU));
}

// Sort particle indices based on their cell IDs and find the start and end indices of each cell.
void sortParticlesAndFindCells(ParticlesSoA &particles, int *cellStart, int *cellEnd, float cellSize, float xLen,
                               float yLen, float zLen, int gridDimX, int gridDimY, int gridDimZ) {
  int particleCount = particles.particleCount;
  if (particleCount == 0) return;

  int totalCells = gridDimX * gridDimY * gridDimZ;
  int threads = 128;
  int blocksParticles = (particleCount + threads - 1) / threads;
  int blocksCells = (totalCells + threads - 1) / threads;

  computeCellIdSoA<<<blocksParticles, threads>>>(particleCount, particles.posX, particles.posY, particles.posZ,
                                                 particles.cellId, cellSize, xLen, yLen, zLen, gridDimX, gridDimY,
                                                 gridDimZ);
  CUDA_CHECK(cudaGetLastError());

  thrust::device_ptr<int> dev_indices_ptr(particles.particleIndices);
  thrust::sort(thrust::device, dev_indices_ptr, dev_indices_ptr + particleCount,
               CompareParticlesByCellId(particles.cellId));
  CUDA_CHECK(cudaGetLastError());

  initCells<<<blocksCells, threads>>>(cellStart, cellEnd, totalCells);
  CUDA_CHECK(cudaGetLastError());

  findCellStartEndSoA<<<blocksParticles, threads>>>(particleCount, particles.particleIndices, particles.cellId,
                                                    cellStart, cellEnd, totalCells);
  CUDA_CHECK(cudaGetLastError());
}

// Place particles initially on the host, then copy to SoA on GPU
void placeParticlesSoA(ParticlesSoA &particles, int &particleCount, int &droppingparticleCount, Sink &sink,
                       Trough &trough) {
  float h = KERNEL_RADIUS;
  float hh = h / 2.0f;

  std::vector<float> posX_h, posY_h, posZ_h;
  std::vector<float> velX_h, velY_h, velZ_h;
  std::vector<char> inSink_h;

  int particlesInSinkEst = 0;
  for (float x = -sink.xLen / 2.0f; x <= sink.xLen / 2.0f; x += hh)
    for (float z = -sink.zLen / 2.0f; z <= sink.zLen / 2.0f; z += hh)
      for (float y = -sink.yLen / 2.0f; y <= sink.yLen / 4.0f; y += hh) particlesInSinkEst++;

  int droppingEst = 0;
  for (float x = trough.vertices[0].x + 0.001f; x <= trough.vertices[1].x - 0.001f; x += hh) {
    float y0 = trough.slope * x + trough.intercept;
    for (float z = -trough.zLen / 2.0f + 0.001f; z <= trough.zLen / 2.0f - 0.001f; z += hh)
      for (float y = y0; y <= y0 + trough.yLen; y += hh) droppingEst++;
  }

  posX_h.reserve(particlesInSinkEst + droppingEst);

  // Particles generated in sink
  for (float x = -sink.xLen / 2.0f; x <= sink.xLen / 2.0f; x += hh) {
    for (float z = -sink.zLen / 2.0f; z <= sink.zLen / 2.0f; z += hh) {
      for (float y = -sink.yLen / 2.0f; y <= sink.yLen / 4.0f; y += hh) {
        posX_h.push_back(x);
        posY_h.push_back(y);
        posZ_h.push_back(z);
        velX_h.push_back(0.0f);
        velY_h.push_back(0.0f);
        velZ_h.push_back(0.0f);
        inSink_h.push_back(1);
      }
    }
  }
  int particlesInSinkActual = posX_h.size();

  float vx_drop = 2.5f;
  float vy_drop = vx_drop * trough.slope;
  for (float x = trough.vertices[0].x + 0.001f; x <= trough.vertices[1].x - 0.001f; x += hh) {
    float y0 = trough.slope * x + trough.intercept;
    for (float z = -trough.zLen / 2.0f + 0.001f; z <= trough.zLen / 2.0f - 0.001f; z += hh) {
      for (float y = y0; y <= y0 + trough.yLen; y += hh) {
        posX_h.push_back(x);
        posY_h.push_back(y);
        posZ_h.push_back(z);
        velX_h.push_back(vx_drop);
        velY_h.push_back(vy_drop);
        velZ_h.push_back(0.0f);
        inSink_h.push_back(0);
      }
    }
  }

  particleCount = posX_h.size();
  droppingparticleCount = particleCount - particlesInSinkActual;
  particles.particleCount = particleCount;
  std::cout << "Particle Count: " << particleCount << " (Sink: " << particlesInSinkActual
            << ", Dropping: " << droppingparticleCount << ")" << std::endl;

  if (particleCount == 0) return;

  size_t floatSize = particleCount * sizeof(float);
  size_t intSize = particleCount * sizeof(int);
  size_t charSize = particleCount * sizeof(char);

  CUDA_CHECK(cudaMalloc(&particles.posX, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.posY, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.posZ, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.velX, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.velY, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.velZ, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.avgVelX, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.avgVelY, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.avgVelZ, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.accX, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.accY, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.accZ, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.density, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.pressure, floatSize));
  CUDA_CHECK(cudaMalloc(&particles.inSink, charSize));
  CUDA_CHECK(cudaMalloc(&particles.cellId, intSize));
  CUDA_CHECK(cudaMalloc(&particles.particleIndices, intSize));

  CUDA_CHECK(cudaMemcpy(particles.posX, posX_h.data(), floatSize, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(particles.posY, posY_h.data(), floatSize, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(particles.posZ, posZ_h.data(), floatSize, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(particles.velX, velX_h.data(), floatSize, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(particles.velY, velY_h.data(), floatSize, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(particles.velZ, velZ_h.data(), floatSize, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(particles.inSink, inSink_h.data(), charSize, cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMemset(particles.avgVelX, 0, floatSize));
  CUDA_CHECK(cudaMemset(particles.avgVelY, 0, floatSize));
  CUDA_CHECK(cudaMemset(particles.avgVelZ, 0, floatSize));
  CUDA_CHECK(cudaMemcpy(particles.avgVelX, particles.velX, floatSize, cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy(particles.avgVelY, particles.velY, floatSize, cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy(particles.avgVelZ, particles.velZ, floatSize, cudaMemcpyDeviceToDevice));

  CUDA_CHECK(cudaMemset(particles.accX, 0, floatSize));
  CUDA_CHECK(cudaMemset(particles.accY, 0, floatSize));
  CUDA_CHECK(cudaMemset(particles.accZ, 0, floatSize));
  CUDA_CHECK(cudaMemset(particles.density, 0, floatSize));
  CUDA_CHECK(cudaMemset(particles.pressure, 0, floatSize));
  CUDA_CHECK(cudaMemset(particles.cellId, 0, intSize));

  thrust::device_ptr<int> dev_indices_ptr(particles.particleIndices);
  thrust::sequence(thrust::device, dev_indices_ptr, dev_indices_ptr + particleCount);
}

int *initCellStart(int totalCells) {
  int *cellStart;
  CUDA_CHECK(cudaMalloc(&cellStart, totalCells * sizeof(int)));
  return cellStart;
}

void cleanupCellStart(int *cellStart) {
  if (cellStart) CUDA_CHECK(cudaFree(cellStart));
}

int *initCellEnd(int totalCells) {
  int *cellEnd;
  CUDA_CHECK(cudaMalloc(&cellEnd, totalCells * sizeof(int)));
  return cellEnd;
}

void cleanupCellEnd(int *cellEnd) {
  if (cellEnd) CUDA_CHECK(cudaFree(cellEnd));
}

float normalizeMassSoA(ParticlesSoA &particles, int particlesInSinkCount, const Sink &sink, int *cellStart,
                       int *cellEnd, float POLY6, float WEIGHT_AT_0) {
  if (particlesInSinkCount == 0) {
    std::cerr << "Warning: No particles in sink for mass normalization. Using default mass 1.0." << std::endl;
    return 1.0f;
  }

  float mass = 1.0f;
  int particleCount = particles.particleCount;
  int blockDim = 128;
  int gridDim = (particleCount + (blockDim - 1)) / blockDim;

  float cellSize = KERNEL_RADIUS;
  int gridDimX = (int)ceilf(sink.xLen / cellSize);
  int gridDimY = (int)ceilf(sink.yLen / cellSize);
  int gridDimZ = (int)ceilf(sink.zLen / cellSize);

  sortParticlesAndFindCells(particles, cellStart, cellEnd, cellSize, sink.xLen, sink.yLen, sink.zLen, gridDimX,
                            gridDimY, gridDimZ);

  computeDensityPressureSoA<<<gridDim, blockDim>>>(
      particleCount, particles.particleIndices, particles.posX, particles.posY, particles.posZ, particles.density,
      particles.pressure, mass, cellStart, cellEnd, cellSize, gridDimX, gridDimY, gridDimZ, sink.xLen, sink.yLen,
      sink.zLen, POLY6, WEIGHT_AT_0);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

  std::vector<float> density_h(particleCount);
  std::vector<char> inSink_h(particleCount);
  CUDA_CHECK(cudaMemcpy(density_h.data(), particles.density, particleCount * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(inSink_h.data(), particles.inSink, particleCount * sizeof(char), cudaMemcpyDeviceToHost));

  float rho0 = REST_DENSITY;
  double rho2s_sum = 0.0;
  double rhos_sum = 0.0;
  int sink_count_for_avg = 0;

  for (int i = 0; i < particlesInSinkCount; ++i) {
    if (density_h[i] > 1e-6) {
      rho2s_sum += (double)density_h[i] * density_h[i];
      rhos_sum += (double)density_h[i];
      sink_count_for_avg++;
    }
  }

  if (sink_count_for_avg > 0 && rho2s_sum > 1e-12) {
    mass *= (rho0 * rhos_sum / rho2s_sum);
    std::cout << "Normalized Mass: " << mass << std::endl;
  } else {
    std::cerr << "Warning: Could not normalize mass (zero density sum or count). Using default mass 1.0." << std::endl;
    mass = 1.0f;
  }

  return mass;
}

void initParticlesSoA(ParticlesSoA &particles, int &particleCount, float &mass, Sink &sink, Trough &trough,
                      int *&cellStart, int *&cellEnd, float POLY6, float WEIGHT_AT_0) {
  int droppingparticleCount = 0;
  int particlesInSinkCount = 0;

  placeParticlesSoA(particles, particleCount, droppingparticleCount, sink, trough);
  particlesInSinkCount = particleCount - droppingparticleCount;

  if (particleCount > 0) {
    mass = normalizeMassSoA(particles, particlesInSinkCount, sink, cellStart, cellEnd, POLY6, WEIGHT_AT_0);
  } else {
    mass = 1.0f;
  }
}

void cleanupParticlesSoA(ParticlesSoA &particles) {
  CUDA_CHECK(cudaFree(particles.posX));
  CUDA_CHECK(cudaFree(particles.posY));
  CUDA_CHECK(cudaFree(particles.posZ));
  CUDA_CHECK(cudaFree(particles.velX));
  CUDA_CHECK(cudaFree(particles.velY));
  CUDA_CHECK(cudaFree(particles.velZ));
  CUDA_CHECK(cudaFree(particles.avgVelX));
  CUDA_CHECK(cudaFree(particles.avgVelY));
  CUDA_CHECK(cudaFree(particles.avgVelZ));
  CUDA_CHECK(cudaFree(particles.accX));
  CUDA_CHECK(cudaFree(particles.accY));
  CUDA_CHECK(cudaFree(particles.accZ));
  CUDA_CHECK(cudaFree(particles.density));
  CUDA_CHECK(cudaFree(particles.pressure));
  CUDA_CHECK(cudaFree(particles.inSink));
  CUDA_CHECK(cudaFree(particles.cellId));
  CUDA_CHECK(cudaFree(particles.particleIndices));
  particles.particleCount = 0;
}

Vec2 *initScreenPos(int particleCount) {
  Vec2 *screenPos;
  CUDA_CHECK(cudaMalloc(&screenPos, particleCount * sizeof(Vec2)));
  return screenPos;
}

void cleanupScreenPos(Vec2 *screenPosOnGPU) {
  if (screenPosOnGPU) CUDA_CHECK(cudaFree(screenPosOnGPU));
}

Vec2 *initAllFrames(int particleCount, int frameCount) {
  Vec2 *allFrames;
  size_t totalSize = (size_t)particleCount * frameCount * sizeof(Vec2);
  CUDA_CHECK(cudaMalloc(&allFrames, totalSize));
  return allFrames;
}

void cleanupAllFrames(Vec2 *allFramesOnGPU) {
  if (allFramesOnGPU) CUDA_CHECK(cudaFree(allFramesOnGPU));
}

void copyAllFramesToCPU(Vec2 *allFramesOnGPU, Vec2 *allFramesOnCPU, int particleCount, int frameCount) {
  size_t totalSize = (size_t)particleCount * frameCount * sizeof(Vec2);
  CUDA_CHECK(cudaMemcpy(allFramesOnCPU, allFramesOnGPU, totalSize, cudaMemcpyDeviceToHost));
}

void updateSimulationSoA(ParticlesSoA &particles, const Sink &sink, const Trough &trough, float mass,
                         float *transformMatOnGPU, int *cellStart, int *cellEnd, Vec2 *screenPosOnGPU, float POLY6,
                         float VISCOSITY_LAPLACIAN, float WEIGHT_AT_0, int frameCount, Vec2 *allFramesOnGPU) {
  int particleCount = particles.particleCount;
  if (particleCount == 0) return;

  int blockDim = 128;
  int gridDim = (particleCount + (blockDim - 1)) / blockDim;

  float cellSize = KERNEL_RADIUS;
  float xLen = sink.xLen;
  float yLen = sink.yLen;
  float zLen = sink.zLen;
  int gridDimX = (int)ceilf(xLen / cellSize);
  int gridDimY = (int)ceilf(yLen / cellSize);
  int gridDimZ = (int)ceilf(zLen / cellSize);

  sortParticlesAndFindCells(particles, cellStart, cellEnd, cellSize, xLen, yLen, zLen, gridDimX, gridDimY, gridDimZ);

  computeDensityPressureSoA<<<gridDim, blockDim>>>(particleCount, particles.particleIndices, particles.posX,
                                                   particles.posY, particles.posZ, particles.density,
                                                   particles.pressure, mass, cellStart, cellEnd, cellSize, gridDimX,
                                                   gridDimY, gridDimZ, xLen, yLen, zLen, POLY6, WEIGHT_AT_0);
  CUDA_CHECK(cudaGetLastError());

  computeAccelSoA<<<gridDim, blockDim>>>(
      particleCount, particles.particleIndices, particles.posX, particles.posY, particles.posZ, particles.velX,
      particles.velY, particles.velZ, particles.avgVelX, particles.avgVelY, particles.avgVelZ, particles.accX,
      particles.accY, particles.accZ, particles.density, particles.pressure, mass, cellStart, cellEnd, cellSize,
      gridDimX, gridDimY, gridDimZ, xLen, yLen, zLen, VISCOSITY_LAPLACIAN);
  CUDA_CHECK(cudaGetLastError());

  integrationSoA<<<gridDim, blockDim>>>(particleCount, particles.particleIndices, particles.posX, particles.posY,
                                        particles.posZ, particles.velX, particles.velY, particles.velZ,
                                        particles.avgVelX, particles.avgVelY, particles.avgVelZ, particles.accX,
                                        particles.accY, particles.accZ, particles.inSink, sink.xLen, sink.yLen,
                                        sink.zLen, trough.zLen, trough.slope, trough.intercept, trough.normal);
  CUDA_CHECK(cudaGetLastError());

  coordTransformSoA<<<gridDim, blockDim>>>(particleCount, particles.posX, particles.posY, particles.posZ,
                                           transformMatOnGPU, screenPosOnGPU);
  CUDA_CHECK(cudaGetLastError());

  size_t frameOffset = (size_t)frameCount * particleCount;
  size_t frameSizeBytes = particleCount * sizeof(Vec2);
  CUDA_CHECK(cudaMemcpy(allFramesOnGPU + frameOffset, screenPosOnGPU, frameSizeBytes, cudaMemcpyDeviceToDevice));
}