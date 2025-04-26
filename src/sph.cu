#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "const.h"
#include "kernel.h"
#include "sph.h"

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
float normalizeMass(Particle *particles, int particleCount, const Sink &sink, int *cellStart, int *cellEnd, float POLY6,
                    float WEIGHT_AT_0) {
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
                                                      gridDimX, gridDimY, gridDimZ, xLen, yLen, zLen, POLY6,
                                                      WEIGHT_AT_0);
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
  mass = normalizeMass(particlesOnGPU, particleCount - droppingparticleCount, sink, cellStart, cellEnd, POLY6,
                       WEIGHT_AT_0);
  return particlesOnGPU;
}

Vec2 *initScreenPos(int particleCount) {
  Vec2 *screenPos;
  cudaMalloc(&screenPos, particleCount * sizeof(Vec2));
  return screenPos;
}

Vec2 *initAllFrames(int particleCount, int frameCount) {
  Vec2 *allFrames;
  cudaMalloc(&allFrames, sizeof(Vec2) * particleCount * frameCount);
  return allFrames;
}

void copyAllFramesToCPU(Vec2 *allFramesOnGPU, Vec2 *allFramesOnCPU, int particleCount, int frameCount) {
  cudaError_t copyErr =
      cudaMemcpy(allFramesOnCPU, allFramesOnGPU, sizeof(Vec2) * particleCount * frameCount, cudaMemcpyDeviceToHost);
  if (copyErr != cudaSuccess) {
    std::cerr << "cudaMemcpy Error: " << cudaGetErrorString(copyErr) << std::endl;
  }
}

void updateSimulation(Particle *particles, int particleCount, const Sink &sink, const Trough &trough, float mass,
                      float *transformMat, int *cellStart, int *cellEnd, Vec2 *screenPosOnGPU, Vec2 *screenPosOnCPU,
                      float POLY6, float VISCOSITY_LAPLACIAN, float WEIGHT_AT_0, int frameCount, Vec2 *allFramesOnGPU) {
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
                                                      gridDimX, gridDimY, gridDimZ, xLen, yLen, zLen, POLY6,
                                                      WEIGHT_AT_0);
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
  cudaError_t copyErr = cudaMemcpy(allFramesOnGPU + frameCount * particleCount, screenPosOnGPU,
                                   sizeof(Vec2) * particleCount, cudaMemcpyDeviceToDevice);
  if (copyErr != cudaSuccess) {
    std::cerr << "cudaMemcpy Error: " << cudaGetErrorString(copyErr) << std::endl;
  }
}
