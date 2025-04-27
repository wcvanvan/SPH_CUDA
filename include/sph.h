#pragma once

#include <cuda_runtime.h>
#include <stdint.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include "const.h"
#include "scene.h"
#include "vec.h"

#ifndef __CUDACC__
#define __host__
#define __device__
#endif

/**
 * Particle data structure
 * This structure is used to store the particle's position, velocity, density, pressure, and other properties.
 */
class Particle {
 public:
  __host__ __device__ Particle()
      : density(0.0f),
        pressure(0.0f),
        inSink(true),
        position(),
        velocity(),
        averageVelocity(),
        acceleration(),
        cellId(0) {}
  __host__ __device__ Particle(Vec3 position, int id);
  float density = 0.0f;
  float pressure = 0.0f;
  bool inSink;

  Vec3 position;
  Vec3 velocity;
  Vec3 averageVelocity;
  Vec3 acceleration;
  int cellId;
};

/**
 * Create particles
 */
Particle *initParticles(int &particleCount, float &mass, Sink &sink, Trough &trough, int *cellStart, int *cellEnd,
                        float POLY6, float WEIGHT_AT_0);

/**
 * Create a buffer to store one frame particle positions on GPU side.
 * This buffer will be copied DtoD to the allFrames buffer.
 */
Vec2 *initScreenPos(int particleCount);

/**
 * Create a buffer to store all frames particle positions on GPU side.
 * Once the simulation is done, we will copy this buffer to CPU side at once.
 * This is more efficient than copying each frame one by one.
 */
Vec2 *initAllFrames(int particleCount, int frameCount);

/**
 * Compute interaction and update the particles
 */
void updateSimulation(Particle *particles, int particleCount, const Sink &sink, const Trough &trough, float mass,
                      float *transformMat, int *cellStart, int *cellEnd, Vec2 *screenPosOnGPU, float POLY6,
                      float VISCOSITY_LAPLACIAN, float WEIGHT_AT_0, int frameCount, Vec2 *allFramesOnGPU,
                      cudaStream_t stream);

float *allocateMatOnGPU(Mat4 &mat);

/*
 * This function will return a pointer to an array of integers
 * that will be used to store the end index of each cell (exclusive) in the cellStart array.
 * cellEnd[i] = the first particle in cell i+1
 */
int *initCellEnd(int totalCells);

/*
 * This function will return a pointer to an array of integers
 * that will be used to store the start index of each cell (inclusive) in the cellEnd array.
 * cellStart[i] = the first particle in cell i
 */
int *initCellStart(int totalCells);

/**
 * This function will copy all frames from GPU to CPU.
 */
void copyAllFramesToCPU(Vec2 *allFramesOnGPU, Vec2 *allFramesOnCPU, int particleCount, int frameCount);

void sortParticles(Particle *particles, int particleCount, int *&cellStart, int *&cellEnd, float cellSize, float xLen,
                   float yLen, float zLen, int gridDimX, int gridDimY, int gridDimZ, cudaStream_t stream);
