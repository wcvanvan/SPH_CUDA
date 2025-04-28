#pragma once

#include <stdint.h>
#include "const.h"
#include "scene.h"
#include "vec.h"

#ifndef __CUDACC__
#define __host__
#define __device__
#endif

/**
 * Structure of Arrays (SoA) for particle data on the GPU.
 */
struct ParticlesSoA {
  int particleCount = 0;
  // Position
  float *posX = nullptr;
  float *posY = nullptr;
  float *posZ = nullptr;
  // Velocity
  float *velX = nullptr;
  float *velY = nullptr;
  float *velZ = nullptr;
  // Average Velocity (for viscosity)
  float *avgVelX = nullptr;
  float *avgVelY = nullptr;
  float *avgVelZ = nullptr;
  // Acceleration
  float *accX = nullptr;
  float *accY = nullptr;
  float *accZ = nullptr;
  // SPH properties
  float *density = nullptr;
  float *pressure = nullptr;
  // State and Grid
  char *inSink = nullptr;  // Changed from bool to char
  int *cellId = nullptr;
  // Sorting index map
  int *particleIndices =
      nullptr;  // Stores sorted order: particleIndices[i] is the original index of the i-th sorted particle
};

/**
 * Allocate and initialize particle data structures on the GPU using SoA layout.
 * Calculates normalized particle mass.
 */
void initParticlesSoA(ParticlesSoA &particles, int &particleCount, float &mass, Sink &sink, Trough &trough,
                      int *&cellStart, int *&cellEnd, float POLY6, float WEIGHT_AT_0);

/**
 * Free all device memory associated with ParticlesSoA.
 */
void cleanupParticlesSoA(ParticlesSoA &particles);

/**
 * Create a buffer to store one frame particle positions on GPU side.
 */
Vec2 *initScreenPos(int particleCount);

/**
 * Free the screen position buffer on GPU.
 */
void cleanupScreenPos(Vec2 *screenPosOnGPU);

/**
 * Create a buffer to store all frames particle positions on GPU side.
 */
Vec2 *initAllFrames(int particleCount, int frameCount);

/**
 * Free the all-frames buffer on GPU.
 */
void cleanupAllFrames(Vec2 *allFramesOnGPU);

/**
 * Compute interaction and update the particles for one simulation step using SoA.
 */
void updateSimulationSoA(ParticlesSoA &particles, const Sink &sink, const Trough &trough, float mass,
                         float *transformMatOnGPU, int *cellStart, int *cellEnd, Vec2 *screenPosOnGPU, float POLY6,
                         float VISCOSITY_LAPLACIAN, float WEIGHT_AT_0, int frameCount, Vec2 *allFramesOnGPU);

/**
 * Allocate a 4x4 matrix on the GPU.
 */
float *allocateMatOnGPU(Mat4 &mat);

/**
 * Free the 4x4 matrix on the GPU.
 */
void cleanupMatOnGPU(float *matOnGPU);

/**
 * Allocate cell start index array on GPU.
 */
int *initCellStart(int totalCells);

/**
 * Free cell start index array on GPU.
 */
void cleanupCellStart(int *cellStart);

/**
 * Allocate cell end index array on GPU.
 */
int *initCellEnd(int totalCells);

/**
 * Free cell end index array on GPU.
 */
void cleanupCellEnd(int *cellEnd);

/**
 * Copy all frames (Vec2 screen positions) from GPU to CPU.
 */
void copyAllFramesToCPU(Vec2 *allFramesOnGPU, Vec2 *allFramesOnCPU, int particleCount, int frameCount);