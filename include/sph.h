#pragma once

#include <stdint.h>
#include "const.h"
#include "scene.h"
#include "vec.h"

#ifndef __CUDACC__
#define __host__
#define __device__
#endif

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
Particle *initParticles(int &particleCount, float &mass, Sink &sink, Trough &trough, int *cellStart, int *cellEnd);

Vec2 *initScreenPos(int particleCount);
/**
 * Compute interaction and update the particles
 */
void updateSimulation(Particle *particles, int particleCount, const Sink &sink, const Trough &trough, float mass,
                      float *transformMat, int *cellStart, int *cellEnd, Vec2 *screenPosOnGPU, Vec2 *screenPosOnCPU);
float *allocateMatOnGPU(Mat4 &mat);
int *initCellEnd(int totalCells);
int *initCellStart(int totalCells);
