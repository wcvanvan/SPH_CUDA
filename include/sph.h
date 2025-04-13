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
      : density(0.0f), position(), screenPos(), velocity(), velocityHalf(), acceleration(), cellId(0) {}
  __host__ __device__ Particle(Vec3 position, int id);
  float density = 0.0f;

  Vec3 position;
  Vec2 screenPos;
  Vec3 velocity;
  Vec3 velocityHalf;
  Vec3 acceleration;
  int cellId;
};

/**
 * Create particles
 */
Particle *initParticles(int &particleCount, float &mass, Sink &sink);
void initSimulation(Particle *particles, int particleCount, const Sink &sink, float mass, float *transformMat);
/**
 * Compute interaction and update the particles
 */
void updateSimulation(Particle *particles, int particleCount, const Sink &sink, float mass, float *transformMat);
float *allocateMatOnGPU(Mat4 &mat);
