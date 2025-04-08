#pragma once

#include <stdint.h>
#include "const.h"
#include "vec.h"
#include "sink.h"

class Particle
{
public:
    Particle(Vec3 position, int id);
    float density = 0.0f;

    Vec3 position;
    Vec3 velocity;
    Vec3 velocityHalf;
    Vec3 acceleration;
};

/**
 * Create particles
 */
 Particle *initParticles(int particleCount, float *mass, Sink &sink);
 void initSimulation(Particle *particles, int particleCount, const Sink &sink, float mass);
/**
 * Compute interaction and update the particles
 */
 void updateSimulation(Particle *particles, int particleCount, const Sink &sink, float mass);

