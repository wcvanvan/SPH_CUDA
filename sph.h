#pragma once

#include <stdint.h>
#include "const.h"
#include "vec.h"

class Particle
{
public:
    Particle(Vec3 position, int id);
    float mass = 1.0f;
    float density = 0.0f;
    float pressure = 0.0f;

    Vec3 position;
    Vec3 velocity;
    Vec3 evelocity;
    Vec3 acceleration;
};

/**
 * Create particles
 */
 Particle *initParticles(int particleCount);

/**
 * Compute interaction and update the particles
 */
void updateSimulation(Particle *particles, int particleCount);


