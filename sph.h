#pragma once

#include <stdint.h>
#include "const.h"

class vec2
{
public:
    float x = 0, y = 0;
};

class Particle
{
public:
    Particle(vec2 position, int id);
    float mass = 1.0f;
    float density = 0.0f;
    float pressure = 0.0f;

    vec2 position;
    vec2 velocity;
    vec2 evelocity;
    vec2 acceleration;
};

/**
 * Create particles
 */
 Particle *SPHInit();

/**
 * Compute interaction and update the particles
 */
void updateSimulation(Particle *particles, int particleCount);


