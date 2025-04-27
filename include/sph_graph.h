// sph_graph.h
#pragma once
#include "sph.h"

void simulateAllFramesCuda(Particle *particlesOnGPU, int particleCount, const Sink &sink, const Trough &trough,
                           float mass, float *transformMatOnGPU, int *cellStart, int *cellEnd, Vec2 *screenPosOnGPU,
                           Vec2 *screenPosOnCPU, float POLY6, float VISCOSITY_LAPLACIAN, float WEIGHT_AT_0,
                           Vec2 *allFramesOnGPU);
