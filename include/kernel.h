#pragma once

#include "sph.h"

// kernel declarations

/**
 * Assign a spatial cell ID to each particle based on its 3D position.
 * Used for spatial partitioning to speed up neighbor search.
 */
__global__ void computeCellId(Particle *particles, int particleCount, float cellSize, float xLen, float yLen,
                              float zLen, int gridDimX, int gridDimY, int gridDimZ);

/**
 * Initialize cell start and end arrays to -1, indicating empty cells.
 */
__global__ void initCells(int *cellStart, int *cellEnd, int totalCells);

/**
 * Find the start and end indices for each cell after particles have been sorted by cell ID.
 * Assumes particles are sorted by cellId.
 */
__global__ void findCellStartEnd(Particle *particles, int particleCount, int *cellStart, int *cellEnd, int totalCells);

/**
 * Compute particle densities and pressures based on neighboring particles within the smoothing kernel radius.
 * Each thread processes one particle by looping over neighboring cells.
 */
__global__ void computeDensityPressureSorted(Particle *particles, int particleCount, float mass, int *cellStart,
                                             int *cellEnd, float cellSize, int gridDimX, int gridDimY, int gridDimZ,
                                             float xLen, float yLen, float zLen, float POLY6, float WEIGHT_AT_0);

/**
 * Compute particle accelerations due to pressure and viscosity forces.
 * Each thread processes one particle by iterating over neighboring cells.
 */
__global__ void computeAccelSorted(Particle *particles, int particleCount, float mass, int *cellStart, int *cellEnd,
                                   float cellSize, int gridDimX, int gridDimY, int gridDimZ, float xLen, float yLen,
                                   float zLen, float VISCOSITY_LAPLACIAN);

/**
 * Update particle velocities and positions using explicit time integration.
 * Handle collision reflections against sink and trough surfaces.
 */
__global__ void integration(Particle *particles, int particleCount, float sinkXLen, float sinkYLen, float sinkZLen,
                            float troughZLen, float slope, float intercept, Vec3 normal);

/**
 * Transform 3D particle positions into 2D screen space coordinates.
 * Apply a 4x4 transformation matrix and perform homogeneous division.
 */
__global__ void coordTransform(Particle *particles, int particleCount, float *transformMat, Vec2 *screenPosOnGPU);

// device functions

/**
 * Reflect a particle inside the sink when it hits the sink walls.
 * Applies damping to the reflected velocity.
 */
__device__ void reflectInSink(Particle &particle, float xLen, float yLen, float zLen);

/**
 * Reflect a particle inside the trough when it hits the trough surface or side walls.
 * Applies damping to the reflected velocity.
 */
__device__ void reflectInTrough(Particle &particle, float zLen, float slope, float intercept, Vec3 normal);

// comparator
/**
 * Comparator used by thrust::sort to sort particles based on their cell ID.
 */
struct ParticleComparator {
  __host__ __device__ bool operator()(const Particle &a, const Particle &b) const { return a.cellId < b.cellId; }
};
