#pragma once

#include "sph.h"  // Includes ParticlesSoA definition

// kernel declarations

/**
 * Assign a spatial cell ID to each particle based on its 3D position (SoA version).
 */
__global__ void computeCellIdSoA(int particleCount, float *posX, float *posY, float *posZ, int *cellId, float cellSize,
                                 float xLen, float yLen, float zLen, int gridDimX, int gridDimY, int gridDimZ);

/**
 * Initialize cell start and end arrays to -1, indicating empty cells.
 */
__global__ void initCells(int *cellStart, int *cellEnd, int totalCells);

/**
 * Find the start and end indices for each cell after particle indices have been sorted by cell ID.
 * Uses the sorted particleIndices array.
 */
__global__ void findCellStartEndSoA(int particleCount, const int *particleIndices, const int *cellId, int *cellStart,
                                    int *cellEnd, int totalCells);

/**
 * Compute particle densities and pressures using SoA and sorted indices.
 */
__global__ void computeDensityPressureSoAGlobal(int particleCount, const int *particleIndices, float *posX, float *posY,
                                                float *posZ, float *density, float *pressure, float mass,
                                                const int *cellStart, const int *cellEnd, float cellSize, int gridDimX,
                                                int gridDimY, int gridDimZ, float xLen, float yLen, float zLen,
                                                float POLY6, float WEIGHT_AT_0);

// device functions

/**
 * Reflect a particle inside the sink (operates on individual attributes by reference).
 */
__device__ void reflectInSinkSoA(float &px, float &py, float &pz, float &vx, float &vy, float &vz, float xLen,
                                 float yLen, float zLen);

/**
 * Reflect a particle inside the trough (operates on individual attributes by reference).
 */
__device__ void reflectInTroughSoA(float &px, float &py, float &pz, float &vx, float &vy, float &vz, float zLen,
                                   float slope, float intercept, Vec3 normal);

/**
 * Calling density and pressure calculation, force calculation, time integration and coord transforming
 */
__global__ void computeParticlePosition(int particleCount, const int *particleIndices, float *posX, float *posY,
                                        float *posZ, float *velX, float *velY, float *velZ, float *avgVelX,
                                        float *avgVelY, float *avgVelZ, float *accX, float *accY, float *accZ,
                                        float *density, float *pressure, float mass, const int *cellStart,
                                        const int *cellEnd, float cellSize, int gridDimX, int gridDimY, int gridDimZ,
                                        float xLen, float yLen, float zLen, float POLY6, float WEIGHT_AT_0,
                                        float VISCOSITY_LAPLACIAN, char *inSink, float sinkXLen, float sinkYLen,
                                        float sinkZLen, float troughZLen, float slope, float intercept, Vec3 normal,
                                        const float *transformMat, Vec2 *screenPosOnGPU);

// comparator for thrust::sort based on cellId using an index array
struct CompareParticlesByCellId {
  const int *cellId_d;  // Pointer to device memory for cell IDs

  CompareParticlesByCellId(const int *cid_d) : cellId_d(cid_d) {}

  __host__ __device__ bool operator()(const int &a_idx, const int &b_idx) const {
    // Compare the cell IDs of the particles pointed to by the indices
    return cellId_d[a_idx] < cellId_d[b_idx];
  }
};