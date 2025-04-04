#include <iostream>
#include "sph.h"

Particle *SPHInit()
{
    Particle *particles;
    cudaMallocManaged(&particles, MAX_PARTICLES * sizeof(Particle));
    return particles;
}

void updateSimulation(Particle *particles, int particleCount)
{
    int blockDim = 32;
    int gridDim = (particleCount + (blockDim - 1)) / blockDim;
    cudaError_t err;
    // TODO: call computation kernels here
	// cudaDeviceSynchronize();

    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error: " << cudaGetErrorString(err) << std::endl;
    }
}
