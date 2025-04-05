#include <iostream>
#include "sph.h"

Particle *initParticles(int particleCount)
{
    Particle *particles;
    cudaMallocManaged(&particles, particleCount * sizeof(Particle));
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
