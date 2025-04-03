#include <iostream>
#include "sph.h"

uint32_t *gpuAlloc()
{
    uint32_t *host_ptr;

    cudaError_t err = cudaHostAlloc(&host_ptr, SCREEN_SIZE * sizeof(uint32_t), cudaHostAllocMapped);
    if (err != cudaSuccess)
    {
        std::cerr << "cudaHostAlloc Error " << cudaGetErrorString(err) << std::endl;
    }
    return host_ptr;
};

void gpuFree(void *host_ptr)
{
    cudaFreeHost(host_ptr);
}

__global__ void fillScreenKernel(uint32_t *buf)
{
    const int pixelX = blockDim.x * blockIdx.x + threadIdx.x;
    const int pixelY = blockDim.y * blockIdx.y + threadIdx.y;
    if (pixelX >= SCREEN_WIDTH || pixelY >= SCREEN_HEIGHT)
    {
        return;
    }
    unsigned int pos = SCREEN_WIDTH * pixelY + pixelX;
    buf[pos] = 0xFFFF0000;
}

void SPHSimulation(uint32_t *host_ptr)
{
    uint32_t *device_ptr;
    cudaError_t err = cudaHostGetDevicePointer(&device_ptr, host_ptr, 0);
    if (err != cudaSuccess)
    {
        std::cerr << "cudaHostGetDevicePointer error: " << cudaGetErrorString(err) << std::endl;
        return;
    }
    const dim3 gridDim(H_TILES, V_TILES);
    const dim3 blockDim(TILE_WIDTH, TILE_HEIGHT);
    fillScreenKernel<<<gridDim, blockDim>>>(device_ptr);

    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error: " << cudaGetErrorString(err) << std::endl;
    }
    cudaDeviceSynchronize();
}
