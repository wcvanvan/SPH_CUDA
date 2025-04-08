#include <iostream>
#include "sph.h"

__global__ void computeDensity(Particle *particles, int particleCount, float mass)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= particleCount)
    {
        return;
    }
    particles[i].density = 0.0f;
    float h2 = KERNEL_RADIUS * KERNEL_RADIUS;
    float h4 = h2 * h2;
    float h8 = h4 * h4;
    float C = 4 * mass / M_PI / h8;
    // naive method. iterate over all the other particles
    particles[i].density += 4 * mass / M_PI / h2;
    for (int j = 0; j < particleCount; j++)
    {
        if (i == j)
        {
            continue;
        }
        float dx = particles[i].position.x - particles[j].position.x;
        float dy = particles[i].position.y - particles[j].position.y;
        float dz = particles[i].position.z - particles[j].position.z;
        float r2 = dx * dx + dy * dy + dz * dz;
        float z = h2 - r2;
        if (z <= 0) continue;
        float rho = C * z * z * z;
        particles[i].density += rho;
    }
}

__global__ void computeAccel(Particle *particles, int particleCount, float mass)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= particleCount)
    {
        return;
    }
    particles[i].acceleration.x = 0.0f;
    particles[i].acceleration.y = GRAVITY;
    particles[i].acceleration.z = 0.0f;

    float h2 = KERNEL_RADIUS * KERNEL_RADIUS;
    float C0 = mass / M_PI / (h2 * h2);
    float Cp = 15 * STIFFNESS;
    float Cv = -40.0f * VISCOSITY;

    // naive method. iterate over all the other particles
    for (int j = 0; j < particleCount; j++)
    {
        if (i == j)
        {
            continue;
        }
        float dx = particles[i].position.x - particles[j].position.x;
        float dy = particles[i].position.y - particles[j].position.y;
        float dz = particles[i].position.z - particles[j].position.z;
        float r2 = dx * dx + dy * dy + dz * dz;
        if (r2 >= h2) continue;
        float q = sqrt(r2) / KERNEL_RADIUS;
        float u = 1 - q;
        float w0 = C0 * u / particles[i].density / particles[j].density;
        float wp = w0 * Cp * (particles[i].density + particles[j].density - 2 * REST_DENSITY) * u / q;
        float wv = w0 * Cv;
        float dvx = particles[i].velocity.x - particles[j].velocity.x;
        float dvy = particles[i].velocity.y - particles[j].velocity.y;
        float dvz = particles[i].velocity.z - particles[j].velocity.z;
        particles[i].acceleration.x += (wp * dx + wv * dvx);
        particles[i].acceleration.y += (wp * dy + wv * dvy);
        particles[i].acceleration.z += (wp * dz + wv * dvz);
    }
}

float normalizeMass(Particle *particles, int particleCount)
{
    float mass = 1.0f;
    int blockDim = 32;
    int gridDim = (particleCount + (blockDim - 1)) / blockDim;
    cudaError_t err;
    computeDensity<<<gridDim, blockDim>>>(particles, particleCount, mass);
    cudaDeviceSynchronize();
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (computeDensity): " << cudaGetErrorString(err) << std::endl;
    }
    float rho0 = REST_DENSITY;
    float rho2s = 0.0f;
    float rhos = 0.0f;
    for (int i = 0; i < particleCount; i++)
    {
        rho2s += particles[i].density * particles[i].density;
        rhos += particles[i].density;
    }
    mass *= (rho0 * rhos / rho2s);
    std::cout << "Mass: " << mass << std::endl;
    return mass;
}

Particle *initParticles(int particleCount, float *mass, Sink &sink)
{
    Particle *particles;
    cudaMallocManaged(&particles, particleCount * sizeof(Particle));
    srand(time(0));
    for (int i = 0; i < particleCount; i++)
    {
        particles[i].velocity = Vec3(0.0f, 0.0f, 0.0f);
        particles[i].acceleration = Vec3(0.0f, 0.0f, 0.0f);
        particles[i].position.x = ((float)rand() / RAND_MAX) * sink.xLen - sink.xLen / 2;
        particles[i].position.y = ((float)rand() / RAND_MAX) * sink.yLen - sink.yLen / 2;
        particles[i].position.z = ((float)rand() / RAND_MAX) * sink.zLen - sink.zLen / 2;
    }
    *mass = normalizeMass(particles, particleCount);
    return particles;
}

__device__ void reflect(Particle &particle, float xLen, float yLen, float zLen)
{
    float tbounce = 0.0f;
    if (particle.velocity.x != 0 && (particle.position.x > xLen / 2 || particle.position.x < -xLen / 2))
    {
        if (particle.position.x > xLen / 2)
        {
            tbounce = (particle.position.x - xLen / 2) / particle.velocity.x;
            particle.position.x = xLen - particle.position.x;
        }
        else
        {
            tbounce = (particle.position.x + xLen / 2) / particle.velocity.x;
            particle.position.x = -xLen - particle.position.x;
        }
        // revert the movement for the period
        particle.position.y -= particle.velocity.y * (1 - REFLECT_DAMP) * tbounce;
        particle.position.z -= particle.velocity.z * (1 - REFLECT_DAMP) * tbounce;
        particle.velocity.x = -particle.velocity.x;
        particle.velocityHalf.x = -particle.velocityHalf.x;
        particle.velocity.x *= REFLECT_DAMP;
        particle.velocity.y *= REFLECT_DAMP;
        particle.velocity.z *= REFLECT_DAMP;
        particle.velocityHalf.x *= REFLECT_DAMP;
        particle.velocityHalf.y *= REFLECT_DAMP;
        particle.velocityHalf.z *= REFLECT_DAMP;
    }
    if (particle.velocity.y != 0 && (particle.position.y > yLen / 2 || particle.position.y < -yLen / 2))
    {
        // bounce back
        if (particle.position.y > yLen / 2)
        {
            tbounce = (particle.position.y - yLen / 2) / particle.velocity.y;
            particle.position.y = yLen - particle.position.y;
        }
        else
        {
            tbounce = (particle.position.y + yLen / 2) / particle.velocity.y;
            particle.position.y = -yLen - particle.position.y;
        }
        // revert the movement for the period
        particle.position.x -= particle.velocity.x * (1 - REFLECT_DAMP) * tbounce;
        particle.position.z -= particle.velocity.z * (1 - REFLECT_DAMP) * tbounce;
        particle.velocity.y = -particle.velocity.y;
        particle.velocityHalf.y = -particle.velocityHalf.y;
        particle.velocity.x *= REFLECT_DAMP;
        particle.velocity.y *= REFLECT_DAMP;
        particle.velocity.z *= REFLECT_DAMP;
        particle.velocityHalf.x *= REFLECT_DAMP;
        particle.velocityHalf.y *= REFLECT_DAMP;
        particle.velocityHalf.z *= REFLECT_DAMP;
    }
    if (particle.velocity.z != 0 && (particle.position.z > zLen / 2 || particle.position.z < -zLen / 2))
    {
        // bounce back
        if (particle.position.z > zLen / 2)
        {
            tbounce = (particle.position.z - zLen / 2) / particle.velocity.z;
            particle.position.z = zLen - particle.position.z;
        }
        else
        {
            tbounce = (particle.position.z + zLen / 2) / particle.velocity.z;
            particle.position.z = -zLen - particle.position.z;
        }
        // revert the movement for the period
        particle.position.x -= particle.velocity.x * (1 - REFLECT_DAMP) * tbounce;
        particle.position.y -= particle.velocity.y * (1 - REFLECT_DAMP) * tbounce;
        particle.velocity.z = -particle.velocity.z;
        particle.velocityHalf.z = -particle.velocityHalf.z;
        particle.velocity.x *= REFLECT_DAMP;
        particle.velocity.y *= REFLECT_DAMP;
        particle.velocity.z *= REFLECT_DAMP;
        particle.velocityHalf.x *= REFLECT_DAMP;
        particle.velocityHalf.y *= REFLECT_DAMP;
        particle.velocityHalf.z *= REFLECT_DAMP;
    }
}

__global__ void leapfrogStart(Particle *particles, int particleCount, float xLen, float yLen, float zLen)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= particleCount)
    {
        return;
    }
    particles[i].velocityHalf.x = particles[i].velocity.x + particles[i].acceleration.x * DELTA_T / 2;
    particles[i].velocityHalf.y = particles[i].velocity.y + particles[i].acceleration.y * DELTA_T / 2;
    particles[i].velocityHalf.z = particles[i].velocity.z + particles[i].acceleration.z * DELTA_T / 2;
    particles[i].velocity.x = particles[i].velocity.x + particles[i].acceleration.x * DELTA_T;
    particles[i].velocity.y = particles[i].velocity.y + particles[i].acceleration.y * DELTA_T;
    particles[i].velocity.z = particles[i].velocity.z + particles[i].acceleration.z * DELTA_T;
    particles[i].position.x = particles[i].position.x + particles[i].velocityHalf.x * DELTA_T;
    particles[i].position.y = particles[i].position.y + particles[i].velocityHalf.y * DELTA_T;
    particles[i].position.z = particles[i].position.z + particles[i].velocityHalf.z * DELTA_T;
    reflect(particles[i], xLen, yLen, zLen);
}

__global__ void leapfrogStep(Particle *particles, int particleCount, float xLen, float yLen, float zLen)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= particleCount)
    {
        return;
    }
    particles[i].velocityHalf.x = particles[i].velocityHalf.x + particles[i].acceleration.x * DELTA_T;
    particles[i].velocityHalf.y = particles[i].velocityHalf.y + particles[i].acceleration.y * DELTA_T;
    particles[i].velocityHalf.z = particles[i].velocityHalf.z + particles[i].acceleration.z * DELTA_T;
    particles[i].velocity.x = particles[i].velocityHalf.x + particles[i].acceleration.x * DELTA_T / 2;
    particles[i].velocity.y = particles[i].velocityHalf.y + particles[i].acceleration.y * DELTA_T / 2;
    particles[i].velocity.z = particles[i].velocityHalf.z + particles[i].acceleration.z * DELTA_T / 2;
    particles[i].position.x = particles[i].position.x + particles[i].velocityHalf.x * DELTA_T;
    particles[i].position.y = particles[i].position.y + particles[i].velocityHalf.y * DELTA_T;
    particles[i].position.z = particles[i].position.z + particles[i].velocityHalf.z * DELTA_T;
    reflect(particles[i], xLen, yLen, zLen);
}

void initSimulation(Particle *particles, int particleCount, const Sink &sink, float mass)
{
    int blockDim = 32;
    int gridDim = (particleCount + (blockDim - 1)) / blockDim;
    cudaError_t err;
    computeDensity<<<gridDim, blockDim>>>(particles, particleCount, mass);
    cudaDeviceSynchronize(); // computing acceleration needs all particles' density
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (computeDensity): " << cudaGetErrorString(err) << std::endl;
    }
    computeAccel<<<gridDim, blockDim>>>(particles, particleCount, mass);
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (ComputeAccel): " << cudaGetErrorString(err) << std::endl;
    }
    // no need to synchronize between computing acceleration and integration
    leapfrogStart<<<gridDim, blockDim>>>(particles, particleCount, sink.xLen, sink.yLen, sink.zLen);
    cudaDeviceSynchronize();
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (leapfrogStart): " << cudaGetErrorString(err) << std::endl;
    }
}

void updateSimulation(Particle *particles, int particleCount, const Sink &sink, float mass)
{
    int blockDim = 32;
    int gridDim = (particleCount + (blockDim - 1)) / blockDim;
    cudaError_t err;
    computeDensity<<<gridDim, blockDim>>>(particles, particleCount, mass);
    cudaDeviceSynchronize(); // computing acceleration needs all particles' density
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (computeDensity): " << cudaGetErrorString(err) << std::endl;
    }
    computeAccel<<<gridDim, blockDim>>>(particles, particleCount, mass);
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (ComputeAccel): " << cudaGetErrorString(err) << std::endl;
    }
    // no need to synchronize between computing acceleration and integration
    leapfrogStep<<<gridDim, blockDim>>>(particles, particleCount, sink.xLen, sink.yLen, sink.zLen);
    cudaDeviceSynchronize();
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (leapfrogStep): " << cudaGetErrorString(err) << std::endl;
    }
}
