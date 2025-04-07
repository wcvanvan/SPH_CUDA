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
        if (z > 0)
        {
            float rho = C * z * z * z;
            particles[i].density += rho;
        }
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
        if (r2 >= h2 || r2 < 1e-12)
        {
            continue;
        }
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
        std::cerr << "Kernel error: " << cudaGetErrorString(err) << std::endl;
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

void reflect(Particle &particle, const Sink &sink)
{
    const float DAMP = 0.75f;
    float tbounce = 0.0f;
    if (particle.position.x > sink.xLen / 2 || particle.position.x < -sink.xLen / 2)
    {
        if (particle.position.x > sink.xLen / 2)
        {
            tbounce = (particle.position.x - sink.xLen / 2) / particle.velocity.x;
            particle.position.x = sink.xLen - particle.position.x;
        }
        else
        {
            tbounce = (particle.position.x + sink.xLen / 2) / particle.velocity.x;
            particle.position.x = -sink.xLen - particle.position.x;
        }
        particle.position.y -= particle.velocity.y * (1 - DAMP) * tbounce;
        particle.position.z -= particle.velocity.z * (1 - DAMP) * tbounce;
        particle.velocity.x = -particle.velocity.x;
        particle.velocityHalf.x = -particle.velocityHalf.x;
        particle.velocity.x *= DAMP;
        particle.velocity.y *= DAMP;
        particle.velocity.z *= DAMP;
        particle.velocityHalf.x *= DAMP;
        particle.velocityHalf.y *= DAMP;
        particle.velocityHalf.z *= DAMP;
    }
    if (particle.position.y > sink.yLen / 2 || particle.position.y < -sink.yLen / 2)
    {
        if (particle.position.y > sink.yLen / 2)
        {
            tbounce = (particle.position.y - sink.yLen / 2) / particle.velocity.y;
            particle.position.y = sink.yLen - particle.position.y;
        }
        else
        {
            tbounce = (particle.position.y + sink.yLen / 2) / particle.velocity.y;
            particle.position.y = -sink.yLen - particle.position.y;
        }
        particle.position.x -= particle.velocity.x * (1 - DAMP) * tbounce;
        particle.position.z -= particle.velocity.z * (1 - DAMP) * tbounce;
        particle.velocity.y = -particle.velocity.y;
        particle.velocityHalf.y = -particle.velocityHalf.y;
        particle.velocity.x *= DAMP;
        particle.velocity.y *= DAMP;
        particle.velocity.z *= DAMP;
        particle.velocityHalf.x *= DAMP;
        particle.velocityHalf.y *= DAMP;
        particle.velocityHalf.z *= DAMP;
    }
    if (particle.position.z > sink.zLen / 2 || particle.position.z < -sink.zLen / 2)
    {
        if (particle.position.z > sink.zLen / 2)
        {
            tbounce = (particle.position.z - sink.zLen / 2) / particle.velocity.z;
            particle.position.z = sink.zLen - particle.position.z;
        }
        else
        {
            tbounce = (particle.position.z + sink.zLen / 2) / particle.velocity.z;
            particle.position.z = -sink.zLen - particle.position.z;
        }
        particle.position.x -= particle.velocity.x * (1 - DAMP) * tbounce;
        particle.position.y -= particle.velocity.y * (1 - DAMP) * tbounce;
        particle.velocity.z = -particle.velocity.z;
        particle.velocityHalf.z = -particle.velocityHalf.z;
        particle.velocity.x *= DAMP;
        particle.velocity.y *= DAMP;
        particle.velocity.z *= DAMP;
        particle.velocityHalf.x *= DAMP;
        particle.velocityHalf.y *= DAMP;
        particle.velocityHalf.z *= DAMP;
    }
}

void leapfrogStart(Particle *particles, int particleCount, const Sink &sink)
{
    for (int i = 0; i < particleCount; i++)
    {
        particles[i].velocityHalf = particles[i].velocity + particles[i].acceleration * DELTA_T / 2;
        particles[i].velocity = particles[i].velocity + particles[i].acceleration * DELTA_T;
        particles[i].position = particles[i].position + particles[i].velocityHalf * DELTA_T;
        reflect(particles[i], sink);
    }
}

void leapfrogStep(Particle *particles, int particleCount, const Sink &sink)
{
    for (int i = 0; i < particleCount; i++)
    {
        particles[i].velocityHalf = particles[i].velocityHalf + particles[i].acceleration * DELTA_T;
        particles[i].velocity = particles[i].velocityHalf + particles[i].acceleration * DELTA_T / 2;
        particles[i].position = particles[i].position + particles[i].velocityHalf * DELTA_T;
        reflect(particles[i], sink);
    }
}

void initSimulation(Particle *particles, int particleCount, const Sink &sink, float mass)
{
    int blockDim = 32;
    int gridDim = (particleCount + (blockDim - 1)) / blockDim;
    cudaError_t err;
    computeDensity<<<gridDim, blockDim>>>(particles, particleCount, mass);
    cudaDeviceSynchronize();
    computeAccel<<<gridDim, blockDim>>>(particles, particleCount, mass);
    cudaDeviceSynchronize();
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error: " << cudaGetErrorString(err) << std::endl;
    }
    leapfrogStart(particles, particleCount, sink);
}

void updateSimulation(Particle *particles, int particleCount, const Sink &sink, float mass)
{
    int blockDim = 32;
    int gridDim = (particleCount + (blockDim - 1)) / blockDim;
    cudaError_t err;
    computeDensity<<<gridDim, blockDim>>>(particles, particleCount, mass);
    cudaDeviceSynchronize();
    computeAccel<<<gridDim, blockDim>>>(particles, particleCount, mass);
    cudaDeviceSynchronize();
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error: " << cudaGetErrorString(err) << std::endl;
    }
    leapfrogStep(particles, particleCount, sink);
}
