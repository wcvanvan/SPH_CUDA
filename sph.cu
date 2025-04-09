#include <iostream>
#include "sph.h"
#include "const.h"

float *allocateMatOnGPU(Mat4 &mat) {
    float data[16];
    int count = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            data[count++] = mat[i][j];
        }
    }
    float* matOnGPU;
    cudaMalloc((void**)&matOnGPU, 16 * sizeof(float));
    cudaMemcpy(matOnGPU, data, 16 * sizeof(float), cudaMemcpyHostToDevice);
    return matOnGPU;
}

__global__ void computeDensity(Particle *particles, int particlesCount, float mass)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= particlesCount)
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
    for (int j = 0; j < particlesCount; j++)
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

__global__ void computeAccel(Particle *particles, int particlesCount, float mass)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= particlesCount)
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
    for (int j = 0; j < particlesCount; j++)
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

float normalizeMass(Particle *particles, int particlesCount)
{
    float mass = 1.0f;
    int blockDim = 32;
    int gridDim = (particlesCount + (blockDim - 1)) / blockDim;
    cudaError_t err;
    computeDensity<<<gridDim, blockDim>>>(particles, particlesCount, mass);
    cudaDeviceSynchronize();
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (computeDensity): " << cudaGetErrorString(err) << std::endl;
    }
    float rho0 = REST_DENSITY;
    float rho2s = 0.0f;
    float rhos = 0.0f;
    for (int i = 0; i < particlesCount; i++)
    {
        rho2s += particles[i].density * particles[i].density;
        rhos += particles[i].density;
    }
    mass *= (rho0 * rhos / rho2s);
    std::cout << "Mass: " << mass << std::endl;
    return mass;
}

Particle * placeParticles(int &particlesCount, int droppingParticlesCount, Sink &sink) {
    float h = KERNEL_RADIUS;
    float hh = h / 1.3f;
    int particlesInSink = 0;
    for (float x = -sink.xLen / 2.0f; x <= sink.xLen / 2.0f; x += hh) {
        for (float z = -sink.zLen / 2.0f; z <= sink.zLen / 2.0f; z += hh) {
            for (float y = -sink.yLen / 2.0f; y <= 0; y += hh) {
                particlesInSink++;
            }
        }
    }
    particlesCount = particlesInSink + droppingParticlesCount;
    std::cout << "Particle Count: " << particlesCount << std::endl;
    Particle *particles;
    cudaMallocManaged(&particles, particlesCount * sizeof(Particle));
    int count = 0;
    for (float x = -sink.xLen / 2.0f; x <= sink.xLen / 2.0f; x += hh) {
        for (float z = -sink.zLen / 2.0f; z <= sink.zLen / 2.0f; z += hh) {
            for (float y = -sink.yLen / 2.0f; y <= 0; y += hh) {
                particles[count].position = {x, y, z};
                particles[count].density = 0.0f;
                particles[count].velocity = {0.0f, 0.0f, 0.0f};
                particles[count].velocityHalf = {0.0f, 0.0f, 0.0f};
                particles[count].acceleration = {0.0f, 0.0f, 0.0f};
                count++;
            }
        }
    }
    srand(time(0));
    for (int i = particlesInSink; i < particlesCount; i++) {
        particles[i].position.x = ((float)rand() / RAND_MAX) * sink.xLen / 5.0f - sink.xLen / 10.0f;
        particles[i].position.y = ((float)rand() / RAND_MAX) * sink.yLen;
        particles[i].position.z = ((float)rand() / RAND_MAX) * sink.zLen / 5.0f - sink.zLen / 10.0f;
        particles[i].density = 0.0f;
        particles[i].velocity = {0.0f, 0.0f, 0.0f};
        particles[i].velocityHalf = {0.0f, 0.0f, 0.0f};
        particles[i].acceleration = {0.0f, 0.0f, 0.0f};
    }
    return particles;
}

Particle *initParticles(int &particlesCount, float &mass, Sink &sink)
{
    int droppingParticlesCount = 0;
    Particle *particles = placeParticles(particlesCount, droppingParticlesCount, sink);
    mass = normalizeMass(particles, particlesCount - droppingParticlesCount);
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

__global__ void leapfrogStart(Particle *particles, int particlesCount, float xLen, float yLen, float zLen)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= particlesCount)
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

__global__ void leapfrogStep(Particle *particles, int particlesCount, float xLen, float yLen, float zLen)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= particlesCount)
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

__global__ void coordTransform(Particle *particles, int particlesCount, float *transformMat) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= particlesCount)
    {
        return;
    }
    float worldPos[4], result[4];
    worldPos[0] = particles[i].position.x;
    worldPos[1] = particles[i].position.y;
    worldPos[2] = particles[i].position.z;
    worldPos[3] = 1.0f;
    for (int i = 0; i < 4; i++)
    {
        result[i] = 0.0f;
        for (int j = 0; j < 4; j++)
        {
            result[i] += transformMat[i * 4 + j] * worldPos[j];
        }
    }
    float x = result[0] / result[3];
    float y = result[1] / result[3];

    if (x < -1.0f || x > 1.0f || y < -1.0f || y > 1.0f) {
        particles[i].screenPos.x = -1.0f;
        particles[i].screenPos.y = -1.0f;
    } else {
        float screenX = fmaxf(0.0f, fminf(1.0f, (x + 1.0f) * 0.5f)) * SCREEN_WIDTH;
        float screenY = fmaxf(0.0f, fminf(1.0f, (1.0f - y) * 0.5f)) * SCREEN_HEIGHT;
        particles[i].screenPos.x = screenX;
        particles[i].screenPos.y = screenY;
    }
}

void initSimulation(Particle *particles, int particlesCount, const Sink &sink, float mass, float *transformMat)
{
    int blockDim = 32;
    int gridDim = (particlesCount + (blockDim - 1)) / blockDim;
    cudaError_t err;
    computeDensity<<<gridDim, blockDim>>>(particles, particlesCount, mass);
    cudaDeviceSynchronize(); // computing acceleration needs all particles' density
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (computeDensity): " << cudaGetErrorString(err) << std::endl;
    }
    computeAccel<<<gridDim, blockDim>>>(particles, particlesCount, mass);
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (ComputeAccel): " << cudaGetErrorString(err) << std::endl;
    }
    // no need to synchronize between computing acceleration and integration
    leapfrogStart<<<gridDim, blockDim>>>(particles, particlesCount, sink.xLen, sink.yLen, sink.zLen);
    coordTransform<<<gridDim, blockDim>>>(particles, particlesCount,transformMat);
    cudaDeviceSynchronize();
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (leapfrogStart or coordTransform): " << cudaGetErrorString(err) << std::endl;
    }
}

void updateSimulation(Particle *particles, int particlesCount, const Sink &sink, float mass, float *transformMat)
{
    int blockDim = 32;
    int gridDim = (particlesCount + (blockDim - 1)) / blockDim;
    cudaError_t err;
    computeDensity<<<gridDim, blockDim>>>(particles, particlesCount, mass);
    cudaDeviceSynchronize(); // computing acceleration needs all particles' density
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (computeDensity): " << cudaGetErrorString(err) << std::endl;
    }
    computeAccel<<<gridDim, blockDim>>>(particles, particlesCount, mass);
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (ComputeAccel): " << cudaGetErrorString(err) << std::endl;
    }
    // no need to synchronize between computing acceleration and integration
    leapfrogStep<<<gridDim, blockDim>>>(particles, particlesCount, sink.xLen, sink.yLen, sink.zLen);
    coordTransform<<<gridDim, blockDim>>>(particles, particlesCount,transformMat);
    cudaDeviceSynchronize();
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
        std::cerr << "Kernel error (leapfrogStart or coordTransform): " << cudaGetErrorString(err) << std::endl;
    }
}
