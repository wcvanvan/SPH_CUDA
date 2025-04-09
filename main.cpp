#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <string.h>
#include <assert.h>
#include "const.h"
#include "sph.h"
#include "sink.h"
#include "vec.h"

const char *filename = "./particles.dat";
struct FrameData {
    std::vector<std::vector<Vec2>*> frames;
    ~FrameData() {
        for (auto frame : frames) delete frame;
    }
};

void writeDataToFile(FILE *file, const FrameData *frameData, int particleCount)
{
    fprintf(file, "%d\n", particleCount);
    int count = 0;
    for (auto &&frame : frameData->frames)
    {
        std::cout << "frame " << count++ << std::endl;
        fprintf(file, "FRAMESTART\n");
        for (const auto &pos : *frame)
        {
            fprintf(file, "%.2f %.2f\n", pos.x, pos.y);
        }
        fprintf(file, "FRAMEEND\n");
    }
    fflush(file);
}

int main(int argc, char *args[])
{
    Vec3 camera(1.25, 1.5, 1.8);
    Vec3 center(0, 0, 0);
    Vec3 up(0, 1, 0);

    float fov = 45.0f * M_PI / 180.0f;
    float aspectRatio = (float)SCREEN_WIDTH / SCREEN_HEIGHT;
    float near = 0.1f;
    float far = 10.0f;
    Mat4 projMat = proj(fov, near, far, aspectRatio);
    Mat4 viewMat = lookat(camera, center, up);
    Mat4 transformMat = projMat * viewMat;
    float *transformMatOnGPU = allocateMatOnGPU(transformMat);

    FILE *file = fopen(filename, "w");
    if (!file)
    {
        std::cerr << "Error opening file" << std::endl;
        return 1;
    }
    // Init scene
    Sink sink;
    int particleCount = 0;
    float mass = 1.0f;
    Particle *particles = initParticles(particleCount, mass, sink);
    initSimulation(particles, particleCount, sink, mass, transformMatOnGPU);
    FrameData *frameData = new FrameData();
    int count = 0;
    while (count < FRAMES)
    {
        updateSimulation(particles, particleCount, sink, mass, transformMatOnGPU);
        std::vector<Vec2> *framePositions = new std::vector<Vec2>();
        framePositions->reserve(particleCount);
        for (int i = 0; i < particleCount; i++)
        {
            framePositions->push_back(particles[i].screenPos);
        }
        frameData->frames.push_back(framePositions);
        count++;
    }
    writeDataToFile(file, frameData, particleCount);
    delete frameData;
    return 0;
}
