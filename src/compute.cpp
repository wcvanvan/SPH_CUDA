#include "compute.h"
#include <chrono>
#include <cmath>
#include <iostream>
#include "const.h"
#include "sink.h"
#include "sph.h"
#include "vec.h"

const char *filename = "./particles.dat";

void writeDataToFile(FILE *file, const FrameData *frameData, int particleCount) {
  fprintf(file, "%d\n", particleCount);
  int count = 0;
  for (auto &&frame : frameData->frames) {
    fprintf(file, "FRAMESTART\n");
    for (const auto &pos : *frame) {
      fprintf(file, "%.2f %.2f\n", pos.x, pos.y);
    }
    fprintf(file, "FRAMEEND\n");
  }
  fflush(file);
}

int compute() {
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
  if (!file) {
    std::cerr << "Error opening file" << std::endl;
    return 1;
  }

  // [Optional] Timer
  auto start = std::chrono::high_resolution_clock::now();

  // Init scene
  Sink sink;
  int particleCount = 0;
  float mass = 1.0f;
  Particle *particles = initParticles(particleCount, mass, sink);
  initSimulation(particles, particleCount, sink, mass, transformMatOnGPU);

  // [Optional] Timer
  auto init_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> init_elapsed = init_end - start;
  std::cout << "Initialization time: " << init_elapsed.count() << " seconds" << std::endl;

  FrameData *frameData = new FrameData();
  int count = 0;
  while (count < FRAMES) {
    updateSimulation(particles, particleCount, sink, mass, transformMatOnGPU);
    std::vector<Vec2> *framePositions = new std::vector<Vec2>();
    framePositions->reserve(particleCount);
    for (int i = 0; i < particleCount; i++) {
      framePositions->push_back(particles[i].screenPos);
    }
    frameData->frames.push_back(framePositions);
    count++;
  }

  // [Optional] Timer
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - init_end;
  std::cout << "Frame count: " << frameData->frames.size() << std::endl;
  std::cout << "Computation time: " << elapsed.count() << " seconds" << std::endl;

  writeDataToFile(file, frameData, particleCount);
  delete frameData;
  return 0;
}
