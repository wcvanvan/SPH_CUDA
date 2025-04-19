#include "compute.h"
#include <chrono>
#include <cmath>
#include <iostream>
#include "const.h"
#include "scene.h"
#include "sph.h"
#include "vec.h"

const char *filename = "./particles.dat";

void getTroughPosition(Trough &trough, Sink &sink) {
  // translate trough to the side and above the sink
  float translationX = trough.xLen / 2.0f + sink.xLen / 2.0f;
  float translationY = sink.yLen / 2.0f + trough.yLen / 2.0f;
  for (int i = 0; i < 8; i++) {
    trough.vertices[i].x -= translationX;
    trough.vertices[i].y += translationY;
  }
  // calculate the bottom's normal vector
  Vec3 v1 = trough.vertices[1] - trough.vertices[0];
  Vec3 v2 = trough.vertices[3] - trough.vertices[0];
  trough.normal = cross(v1, v2).normalize();
  // std::cout << "trough's normal: " << trough.normal.x << " " << trough.normal.y << " " << trough.normal.z <<
  // std::endl;

  // calculate the slope and intercept of the x-y plane of the bottom
  float slope = (trough.vertices[0].y - trough.vertices[1].y) / (trough.vertices[0].x - trough.vertices[1].x);
  float intercept = trough.vertices[1].y - slope * trough.vertices[1].x;
  trough.slope = slope;
  trough.intercept = intercept;
  // std::cout << "slope: " << slope << " intercept: " << intercept << std::endl;
}

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

int main() {
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
  Trough trough;
  getTroughPosition(trough, sink);
  int particleCount = 0;
  float mass = 1.0f;
  float cellSize = KERNEL_RADIUS;
  int gridDimX = (int)ceil(sink.xLen / cellSize);
  int gridDimY = (int)ceil(sink.yLen / cellSize);
  int gridDimZ = (int)ceil(sink.zLen / cellSize);
  int totalCells = gridDimX * gridDimY * gridDimZ;
  int *cellStart = initCellStart(totalCells);
  int *cellEnd = initCellEnd(totalCells);
  Particle *particles = initParticles(particleCount, mass, sink, trough, cellStart, cellEnd);

  // [Optional] Timer
  auto init_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> init_elapsed = init_end - start;
  std::cout << "Initialization time: " << init_elapsed.count() << " seconds" << std::endl;

  FrameData *frameData = new FrameData();
  int count = 0;
  while (count < FRAMES) {
    updateSimulation(particles, particleCount, sink, trough, mass, transformMatOnGPU, cellStart, cellEnd);
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
