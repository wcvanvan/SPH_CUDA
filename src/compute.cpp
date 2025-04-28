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

  // calculate the slope and intercept of the x-y plane of the bottom
  if (fabs(trough.vertices[0].x - trough.vertices[1].x) > 1e-6) {
    trough.slope = (trough.vertices[0].y - trough.vertices[1].y) / (trough.vertices[0].x - trough.vertices[1].x);
    trough.intercept = trough.vertices[1].y - trough.slope * trough.vertices[1].x;
  } else {
    trough.slope = 0.0f;
    trough.intercept = 0.0f;
    std::cerr << "Warning: Trough bottom appears vertical, slope calculation may be inaccurate." << std::endl;
  }
}

void writeDataToFile(FILE *file, const Vec2 *screenPosOnCPU, int particleCount, int totalFrameCount) {
  fprintf(file, "%d\n", particleCount);
  for (int frame = 0; frame < totalFrameCount; frame++) {
    fprintf(file, "FRAMESTART\n");
    if (particleCount > 0) {
      for (int i = 0; i < particleCount; i++) {
        const Vec2 &pos = screenPosOnCPU[(size_t)frame * particleCount + i];
        fprintf(file, "%.2f %.2f\n", pos.x, pos.y);
      }
    }
    fprintf(file, "FRAMEEND\n");
  }
  fflush(file);  // Ensure data is written
  fclose(file);
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
    std::cerr << "Error opening file: " << filename << std::endl;
    cleanupMatOnGPU(transformMatOnGPU);
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();

  float h = KERNEL_RADIUS;
  float h_pow_9 = powf(h, 9.0f);
  float h_pow_6 = powf(h, 6.0f);
  float POLY6 = (315.0f / (64.0f * M_PI * h_pow_9));
  float WEIGHT_AT_0 = POLY6 * h_pow_6;
  float VISCOSITY_LAPLACIAN = 45.0f / (M_PI * h_pow_6);

  Sink sink;
  Trough trough;
  getTroughPosition(trough, sink);

  ParticlesSoA particles;
  int particleCount = 0;
  float mass = 1.0f;

  float cellSize = KERNEL_RADIUS;
  int gridDimX = (int)ceilf(sink.xLen / cellSize);
  int gridDimY = (int)ceilf(sink.yLen / cellSize);
  int gridDimZ = (int)ceilf(sink.zLen / cellSize);
  int totalCells = gridDimX * gridDimY * gridDimZ;
  if (totalCells <= 0) {
    std::cerr << "Error: Invalid grid dimensions calculated. Check scene dimensions and kernel radius." << std::endl;
    fclose(file);
    cleanupMatOnGPU(transformMatOnGPU);
    return 1;
  }

  int *cellStart = initCellStart(totalCells);
  int *cellEnd = initCellEnd(totalCells);

  initParticlesSoA(particles, particleCount, mass, sink, trough, cellStart, cellEnd, POLY6, WEIGHT_AT_0);

  Vec2 *screenPosOnGPU = nullptr;
  Vec2 *allFramesOnGPU = nullptr;
  Vec2 *screenPosOnCPU = nullptr;

  if (particleCount > 0) {
    screenPosOnGPU = initScreenPos(particleCount);
    allFramesOnGPU = initAllFrames(particleCount, FRAMES);
    try {
      screenPosOnCPU = new Vec2[(size_t)particleCount * FRAMES];
    } catch (const std::bad_alloc &e) {
      std::cerr << "Error: Failed to allocate host memory for screen positions: " << e.what() << std::endl;
      cleanupParticlesSoA(particles);
      cleanupCellStart(cellStart);
      cleanupCellEnd(cellEnd);
      cleanupMatOnGPU(transformMatOnGPU);
      cleanupScreenPos(screenPosOnGPU);
      cleanupAllFrames(allFramesOnGPU);
      fclose(file);
      return 1;
    }
  } else {
    std::cout << "Warning: 0 particles initialized. Simulation will not run." << std::endl;
  }

  auto init_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> init_elapsed = init_end - start;
  std::cout << "Initialization time: " << init_elapsed.count() << " seconds" << std::endl;

  int frameCount = 0;
  if (particleCount > 0) {
    while (frameCount < FRAMES) {
      updateSimulationSoA(particles, sink, trough, mass, transformMatOnGPU, cellStart, cellEnd, screenPosOnGPU, POLY6,
                          VISCOSITY_LAPLACIAN, WEIGHT_AT_0, frameCount, allFramesOnGPU);
      frameCount++;
    }

    copyAllFramesToCPU(allFramesOnGPU, screenPosOnCPU, particleCount, FRAMES);
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - init_end;
  std::cout << "Frame count: " << frameCount << std::endl;
  std::cout << "Computation time: " << elapsed.count() << " seconds" << std::endl;

  if (particleCount > 0) {
    writeDataToFile(file, screenPosOnCPU, particleCount, frameCount);
  } else {
    fprintf(file, "0\n");
    fclose(file);
  }

  delete[] screenPosOnCPU;

  cleanupParticlesSoA(particles);
  cleanupCellStart(cellStart);
  cleanupCellEnd(cellEnd);
  cleanupMatOnGPU(transformMatOnGPU);
  cleanupScreenPos(screenPosOnGPU);
  cleanupAllFrames(allFramesOnGPU);

  std::cout << "Simulation finished and data written to " << filename << std::endl;
  return 0;
}