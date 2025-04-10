#ifndef COMPUTE_H
#define COMPUTE_H

#include <cstdio>
#include <vector>
#include "const.h"

// Forward declarations for types and functions defined elsewhere
struct Vec2;
struct Vec3;
struct Mat4;
struct Sink;
struct Particle;

// Frame data structure storing simulation results
struct FrameData {
  std::vector<std::vector<Vec2> *> frames;

  ~FrameData() {
    for (auto frame : frames) delete frame;
  }
};

// Global filename for particle output
extern const char *filename;

// Main compute entry point
int compute();

// Utility function to serialize frame data
void writeDataToFile(FILE *file, const FrameData *frameData, int particleCount);

#endif