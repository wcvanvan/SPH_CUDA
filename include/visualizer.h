#ifndef VISUALIZER_H
#define VISUALIZER_H

#include <SDL2/SDL.h>
#include <vector>

#include "const.h"
#include "scene.h"
#include "vec.h"

// External global transform matrix used for coordinate transforms
extern Mat4 transformMat;

// Path to the particle data file
extern const char *filename;

// Struct to store all frame data
struct FrameData2 {
  std::vector<std::vector<Vec2>> frames;
};

// Core functions for visualization
int visualize();
FrameData2 *readAllFramesFromFile(int &particleCount);

// Utility functions
Vec2 worldToScreen(Vec3 &worldCoord);
void transformSinkPoints(Sink &sink);
void renderSinkBackFaces(SDL_Renderer *sdlRenderer, Sink &sink);
void renderSinkFrontFaces(SDL_Renderer *sdlRenderer, Sink &sink);
void drawPoints(std::vector<Vec2> &points, SDL_Renderer *sdlRenderer);

#endif  // VISUALIZER_H
