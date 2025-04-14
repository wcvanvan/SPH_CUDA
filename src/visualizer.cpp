#include "visualizer.h"
#include <SDL2/SDL.h>
#include <SDL2_gfxPrimitives.h>
#include <cmath>
#include <iostream>

const char *filename = "./particles.dat";

Mat4 transformMat;

void drawPoints(std::vector<Vec2> &points, SDL_Renderer *sdlRenderer) {
  for (auto &&v : points) {
    SDL_RenderDrawPoint(sdlRenderer, static_cast<int>(v.x), static_cast<int>(v.y));
  }
}

Vec2 worldToScreen(Vec3 &worldCoord) {
  Vec4 p(worldCoord);
  Vec4 ndcCoord = transformMat * p;
  ndcCoord.x /= ndcCoord.w;
  ndcCoord.y /= ndcCoord.w;

  if (ndcCoord.x < -1.0f || ndcCoord.x > 1.0f || ndcCoord.y < -1.0f || ndcCoord.y > 1.0f) return {-1, -1};

  float screenX = std::max(0.0f, std::min(1.0f, (ndcCoord.x + 1.0f) * 0.5f)) * SCREEN_WIDTH;
  float screenY = std::max(0.0f, std::min(1.0f, (1.0f - ndcCoord.y) * 0.5f)) * SCREEN_HEIGHT;
  return {screenX, screenY};
}

void transformSinkPoints(Sink &sink) {
  // change from world coordinates to screen coordinates
  for (int i = 0; i < 8; i++) {
    Vec2 screenCoord = worldToScreen(sink.vertices[i]);
    sink.screenPoints[i] = {(int)screenCoord.x, (int)screenCoord.y};
  }
}

void transformTroughPoints(Trough &trough, Sink &sink) {
  // translate trough to the side and above the sink
  float translationX = trough.xLen / 2.0f + sink.xLen / 2.0f;
  float translationY = sink.yLen / 2.0f + trough.yLen / 2.0f;
  for (int i = 0; i < 8; i++) {
    trough.vertices[i].x -= translationX;
    trough.vertices[i].y += translationY;
  }
  // change from world coordinates to screen coordinates
  for (int i = 0; i < 8; i++) {
    Vec2 screenCoord = worldToScreen(trough.vertices[i]);
    trough.screenPoints[i] = {(int)screenCoord.x, (int)screenCoord.y};
  }
  // std::cout << "trough" << std::endl;
  // for (int i = 0; i < 8; i++) {
  //   std::cout << trough.vertices[i].x << " " << trough.vertices[i].y << " " << trough.vertices[i].z << std::endl;
  // }
  // std::cout << "sink" << std::endl;
  // for (int i = 0; i < 8; i++) {
  //   std::cout << sink.vertices[i].x << " " << sink.vertices[i].y << " " << sink.vertices[i].z << std::endl;
  // }
}

void renderSinkBackFaces(SDL_Renderer *sdlRenderer, Sink &sink) {
  for (int i = 0; i < 3; i++) {
    Sint16 vx[4], vy[4];
    for (int j = 0; j < 4; j++) {
      vx[j] = (Sint16)sink.screenPoints[sink.backFaces[i][j]].x;
      vy[j] = (Sint16)sink.screenPoints[sink.backFaces[i][j]].y;
    }
    filledPolygonRGBA(sdlRenderer, vx, vy, 4, 220, 220, 220, 255);
    polygonRGBA(sdlRenderer, vx, vy, 4, 0, 0, 0, 255);
  }
}

void renderSinkFrontFaces(SDL_Renderer *sdlRenderer, Sink &sink) {
  for (int i = 0; i < 2; i++) {
    Sint16 vx[4], vy[4];
    for (int j = 0; j < 4; j++) {
      vx[j] = (Sint16)sink.screenPoints[sink.frontFaces[i][j]].x;
      vy[j] = (Sint16)sink.screenPoints[sink.frontFaces[i][j]].y;
    }
    filledPolygonRGBA(sdlRenderer, vx, vy, 4, 220, 220, 220, 255);
    polygonRGBA(sdlRenderer, vx, vy, 4, 0, 0, 0, 255);
  }
}

void renderTroughBackFaces(SDL_Renderer *sdlRenderer, Trough &trough) {
  for (int i = 0; i < 3; i++) {
    Sint16 vx[4], vy[4];
    for (int j = 0; j < 4; j++) {
      vx[j] = (Sint16)trough.screenPoints[trough.backFaces[i][j]].x;
      vy[j] = (Sint16)trough.screenPoints[trough.backFaces[i][j]].y;
    }
    filledPolygonRGBA(sdlRenderer, vx, vy, 4, 220, 220, 220, 255);
    polygonRGBA(sdlRenderer, vx, vy, 4, 0, 0, 0, 255);
  }
}

void renderTroughFrontFaces(SDL_Renderer *sdlRenderer, Trough &trough) {
  for (int i = 0; i < 2; i++) {
    Sint16 vx[4], vy[4];
    for (int j = 0; j < 4; j++) {
      vx[j] = (Sint16)trough.screenPoints[trough.frontFaces[i][j]].x;
      vy[j] = (Sint16)trough.screenPoints[trough.frontFaces[i][j]].y;
    }
    filledPolygonRGBA(sdlRenderer, vx, vy, 4, 220, 220, 220, 255);
    polygonRGBA(sdlRenderer, vx, vy, 4, 0, 0, 0, 255);
  }
}

FrameData2 *readAllFramesFromFile(int &particleCount) {
  FrameData2 *frameData2 = new FrameData2();
  FILE *file = fopen(filename, "r");
  if (!file) {
    std::cerr << "Error opening file" << std::endl;
    return frameData2;
  }

  // Read particle count
  char line[256];
  if (fgets(line, sizeof(line), file) && sscanf(line, "%d", &particleCount) == 1) {
    std::cout << "Particle Count: " << particleCount << std::endl;
  } else {
    std::cout << "Fail to read particle count" << std::endl;
    fclose(file);
    return frameData2;
  }

  // Read frames
  while (!feof(file)) {
    while (fgets(line, sizeof(line), file)) {
      if (strncmp(line, "FRAMESTART", 10) == 0) {
        break;
      }
    }
    std::vector<Vec2> framePositions;
    while (fgets(line, sizeof(line), file)) {
      if (strncmp(line, "FRAMEEND", 8) == 0) {
        break;
      }

      float x, y;
      if (sscanf(line, "%f %f", &x, &y) == 2) {
        framePositions.push_back({x, y});
      }
    }
    frameData2->frames.push_back(framePositions);
  }

  fclose(file);
  return frameData2;
}

int main() {
  int particleCount = 0;
  std::cout << "Reading frame data" << std::endl;
  FrameData2 *frameData2 = readAllFramesFromFile(particleCount);
  std::cout << "Finished reading frame data" << std::endl;
  if (SDL_Init(SDL_INIT_EVERYTHING) < 0) {
    std::cerr << "SDL init failed: " << SDL_GetError() << std::endl;
    return 1;
  }

  SDL_Window *window = SDL_CreateWindow("Fluid Simulation", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                        SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
  if (window == NULL) {
    std::cerr << "could not create window: " << SDL_GetError() << std::endl;
    return 1;
  }

  SDL_Renderer *sdlRenderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
  SDL_SetRenderDrawColor(sdlRenderer, 0, 0, 0, 255);
  SDL_RenderClear(sdlRenderer);

  // Init scene
  Sink sink;
  Trough trough;
  Vec3 camera(1.25, 1.5, 1.8);
  Vec3 center(0, 0, 0);
  Vec3 up(0, 1, 0);

  float fov = 45.0f * M_PI / 180.0f;
  float aspectRatio = (float)SCREEN_WIDTH / SCREEN_HEIGHT;
  float near = 0.1f;
  float far = 10.0f;
  Mat4 projMat = proj(fov, near, far, aspectRatio);
  Mat4 viewMat = lookat(camera, center, up);
  transformMat = projMat * viewMat;
  transformSinkPoints(sink);
  transformTroughPoints(trough, sink);

  Uint32 frameStart = 0;
  Uint32 frameTime = 0;
  int count = 0;
  while (count < FRAMES) {
    Uint32 currentTime = SDL_GetTicks();
    frameTime = frameStart > 0 ? currentTime - frameStart : 0;
    frameStart = currentTime;
    float fps = frameTime > 0 ? 1000.0f / frameTime : 0;
    printf("Frame time: %u ms (%.1f FPS), count = %d\n", frameTime, fps, count);

    SDL_Event e;
    if (SDL_PollEvent(&e)) {
      if (e.type == SDL_QUIT) {
        break;
      }
    }
    SDL_SetRenderDrawColor(sdlRenderer, 0, 0, 0, 255);
    SDL_RenderClear(sdlRenderer);
    renderSinkBackFaces(sdlRenderer, sink);
    SDL_SetRenderDrawColor(sdlRenderer, 18, 149, 217, 255);
    drawPoints(frameData2->frames[count], sdlRenderer);
    renderSinkFrontFaces(sdlRenderer, sink);
    renderTroughBackFaces(sdlRenderer, trough);
    renderTroughFrontFaces(sdlRenderer, trough);
    SDL_RenderPresent(sdlRenderer);
    count++;
  }
  delete frameData2;

  SDL_DestroyRenderer(sdlRenderer);
  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}
