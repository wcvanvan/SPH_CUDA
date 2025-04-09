#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <string.h>
#include <assert.h>
#include <SDL2/SDL.h>
#include <SDL2_gfxPrimitives.h>
#include "const.h"
#include "sph.h"
#include "sink.h"
#include "vec.h"

Mat4 transformMat;
const char *filename = "./particles.dat";
struct FrameData {
    std::vector<std::vector<Vec2>> frames;
};

void drawPoints(std::vector<Vec2> &points, SDL_Renderer *sdlRenderer) {
    for (auto && v : points) {
        SDL_RenderDrawPoint(sdlRenderer, static_cast<int>(v.x), static_cast<int>(v.y));
    }
}

Vec2 worldToScreen(Vec3 &worldCoord)
{
    Vec4 p(worldCoord);
    Vec4 ndcCoord = transformMat * p;
    ndcCoord.x /= ndcCoord.w;
    ndcCoord.y /= ndcCoord.w;

    if (ndcCoord.x < -1.0f || ndcCoord.x > 1.0f || ndcCoord.y < -1.0f || ndcCoord.y > 1.0f)
        return {-1, -1};

    float screenX = std::max(0.0f, std::min(1.0f, (ndcCoord.x + 1.0f) * 0.5f)) * SCREEN_WIDTH;
    float screenY = std::max(0.0f, std::min(1.0f, (1.0f - ndcCoord.y) * 0.5f)) * SCREEN_HEIGHT;
    return {screenX, screenY};
}

void transformSinkPoints(Sink &sink)
{
    // change from world coordinates to screen coordinates
    for (int i = 0; i < 8; i++)
    {
        Vec2 screenCoord = worldToScreen(sink.sinkVertices[i]);
        sink.screenPoints[i] = {(int)screenCoord.x, (int)screenCoord.y};
    }
}

void renderSinkBackFaces(SDL_Renderer *sdlRenderer, Sink &sink)
{
    for (int i = 0; i < 3; i++)
    {
        Sint16 vx[4], vy[4];
        for (int j = 0; j < 4; j++)
        {
            vx[j] = (Sint16)sink.screenPoints[sink.backFaces[i][j]].x;
            vy[j] = (Sint16)sink.screenPoints[sink.backFaces[i][j]].y;
        }
        filledPolygonRGBA(sdlRenderer, vx, vy, 4, 220, 220, 220, 255);
        polygonRGBA(sdlRenderer, vx, vy, 4, 0, 0, 0, 255);
    }
}

void renderSinkFrontFaces(SDL_Renderer *sdlRenderer, Sink &sink)
{
    for (int i = 0; i < 2; i++)
    {
        Sint16 vx[4], vy[4];
        for (int j = 0; j < 4; j++)
        {
            vx[j] = (Sint16)sink.screenPoints[sink.frontFaces[i][j]].x;
            vy[j] = (Sint16)sink.screenPoints[sink.frontFaces[i][j]].y;
        }
        filledPolygonRGBA(sdlRenderer, vx, vy, 4, 220, 220, 220, 255);
        polygonRGBA(sdlRenderer, vx, vy, 4, 0, 0, 0, 255);
    }
}

FrameData *readAllFramesFromFile(int &particleCount) {
    FrameData *frameData = new FrameData();
    FILE* file = fopen(filename, "r");
    if (!file) {
        std::cerr << "Error opening file" << std::endl;
        return frameData;
    }
    
    // Read particle count
    char line[256];
    if (fgets(line, sizeof(line), file) && sscanf(line, "%d", &particleCount) == 1) {
        std::cout << "Particle Count: " << particleCount << std::endl;
    } else {
        std::cout << "Fail to read particle count" << std::endl;
        fclose(file);
        return frameData;
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
        frameData->frames.push_back(framePositions);
    }
    
    fclose(file);
    return frameData;
}

int main()
{
    int particleCount = 0;
    std::cout << "Reading frame data" << std::endl;
    FrameData *frameData = readAllFramesFromFile(particleCount);
    std::cout << "Finished reading frame data" << std::endl;
    if (SDL_Init(SDL_INIT_EVERYTHING) < 0)
    {
        std::cerr << "SDL init failed: " << SDL_GetError() << std::endl;
        return 1;
    }

    SDL_Window *window = SDL_CreateWindow(
        "Fluid Simulation",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        SCREEN_WIDTH, SCREEN_HEIGHT,
        SDL_WINDOW_SHOWN);
    if (window == NULL)
    {
        std::cerr << "could not create window: " << SDL_GetError() << std::endl;
        return 1;
    }

    SDL_Renderer *sdlRenderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    SDL_SetRenderDrawColor(sdlRenderer, 0, 0, 0, 255);
    SDL_RenderClear(sdlRenderer);

    // Init scene
    Sink sink;
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

    Uint32 frameStart = 0;
    Uint32 frameTime = 0;
    int count = 0;
    while (count < FRAMES)
    {
        Uint32 currentTime = SDL_GetTicks();
        frameTime = frameStart > 0 ? currentTime - frameStart : 0;
        frameStart = currentTime;
        float fps = frameTime > 0 ? 1000.0f / frameTime : 0;
        printf("Frame time: %u ms (%.1f FPS)\n", frameTime, fps);

        SDL_Event e;
        if (SDL_PollEvent(&e))
        {
            if (e.type == SDL_QUIT)
            {
                break;
            }
        }
        SDL_SetRenderDrawColor(sdlRenderer, 0, 0, 0, 255);
        SDL_RenderClear(sdlRenderer);
        renderSinkBackFaces(sdlRenderer, sink);
        SDL_SetRenderDrawColor(sdlRenderer, 18, 149, 217, 255);
        drawPoints(frameData->frames[count], sdlRenderer);
        renderSinkFrontFaces(sdlRenderer, sink);
        SDL_RenderPresent(sdlRenderer);
        count++;
    }
    delete frameData;

    SDL_DestroyRenderer(sdlRenderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}