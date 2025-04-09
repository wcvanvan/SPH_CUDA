#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <string.h>
#include <SDL2/SDL.h>
#include "SDL2_gfxPrimitives.h"
#include "const.h"
#include "sph.h"
#include "sink.h"
#include "vec.h"

Mat4 viewMat;
Mat4 projMat;

Vec2 worldToScreen(Vec3 &worldCoord)
{
    Vec4 p(worldCoord);
    Vec4 ndcCoord = viewMat * p;
    ndcCoord = projMat * ndcCoord;
    ndcCoord.x /= ndcCoord.w;
    ndcCoord.y /= ndcCoord.w;

    if (ndcCoord.x < -1.0f || ndcCoord.x > 1.0f || ndcCoord.y < -1.0f || ndcCoord.y > 1.0f)
        return {-1, -1};

    float screenX = std::max(0.0f, std::min(1.0f, (ndcCoord.x + 1.0f) * 0.5f)) * SCREEN_WIDTH;
    float screenY = std::max(0.0f, std::min(1.0f, (1.0f - ndcCoord.y) * 0.5f)) * SCREEN_HEIGHT;
    return {screenX, screenY};
}

Sink::Sink()
{
    // change from world coordinates to screen coordinates
    for (int i = 0; i < 8; i++)
    {
        Vec2 screenCoord = worldToScreen(sinkVertices[i]);
        screenPoints[i] = {(int)screenCoord.x, (int)screenCoord.y};
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

int main(int argc, char *args[])
{

    uint32_t time_step = 1000. / FPS;

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

    Vec3 camera(1.25, 1.5, 1.8);
    Vec3 center(0, 0, 0);
    Vec3 up(0, 1, 0);

    float fov = 45.0f * M_PI / 180.0f;
    float aspectRatio = (float)SCREEN_WIDTH / SCREEN_HEIGHT;
    float near = 0.1f;
    float far = 10.0f;
    projMat = proj(fov, near, far, aspectRatio);
    viewMat = lookat(camera, center, up);
    Mat4 transformMat = projMat * viewMat;
    float *transformMatOnGPU = allocateMatOnGPU(transformMat);

    // Init scene
    Sink sink;
    int particleCount = 0;
    float mass = 1.0f;
    Particle *particles = initParticles(particleCount, mass, sink );
    initSimulation(particles, particleCount, sink, mass, transformMatOnGPU);

    // Uint32 frameStart = 0;
    // Uint32 frameTime = 0;

    while (1)
    {
        // Uint32 currentTime = SDL_GetTicks();
        // frameTime = frameStart > 0 ? currentTime - frameStart : 0;
        // frameStart = currentTime;
        // float fps = frameTime > 0 ? 1000.0f / frameTime : 0;
        // printf("Frame time: %u ms (%.1f FPS)\n", frameTime, fps);
        
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
        updateSimulation(particles, particleCount, sink, mass, transformMatOnGPU);
        for (int i = 0; i < particleCount; i++)
        {
            if (particles[i].screenPos.x == -1)
            {
                continue;
            }
            SDL_SetRenderDrawColor(sdlRenderer, 18, 149, 217, 255);
            SDL_RenderDrawPoint(sdlRenderer, particles[i].screenPos.x, particles[i].screenPos.y);
        }
        renderSinkFrontFaces(sdlRenderer, sink);
        SDL_RenderPresent(sdlRenderer);
    }

    SDL_DestroyRenderer(sdlRenderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
