#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <string.h>
#include <SDL2/SDL.h>
#include "SDL2_gfxPrimitives.h"
#include "const.h"
#include "sph.h"
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

class Sink
{
public:
    float xLen = 1;
    float yLen = 0.2;
    float zLen = 1;
    Vec3 sinkVertices[8] = {
        {-xLen / 2, -yLen / 2, zLen / 2},  // Bottom front left
        {xLen / 2, -yLen / 2, zLen / 2},   // Bottom front right
        {xLen / 2, -yLen / 2, -zLen / 2},  // Bottom back right
        {-xLen / 2, -yLen / 2, -zLen / 2}, // Bottom back left
        {-xLen / 2, yLen / 2, zLen / 2},   // Top front left
        {xLen / 2, yLen / 2, zLen / 2},    // Top front right
        {xLen / 2, yLen / 2, -zLen / 2},   // Top back right
        {-xLen / 2, yLen / 2, -zLen / 2}   // Top back left
    };
    int backFaces[3][4] = {
        // In order. Don't change when not necessary
        {2, 3, 7, 6}, // Back
        {0, 1, 2, 3}, // Bottom
        {3, 0, 4, 7}, // Left
    };
    int frontFaces[2][4] = {
        // In order. Don't change when not necessary
        {1, 2, 6, 5}, // Right
        {0, 1, 5, 4}, // Front
    };
    SDL_Point screenPoints[8];
    Sink()
    {
        // change from world coordinates to screen coordinates
        for (int i = 0; i < 8; i++)
        {
            Vec2 screenCoord = worldToScreen(sinkVertices[i]);
            screenPoints[i] = {(int)screenCoord.x, (int)screenCoord.y};
        }
    }
};

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

    Vec3 camera(1, 1.25, 1.5);
    Vec3 center(0, 0, 0);
    Vec3 up(0, 1, 0);

    float fov = 45.0f * PI / 180.0f;
    float aspectRatio = (float)SCREEN_WIDTH / SCREEN_HEIGHT;
    float near = 0.1f;
    float far = 5.0f;
    projMat = proj(fov, near, far, aspectRatio);
    viewMat = lookat(camera, center, up);
    
    // Init scene
    Sink sink;
    int particleCount = 5000;
    Particle *particles = initParticles(particleCount);
    srand(time(0));
    for (int i = 0; i < particleCount; i++) {
        particles[i].position.x = ((float)rand() / RAND_MAX) * sink.xLen - sink.xLen / 2;
        particles[i].position.y = ((float)rand() / RAND_MAX) * sink.yLen - sink.yLen / 2;
        particles[i].position.z = ((float)rand() / RAND_MAX) * sink.zLen - sink.zLen / 2;
    }

    uint32_t next_frame_time = SDL_GetTicks();

    while (1)
    {
        SDL_Event e;
        if (SDL_PollEvent(&e))
        {
            if (e.type == SDL_QUIT)
            {
                break;
            }
        }

        uint32_t now = SDL_GetTicks();
        if (next_frame_time <= now)
        {
            SDL_SetRenderDrawColor(sdlRenderer, 0, 0, 0, 255);
            SDL_RenderClear(sdlRenderer);
            renderSinkBackFaces(sdlRenderer, sink);
            updateSimulation(particles, particleCount);
            for (int i = 0; i < particleCount; i++)
            {
                Vec2 screenCoord = worldToScreen(particles[i].position);
                if (screenCoord.x == -1) {
                    continue;
                }
                SDL_SetRenderDrawColor(sdlRenderer, 18, 149, 217, 255);
                SDL_RenderDrawPoint(sdlRenderer, screenCoord.x, screenCoord.y);
            }
            renderSinkFrontFaces(sdlRenderer, sink);
            next_frame_time += time_step;
            SDL_RenderPresent(sdlRenderer);
        }
        else
        {
            SDL_Delay(next_frame_time - now);
        }
    }

    SDL_DestroyRenderer(sdlRenderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
