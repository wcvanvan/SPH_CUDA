#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <string.h>
#include <SDL2/SDL.h>
#include "const.h"
#include "sph.h"
#include "vec.h"

Mat4 viewMat;
Mat4 projMat;

Vec2 getScreenCoord(Vec3 worldCoord)
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

    Vec3 eye(5, 10, 10);
    Vec3 center(0, 0, 0);
    Vec3 up(0, 1, 0);

    float fov = 45.0f * PI / 180.0f;
    float aspectRatio = (float)SCREEN_WIDTH / SCREEN_HEIGHT;
    float near = 0.1f;
    float far = 10.0f;
    projMat = proj(fov, near, far, aspectRatio);
    viewMat = lookat(eye, center, up);

    // Init particles randomly
    int particleCount = 10;
    Particle *particles = SPHInit();
    for (int i = 0; i < particleCount; i++)
    {
        particles[i].position.x = 0;
        particles[i].position.y = 0;
        particles[i].position.z = 5 * (float)i / particleCount;
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
            updateSimulation(particles, particleCount);
            for (int i = 0; i < particleCount; i++)
            {
                Vec2 screenCoord = getScreenCoord(particles[i].position);
                if (screenCoord.x == -1) {
                    continue;
                }
                SDL_SetRenderDrawColor(sdlRenderer, 18, 149, 217, 255);
                SDL_RenderDrawPoint(sdlRenderer, screenCoord.x, screenCoord.y);
            }
            next_frame_time += time_step;
            SDL_RenderPresent(sdlRenderer);
            SDL_SetRenderDrawColor(sdlRenderer, 0, 0, 0, 255);
            SDL_RenderClear(sdlRenderer);
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
