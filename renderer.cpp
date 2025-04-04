#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <SDL2/SDL.h>
#include "const.h"
#include "sph.h"

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

    // Init particles randomly
    int particleCount = 1000;
    Particle *particles = SPHInit();
    for (int i = 0; i < particleCount; i++)
    {
        particles[i].position.x = (float)rand() / RAND_MAX;
        particles[i].position.y = (float)rand() / RAND_MAX;
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
                SDL_SetRenderDrawColor(sdlRenderer, 18, 149, 217, 255);
                SDL_RenderDrawPoint(sdlRenderer, particles[i].position.x * SCREEN_WIDTH, particles[i].position.y * SCREEN_HEIGHT);
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
