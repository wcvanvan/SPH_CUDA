#include <stdio.h>
#include <iostream>
#include <string.h>
#include <SDL2/SDL.h>
#include "const.h"
#include "sph.h"

void render(void *buf)
{
    SPHSimulation((uint32_t *)buf);
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
        "main",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        SCREEN_WIDTH, SCREEN_HEIGHT,
        SDL_WINDOW_SHOWN);
    if (window == NULL)
    {
        std::cerr << "could not create window: " << SDL_GetError() << std::endl;
        return 1;
    }

    uint32_t *buf = gpuAlloc();

    SDL_Surface *screen = SDL_CreateRGBSurfaceFrom(buf, SCREEN_WIDTH, SCREEN_HEIGHT, 32, SCREEN_WIDTH * sizeof(uint32_t),
                                                   0x00FF0000,
                                                   0x0000FF00,
                                                   0x000000FF,
                                                   0xFF000000);

    if (screen == NULL)
    {
        SDL_Log("SDL_CreateRGBSurfaceFrom() failed: %s", SDL_GetError());
        return 1;
    }

    SDL_Renderer *sdlRenderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);

    SDL_Texture *sdlTexture = SDL_CreateTexture(sdlRenderer,
                                                SDL_PIXELFORMAT_ARGB8888,
                                                SDL_TEXTUREACCESS_STREAMING | SDL_TEXTUREACCESS_TARGET,
                                                SCREEN_WIDTH, SCREEN_HEIGHT);

    if (sdlTexture == NULL)
    {
        SDL_Log("SDL_Error failed: %s", SDL_GetError());
        return 1;
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
            render(buf);
            SDL_UpdateTexture(sdlTexture, NULL, screen->pixels, screen->pitch);
            SDL_RenderClear(sdlRenderer);
            SDL_RenderCopy(sdlRenderer, sdlTexture, NULL, NULL);
            SDL_RenderPresent(sdlRenderer);

            next_frame_time += time_step;
        }
        else
        {
            SDL_Delay(next_frame_time - now);
        }
    }

    gpuFree(buf);
    SDL_DestroyTexture(sdlTexture);
    SDL_DestroyRenderer(sdlRenderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
