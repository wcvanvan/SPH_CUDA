#pragma once
#include <SDL2/SDL.h>
#include "vec.h"

class Sink
{
public:
    float xLen = 1.0f;
    float yLen = 0.2f;
    float zLen = 1.0f;
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
};