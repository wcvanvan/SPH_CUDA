#pragma once
#define SCREEN_WIDTH 1000
#define SCREEN_HEIGHT 800
#define SCREEN_SIZE (SCREEN_WIDTH * SCREEN_HEIGHT)
#define FPS 60
#define KERNEL_RADIUS 0.06f
#define REST_DENSITY 1000.0f
#define STIFFNESS 1.0f
#define VISCOSITY 6.5f
#define DELTA_T 0.003f
#define GRAVITY -0.5f
#define REFLECT_DAMP 0.75f
#define FRAMES 2000
#define POLY6                                                                                                 \
  (315.0f / (64.0f * M_PI *                                                                                   \
             (KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS * \
              KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS)))
#define WEIGHT_AT_0 \
  POLY6 *(KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS)
#define VISCOSITY_LAPLACIAN \
  (45.0 / (M_PI * (KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS)))
#define SPIKY_GRAD (-VISCOSITY_LAPLACIAN)