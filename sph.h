#pragma once

#include <stdint.h>
#include "const.h"

uint32_t * gpuAlloc();
void gpuFree(void* host_ptr);
void SPHSimulation(uint32_t* buf);
