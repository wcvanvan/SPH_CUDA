CXX=g++
NVCC=nvcc

CUDA_VERSION=12.4
CUDA_PATH=/usr/local/cuda
SDL2_PATH=$(HOME)/private/SDL-release-2.32.4
SDL2_GFX_PATH=$(HOME)/private/SDL2_gfx-1.0.4

INC_DIRS=$(CUDA_PATH)/include $(SDL2_PATH)/include $(SDL2_GFX_PATH)
INC=$(foreach d, $(INC_DIRS), -I$d)

LIB_DIRS=$(CUDA_PATH)/lib64 $(CUDA_PATH)/lib $(SDL2_PATH)/build $(SDL2_GFX_PATH)/.libs
LIBS=$(foreach d, $(LIB_DIRS), -L$d)

LIBFLAGS=-lSDL2 -lSDL2_gfx -lcudart

CPP_SRCS=main.cpp compute.cpp visualizer.cpp vec.cpp
CU_SRCS=sph.cu
OBJS=$(CPP_SRCS:.cpp=.o) $(CU_SRCS:.cu=.o)

# We do not specify compute mode or visual mode in compile stage.
# Instead, we use a flag to determine the mode at runtime.
# Usage: 
# Step 1: make -> generate SPH_CUDA file.
# Step 2: ./SPH_CUDA -> First generate the temp file, then visualize immediately.
# or
# Step 2: ./SPH_CUDA -c -> Only generate the temp file, then use ./SPH_CUDA -v to visualize separately.

all: SPH_CUDA

SPH_CUDA: $(OBJS)
	$(CXX) -o $@ $^ $(LIBS) $(LIBFLAGS)

%.o: %.cpp
	$(CXX) $(INC) -c -o $@ $<

%.o: %.cu
	$(NVCC) $(INC) -c -o $@ $<

clean:
	rm -f *.o particles.dat SPH_CUDA
