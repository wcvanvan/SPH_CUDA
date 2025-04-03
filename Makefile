CXX=g++
NVCC=nvcc
CUDA_VERSION=12.4  # Updated to 12.4 as you mentioned
# Find CUDA installation
CUDA_PATH=/usr/local/cuda
# Find SDL2 path
SDL2_PATH=$(HOME)/private/SDL-release-2.32.4

# Include directories
INC_DIRS=$(CUDA_PATH)/include $(SDL2_PATH)/include
INC=$(foreach d, $(INC_DIRS), -I$d)

# Library directories
LIB_DIRS=$(CUDA_PATH)/lib64 $(CUDA_PATH)/lib $(SDL2_PATH)/build
LIBS=$(foreach d, $(LIB_DIRS), -L$d)

all: SPH_CUDA

SPH_CUDA: renderer.o sph.o
	${CXX} -o SPH_CUDA renderer.o sph.o $(LIBS) -lSDL2 -lcudart

renderer.o: renderer.cpp
	${CXX} $(INC) -c -o renderer.o renderer.cpp

sph.o: sph.cu sph.h
	$(NVCC) $(INC) -c -o sph.o sph.cu

run: SPH_CUDA
	LD_LIBRARY_PATH=$(SDL2_PATH)/build:$(CUDA_PATH)/lib64:$$LD_LIBRARY_PATH ./SPH_CUDA

clean:
	rm -f *.o SPH_CUDA