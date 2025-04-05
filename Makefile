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

all: SPH_CUDA

SPH_CUDA: renderer.o sph.o vec.o
	${CXX} -o SPH_CUDA renderer.o sph.o vec.o $(LIBS) -lSDL2 -lSDL2_gfx -lcudart

renderer.o: renderer.cpp
	${CXX} $(INC) -c -o renderer.o renderer.cpp

sph.o: sph.cu sph.h vec.h
	$(NVCC) $(INC) -c -o sph.o sph.cu

vec.o: vec.cpp vec.h
	${CXX} $(INC) -c -o vec.o vec.cpp

run: SPH_CUDA
	LD_LIBRARY_PATH=$(SDL2_PATH)/build:$(SDL2_GFX_PATH)/.libs:$(CUDA_PATH)/lib64:$$LD_LIBRARY_PATH ./SPH_CUDA

clean:
	rm -f *.o SPH_CUDA