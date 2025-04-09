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

# Clean target (doesn't require mode)
clean:
	rm -f *.o particles.dat SPH_CUDA Visualizer

# Check mode parameter for non-clean targets
ifeq ($(filter clean,$(MAKECMDGOALS)),)
  ifndef mode
    $(error Usage: "make run mode=c" (computation) or "make run mode=v" (visualization))
  endif
  ifneq ($(mode),c)
    ifneq ($(mode),v)
      $(error Usage: "make run mode=c" (computation) or "make run mode=v" (visualization))
    endif
  endif
endif

ifeq ($(mode),v)
# Visualizer mode
all: Visualizer
Visualizer: visualizer.o vec.o
	${CXX} -o Visualizer visualizer.o vec.o $(LIBS) -lSDL2 -lSDL2_gfx
visualizer.o: visualizer.cpp
	${CXX} $(INC) -c -o visualizer.o visualizer.cpp
vec.o: vec.cpp vec.h
	${CXX} $(INC) -c -o vec.o vec.cpp
run: Visualizer
	LD_LIBRARY_PATH=$(SDL2_PATH)/build:$(SDL2_GFX_PATH)/.libs:$$LD_LIBRARY_PATH ./Visualizer
endif
ifeq ($(mode),c)
# Compute mode
all: SPH_CUDA
SPH_CUDA: main.o sph.o vec.o
	${CXX} -o SPH_CUDA main.o sph.o vec.o $(LIBS) -lSDL2 -lSDL2_gfx -lcudart
main.o: main.cpp
	${CXX} $(INC) -c -o main.o main.cpp
sph.o: sph.cu sph.h vec.h
	$(NVCC) $(INC) -c -o sph.o sph.cu
vec.o: vec.cpp vec.h
	${CXX} $(INC) -c -o vec.o vec.cpp
run: SPH_CUDA
	LD_LIBRARY_PATH=$(SDL2_PATH)/build:$(SDL2_GFX_PATH)/.libs:$(CUDA_PATH)/lib64:$$LD_LIBRARY_PATH ./SPH_CUDA $(mode)
endif