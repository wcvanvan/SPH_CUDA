CXX=g++
NVCC=nvcc
CXXFLAGS += -Iinclude

CUDA_VERSION=12.4
CUDA_PATH=/usr/local/cuda
SDL2_PATH=$(HOME)/private/SDL-release-2.32.4
SDL2_GFX_PATH=$(HOME)/private/SDL2_gfx-1.0.4

INC_DIRS=$(CUDA_PATH)/include $(SDL2_PATH)/include $(SDL2_GFX_PATH) include
INC=$(foreach d, $(INC_DIRS), -I$d)

LIB_DIRS=$(CUDA_PATH)/lib64 $(CUDA_PATH)/lib $(SDL2_PATH)/build $(SDL2_GFX_PATH)/.libs
LIBS=$(foreach d, $(LIB_DIRS), -L$d)

LIBFLAGS=-lSDL2 -lSDL2_gfx -lcudart

SRC_DIR=src
INC_DIR=include

CPP_SRCS=$(wildcard $(SRC_DIR)/*.cpp) main.cpp
CU_SRCS=$(wildcard $(SRC_DIR)/*.cu)
OBJS=$(CPP_SRCS:.cpp=.o) $(CU_SRCS:.cu=.o)

OBJ_NAMES=$(notdir $(OBJS))

FORMAT_SRCS=$(CPP_SRCS) $(CU_SRCS) $(wildcard $(INC_DIR)/*.h) main.cpp

# We do not specify compute mode or visual mode in compile stage.
# Instead, we use a flag to determine the mode at runtime.
# Usage: 
# Step 1: make -> generate SPH_CUDA file.
# Step 2: ./SPH_CUDA -> First generate the temp file, then visualize immediately.
# or
# Step 2: ./SPH_CUDA -c -> Only generate the temp file, then use ./SPH_CUDA -v to visualize separately.

.PHONY: all
all: format SPH_CUDA

SPH_CUDA: $(OBJ_NAMES)
	$(CXX) -o $@ $^ $(LIBS) $(LIBFLAGS)

%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c -o $@ $<

%.o: $(SRC_DIR)/%.cu
	$(NVCC) $(INC) -c -o $@ $<

main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c -o $@ $<

.PHONY: format
format:
	@echo "Formatting source files..."
	clang-format -i $(FORMAT_SRCS)

.PHONY: clean
clean:
	rm -f *.o particles.dat SPH_CUDA