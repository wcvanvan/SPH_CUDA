// sph_graph.cu: Implements simulateAllFramesCuda using CUDA Graphs

#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include "const.h"  // FRAMES
#include "kernel.h"
#include "sph.h"
#include "sph_graph.h"

void simulateAllFramesCuda(Particle *particlesOnGPU, int particleCount, const Sink &sink, const Trough &trough,
                           float mass, float *transformMatOnGPU, int *cellStart, int *cellEnd, Vec2 *screenPosOnGPU,
                           Vec2 *screenPosOnCPU, float POLY6, float VISCOSITY_LAPLACIAN, float WEIGHT_AT_0,
                           Vec2 *allFramesOnGPU) {
  cudaStream_t stream;
  cudaStreamCreate(&stream);

  cudaGraph_t graph;
  cudaGraphExec_t graphExec;
  cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);

  int blockDim = 32;
  int gridDim = (particleCount + blockDim - 1) / blockDim;
  float cellSize = KERNEL_RADIUS;
  int gridDimX = (int)ceil(sink.xLen / cellSize);
  int gridDimY = (int)ceil(sink.yLen / cellSize);
  int gridDimZ = (int)ceil(sink.zLen / cellSize);

  computeDensityPressureSorted<<<gridDim, blockDim, 0, stream>>>(particlesOnGPU, particleCount, mass, cellStart,
                                                                 cellEnd, cellSize, gridDimX, gridDimY, gridDimZ,
                                                                 sink.xLen, sink.yLen, sink.zLen, POLY6, WEIGHT_AT_0);
  computeAccelSorted<<<gridDim, blockDim, 0, stream>>>(particlesOnGPU, particleCount, mass, cellStart, cellEnd,
                                                       cellSize, gridDimX, gridDimY, gridDimZ, sink.xLen, sink.yLen,
                                                       sink.zLen, VISCOSITY_LAPLACIAN);
  integration<<<gridDim, blockDim, 0, stream>>>(particlesOnGPU, particleCount, sink.xLen, sink.yLen, sink.zLen,
                                                trough.zLen, trough.slope, trough.intercept, trough.normal);
  coordTransform<<<gridDim, blockDim, 0, stream>>>(particlesOnGPU, particleCount, transformMatOnGPU, screenPosOnGPU);
  cudaMemcpyAsync(allFramesOnGPU, screenPosOnGPU, sizeof(Vec2) * particleCount, cudaMemcpyDeviceToDevice, stream);

  cudaStreamEndCapture(stream, &graph);
  cudaGraphInstantiate(&graphExec, graph, nullptr, nullptr, 0);

  size_t numNodes = 0;
  cudaGraphGetNodes(graph, nullptr, &numNodes);
  std::vector<cudaGraphNode_t> nodes(numNodes);
  cudaGraphGetNodes(graph, nodes.data(), &numNodes);
  cudaGraphNode_t memcpyNode = nullptr;
  for (auto n : nodes) {
    cudaGraphNodeType t;
    cudaGraphNodeGetType(n, &t);
    if (t == cudaGraphNodeTypeMemcpy) {
      memcpyNode = n;
      break;
    }
  }
  if (!memcpyNode) {
    std::cerr << "Memcpy node not found" << std::endl;
    return;
  }

  size_t frameBytes = sizeof(Vec2) * particleCount;
  for (int f = 0; f < FRAMES; ++f) {
    sortParticles(particlesOnGPU, particleCount, cellStart, cellEnd, cellSize, sink.xLen, sink.yLen, sink.zLen,
                  gridDimX, gridDimY, gridDimZ,
                  /*stream=*/0);
    cudaMemcpy3DParms copyParams = {};
    copyParams.srcPtr = make_cudaPitchedPtr(screenPosOnGPU, frameBytes, frameBytes, 1);
    Vec2 *dstPtr = allFramesOnGPU + size_t(f) * particleCount;
    copyParams.dstPtr = make_cudaPitchedPtr(dstPtr, frameBytes, frameBytes, 1);
    copyParams.extent = make_cudaExtent(frameBytes, 1, 1);
    copyParams.kind = cudaMemcpyDeviceToDevice;
    cudaGraphExecMemcpyNodeSetParams(graphExec, memcpyNode, &copyParams);

    cudaGraphLaunch(graphExec, stream);
  }

  cudaStreamSynchronize(stream);
  cudaMemcpy(screenPosOnCPU, allFramesOnGPU, frameBytes * FRAMES, cudaMemcpyDeviceToHost);

  cudaGraphExecDestroy(graphExec);
  cudaGraphDestroy(graph);
  cudaStreamDestroy(stream);
}
