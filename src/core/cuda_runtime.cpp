#include "cosmosim/core/cuda_runtime.hpp"

#include <stdexcept>
#include <string>

#if defined(COSMOSIM_ENABLE_CUDA) && COSMOSIM_ENABLE_CUDA
#include <cuda_runtime.h>
#endif

namespace cosmosim::core {

CudaRuntimeInfo queryCudaRuntime() {
  CudaRuntimeInfo info;
#if defined(COSMOSIM_ENABLE_CUDA) && COSMOSIM_ENABLE_CUDA
  info.build_enabled = true;
  int device_count = 0;
  const cudaError_t status = cudaGetDeviceCount(&device_count);
  if (status == cudaSuccess && device_count > 0) {
    info.runtime_available = true;
    info.visible_device_count = device_count;
  } else if (status != cudaSuccess) {
    info.diagnostic = std::string("cudaGetDeviceCount failed: ") + cudaGetErrorString(status);
  } else {
    info.diagnostic = "CUDA runtime reported zero visible devices";
  }
#else
  info.build_enabled = false;
#endif
  return info;
}

int selectCudaDeviceRoundRobin(int local_rank, int visible_device_count) {
  if (local_rank < 0) {
    throw std::invalid_argument("local_rank must be non-negative for CUDA device selection");
  }
  if (visible_device_count <= 0) {
    throw std::invalid_argument("visible_device_count must be positive for CUDA device selection");
  }
  return local_rank % visible_device_count;
}

void setCudaDeviceOrThrow(int device_index) {
#if defined(COSMOSIM_ENABLE_CUDA) && COSMOSIM_ENABLE_CUDA
  if (device_index < 0) {
    throw std::invalid_argument("CUDA device index must be non-negative");
  }
  const cudaError_t status = cudaSetDevice(device_index);
  if (status != cudaSuccess) {
    throw std::runtime_error(
        std::string("cudaSetDevice failed for device ") + std::to_string(device_index) +
        ": " + cudaGetErrorString(status));
  }
#else
  (void)device_index;
  throw std::runtime_error("setCudaDeviceOrThrow requires COSMOSIM_ENABLE_CUDA=ON in this build");
#endif
}

}  // namespace cosmosim::core
