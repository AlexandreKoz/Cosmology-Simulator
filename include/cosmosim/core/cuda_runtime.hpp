#pragma once

#include <stdexcept>

namespace cosmosim::core {

struct CudaRuntimeInfo {
  bool build_enabled = false;
  bool runtime_available = false;
  int visible_device_count = 0;
};

[[nodiscard]] CudaRuntimeInfo queryCudaRuntime();
[[nodiscard]] int selectCudaDeviceRoundRobin(int world_rank, int visible_device_count);
void setCudaDeviceOrThrow(int device_index);

}  // namespace cosmosim::core
