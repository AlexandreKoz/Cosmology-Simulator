#pragma once

#include <cstdint>

namespace cosmosim::core {

struct OpenMpRuntimeInfo {
  bool compiled = false;
  int requested_threads = 1;
  int configured_threads = 1;
  int maximum_threads = 1;
};

[[nodiscard]] bool openMpBackendAvailable() noexcept;
[[nodiscard]] int openMpMaximumThreads() noexcept;
[[nodiscard]] int configuredOpenMpThreadCount() noexcept;

// requested_threads == 0 selects the OpenMP runtime default/max thread count.
// Values greater than one fail when the binary was built without OpenMP.
void configureOpenMpThreads(int requested_threads);

[[nodiscard]] OpenMpRuntimeInfo openMpRuntimeInfo() noexcept;

}  // namespace cosmosim::core
