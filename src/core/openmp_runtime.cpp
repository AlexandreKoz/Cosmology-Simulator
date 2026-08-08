#include "cosmosim/core/openmp_runtime.hpp"

#include <algorithm>
#include <atomic>
#include <stdexcept>

#include "cosmosim/core/build_config.hpp"

#if COSMOSIM_HAVE_OPENMP
#include <omp.h>
#endif

namespace cosmosim::core {
namespace {
std::atomic<int> g_requested_threads{1};
std::atomic<int> g_configured_threads{1};
}

bool openMpBackendAvailable() noexcept {
#if COSMOSIM_HAVE_OPENMP
  return true;
#else
  return false;
#endif
}

int openMpMaximumThreads() noexcept {
#if COSMOSIM_HAVE_OPENMP
  return std::max(1, omp_get_max_threads());
#else
  return 1;
#endif
}

int configuredOpenMpThreadCount() noexcept {
  return g_configured_threads.load(std::memory_order_relaxed);
}

void configureOpenMpThreads(int requested_threads) {
  if (requested_threads < 0) {
    throw std::invalid_argument("OpenMP thread request must be >= 0 (0 means runtime default)");
  }
#if COSMOSIM_HAVE_OPENMP
  omp_set_dynamic(0);
  const int selected_threads = requested_threads == 0
      ? std::max(1, omp_get_max_threads())
      : requested_threads;
  if (selected_threads < 1) {
    throw std::invalid_argument("OpenMP selected thread count must be positive");
  }
  omp_set_num_threads(selected_threads);
  int observed_threads = 1;
#pragma omp parallel
  {
#pragma omp single
    observed_threads = omp_get_num_threads();
  }
  g_requested_threads.store(requested_threads, std::memory_order_relaxed);
  g_configured_threads.store(observed_threads, std::memory_order_relaxed);
#else
  if (requested_threads > 1) {
    throw std::invalid_argument(
        "OpenMP thread request exceeds one but this binary was built without OpenMP");
  }
  g_requested_threads.store(requested_threads, std::memory_order_relaxed);
  g_configured_threads.store(1, std::memory_order_relaxed);
#endif
}

OpenMpRuntimeInfo openMpRuntimeInfo() noexcept {
  return OpenMpRuntimeInfo{
      .compiled = openMpBackendAvailable(),
      .requested_threads = g_requested_threads.load(std::memory_order_relaxed),
      .configured_threads = g_configured_threads.load(std::memory_order_relaxed),
      .maximum_threads = openMpMaximumThreads(),
  };
}

}  // namespace cosmosim::core
