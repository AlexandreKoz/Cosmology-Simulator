#include "cosmosim/core/process_memory.hpp"

#include <cassert>

int main() {
  const cosmosim::core::ProcessMemoryObservation observation =
      cosmosim::core::observeProcessMemory();
#if defined(__linux__)
  assert(observation.current_rss_bytes.has_value());
  assert(*observation.current_rss_bytes > 0U);
  assert(observation.peak_rss_bytes.has_value());
  assert(*observation.peak_rss_bytes >= *observation.current_rss_bytes);
  if (observation.pss_bytes.has_value()) {
    assert(*observation.pss_bytes > 0U);
  }
#else
  assert(!observation.current_rss_bytes.has_value());
  assert(!observation.peak_rss_bytes.has_value());
  assert(!observation.pss_bytes.has_value());
#endif
  return 0;
}
