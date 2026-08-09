#pragma once

#include <cstdint>
#include <optional>

#if defined(__linux__) || defined(__APPLE__)
#include <sys/resource.h>
#include <unistd.h>
#endif

#if defined(__linux__)
#include <fstream>
#endif

namespace cosmosim::bench {

struct ProcessMemorySample {
  std::optional<std::uint64_t> current_rss_bytes;
  std::optional<std::uint64_t> peak_rss_bytes;
};

[[nodiscard]] inline ProcessMemorySample sampleProcessMemory() {
  ProcessMemorySample sample;
#if defined(__linux__)
  std::ifstream statm("/proc/self/statm");
  std::uint64_t total_pages = 0;
  std::uint64_t resident_pages = 0;
  if (statm >> total_pages >> resident_pages) {
    (void)total_pages;
    const long page_size = ::sysconf(_SC_PAGESIZE);
    if (page_size > 0) {
      const auto page_bytes = static_cast<std::uint64_t>(page_size);
      if (resident_pages <= UINT64_MAX / page_bytes) {
        sample.current_rss_bytes = resident_pages * page_bytes;
      }
    }
  }
#endif
#if defined(__linux__) || defined(__APPLE__)
  struct rusage usage {};
  if (::getrusage(RUSAGE_SELF, &usage) == 0) {
#if defined(__APPLE__)
    sample.peak_rss_bytes = static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    const auto peak_kib = static_cast<std::uint64_t>(usage.ru_maxrss);
    if (peak_kib <= UINT64_MAX / 1024U) {
      sample.peak_rss_bytes = peak_kib * 1024U;
    }
#endif
  }
#endif
  return sample;
}

}  // namespace cosmosim::bench
