#pragma once

#include <cstdint>
#include <optional>

namespace cosmosim::core {

// OS-observed process residency. These measurements are diagnostic evidence,
// not allocation authority. Unsupported/unavailable fields remain std::nullopt
// rather than being represented as zero bytes.
struct ProcessMemoryObservation {
  std::optional<std::uint64_t> current_rss_bytes;
  std::optional<std::uint64_t> peak_rss_bytes;
  std::optional<std::uint64_t> pss_bytes;
};

// Samples the current process at an explicit workflow/diagnostic boundary.
// The function is noexcept by contract: unavailable or malformed platform
// sources degrade to std::nullopt fields and never affect simulation science.
[[nodiscard]] ProcessMemoryObservation observeProcessMemory() noexcept;

}  // namespace cosmosim::core
