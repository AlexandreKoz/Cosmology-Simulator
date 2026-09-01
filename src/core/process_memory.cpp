#include "cosmosim/core/process_memory.hpp"

#include <charconv>
#include <cstdint>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

namespace cosmosim::core {
namespace {

[[nodiscard]] std::optional<std::uint64_t> checkedKibToBytes(
    std::uint64_t kib) noexcept {
  constexpr std::uint64_t k_bytes_per_kib = 1024U;
  if (kib > std::numeric_limits<std::uint64_t>::max() / k_bytes_per_kib) {
    return std::nullopt;
  }
  return kib * k_bytes_per_kib;
}

[[nodiscard]] std::optional<std::uint64_t> parseProcStatusKibField(
    std::string_view line,
    std::string_view key) noexcept {
  if (!line.starts_with(key)) {
    return std::nullopt;
  }
  std::string_view value = line.substr(key.size());
  while (!value.empty() && (value.front() == ' ' || value.front() == '\t')) {
    value.remove_prefix(1U);
  }
  const std::size_t end = value.find_first_of(" \t");
  const std::string_view number = value.substr(0U, end);
  if (number.empty()) {
    return std::nullopt;
  }
  std::uint64_t kib = 0U;
  const auto [ptr, ec] = std::from_chars(
      number.data(), number.data() + number.size(), kib);
  if (ec != std::errc{} || ptr != number.data() + number.size()) {
    return std::nullopt;
  }
  return checkedKibToBytes(kib);
}

#if defined(__linux__)
[[nodiscard]] std::pair<std::optional<std::uint64_t>, std::optional<std::uint64_t>>
readLinuxStatusRss() noexcept {
  try {
    std::ifstream status("/proc/self/status");
    if (!status) {
      return {};
    }
    std::optional<std::uint64_t> rss;
    std::optional<std::uint64_t> hwm;
    std::string line;
    while (std::getline(status, line)) {
      if (!rss.has_value()) {
        rss = parseProcStatusKibField(line, "VmRSS:");
      }
      if (!hwm.has_value()) {
        hwm = parseProcStatusKibField(line, "VmHWM:");
      }
      if (rss.has_value() && hwm.has_value()) {
        break;
      }
    }
    return {rss, hwm};
  } catch (...) {
    return {};
  }
}

[[nodiscard]] std::optional<std::uint64_t> readLinuxPss() noexcept {
  try {
    std::ifstream rollup("/proc/self/smaps_rollup");
    if (!rollup) {
      return std::nullopt;
    }
    std::string line;
    while (std::getline(rollup, line)) {
      if (line.starts_with("Pss:")) {
        return parseProcStatusKibField(line, "Pss:");
      }
    }
  } catch (...) {
  }
  return std::nullopt;
}
#endif

}  // namespace

ProcessMemoryObservation observeProcessMemory() noexcept {
  ProcessMemoryObservation observation;
#if defined(__linux__)
  const auto [rss, hwm] = readLinuxStatusRss();
  observation.current_rss_bytes = rss;
  observation.peak_rss_bytes = hwm;
  observation.pss_bytes = readLinuxPss();
#endif
  return observation;
}

}  // namespace cosmosim::core
