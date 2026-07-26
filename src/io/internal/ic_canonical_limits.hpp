#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

namespace cosmosim::io::internal {

inline void validateCanonicalSingleFileCounts(
    const std::array<std::uint64_t, 6>& counts) {
  for (std::size_t type = 0U; type < counts.size(); ++type) {
    if (counts[type] > std::numeric_limits<std::uint32_t>::max()) {
      throw std::length_error(
          "canonical single-file output cannot represent NumPart_ThisFile "
          "for PartType" + std::to_string(type) + ": count " +
          std::to_string(counts[type]) + " exceeds UINT32_MAX; use a "
          "canonical multifile writer");
    }
  }
}

}  // namespace cosmosim::io::internal
