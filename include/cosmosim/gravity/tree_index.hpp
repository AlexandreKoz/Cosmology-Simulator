#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace cosmosim::gravity {

// Compact local tree indices are intentionally 32-bit for cache/bandwidth
// efficiency. Global particle IDs and MPI ownership identities are separate.
// This alias centralizes the policy so a future wide-index build can change one
// contract rather than silently widening every hot structure today.
using TreeLocalIndex = std::uint32_t;
inline constexpr std::uint64_t k_tree_local_index_max =
    static_cast<std::uint64_t>(std::numeric_limits<TreeLocalIndex>::max());

[[nodiscard]] inline TreeLocalIndex checkedTreeLocalIndex(
    std::size_t value,
    const char* context) {
  if (value > static_cast<std::size_t>(k_tree_local_index_max)) {
    throw std::overflow_error(context);
  }
  return static_cast<TreeLocalIndex>(value);
}

}  // namespace cosmosim::gravity
