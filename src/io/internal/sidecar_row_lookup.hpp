#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory_resource>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>
#include <string_view>

#include "cosmosim/core/checked_arithmetic.hpp"

namespace cosmosim::io::internal {

// Compact, non-authoritative index for arbitrary sidecar row ordering.
// The caller owns the bounded arena and keeps it alive through serialization.
class SidecarRowLookup final {
 public:
  SidecarRowLookup(std::span<const std::uint32_t> particle_indices,
                   std::pmr::memory_resource* resource,
                   std::string_view owner)
      : m_keys(resource), m_owner(owner) {
    if (particle_indices.size() > std::numeric_limits<std::uint32_t>::max()) {
      throw std::length_error("sidecar row index exceeds uint32 range");
    }
    m_keys.reserve(particle_indices.size());
    for (std::size_t row = 0; row < particle_indices.size(); ++row) {
      m_keys.push_back((static_cast<std::uint64_t>(particle_indices[row]) << 32U) |
                       static_cast<std::uint64_t>(row));
    }
    std::sort(m_keys.begin(), m_keys.end());
    for (std::size_t row = 1; row < m_keys.size(); ++row) {
      if ((m_keys[row] >> 32U) == (m_keys[row - 1U] >> 32U)) {
        throw std::runtime_error("snapshot writer: duplicate " + std::string(m_owner) +
                                 " sidecar particle index");
      }
    }
  }

  [[nodiscard]] std::size_t rowFor(std::uint32_t particle_index) const {
    const auto it = std::lower_bound(m_keys.begin(), m_keys.end(),
                                    static_cast<std::uint64_t>(particle_index) << 32U);
    if (it == m_keys.end() || (*it >> 32U) != particle_index) {
      throw std::runtime_error("snapshot writer: particle lacks authoritative " +
                               std::string(m_owner) + " sidecar row");
    }
    return static_cast<std::uint32_t>(*it);
  }

 private:
  std::pmr::vector<std::uint64_t> m_keys;
  std::string_view m_owner;
};

// Three independent index allocations, aligned in one monotonic arena.
// No hash nodes, buckets, or population-sized staging copies are retained.
[[nodiscard]] inline std::uint64_t sidecarLookupWorkspaceBytes(
    std::uint64_t star_count, std::uint64_t tracer_count, std::uint64_t bh_count) {
  std::uint64_t total = core::checkedMemoryBytesAdd(star_count, tracer_count,
                                                    "snapshot sidecar rows");
  total = core::checkedMemoryBytesAdd(total, bh_count, "snapshot sidecar rows");
  if (total == 0U) return 0U;
  if (total > std::numeric_limits<std::uint64_t>::max() / sizeof(std::uint64_t)) {
    throw std::overflow_error("snapshot sidecar lookup byte count overflow");
  }
  return core::checkedMemoryBytesAdd(total * sizeof(std::uint64_t), 256U,
                                     "snapshot sidecar lookup alignment");
}

}  // namespace cosmosim::io::internal
