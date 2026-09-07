#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>
#include <stdexcept>
#include <span>
#include <vector>

namespace cosmosim::parallel::internal {

struct SfcCutPoint {
  std::uint64_t key = 0U;
  std::uint64_t entity_id = 0U;
};

[[nodiscard]] inline bool lessSfcPoint(
    const SfcCutPoint& lhs, const SfcCutPoint& rhs) noexcept {
  return lhs.key != rhs.key ? lhs.key < rhs.key : lhs.entity_id < rhs.entity_id;
}

// Exact feasibility-first repair of sampled work cuts. The callbacks expose
// inclusive global byte prefixes and actual predecessor points, not a globally
// materialized particle list. A distributed caller supplies bounded collective
// implementations; tests can supply deterministic synthetic shards. Equal
// (key,ID) groups remain indivisible. Nonnegative bytes are required.
//
// Right-greedy suffix packing establishes mandatory boundaries. For each rank
// the work-preferred cut is clamped between its suffix-feasibility boundary and
// the largest prefix that fits the remaining rank memory. The construction is
// exact for ordered indivisible nonnegative-byte groups, including zero-byte
// groups and empty ranks. It never changes the relative order of entities.
template <class GlobalPrefix, class FirstPrefixAtLeast, class GlobalNeighbor>
[[nodiscard]] std::vector<SfcCutPoint> repairMemoryConstrainedSfcCuts(
    std::span<const SfcCutPoint> proposed_cuts,
    std::uint64_t persistent_limit,
    std::uint64_t total_memory,
    SfcCutPoint first_point,
    SfcCutPoint last_point,
    GlobalPrefix&& global_prefix,
    FirstPrefixAtLeast&& first_prefix_at_least,
    GlobalNeighbor&& global_neighbor,
    std::span<SfcCutPoint> mandatory,
    std::vector<SfcCutPoint>& repaired_cuts) {
  const std::size_t rank_count = proposed_cuts.size() + 1U;
  if (mandatory.size() != proposed_cuts.size() ||
      repaired_cuts.capacity() < proposed_cuts.size()) {
    throw std::invalid_argument("SFC repair metadata must be preallocated before collective execution");
  }
  if (persistent_limit == 0U) {
    throw std::invalid_argument("SFC persistent rank allowance must be positive");
  }
  if (lessSfcPoint(last_point, first_point)) {
    throw std::invalid_argument("SFC first/last point order is invalid");
  }
  const auto prefixAtLeast = [&](std::uint64_t target) {
    return target == 0U ? first_point : first_prefix_at_least(target);
  };
  const auto maximumCut = [&](std::uint64_t allowed_prefix) -> std::optional<SfcCutPoint> {
    if (allowed_prefix >= total_memory) { return last_point; }
    const SfcCutPoint first_excess = first_prefix_at_least(allowed_prefix + 1U);
    return global_neighbor(first_excess, true);
  };
  std::fill(mandatory.begin(), mandatory.end(), first_point);
  std::uint64_t right_prefix = total_memory;
  for (std::size_t remaining = 1U; remaining < rank_count; ++remaining) {
    const std::uint64_t target = right_prefix > persistent_limit
        ? right_prefix - persistent_limit : 0U;
    const SfcCutPoint boundary = prefixAtLeast(target);
    mandatory[rank_count - remaining - 1U] = boundary;
    right_prefix = global_prefix(boundary);
  }
  repaired_cuts.clear();
  SfcCutPoint previous{};
  std::uint64_t previous_prefix = 0U;
  for (std::size_t rank = 0U; rank + 1U < rank_count; ++rank) {
    const std::uint64_t allowed_prefix = previous_prefix >
        std::numeric_limits<std::uint64_t>::max() - persistent_limit
        ? std::numeric_limits<std::uint64_t>::max()
        : previous_prefix + persistent_limit;
    const auto maximum = maximumCut(allowed_prefix);
    if (!maximum.has_value() || lessSfcPoint(*maximum, mandatory[rank])) {
      throw std::runtime_error(
          "SFC decomposition cannot satisfy hard rank memory ceiling with the available ranks");
    }
    SfcCutPoint cut = proposed_cuts[rank];
    if (lessSfcPoint(cut, mandatory[rank])) { cut = mandatory[rank]; }
    if (lessSfcPoint(*maximum, cut)) { cut = *maximum; }
    if (rank != 0U && lessSfcPoint(cut, previous)) { cut = previous; }
    repaired_cuts.push_back(cut);
    previous = cut;
    previous_prefix = global_prefix(cut);
  }
  if (total_memory - previous_prefix > persistent_limit) {
    throw std::runtime_error(
        "SFC decomposition cannot satisfy hard rank memory ceiling with the available ranks");
  }
  return std::move(repaired_cuts);
}

}  // namespace cosmosim::parallel::internal
