#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "cosmosim/core/checked_arithmetic.hpp"

namespace cosmosim::gravity::internal {

struct SparseTreePmRoundPlan {
  std::size_t targets_per_peer_per_round = 1U;
  std::uint64_t aggregate_bytes_per_target = 0U;
};

[[nodiscard]] inline SparseTreePmRoundPlan planSparseTreePmRound(
    std::size_t configured_targets_per_peer,
    std::size_t peer_degree,
    std::size_t request_wire_bytes,
    std::size_t response_wire_bytes,
    std::uint64_t classic_mpi_byte_limit =
        static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
  if (configured_targets_per_peer == 0U) {
    throw std::invalid_argument("TreePM sparse round requires at least one configured target");
  }
  if (request_wire_bytes == 0U || response_wire_bytes == 0U ||
      classic_mpi_byte_limit == 0U) {
    throw std::invalid_argument("TreePM sparse round byte widths/limit must be positive");
  }
  if (peer_degree == 0U) {
    return SparseTreePmRoundPlan{
        .targets_per_peer_per_round = configured_targets_per_peer,
        .aggregate_bytes_per_target = 0U,
    };
  }

  const std::uint64_t aggregate_bytes_per_target =
      core::checkedIntegralNarrow<std::uint64_t>(
          core::checkedSizeMultiply(
              peer_degree,
              std::max(request_wire_bytes, response_wire_bytes),
              "TreePM sparse round aggregate bytes per target"),
          "TreePM sparse round aggregate byte width");
  if (aggregate_bytes_per_target > classic_mpi_byte_limit) {
    throw std::overflow_error(
        "TreePM sparse peer degree cannot fit one target request/response round "
        "within the configured classic-MPI byte limit");
  }
  const std::uint64_t representable_targets =
      classic_mpi_byte_limit / aggregate_bytes_per_target;
  if (representable_targets == 0U) {
    throw std::overflow_error("TreePM sparse round has zero representable targets");
  }
  return SparseTreePmRoundPlan{
      .targets_per_peer_per_round = std::min(
          configured_targets_per_peer,
          core::checkedIntegralNarrow<std::size_t>(
              representable_targets,
              "TreePM sparse representable target count")),
      .aggregate_bytes_per_target = aggregate_bytes_per_target,
  };
}

[[nodiscard]] inline std::uint64_t sparseTreePmPhysicalRoundCount(
    std::uint64_t logical_target_count,
    std::size_t targets_per_round) {
  if (targets_per_round == 0U) {
    throw std::invalid_argument("TreePM sparse physical round size must be positive");
  }
  if (logical_target_count == 0U) return 0U;
  const std::uint64_t round = core::checkedIntegralNarrow<std::uint64_t>(
      targets_per_round, "TreePM sparse physical round target width");
  return 1U + (logical_target_count - 1U) / round;
}

}  // namespace cosmosim::gravity::internal
