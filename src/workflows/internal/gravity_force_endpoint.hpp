#pragma once

// Synchronized hierarchical force-endpoint lane transfer (2026-10-27).
// Two deliberately separate paths:
// * Ordinary restart import (remapForceEndpointById): a checkpoint may come
//   from a run whose particle rows were ordered differently, so rows are
//   matched by stable particle ID through an O(N) hash map.
// * Transactional rollback (capture/restoreForceEndpointRowAligned): inside
//   a multirate PM cycle the particle layout is frozen, so rows are copied
//   directly by index. No ID lookup, no hash map, no allocation on restore.
//   Every precondition is validated before the first write; a refused
//   restore leaves the destination untouched.

#include <cstddef>
#include <cstdint>
#include <span>

#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/workflows/gravity_runtime.hpp"

namespace cosmosim::workflows::internal {

// Runtime-owned destination lanes of the dense particle force endpoint.
// total_* is the particle force cache; tree_*, pm_* and nearest_pair_r2 are
// the hierarchical split endpoint (may be empty for a non-hierarchical
// restart import).
struct ForceEndpointLanes {
  std::span<double> total_x;
  std::span<double> total_y;
  std::span<double> total_z;
  std::span<double> tree_x;
  std::span<double> tree_y;
  std::span<double> tree_z;
  std::span<double> pm_x;
  std::span<double> pm_y;
  std::span<double> pm_z;
  std::span<double> nearest_pair_r2;
};

struct ConstForceEndpointLanes {
  std::span<const double> total_x;
  std::span<const double> total_y;
  std::span<const double> total_z;
  std::span<const double> tree_x;
  std::span<const double> tree_y;
  std::span<const double> tree_z;
  std::span<const double> pm_x;
  std::span<const double> pm_y;
  std::span<const double> pm_z;
  std::span<const double> nearest_pair_r2;
};

// Read-only view of mutable lanes (extent checks before a write).
[[nodiscard]] inline ConstForceEndpointLanes asConst(const ForceEndpointLanes& lanes) noexcept {
  return ConstForceEndpointLanes{
      .total_x = lanes.total_x, .total_y = lanes.total_y, .total_z = lanes.total_z,
      .tree_x = lanes.tree_x, .tree_y = lanes.tree_y, .tree_z = lanes.tree_z,
      .pm_x = lanes.pm_x, .pm_y = lanes.pm_y, .pm_z = lanes.pm_z,
      .nearest_pair_r2 = lanes.nearest_pair_r2};
}

// Frozen-layout stamp shared by capture and restore.
struct RowAlignedEndpointStamp {
  std::uint64_t particle_count = 0U;
  std::uint64_t particle_index_generation = 0U;
  std::uint64_t layout_certificate = 0U;
};

// Copies the complete endpoint into `endpoint` (lane capacity is reused; it
// allocates only when a lane's capacity is below the row count, which the
// governed caller prevents by reserving first). All source lanes must have
// extent stamp.particle_count. On any exception `endpoint` is left invalid.
void captureForceEndpointRowAligned(
    const ConstForceEndpointLanes& source, const RowAlignedEndpointStamp& stamp,
    std::uint64_t pm_field_version, TransactionalForceEndpoint& endpoint);

// Throws std::logic_error naming the first violated invariant (endpoint not
// valid, particle count, particle-index generation, layout certificate, lane
// extents of the endpoint or of the destination). O(1): no per-row scan (the
// endpoint is the runtime's own validated synchronized endpoint). Writes
// nothing.
void validateRowAlignedForceEndpoint(
    const TransactionalForceEndpoint& endpoint, const RowAlignedEndpointStamp& expected,
    const ConstForceEndpointLanes& destination);

// validateRowAlignedForceEndpoint, then a row-by-row copy of all ten lanes.
// Strong guarantee: either every lane is restored or nothing is written.
void restoreForceEndpointRowAligned(
    const TransactionalForceEndpoint& endpoint, const RowAlignedEndpointStamp& expected,
    const ForceEndpointLanes& destination);

// Ordinary restart path. Matches cache rows to `state_particle_ids` by stable
// ID and writes the total lanes and, when `include_hierarchical_endpoint`,
// the split endpoint lanes. Destination lanes must already have extent
// state_particle_ids.size(). Rejects extent mismatches, non-finite totals,
// zero/duplicate cache IDs and IDs missing from the cache (std::runtime_error,
// messages unchanged from the historical importer). Validation (extents,
// finiteness, ID map) completes before the first write.
void remapForceEndpointById(
    const io::GravityForceCachePersistentState& cache,
    std::span<const std::uint64_t> state_particle_ids,
    const ForceEndpointLanes& destination, bool include_hierarchical_endpoint);

}  // namespace cosmosim::workflows::internal
