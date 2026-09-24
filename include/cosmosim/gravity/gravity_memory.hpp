#pragma once

#include <cstddef>
#include <cstdint>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"

namespace cosmosim::gravity {

// Single auditable arithmetic source for the TreePM short-range exchange
// workspace. Preflight uses the complete-graph peer degree d = rank_count - 1
// so the MemoryGovernor covers the worst-case simultaneous live set; runtime
// telemetry may report a smaller observed degree. Host packet capacity terms
// use the shared wire widths; that host/wire identity is compile-time
// enforced by static_assert on the packet types in tree_pm_coupling.cpp, so
// a padding change on a supported ABI fails the build instead of silently
// substituting wire bytes for host object bytes. This preflight answers a
// conservative "might require" question and need not equal the runtime
// retained-capacity high-water.
struct TreePmExchangeMemoryInput {
  std::uint32_t rank_count = 1U;
  std::uint64_t tree_exchange_batch_bytes = 4ULL * 1024ULL * 1024ULL;
  // Optional runtime peer degree override for telemetry-only re-estimation.
  // Zero means "use the preflight complete-graph degree".
  std::uint32_t runtime_peer_degree = 0U;
};

struct TreePmExchangeMemoryEstimate {
  std::uint64_t peer_degree = 0U;
  std::uint64_t batch_targets_per_peer = 0U;
  std::uint64_t wire_request_bytes_per_peer = 0U;
  std::uint64_t wire_response_bytes_per_peer = 0U;
  std::uint64_t wire_send_bytes = 0U;
  std::uint64_t wire_recv_bytes = 0U;
  std::uint64_t wire_response_send_bytes = 0U;
  std::uint64_t wire_response_recv_bytes = 0U;
  std::uint64_t wire_buffer_total_bytes = 0U;
  std::uint64_t structured_request_capacity_bytes = 0U;
  std::uint64_t response_mask_bytes = 0U;
  std::uint64_t response_count_bytes = 0U;
  std::uint64_t remote_accumulator_bytes = 0U;
  std::uint64_t rank_metadata_bytes = 0U;
  std::uint64_t transient_codec_bytes = 0U;
  std::uint64_t known_workspace_peak_bytes = 0U;
};

struct GravityMemoryEstimateInput {
  std::uint64_t local_source_count = 0U;
  std::uint64_t local_target_count = 0U;
  std::uint64_t local_particle_count = 0U;
  std::uint64_t local_cell_count = 0U;
  std::size_t tree_leaf_size = 16U;
  TreeMultipoleOrder multipole_order = TreeMultipoleOrder::kQuadrupole;
  PmGridShape pm_shape{};
  PmAssignmentScheme assignment_scheme = PmAssignmentScheme::kCic;
  core::PmDecompositionMode decomposition_mode = core::PmDecompositionMode::kSlab;
  std::uint32_t mpi_rank_count = 1U;
  int mpi_world_rank = 0;
  bool zoom_enabled = false;
  PmGridShape zoom_pm_shape{};
  bool periodic_tree_coordinates = true;
  bool indexed_target_coordinates = true;
  bool cuda_resident = false;
  std::uint64_t tree_exchange_batch_bytes = 4ULL * 1024ULL * 1024ULL;
  std::uint64_t pm_exchange_batch_bytes = 16ULL * 1024ULL * 1024ULL;
  std::uint64_t backend_unknown_reserve_bytes = 0U;
  double safety_margin_fraction = 0.0;
};

struct GravityMemoryEstimate {
  core::MemoryReport report;
  std::uint64_t known_peak_bytes = 0U;
  std::uint64_t external_backend_unknown_bytes = 0U;
  std::uint64_t budget_required_bytes = 0U;
  std::uint64_t estimated_tree_nodes = 0U;
  std::uint64_t pm_plan_owned_bytes = 0U;
};

struct DmoProcessMemoryPolicy {
  std::uint32_t mpi_rank_count = 1U;
  // scheduler_owned_bytes remains the compatibility/budget field. The
  // adjacent counters make live/capacity/high-water reporting truthful; a zero
  // high-water falls back to scheduler_owned_bytes for older callers.
  std::uint64_t scheduler_current_size_bytes = 0U;
  std::uint64_t scheduler_owned_bytes = 0U;
  std::uint64_t scheduler_high_water_bytes = 0U;
  std::uint64_t output_restart_overlap_bytes = 0U;
  std::uint64_t mpi_external_reserve_bytes = 0U;
  std::uint64_t fftw_external_reserve_bytes = 0U;
  std::uint64_t hdf5_external_reserve_bytes = 0U;
  std::uint64_t allocator_external_reserve_bytes = 0U;
  double safety_margin_fraction = 0.0;
};

struct DmoProcessMemoryEstimate {
  core::MemoryReport report;
  std::uint64_t known_owned_peak_bytes = 0U;
  std::uint64_t external_unknown_reserve_bytes = 0U;
  std::uint64_t modeled_subtotal_bytes = 0U;
  std::uint64_t safety_margin_bytes = 0U;
  std::uint64_t budget_required_bytes = 0U;
  std::uint64_t aggregate_required_bytes = 0U;
};

[[nodiscard]] GravityMemoryEstimate estimateGravityMemory(
    const GravityMemoryEstimateInput& input);

// Pure arithmetic; throws on overflow or invalid rank/batch inputs.
[[nodiscard]] TreePmExchangeMemoryEstimate estimateTreePmExchangeMemory(
    const TreePmExchangeMemoryInput& input);

[[nodiscard]] DmoProcessMemoryEstimate estimateDmoProcessMemory(
    const core::MemoryReport& canonical_runtime_report,
    const GravityMemoryEstimate& gravity_estimate,
    const DmoProcessMemoryPolicy& policy);

void enforceGravityMemoryBudget(
    const GravityMemoryEstimate& estimate,
    std::uint64_t budget_bytes);

void enforceDmoProcessMemoryBudget(
    const DmoProcessMemoryEstimate& estimate,
    std::uint64_t budget_bytes);

}  // namespace cosmosim::gravity
