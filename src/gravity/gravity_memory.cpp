#include "cosmosim/gravity/gravity_memory.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/core/openmp_runtime.hpp"
#include "internal/tree_pm_transport_planner.hpp"
#include "cosmosim/parallel/distributed_mesh.hpp"

namespace cosmosim::gravity {
namespace {

[[nodiscard]] std::uint64_t checkedAdd(std::uint64_t lhs, std::uint64_t rhs, const char* context) {
  if (lhs > std::numeric_limits<std::uint64_t>::max() - rhs) {
    throw std::overflow_error(context);
  }
  return lhs + rhs;
}

[[nodiscard]] std::uint64_t checkedMul(std::uint64_t lhs, std::uint64_t rhs, const char* context) {
  if (lhs != 0U && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
    throw std::overflow_error(context);
  }
  return lhs * rhs;
}

[[nodiscard]] std::uint64_t gridCells(const PmGridShape& shape, const char* context) {
  if (!shape.isValid()) {
    return 0U;
  }
  const std::uint64_t nx = static_cast<std::uint64_t>(shape.nx);
  const std::uint64_t ny = static_cast<std::uint64_t>(shape.ny);
  const std::uint64_t nz = static_cast<std::uint64_t>(shape.nz);
  return checkedMul(checkedMul(nx, ny, context), nz, context);
}

void addEstimate(
    core::MemoryReportBuilder& builder,
    core::MemorySubsystem subsystem,
    core::MemoryLifetime lifetime,
    std::string label,
    std::uint64_t bytes,
    std::string note = {}) {
  builder.addEntry(core::MemoryEntry{
      .subsystem = subsystem,
      .lifetime = lifetime,
      .label = std::move(label),
      .current_size_bytes = 0U,
      .owned_capacity_bytes = bytes,
      .high_water_bytes = bytes,
      .estimated_next_step_bytes = bytes,
      .estimate_only = true,
      .uncertainty_note = std::move(note),
  });
}

[[nodiscard]] std::uint64_t entryBudgetBytes(const core::MemoryEntry& entry) {
  return entry.estimated_next_step_bytes != 0U
      ? entry.estimated_next_step_bytes
      : entry.owned_capacity_bytes;
}

void copyKnownEntries(
    core::MemoryReportBuilder& builder,
    const core::MemoryReport& report) {
  for (const core::MemoryEntry& entry : report.entries) {
    if (entry.label == "category_present" ||
        entry.lifetime == core::MemoryLifetime::kUnknown) {
      continue;
    }
    builder.addEntry(entry);
  }
}

}  // namespace

TreePmExchangeMemoryEstimate estimateTreePmExchangeMemory(
    const TreePmExchangeMemoryInput& input) {
  if (input.rank_count == 0U) {
    throw std::invalid_argument("TreePM exchange memory estimate requires a non-zero rank count");
  }
  TreePmExchangeMemoryEstimate estimate;
  if (input.rank_count == 1U) {
    return estimate;
  }

  const std::uint64_t remote_rank_count =
      static_cast<std::uint64_t>(input.rank_count) - 1U;
  const std::uint64_t peer_degree = input.runtime_peer_degree != 0U
      ? static_cast<std::uint64_t>(input.runtime_peer_degree)
      : remote_rank_count;

  // Wire batch bound: each short-range request occupies a full wire record, so
  // b = floor(B / request_wire_bytes). Cap by the classic-MPI aggregate-round
  // planner for the complete-graph degree so preflight and runtime share one
  // clamp policy.
  if (input.tree_exchange_batch_bytes < kTreePmShortRangeRequestWireBytes) {
    throw std::invalid_argument(
        "TreePM exchange batch budget must fit at least one short-range request wire record");
  }
  std::uint64_t batch_targets = input.tree_exchange_batch_bytes /
      static_cast<std::uint64_t>(kTreePmShortRangeRequestWireBytes);
  if (batch_targets == 0U) {
    throw std::invalid_argument("TreePM exchange batch budget yields zero batch targets");
  }
  const internal::SparseTreePmRoundPlan round_plan = internal::planSparseTreePmRound(
      core::checkedIntegralNarrow<std::size_t>(
          batch_targets, "TreePM exchange batch target count"),
      core::checkedIntegralNarrow<std::size_t>(
          peer_degree, "TreePM exchange peer degree"),
      kTreePmShortRangeRequestWireBytes,
      kTreePmShortRangeResponseWireBytes);
  batch_targets = static_cast<std::uint64_t>(round_plan.targets_per_peer_per_round);

  estimate.peer_degree = peer_degree;
  estimate.batch_targets_per_peer = batch_targets;
  estimate.wire_request_bytes_per_peer = checkedMul(
      batch_targets, kTreePmShortRangeRequestWireBytes,
      "TreePM exchange request wire bytes per peer overflow");
  estimate.wire_response_bytes_per_peer = checkedMul(
      batch_targets, kTreePmShortRangeResponseWireBytes,
      "TreePM exchange response wire bytes per peer overflow");
  estimate.wire_send_bytes = checkedMul(
      peer_degree, estimate.wire_request_bytes_per_peer,
      "TreePM exchange request send bytes overflow");
  estimate.wire_recv_bytes = estimate.wire_send_bytes;
  estimate.wire_response_send_bytes = checkedMul(
      peer_degree, estimate.wire_response_bytes_per_peer,
      "TreePM exchange response send bytes overflow");
  estimate.wire_response_recv_bytes = estimate.wire_response_send_bytes;
  estimate.wire_buffer_total_bytes = checkedAdd(
      checkedAdd(estimate.wire_send_bytes, estimate.wire_recv_bytes,
                 "TreePM exchange wire total overflow"),
      checkedAdd(estimate.wire_response_send_bytes,
                 estimate.wire_response_recv_bytes,
                 "TreePM exchange wire total overflow"),
      "TreePM exchange wire total overflow");

  // Structured host capacity: one reserved request packet vector per non-self
  // peer (self is never reserved), simultaneous with the four wire buffers.
  // Host sizeof(ShortRangeTargetRequestPacket) == wire width is enforced by
  // static_assert in tree_pm_coupling.cpp; wire bytes are used here as host
  // object bytes only because that proof is compile-time.
  estimate.structured_request_capacity_bytes = checkedMul(
      remote_rank_count,
      checkedMul(batch_targets, kTreePmShortRangeRequestWireBytes,
                 "TreePM exchange structured request capacity overflow"),
      "TreePM exchange structured request capacity overflow");
  // response_expected_by_peer + response_seen_by_peer: R lanes x batch bytes.
  estimate.response_mask_bytes = checkedMul(
      static_cast<std::uint64_t>(input.rank_count),
      checkedMul(batch_targets, 2U, "TreePM exchange response mask capacity overflow"),
      "TreePM exchange response mask capacity overflow");
  // expected_response_count + received_response_count: 2 x batch x uint32.
  estimate.response_count_bytes = checkedMul(
      batch_targets, 2U * sizeof(std::uint32_t),
      "TreePM exchange response count capacity overflow");
  // remote_batch accel x/y/z lanes.
  estimate.remote_accumulator_bytes = checkedMul(
      batch_targets, 3U * sizeof(double),
      "TreePM exchange remote accumulator capacity overflow");
  // Arena-backed rank metadata includes the eight communicator-wide
  // count/displacement vectors, eight sparse-neighbor count/displacement
  // vectors, requested-peer identities, and one communicated-peer byte lane.
  estimate.rank_metadata_bytes = checkedMul(
      static_cast<std::uint64_t>(input.rank_count),
      17U * sizeof(int) + sizeof(std::uint8_t),
      "TreePM exchange rank metadata capacity overflow");
  // One-peer codec/validation capacity retained inside the arena lease:
  // request encode + request decode + response encode + response construction
  // + response decode + sorted uint64 identity scratch. Duplicate detection is
  // exact but node-free, so there is no unbounded unordered_set allocator
  // overhead outside the modeled workspace. Host packet sizes equal wire
  // widths by static_assert in tree_pm_coupling.cpp.
  estimate.transient_codec_bytes = checkedMul(
      batch_targets,
      2U * kTreePmShortRangeRequestWireBytes +
          3U * kTreePmShortRangeResponseWireBytes + sizeof(std::uint64_t),
      "TreePM exchange transient codec capacity overflow");

  std::uint64_t known = estimate.wire_buffer_total_bytes;
  known = checkedAdd(known, estimate.structured_request_capacity_bytes,
                     "TreePM exchange known workspace peak overflow");
  known = checkedAdd(known, estimate.response_mask_bytes,
                     "TreePM exchange known workspace peak overflow");
  known = checkedAdd(known, estimate.response_count_bytes,
                     "TreePM exchange known workspace peak overflow");
  known = checkedAdd(known, estimate.remote_accumulator_bytes,
                     "TreePM exchange known workspace peak overflow");
  known = checkedAdd(known, estimate.rank_metadata_bytes,
                     "TreePM exchange known workspace peak overflow");
  known = checkedAdd(known, estimate.transient_codec_bytes,
                     "TreePM exchange known workspace peak overflow");
  estimate.known_workspace_peak_bytes = known;
  return estimate;
}

GravityCommunicationArenaMemoryEstimate estimateGravityCommunicationArenaMemory(
    const GravityCommunicationArenaMemoryInput& input) {
  if (!input.pm_shape.isValid() || !input.pm_layout.isValid()) {
    throw std::invalid_argument(
        "gravity communication arena estimate requires valid PM shape/layout");
  }
  if (input.pm_layout.global_nx != input.pm_shape.nx ||
      input.pm_layout.global_ny != input.pm_shape.ny ||
      input.pm_layout.global_nz != input.pm_shape.nz) {
    throw std::invalid_argument(
        "gravity communication arena PM shape/layout mismatch");
  }
  GravityCommunicationArenaMemoryEstimate estimate;
  if (input.pm_layout.world_size <= 1) {
    return estimate;
  }

  // Reuse the PM routing capacity model. Metadata widths come directly from
  // the current density/interpolation workspace layouts; the 64 KiB structural
  // headroom covers allocator alignment/bookkeeping inside the bounded PMR
  // resource without turning it into a second payload allowance.
  const std::uint64_t ranks =
      static_cast<std::uint64_t>(input.pm_layout.world_size);
  const std::uint64_t density_metadata_bytes = checkedMul(
      ranks, 8U * sizeof(int) + sizeof(std::size_t),
      "gravity PM density metadata estimate overflow");
  const std::uint64_t interpolation_metadata_bytes = checkedMul(
      ranks, 12U * sizeof(int) + sizeof(std::size_t),
      "gravity PM interpolation metadata estimate overflow");
  const PmRoutingCapacityModel density_model = modelPmRoutingCapacity(
      input.pm_layout.world_size,
      input.pm_exchange_batch_bytes,
      k_pm_routing_modeled_workspace_limit_bytes,
      density_metadata_bytes,
      k_pm_routing_max_wire_record_bytes);
  const PmRoutingCapacityModel interpolation_model = modelPmRoutingCapacity(
      input.pm_layout.world_size,
      input.pm_exchange_batch_bytes,
      k_pm_routing_modeled_workspace_limit_bytes,
      interpolation_metadata_bytes,
      k_pm_routing_max_wire_record_bytes);
  estimate.pm_density_required_bytes = checkedAdd(
      density_model.max_simultaneous_workspace_bytes,
      k_pm_routing_workspace_headroom_bytes,
      "gravity PM density arena estimate overflow");
  estimate.pm_interpolation_required_bytes = checkedAdd(
      interpolation_model.max_simultaneous_workspace_bytes,
      k_pm_routing_workspace_headroom_bytes,
      "gravity PM interpolation arena estimate overflow");

  // TreePM production receives directly into the persistent six-lane force
  // halo cache.  Shared scratch therefore needs only left/right send staging.
  if (input.pm_layout.local_nx() != 0U) {
    const std::uint64_t plane_values = checkedMul(
        static_cast<std::uint64_t>(input.pm_shape.ny),
        static_cast<std::uint64_t>(input.pm_shape.nz),
        "gravity PM halo plane estimate overflow");
    estimate.pm_halo_staging_required_bytes = checkedMul(
        plane_values, 2U * sizeof(double),
        "gravity PM halo staging estimate overflow");
  }

  const TreePmExchangeMemoryEstimate tree_exchange =
      estimateTreePmExchangeMemory(TreePmExchangeMemoryInput{
          .rank_count = static_cast<std::uint32_t>(input.pm_layout.world_size),
          .tree_exchange_batch_bytes = input.tree_exchange_batch_bytes,
      });
  estimate.tree_exchange_required_bytes = checkedAdd(
      tree_exchange.known_workspace_peak_bytes,
      k_pm_routing_workspace_headroom_bytes,
      "gravity Tree exchange arena headroom overflow");
  estimate.required_capacity_bytes = std::max(
      std::max(estimate.pm_density_required_bytes,
               estimate.pm_interpolation_required_bytes),
      std::max(estimate.pm_halo_staging_required_bytes,
               estimate.tree_exchange_required_bytes));
  if (estimate.required_capacity_bytes > estimate.certified_limit_bytes) {
    throw std::runtime_error(
        "TreePM gravity communication arena requires " +
        std::to_string(estimate.required_capacity_bytes) +
        " bytes, exceeding certified 256 MiB/rank limit " +
        std::to_string(estimate.certified_limit_bytes));
  }
  return estimate;
}

GravityMemoryEstimate estimateGravityMemory(const GravityMemoryEstimateInput& input) {
  if (input.tree_leaf_size == 0U || input.mpi_rank_count == 0U) {
    throw std::invalid_argument("gravity memory estimate requires non-zero leaf size and rank count");
  }
  if (input.mpi_world_rank < 0 ||
      input.mpi_world_rank >= static_cast<int>(input.mpi_rank_count)) {
    throw std::invalid_argument("gravity memory estimate MPI rank is outside configured rank count");
  }
  if (!std::isfinite(input.safety_margin_fraction) || input.safety_margin_fraction < 0.0 ||
      input.safety_margin_fraction > 1.0) {
    throw std::invalid_argument("gravity memory safety margin must be finite and within [0,1]");
  }
  const std::uint64_t leaf_size = static_cast<std::uint64_t>(input.tree_leaf_size);
  const std::uint64_t leaf_count = input.local_source_count == 0U
      ? 0U
      : checkedAdd(input.local_source_count, leaf_size - 1U, "gravity leaf estimate overflow") / leaf_size;
  const std::uint64_t estimated_tree_nodes = leaf_count == 0U
      ? 0U
      : checkedAdd(1U, checkedMul(2U, leaf_count, "gravity node estimate overflow"), "gravity node estimate overflow");

  // Hot node lanes: center/half-size, mass/COM, softening envelope = 10 doubles;
  // topology/range lanes = 12 uint32-equivalent bytes + child fanout = 32 bytes.
  // Quadrupole adds seven double lanes only when selected.
  const std::uint64_t hot_node_bytes = 10U * sizeof(double) +
      3U * sizeof(std::uint32_t) + sizeof(std::uint8_t) +
      8U * sizeof(std::uint32_t);
  const std::uint64_t cold_node_bytes = input.multipole_order == TreeMultipoleOrder::kQuadrupole
      ? 7U * sizeof(double)
      : 0U;
  const std::uint64_t tree_nodes_bytes = checkedMul(
      estimated_tree_nodes,
      checkedAdd(hot_node_bytes, cold_node_bytes, "gravity node byte estimate overflow"),
      "gravity node byte estimate overflow");

  // Staging is intentionally limited to fields used by gravity. Targets are
  // source-index views, so there is no second target coordinate triplet here.
  // High-resolution classification is cold unless the zoom long-range
  // correction is active; do not charge its byte lanes to homogeneous runs.
  const bool borrowed_homogeneous_dmo =
      input.source_representation ==
      GravitySourceRepresentation::kBorrowedHomogeneousDmo;
  if (borrowed_homogeneous_dmo &&
      (input.local_cell_count != 0U || input.zoom_enabled)) {
    throw std::invalid_argument(
        "borrowed homogeneous DMO gravity estimate requires zero cells and disabled zoom correction");
  }
  const std::uint64_t zoom_mask_bytes_per_source =
      input.zoom_enabled ? sizeof(std::uint8_t) : 0U;
  const std::uint64_t zoom_mask_bytes_per_target =
      input.zoom_enabled ? sizeof(std::uint8_t) : 0U;
  const std::uint64_t source_staging_bytes = borrowed_homogeneous_dmo
      ? 0U
      : checkedMul(
            input.local_source_count,
            5U * sizeof(double) + 3U * sizeof(std::uint32_t) +
                sizeof(std::uint8_t) + zoom_mask_bytes_per_source,
            "gravity source staging estimate overflow");
  const std::uint64_t borrowed_target_bytes_per_target = checkedAdd(
      sizeof(std::uint32_t),
      input.relative_force_mac_enabled ? sizeof(double) : 0U,
      "gravity borrowed target view byte estimate overflow");
  const std::uint64_t target_view_bytes = borrowed_homogeneous_dmo
      ? checkedMul(
            input.local_target_count, borrowed_target_bytes_per_target,
            "gravity borrowed target-view estimate overflow")
      : checkedMul(
            input.local_target_count,
            5U * sizeof(std::uint32_t) + 2U * sizeof(double) +
                sizeof(std::uint8_t) + zoom_mask_bytes_per_target,
            "gravity target view estimate overflow");
  const bool uniform_source_softening =
      borrowed_homogeneous_dmo || input.source_softening_uniform;
  const std::uint64_t tree_construction_bytes = checkedMul(
      input.local_source_count,
      2U * sizeof(std::uint64_t) + 2U * sizeof(TreeLocalIndex) +
          (uniform_source_softening ? 0U : sizeof(double)),
      "gravity tree construction estimate overflow");
  const std::uint64_t acceleration_bytes = checkedMul(
      input.local_target_count,
      3U * sizeof(double),
      "gravity acceleration estimate overflow");
  const std::uint64_t periodic_tree_coordinate_bytes = input.periodic_tree_coordinates
      ? checkedMul(
          input.local_source_count,
          3U * sizeof(double),
          "gravity periodic tree staging estimate overflow")
      : 0U;
  const std::uint64_t zoom_active_correction_bytes = input.zoom_enabled
      ? checkedMul(
          input.local_target_count,
          3U * sizeof(double),
          "gravity zoom active correction estimate overflow")
      : 0U;
  const std::uint64_t persistent_force_cache_bytes = checkedMul(
      input.local_source_count,
      3U * sizeof(double) + sizeof(std::uint8_t) +
          (input.hierarchical_kdk_enabled ? 6U * sizeof(double) + sizeof(std::uint64_t) : 0U),
      "gravity persistent force cache estimate overflow");
  const std::uint64_t runtime_particle_map_bytes = borrowed_homogeneous_dmo
      ? 0U
      : checkedMul(
            input.local_particle_count,
            3U * sizeof(std::int32_t),
            "gravity runtime particle map estimate overflow");
  const std::uint64_t runtime_cell_map_bytes = borrowed_homogeneous_dmo
      ? 0U
      : checkedMul(
            input.local_cell_count,
            2U * sizeof(std::int32_t) + sizeof(std::uint8_t),
            "gravity runtime cell map estimate overflow");
  const std::uint64_t runtime_refresh_list_bytes = borrowed_homogeneous_dmo
      ? 0U
      : checkedMul(
            input.local_particle_count, sizeof(std::uint32_t),
            "gravity runtime refresh-list estimate overflow");
  const std::uint64_t runtime_mapping_bytes = checkedAdd(
      checkedAdd(runtime_particle_map_bytes, runtime_cell_map_bytes,
                 "gravity runtime mapping estimate overflow"),
      runtime_refresh_list_bytes, "gravity runtime mapping estimate overflow");

  (void)gridCells(input.pm_shape, "gravity PM grid estimate overflow");
  const parallel::PmSlabLayout pm_layout = parallel::makePmSlabLayout(
      input.pm_shape.nx,
      input.pm_shape.ny,
      input.pm_shape.nz,
      static_cast<int>(input.mpi_rank_count),
      input.mpi_world_rank);
  const std::uint64_t local_pm_cells =
      static_cast<std::uint64_t>(pm_layout.localCellCount());
  // Production periodic TreePM can place density directly in the FFT real
  // owner, leaving PmGridStorage with only the three force components. Generic
  // compatibility/isolated paths retain compact density explicitly.
  const std::uint64_t pm_grid_components =
      input.periodic_fft_backed_density ? 3U : 4U;
  const std::uint64_t pm_owned_bytes = checkedMul(
      local_pm_cells, pm_grid_components * sizeof(double),
      "gravity PM owned estimate overflow");
  const PmPlanResourcesMemoryEstimate pm_plan_memory =
      estimatePmPlanResourcesMemory(input.pm_shape, pm_layout, input.decomposition_mode);
  const std::uint64_t zoom_cells = input.zoom_enabled
      ? gridCells(input.zoom_pm_shape, "gravity zoom PM estimate overflow")
      : 0U;
  const std::uint64_t zoom_local_cells = zoom_cells == 0U
      ? 0U
      : checkedAdd(zoom_cells, static_cast<std::uint64_t>(input.mpi_rank_count) - 1U,
                   "gravity zoom local cell estimate overflow") /
          static_cast<std::uint64_t>(input.mpi_rank_count);
  const std::uint64_t zoom_bytes = checkedMul(
      zoom_local_cells, 5U * sizeof(double), "gravity zoom PM owned estimate overflow");

  const GravityCommunicationArenaMemoryEstimate communication_arena =
      estimateGravityCommunicationArenaMemory(GravityCommunicationArenaMemoryInput{
          .pm_shape = input.pm_shape,
          .pm_layout = pm_layout,
          .pm_exchange_batch_bytes = input.pm_exchange_batch_bytes,
          .tree_exchange_batch_bytes = input.tree_exchange_batch_bytes,
      });
  const std::uint64_t communication_arena_bytes =
      communication_arena.required_capacity_bytes;
  const std::uint64_t pm_force_halo_cache_bytes = input.mpi_rank_count > 1U
      ? checkedMul(
            checkedMul(
                static_cast<std::uint64_t>(input.pm_shape.ny),
                static_cast<std::uint64_t>(input.pm_shape.nz),
                "gravity PM force halo cache plane overflow"),
            6U * sizeof(double),
            "gravity PM force halo cache estimate overflow")
      : 0U;
  const core::OpenMpRuntimeInfo openmp_info = core::openMpRuntimeInfo();
  const std::uint64_t residual_worker_count =
      openmp_info.compiled
          ? static_cast<std::uint64_t>(std::max(
                1,
                std::max(openmp_info.configured_threads,
                         openmp_info.maximum_threads)))
          : 1U;
  const std::uint64_t residual_depth_bound = kMaximumTreeDepth;
  const std::uint64_t residual_stack_slots = checkedAdd(
      1U, checkedMul(7U, residual_depth_bound, "gravity residual stack depth overflow"),
      "gravity residual stack slot overflow");
  const std::uint64_t treepm_worker_stack_bytes = checkedMul(
      checkedMul(residual_worker_count, residual_stack_slots,
                 "gravity residual worker stack overflow"),
      static_cast<std::uint64_t>(sizeof(TreeLocalIndex)),
      "gravity residual worker stack byte overflow");
  const std::uint64_t treepm_worker_counter_bytes = checkedMul(
      residual_worker_count,
      checkedAdd(kTreePmResidualCounterBytes, kTreePmResidualSpatialCounterBytes,
                 "gravity residual worker counter bundle overflow"),
      "gravity residual worker counter overflow");
  const std::uint64_t residual_block_count = checkedAdd(
      input.local_target_count, kTreePmResidualBlockSize - 1U,
      "gravity residual diagnostic block count overflow") /
      kTreePmResidualBlockSize;
  const std::uint64_t treepm_block_diagnostic_bytes = checkedMul(
      residual_block_count, sizeof(double),
      "gravity residual block diagnostic overflow");
  const std::uint64_t treepm_worker_scratch_bytes = checkedAdd(
      treepm_worker_stack_bytes, treepm_worker_counter_bytes,
      "gravity residual worker scratch overflow");
  const std::uint64_t cuda_owned_workspace = input.cuda_resident
      ? checkedAdd(
            checkedMul(input.local_source_count, 4U * sizeof(double),
                       "gravity CUDA source workspace estimate overflow"),
            checkedAdd(
                checkedMul(input.local_target_count, 3U * sizeof(double),
                           "gravity CUDA target workspace estimate overflow"),
                checkedMul(local_pm_cells, 4U * sizeof(double),
                           "gravity CUDA mesh workspace estimate overflow"),
                "gravity CUDA workspace estimate overflow"),
            "gravity CUDA workspace estimate overflow")
      : 0U;
  // The solver-owned PlanResources vectors are deterministic and accounted as
  // known memory. Only genuinely backend-owned plan/runtime allocations remain
  // in this legacy gravity reserve.
  const std::uint64_t backend_unknown = input.backend_unknown_reserve_bytes;

  core::MemoryReportBuilder builder;
  addEstimate(builder, core::MemorySubsystem::kActiveSets, core::MemoryLifetime::kTransient,
              "gravity.estimate.source_staging", source_staging_bytes,
              borrowed_homogeneous_dmo
                  ? "borrowed homogeneous DMO uses canonical SimulationState XYZ/mass and implicit species/row identity"
                  : "authoritative source staging; excludes canonical SimulationState");
  addEstimate(builder, core::MemorySubsystem::kActiveSets, core::MemoryLifetime::kTransient,
              "gravity.estimate.target_index_views", target_view_bytes,
              borrowed_homogeneous_dmo
                  ? (input.relative_force_mac_enabled
                         ? "one governed uint32 all-particle source/target index lane plus lazy relative-MAC previous-acceleration magnitude"
                         : "one governed uint32 all-particle source/target index lane")
                  : "targets alias source coordinates by compact local index");
  addEstimate(builder, core::MemorySubsystem::kActiveSets, core::MemoryLifetime::kTransient,
              "gravity.estimate.force_accumulators", acceleration_bytes,
              "authoritative compact active acceleration triplet");
  if (periodic_tree_coordinate_bytes > 0U) {
    addEstimate(builder, core::MemorySubsystem::kTree, core::MemoryLifetime::kTransient,
                "gravity.estimate.periodic_tree_coordinate_staging", periodic_tree_coordinate_bytes,
                "three final unwrapped tree coordinates; periodic axis sorting borrows the shared tree construction key lanes");
  }
  if (zoom_active_correction_bytes > 0U) {
    addEstimate(builder, core::MemorySubsystem::kActiveSets, core::MemoryLifetime::kTransient,
                "gravity.estimate.zoom_active_correction", zoom_active_correction_bytes,
                "allocated only while zoom long-range correction is enabled");
  }
  addEstimate(builder, core::MemorySubsystem::kSidecars, core::MemoryLifetime::kPersistent,
              "gravity.estimate.persistent_force_cache", persistent_force_cache_bytes,
              "three acceleration lanes plus validity; exact particle/cell split is runtime-owned");
  addEstimate(builder, core::MemorySubsystem::kActiveSets, core::MemoryLifetime::kTransient,
              "gravity.estimate.runtime_index_and_selection_maps", runtime_mapping_bytes,
              borrowed_homogeneous_dmo
                  ? "identity particle/source/target/active-slot mappings are implicit; no cell maps or refresh list"
                  : "active-slot/owned-local maps, leaf mask and force-refresh particle list");
  addEstimate(builder, core::MemorySubsystem::kTree, core::MemoryLifetime::kTransient,
              "gravity.estimate.tree_nodes", tree_nodes_bytes,
              "leaf-derived estimate; dynamic growth remains possible for adversarial geometry");
  addEstimate(builder, core::MemorySubsystem::kScratch, core::MemoryLifetime::kTransient,
              "gravity.estimate.tree_construction_ownership", tree_construction_bytes,
              uniform_source_softening
                  ? "final TreeLocalIndex permutation + scalar source epsilon + shared uint64 key primary/key scratch + TreeLocalIndex radix/partition scratch"
                  : "final TreeLocalIndex permutation + materialized double source epsilon + shared uint64 key primary/key scratch + TreeLocalIndex radix/partition scratch");
  addEstimate(builder, core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kTransient,
              "gravity.estimate.pm_owned_fields", pm_owned_bytes,
              input.periodic_fft_backed_density
                  ? "periodic TreePM force-only grid: three force fields; density is owned once by the FFT real plan allocation"
                  : "generic PM grid: compact density plus three force fields; real potential is demand-driven");
  addEstimate(builder, core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kPersistent,
              "gravity.estimate.pm_plan_resources_owned_arrays",
              pm_plan_memory.total_owned_bytes,
              pm_plan_memory.used_backend_allocation_query
                  ? "CHUI-owned FFT real/fourier/potential_k arrays plus O(nx+ny+nz) spectral axis metadata sized from the active FFTW MPI allocation query; backend plan internals excluded"
                  : "CHUI-owned FFT real/fourier/potential_k arrays plus O(nx+ny+nz) spectral axis metadata sized from conservative PM decomposition geometry; backend plan internals excluded");
  if (zoom_bytes > 0U) {
    addEstimate(builder, core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kTransient,
                "gravity.estimate.zoom_pm_owned_fields", zoom_bytes,
                "coarse/focused lifetimes are serialized; focused isolated PM retains real potential and defines the modeled correction-grid peak");
  }
  if (pm_force_halo_cache_bytes > 0U) {
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kPmMesh,
        .lifetime = core::MemoryLifetime::kPersistent,
        .memory_class = core::MemoryClass::kPersistentCache,
        .label = "gravity.estimate.pm_force_halo_cache",
        .current_size_bytes = 0U,
        .owned_capacity_bytes = pm_force_halo_cache_bytes,
        .high_water_bytes = pm_force_halo_cache_bytes,
        .estimated_next_step_bytes = pm_force_halo_cache_bytes,
        .estimate_only = true,
        .uncertainty_note =
            "six retained one-plane force halo lanes survive staging reset through PM interpolation",
    });
  }
  if (communication_arena_bytes > 0U) {
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kMpiBuffers,
        .lifetime = core::MemoryLifetime::kTransient,
        .memory_class = core::MemoryClass::kCommunication,
        .label = "gravity.estimate.communication_arena",
        .current_size_bytes = 0U,
        .owned_capacity_bytes = communication_arena_bytes,
        .high_water_bytes = communication_arena_bytes,
        .estimated_next_step_bytes = communication_arena_bytes,
        .estimate_only = true,
        .uncertainty_note =
            "one physical owner = max(PM density routing, PM halo staging, PM interpolation routing, short-range Tree exchange); sequential logical contributors are not summed",
    });
  }
  const auto add_logical_communication_estimate = [&](
      std::string label, std::uint64_t bytes) {
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kMpiBuffers,
        .lifetime = core::MemoryLifetime::kTransient,
        .memory_class = core::MemoryClass::kCommunication,
        .label = std::move(label),
        .current_size_bytes = 0U,
        .owned_capacity_bytes = 0U,
        .high_water_bytes = bytes,
        .estimated_next_step_bytes = 0U,
        .estimate_only = true,
        .uncertainty_note =
            "logical contributor to gravity.estimate.communication_arena; non-owning and excluded from physical peak sum",
    });
  };
  if (communication_arena_bytes > 0U) {
    add_logical_communication_estimate(
        "gravity.estimate.communication_arena.pm_density_logical",
        communication_arena.pm_density_required_bytes);
    add_logical_communication_estimate(
        "gravity.estimate.communication_arena.pm_halo_logical",
        communication_arena.pm_halo_staging_required_bytes);
    add_logical_communication_estimate(
        "gravity.estimate.communication_arena.pm_interpolation_logical",
        communication_arena.pm_interpolation_required_bytes);
    add_logical_communication_estimate(
        "gravity.estimate.communication_arena.tree_exchange_logical",
        communication_arena.tree_exchange_required_bytes);
  }
   addEstimate(builder, core::MemorySubsystem::kScratch, core::MemoryLifetime::kTransient,
               "gravity.estimate.treepm_residual_worker_scratch", treepm_worker_scratch_bytes,
               "OpenMP residual worker storage: enforced kMaximumTreeDepth stacks, aligned hot counters and a conservative bound for optional spatial scratch per planned worker");
   if (treepm_block_diagnostic_bytes > 0U) {
     addEstimate(builder, core::MemorySubsystem::kScratch, core::MemoryLifetime::kTransient,
                 "gravity.estimate.treepm_residual_block_diagnostics", treepm_block_diagnostic_bytes,
                 "deterministic block floating diagnostic storage; force accumulation order remains per target");
   }

  const std::uint64_t bounded_feedback_bytes = checkedAdd(
      checkedMul(2U, sizeof(TreePmDiagnostics), "TreePM diagnostic object bytes"),
      checkedMul(3U, sizeof(std::array<double, parallel::k_spatial_work_bin_count>),
          "TreePM spatial history bytes"), "TreePM fixed feedback bytes");
  addEstimate(builder, core::MemorySubsystem::kTree, core::MemoryLifetime::kPersistent,
      "gravity.estimate.bounded_force_feedback", bounded_feedback_bytes,
      "coordinator and workflow current-force summaries, two decayed bins and one planner view; fixed size independent of particle and solve counts");

  if (cuda_owned_workspace > 0U) {
    addEstimate(builder, core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kPersistent,
                "gravity.estimate.cuda_owned_persistent_workspace", cuda_owned_workspace,
                "known PmSolver-owned CUDA source/mesh/acceleration buffers; cuFFT/runtime internals excluded");
  }
  addEstimate(builder, core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kUnknown,
              "gravity.estimate.external_fft_or_device_workspace", backend_unknown,
              "backend allocator/plan workspace is not owned by gravity containers");

  GravityMemoryEstimate result;
  result.report = std::move(builder).finish();
  result.report.notes.push_back(
      borrowed_homogeneous_dmo
          ? "Gravity pre-run estimate selected borrowed_homogeneous_dmo: canonical XYZ/mass are aliased, homogeneous DM species and row mappings are implicit, and one governed uint32 target/source index lane remains; canonical SimulationState is reported separately."
          : "Gravity pre-run estimate selected materialized_generic: owned source staging, compact target/force lanes, runtime index/selection maps, PM indexed-target scratch, periodic tree staging, tree workspace, PM fields, explicit PM PlanResources arrays, optional zoom lanes, persistent force cache, persistent PM force-halo cache, one shared bounded gravity communication arena, and known CUDA buffers are modeled; canonical SimulationState is reported separately.");
  result.report.notes.push_back(
      std::string("PM estimate profile assignment=") +
      (input.assignment_scheme == PmAssignmentScheme::kTsc ? "tsc" : "cic") +
      " decomposition=" +
      (input.decomposition_mode == core::PmDecompositionMode::kPencil
           ? "transposed_slab_contract"
           : "x_slab"));
  result.report.notes.push_back(
      "CHUI-owned PM PlanResources arrays are explicit known memory. FFTW/cuFFT plan internals and other library-owned allocations remain uncertain and must be carried by configured external reserves.");
  result.external_backend_unknown_bytes = backend_unknown;
  result.estimated_tree_nodes = estimated_tree_nodes;
  result.pm_plan_owned_bytes = pm_plan_memory.total_owned_bytes;
  result.communication_arena_bytes = communication_arena_bytes;
  std::uint64_t known_peak = 0U;
  for (const core::MemoryEntry& entry : result.report.entries) {
    if (entry.lifetime != core::MemoryLifetime::kUnknown) {
      known_peak = checkedAdd(known_peak, entry.estimated_next_step_bytes, "gravity known peak estimate overflow");
    }
  }
  result.known_peak_bytes = known_peak;
  const std::uint64_t local_owned_estimate = checkedAdd(
      checkedAdd(result.report.totals.persistent_total_bytes, result.report.totals.transient_total_bytes,
                 "gravity distributed local memory summary overflow"),
      result.report.totals.unknown_total_bytes,
      "gravity distributed local memory summary overflow");
  result.report.distributed.valid = true;
  result.report.distributed.rank_count = static_cast<int>(input.mpi_rank_count);
  result.report.distributed.local_owned_bytes = local_owned_estimate;
  result.report.distributed.rank_max_owned_bytes = local_owned_estimate;
  result.report.distributed.global_sum_owned_bytes = checkedMul(
      local_owned_estimate, static_cast<std::uint64_t>(input.mpi_rank_count),
      "gravity distributed aggregate memory summary overflow");
  result.report.distributed.max_to_mean_imbalance_ratio = 1.0;
  const std::uint64_t base_budget_requirement = checkedAdd(
      known_peak, backend_unknown, "gravity budget requirement overflow");
  const long double scaled_requirement = static_cast<long double>(base_budget_requirement) *
      (1.0L + static_cast<long double>(input.safety_margin_fraction));
  if (scaled_requirement > static_cast<long double>(std::numeric_limits<std::uint64_t>::max())) {
    throw std::overflow_error("gravity budget safety-margin estimate overflow");
  }
  result.budget_required_bytes = static_cast<std::uint64_t>(std::ceil(scaled_requirement));
  result.report.notes.push_back(
      "gravity_budget_required_bytes=" + std::to_string(result.budget_required_bytes) +
      " (known + backend reserve, then configured safety margin)");
  return result;
}

DmoProcessMemoryEstimate estimateDmoProcessMemory(
    const core::MemoryReport& canonical_runtime_report,
    const GravityMemoryEstimate& gravity_estimate,
    const DmoProcessMemoryPolicy& policy) {
  if (policy.mpi_rank_count == 0U) {
    throw std::invalid_argument("DMO process memory estimate requires a non-zero MPI rank count");
  }
  if (!std::isfinite(policy.safety_margin_fraction) ||
      policy.safety_margin_fraction < 0.0 ||
      policy.safety_margin_fraction > 1.0) {
    throw std::invalid_argument(
        "DMO process memory safety margin must be finite and within [0,1]");
  }

  core::MemoryReportBuilder builder;
  copyKnownEntries(builder, canonical_runtime_report);
  copyKnownEntries(builder, gravity_estimate.report);

  const std::uint64_t scheduler_high_water_bytes =
      policy.scheduler_high_water_bytes != 0U
      ? policy.scheduler_high_water_bytes
      : policy.scheduler_owned_bytes;
  if (policy.scheduler_current_size_bytes > policy.scheduler_owned_bytes ||
      policy.scheduler_owned_bytes > scheduler_high_water_bytes) {
    throw std::invalid_argument(
        "DMO process scheduler memory accounting requires logical <= capacity <= high-water");
  }
  if (policy.scheduler_owned_bytes > 0U || scheduler_high_water_bytes > 0U) {
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kSidecars,
        .lifetime = core::MemoryLifetime::kPersistent,
        .label = "dmo_process.scheduler_owned_state",
        .current_size_bytes = policy.scheduler_current_size_bytes,
        .owned_capacity_bytes = policy.scheduler_owned_bytes,
        .high_water_bytes = scheduler_high_water_bytes,
        .estimated_next_step_bytes = policy.scheduler_owned_bytes,
        .estimate_only = false,
        .uncertainty_note =
            "authoritative hierarchical scheduler logical bytes, retained vector capacity, and historical retained-capacity high-water",
    });
  }
  if (policy.output_restart_overlap_bytes > 0U) {
    addEstimate(builder,
                core::MemorySubsystem::kOutputBuffers,
                core::MemoryLifetime::kTransient,
                "dmo_process.output_restart_overlap",
                policy.output_restart_overlap_bytes,
                "configured owned output/restart staging allowed to overlap resident DMO/PM state");
  }

  const auto add_external = [&builder](std::string label,
                                       std::uint64_t bytes,
                                       std::string note) {
    if (bytes == 0U) {
      return;
    }
    addEstimate(builder,
                core::MemorySubsystem::kScratch,
                core::MemoryLifetime::kUnknown,
                std::move(label),
                bytes,
                std::move(note));
  };
  add_external(
      "dmo_process.external_legacy_gravity_backend_reserve",
      gravity_estimate.external_backend_unknown_bytes,
      "legacy gravity backend reserve retained for backward-compatible parameter files; additive with process-specific reserves");
  add_external(
      "dmo_process.external_mpi_reserve",
      policy.mpi_external_reserve_bytes,
      "configured reserve for MPI implementation-owned eager/rendezvous/collective buffers");
  add_external(
      "dmo_process.external_fftw_reserve",
      policy.fftw_external_reserve_bytes,
      "configured reserve for FFTW plan/runtime allocations not represented by CHUI-owned PlanResources vectors");
  add_external(
      "dmo_process.external_hdf5_reserve",
      policy.hdf5_external_reserve_bytes,
      "configured reserve for HDF5 chunk/cache/library allocations");
  add_external(
      "dmo_process.external_allocator_reserve",
      policy.allocator_external_reserve_bytes,
      "configured reserve for allocator fragmentation and other process-owned overhead not attributable to a tracked container");

  DmoProcessMemoryEstimate result;
  result.report = std::move(builder).finish();

  std::uint64_t known = 0U;
  std::uint64_t unknown = 0U;
  for (const core::MemoryEntry& entry : result.report.entries) {
    const std::uint64_t bytes = entryBudgetBytes(entry);
    if (entry.lifetime == core::MemoryLifetime::kUnknown) {
      unknown = checkedAdd(unknown, bytes, "DMO process external reserve sum overflow");
    } else {
      known = checkedAdd(known, bytes, "DMO process known memory sum overflow");
    }
  }
  result.known_owned_peak_bytes = known;
  result.external_unknown_reserve_bytes = unknown;
  result.modeled_subtotal_bytes = checkedAdd(
      known, unknown, "DMO process modeled subtotal overflow");

  const long double scaled = static_cast<long double>(result.modeled_subtotal_bytes) *
      (1.0L + static_cast<long double>(policy.safety_margin_fraction));
  if (scaled > static_cast<long double>(std::numeric_limits<std::uint64_t>::max())) {
    throw std::overflow_error("DMO process safety-margin estimate overflow");
  }
  result.budget_required_bytes = static_cast<std::uint64_t>(std::ceil(scaled));
  result.safety_margin_bytes = result.budget_required_bytes - result.modeled_subtotal_bytes;
  result.aggregate_required_bytes = checkedMul(
      result.budget_required_bytes,
      static_cast<std::uint64_t>(policy.mpi_rank_count),
      "DMO process aggregate requirement overflow");

  result.report.distributed.valid = true;
  result.report.distributed.rank_count = static_cast<int>(policy.mpi_rank_count);
  result.report.distributed.local_owned_bytes = result.modeled_subtotal_bytes;
  result.report.distributed.rank_max_owned_bytes = result.modeled_subtotal_bytes;
  result.report.distributed.global_sum_owned_bytes = checkedMul(
      result.modeled_subtotal_bytes,
      static_cast<std::uint64_t>(policy.mpi_rank_count),
      "DMO process aggregate modeled subtotal overflow");
  result.report.distributed.max_to_mean_imbalance_ratio = 1.0;
  result.report.notes.push_back(
      "Authoritative DMO process preflight composes live canonical/runtime ownership with predicted gravity/PM peak entries. Unknown zero-byte placeholders from component reports are intentionally replaced by explicit configured reserves.");
  result.report.notes.push_back(
      "IC import staging is excluded from the gravity-phase peak because the importer lifetime ends before the production integrator owns the SimulationState. Output/restart overlap is included only through its explicit configured process reserve.");
  result.report.notes.push_back(
      "dmo_process_known_owned_peak_bytes=" + std::to_string(result.known_owned_peak_bytes));
  result.report.notes.push_back(
      "dmo_process_external_unknown_reserve_bytes=" +
      std::to_string(result.external_unknown_reserve_bytes));
  result.report.notes.push_back(
      "dmo_process_modeled_subtotal_bytes=" + std::to_string(result.modeled_subtotal_bytes));
  result.report.notes.push_back(
      "dmo_process_safety_margin_bytes=" + std::to_string(result.safety_margin_bytes));
  result.report.notes.push_back(
      "dmo_process_budget_required_bytes=" + std::to_string(result.budget_required_bytes));
  result.report.notes.push_back(
      "dmo_process_projected_equal_rank_aggregate_required_bytes=" +
      std::to_string(result.aggregate_required_bytes));
  return result;
}

void enforceGravityMemoryBudget(
    const GravityMemoryEstimate& estimate,
    std::uint64_t budget_bytes) {
  if (budget_bytes == 0U || estimate.budget_required_bytes <= budget_bytes) {
    return;
  }
  throw std::runtime_error(
      "estimated gravity budget requirement " + std::to_string(estimate.budget_required_bytes) +
      " bytes (known peak=" + std::to_string(estimate.known_peak_bytes) +
      ", backend unknown/reserve=" +
      std::to_string(estimate.external_backend_unknown_bytes) +
      ") exceeds configured per-rank gravity memory budget " +
      std::to_string(budget_bytes) + " bytes before tree/communication allocation");
}

void enforceDmoProcessMemoryBudget(
    const DmoProcessMemoryEstimate& estimate,
    std::uint64_t budget_bytes) {
  if (budget_bytes == 0U || estimate.budget_required_bytes <= budget_bytes) {
    return;
  }

  std::vector<std::pair<std::uint64_t, std::string>> contributors;
  contributors.reserve(estimate.report.entries.size());
  for (const core::MemoryEntry& entry : estimate.report.entries) {
    if (entry.label == "category_present" ||
        entry.lifetime == core::MemoryLifetime::kUnknown) {
      continue;
    }
    const std::uint64_t bytes = entryBudgetBytes(entry);
    if (bytes > 0U) {
      contributors.emplace_back(bytes, entry.label);
    }
  }
  std::sort(contributors.begin(), contributors.end(), [](const auto& lhs, const auto& rhs) {
    return lhs.first > rhs.first;
  });

  std::ostringstream message;
  message << "DMO process memory preflight requires "
          << estimate.budget_required_bytes
          << " bytes/rank (known_owned=" << estimate.known_owned_peak_bytes
          << ", external_reserve=" << estimate.external_unknown_reserve_bytes
          << ", safety_margin=" << estimate.safety_margin_bytes
          << ", aggregate_required=" << estimate.aggregate_required_bytes
          << ") but configured parallel.process_memory_budget_bytes="
          << budget_bytes << ". largest_known_contributors=";
  const std::size_t contributor_count = std::min<std::size_t>(5U, contributors.size());
  for (std::size_t i = 0; i < contributor_count; ++i) {
    if (i != 0U) {
      message << ',';
    }
    message << contributors[i].second << ':' << contributors[i].first;
  }
  throw std::runtime_error(message.str());
}

}  // namespace cosmosim::gravity
