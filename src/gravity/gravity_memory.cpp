#include "cosmosim/gravity/gravity_memory.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

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

}  // namespace

GravityMemoryEstimate estimateGravityMemory(const GravityMemoryEstimateInput& input) {
  if (input.tree_leaf_size == 0U || input.mpi_rank_count == 0U) {
    throw std::invalid_argument("gravity memory estimate requires non-zero leaf size and rank count");
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
  const std::uint64_t source_staging_bytes = checkedMul(
      input.local_source_count,
      5U * sizeof(double) + 3U * sizeof(std::uint32_t) + 2U * sizeof(std::uint8_t),
      "gravity source staging estimate overflow");
  const std::uint64_t target_view_bytes = checkedMul(
      input.local_target_count,
      5U * sizeof(std::uint32_t) + 2U * sizeof(double) + 2U * sizeof(std::uint8_t),
      "gravity target view estimate overflow");
  const std::uint64_t ordering_bytes = checkedMul(
      input.local_source_count,
      sizeof(std::uint64_t) + 2U * sizeof(std::uint32_t) + sizeof(double),
      "gravity ordering estimate overflow");
  const std::uint64_t acceleration_bytes = checkedMul(
      input.local_target_count,
      3U * sizeof(double),
      "gravity acceleration estimate overflow");
  // The legacy PM interpolation backend still requires contiguous target
  // coordinates internally. Source-index targets therefore avoid persistent
  // workflow/coordinator copies, but pay one bounded PM gather triplet while
  // interpolation is active.
  const std::uint64_t pm_target_gather_bytes = input.indexed_target_coordinates
      ? checkedMul(
          input.local_target_count,
          3U * sizeof(double),
          "gravity PM target gather estimate overflow")
      : 0U;
  const std::uint64_t periodic_tree_coordinate_bytes = input.periodic_tree_coordinates
      ? checkedMul(
          input.local_source_count,
          5U * sizeof(double),
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
      3U * sizeof(double) + sizeof(std::uint8_t),
      "gravity persistent force cache estimate overflow");
  const std::uint64_t runtime_particle_map_bytes = checkedMul(
      input.local_particle_count,
      3U * sizeof(std::int32_t),
      "gravity runtime particle map estimate overflow");
  const std::uint64_t runtime_cell_map_bytes = checkedMul(
      input.local_cell_count,
      2U * sizeof(std::int32_t) + sizeof(std::uint8_t),
      "gravity runtime cell map estimate overflow");
  const std::uint64_t runtime_refresh_list_bytes = checkedMul(
      input.local_particle_count, sizeof(std::uint32_t),
      "gravity runtime refresh-list estimate overflow");
  const std::uint64_t runtime_mapping_bytes = checkedAdd(
      checkedAdd(runtime_particle_map_bytes, runtime_cell_map_bytes,
                 "gravity runtime mapping estimate overflow"),
      runtime_refresh_list_bytes, "gravity runtime mapping estimate overflow");

  const std::uint64_t global_pm_cells = gridCells(input.pm_shape, "gravity PM grid estimate overflow");
  const std::uint64_t local_pm_cells = checkedAdd(
      global_pm_cells,
      static_cast<std::uint64_t>(input.mpi_rank_count) - 1U,
      "gravity PM local cell estimate overflow") / static_cast<std::uint64_t>(input.mpi_rank_count);
  // Density, potential and three force components. FFT backend work arrays are
  // explicitly external/unknown below rather than silently guessed as owned.
  const std::uint64_t pm_owned_bytes = checkedMul(
      local_pm_cells, 5U * sizeof(double), "gravity PM owned estimate overflow");
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

  const std::uint64_t mpi_bytes = input.mpi_rank_count > 1U
      ? checkedMul(input.tree_exchange_batch_bytes, 2U, "gravity exchange estimate overflow")
      : 0U;
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
  const std::uint64_t modeled_backend_unknown = checkedMul(
      local_pm_cells, 2U * sizeof(double),
      "gravity backend unknown estimate overflow");
  const std::uint64_t backend_unknown = checkedAdd(
      modeled_backend_unknown, input.backend_unknown_reserve_bytes,
      "gravity backend unknown reserve overflow");

  core::MemoryReportBuilder builder;
  addEstimate(builder, core::MemorySubsystem::kActiveSets, core::MemoryLifetime::kTransient,
              "gravity.estimate.source_staging", source_staging_bytes,
              "authoritative source staging; excludes canonical SimulationState");
  addEstimate(builder, core::MemorySubsystem::kActiveSets, core::MemoryLifetime::kTransient,
              "gravity.estimate.target_index_views", target_view_bytes,
              "targets alias source coordinates by compact local index");
  addEstimate(builder, core::MemorySubsystem::kActiveSets, core::MemoryLifetime::kTransient,
              "gravity.estimate.force_accumulators", acceleration_bytes,
              "authoritative compact active acceleration triplet");
  if (pm_target_gather_bytes > 0U) {
    addEstimate(builder, core::MemorySubsystem::kScratch, core::MemoryLifetime::kTransient,
                "gravity.estimate.pm_indexed_target_coordinate_gather", pm_target_gather_bytes,
                "bounded PM interpolation scratch; no persistent coordinator target triplet");
  }
  if (periodic_tree_coordinate_bytes > 0U) {
    addEstimate(builder, core::MemorySubsystem::kTree, core::MemoryLifetime::kTransient,
                "gravity.estimate.periodic_tree_coordinate_staging", periodic_tree_coordinate_bytes,
                "three unwrapped tree coordinates plus two reused axis-ordering workspaces");
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
              "active-slot/owned-local maps, leaf mask and force-refresh particle list");
  addEstimate(builder, core::MemorySubsystem::kTree, core::MemoryLifetime::kTransient,
              "gravity.estimate.tree_nodes", tree_nodes_bytes,
              "leaf-derived estimate; dynamic growth remains possible for adversarial geometry");
  addEstimate(builder, core::MemorySubsystem::kScratch, core::MemoryLifetime::kTransient,
              "gravity.estimate.tree_ordering_workspace", ordering_bytes);
  addEstimate(builder, core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kTransient,
              "gravity.estimate.pm_owned_fields", pm_owned_bytes);
  if (zoom_bytes > 0U) {
    addEstimate(builder, core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kTransient,
                "gravity.estimate.zoom_pm_owned_fields", zoom_bytes,
                "coarse/focused lifetimes are serialized; estimate reports the focused peak contribution");
  }
  if (mpi_bytes > 0U) {
    addEstimate(builder, core::MemorySubsystem::kMpiBuffers, core::MemoryLifetime::kTransient,
                "gravity.estimate.sparse_tree_exchange", mpi_bytes,
                "request/response high-water bounded by configured exchange batch policy");
  }
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
      "Gravity pre-run estimate includes owned source staging, compact target/force lanes, runtime index/selection maps, PM indexed-target scratch, periodic tree staging, tree workspace, PM fields, optional zoom lanes, persistent force cache, sparse exchange buffers, and known CUDA buffers; canonical SimulationState is reported separately.");
  result.report.notes.push_back(
      std::string("PM estimate profile assignment=") +
      (input.assignment_scheme == PmAssignmentScheme::kTsc ? "tsc" : "cic") +
      " decomposition=" +
      (input.decomposition_mode == core::PmDecompositionMode::kPencil
           ? "transposed_slab_contract"
           : "x_slab"));
  result.report.notes.push_back(
      "External FFT/CUDA/MPI library internal allocations remain explicitly uncertain and are not hidden inside the known peak.");
  result.external_backend_unknown_bytes = backend_unknown;
  result.estimated_tree_nodes = estimated_tree_nodes;
  std::uint64_t known_peak = 0U;
  for (const core::MemoryEntry& entry : result.report.entries) {
    if (entry.lifetime != core::MemoryLifetime::kUnknown) {
      known_peak = checkedAdd(known_peak, entry.estimated_next_step_bytes, "gravity known peak estimate overflow");
    }
  }
  result.known_peak_bytes = known_peak;
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

}  // namespace cosmosim::gravity
