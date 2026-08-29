#include "cosmosim/gravity/gravity_memory.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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

  (void)gridCells(input.pm_shape, "gravity PM grid estimate overflow");
  const parallel::PmSlabLayout pm_layout = parallel::makePmSlabLayout(
      input.pm_shape.nx,
      input.pm_shape.ny,
      input.pm_shape.nz,
      static_cast<int>(input.mpi_rank_count),
      input.mpi_world_rank);
  const std::uint64_t local_pm_cells =
      static_cast<std::uint64_t>(pm_layout.localCellCount());
  // Density, potential and three force components.
  const std::uint64_t pm_owned_bytes = checkedMul(
      local_pm_cells, 5U * sizeof(double), "gravity PM owned estimate overflow");
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

  const std::uint64_t tree_mpi_bytes = input.mpi_rank_count > 1U
      ? checkedMul(input.tree_exchange_batch_bytes, 2U, "gravity tree exchange estimate overflow")
      : 0U;
  // PM density and force routing retain exactly two reusable wire buffers. The
  // configured policy is a per-peer payload ceiling, so a rank can receive one
  // bounded chunk from every remote peer in the same collective round. This is
  // O(rank_count * batch_bytes), independent of local particle count and TSC's
  // 27-cell stencil.
  const std::uint64_t remote_rank_count = input.mpi_rank_count > 1U
      ? static_cast<std::uint64_t>(input.mpi_rank_count - 1U)
      : 0U;
  const std::uint64_t pm_routing_bytes = remote_rank_count == 0U
      ? 0U
      : checkedMul(
          checkedMul(input.pm_exchange_batch_bytes, remote_rank_count,
                     "gravity PM routing estimate overflow"),
          2U,
          "gravity PM routing estimate overflow");
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
  addEstimate(builder, core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kPersistent,
              "gravity.estimate.pm_plan_resources_owned_arrays",
              pm_plan_memory.total_owned_bytes,
              pm_plan_memory.used_backend_allocation_query
                  ? "CHUI-owned FFT/Poisson PlanResources arrays sized from the active FFTW MPI allocation query; backend plan internals excluded"
                  : "CHUI-owned FFT/Poisson PlanResources arrays sized from conservative PM decomposition geometry; backend plan internals excluded");
  if (zoom_bytes > 0U) {
    addEstimate(builder, core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kTransient,
                "gravity.estimate.zoom_pm_owned_fields", zoom_bytes,
                "coarse/focused lifetimes are serialized; estimate reports the focused peak contribution");
  }
  if (tree_mpi_bytes > 0U) {
    addEstimate(builder, core::MemorySubsystem::kMpiBuffers, core::MemoryLifetime::kTransient,
                "gravity.estimate.sparse_tree_exchange", tree_mpi_bytes,
                "request/response high-water bounded by configured tree exchange batch policy");
  }
  if (pm_routing_bytes > 0U) {
    addEstimate(builder, core::MemorySubsystem::kMpiBuffers, core::MemoryLifetime::kTransient,
                "gravity.estimate.pm_routing_exchange", pm_routing_bytes,
                "two reusable PM wire buffers; conservative all-remote-peer high-water from the configured per-peer batch policy, independent of particle count/stencil size");
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
      "Gravity pre-run estimate includes owned source staging, compact target/force lanes, runtime index/selection maps, PM indexed-target scratch, periodic tree staging, tree workspace, PM fields, explicit PM PlanResources arrays, optional zoom lanes, persistent force cache, sparse tree exchange, bounded PM routing buffers, and known CUDA buffers; canonical SimulationState is reported separately.");
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
