#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/workflows/runtime_console_reporter.hpp"

namespace cosmosim::core {
class SimulationState;
}

namespace cosmosim::io {
struct RestartReadResult;
}

namespace cosmosim::workflows {

class RuntimeModuleRegistry;

// Integration-layer ownership note:
// This runner intentionally lives outside core/ because it assembles analysis, physics,
// and I/O callbacks into a concrete end-to-end workflow.
struct ReferenceWorkflowOptions {
  std::uint64_t step_index = 0;
  double dt_time_code = 0.0;
  bool write_outputs = true;

  // Process-local presentation controls. These never enter normalized config,
  // provenance, snapshots, or restart state. Embedders are quiet by default;
  // cosmosim_harness explicitly enables native console reporting.
  RuntimeConsoleOptions console{};

  // Narrow integration-test/benchmark seam.  When present, the runner copies
  // this already-validated state and executes the same production setup,
  // scheduler, hydro, restart, and output pipeline used after IC import.
  // The caller retains ownership for the synchronous duration of run().
  const core::SimulationState* initial_state_override = nullptr;

  // Narrow, explicit restart-resume seam used by integration tests and restart
  // tooling. The restored payload remains the restart-authoritative state:
  // SimulationState, integrator state, and particle/gas-cell scheduler state
  // are imported without reinitializing time bins. The caller retains ownership
  // for the synchronous duration of run().
  const io::RestartReadResult* restart_state_override = nullptr;

  // Narrow MPI regression seam. When supplied for a fresh run, these stable
  // particle-ID records replace only the initialized particle scheduler state
  // before the production step loop. They are rejected for restart-resume so
  // checkpoint scheduler authority remains intact.
  std::vector<core::TimeBinSchedulerIdentityRecord> initial_particle_scheduler_identity_records;

  // A bounded execution segment for tests/benchmarks. Zero retains the
  // configuration's max_global_steps. This is intentionally not persisted.
  std::uint64_t max_steps_override = 0;
  // Test-only seam (multirate_outer_kick): the first N PM cycles that reach an
  // interior boundary are rejected once and retried with the identical outer
  // interval, exercising the transactional rollback; the final state must equal
  // an uninterrupted run. 0 (default) in production; never weakens a check.
  std::uint64_t multirate_test_forced_cycle_rollbacks = 0;
  // Test-only seams (multirate_outer_kick); 0/empty/false in production.
  // * forced_rejections_per_cycle: consecutive forced rejections of one cycle
  //   (each retried with the identical T; must stay below the retry bound 64).
  // * snapshot_refusals: the next N rollback-snapshot admissions are treated
  //   as refused on the listed ranks (all ranks when the list is empty), as a
  //   MemoryGovernor refusal or, with snapshot_allocation_failure, as a
  //   std::bad_alloc while copying. Admission is reduced collectively, so a
  //   refusal on one rank must produce the same fallback on every rank.
  std::uint64_t multirate_test_forced_rejections_per_cycle = 1;
  std::uint64_t multirate_test_snapshot_refusals = 0;
  std::vector<int> multirate_test_snapshot_refusal_ranks;
  bool multirate_test_snapshot_allocation_failure = false;
  // Corrupts the captured force endpoint's row-identity certificate on this
  // world rank (-1: never), so the production row-aligned restore validation
  // refuses the rollback there only;
  // every rank must then fail collectively, with no state restored on any
  // rank and no rank left inside a collective.
  int multirate_test_rollback_failure_rank = -1;
  // Test-only: rebuild the hierarchical identity row list and rescan the
  // local maximum mass every block (the uncached reference for the
  // scheduler-preparation caches; results must be bitwise identical).
  bool hierarchical_test_disable_preparation_caches = false;

  // Test/embedding composition seam. Built-in production descriptors are
  // registered first; callers may add independent descriptors before the
  // registry freezes. This does not alter config or restart authority.
  std::function<void(RuntimeModuleRegistry&)> register_runtime_modules;
};

struct ReferenceWorkflowReport {
  struct TreePmCadenceRecord {
    std::uint64_t step_index = 0;
    std::string stage_name;
    std::string pm_sync_surface;
    std::string pm_refresh_reason;
    std::uint64_t gravity_kick_opportunity = 0;
    std::uint64_t field_version = 0;
    std::uint64_t last_refresh_opportunity = 0;
    std::uint64_t field_built_step_index = 0;
    double field_built_scale_factor = 1.0;
    std::uint64_t field_age_in_kick_opportunities = 0;
    std::uint64_t active_particles_kicked = 0;
    std::uint64_t inactive_particles_skipped = 0;
    bool refreshed_long_range_field = false;
    bool short_range_only = false;
    std::uint64_t active_force_targets = 0U;
  };

  // Derived, versioned migration provenance; original checkpoint remains immutable.
  std::string restart_config_compatibility_record;
  bool config_compatible = false;
  bool schema_compatible = false;
  bool canonical_stage_order = false;
  bool normalized_config_snapshot_written = false;
  bool restart_roundtrip_executed = false;
  bool snapshot_roundtrip_executed = false;
  bool restart_roundtrip_ok = false;
  bool snapshot_roundtrip_ok = false;
  std::uint64_t completed_steps = 0;
  // Committed endpoint observability. These are additive report fields; they
  // do not change snapshot/restart schemas or integration-state ownership.
  double final_time_code = 0.0;
  double final_scale_factor = 1.0;
  // Rank-count and dense-row-order independent digest over authoritative
  // physical entities keyed by stable particle/gas-cell identity.
  std::uint64_t final_state_digest = 0;
  std::uint64_t local_particle_count = 0;
  std::uint64_t global_particle_count = 0;
  std::uint64_t local_cell_count = 0;
  std::uint64_t global_cell_count = 0;
  std::uint64_t local_particle_id_sum = 0;
  std::uint64_t global_particle_id_sum = 0;
  std::uint64_t local_particle_id_xor = 0;
  std::uint64_t global_particle_id_xor = 0;
  bool local_particle_ids_unique = true;
  bool global_particle_partition_identity_match = false;
  int world_size = 1;
  int world_rank = 0;
  std::uint64_t treepm_long_range_refresh_count = 0;
  std::uint64_t treepm_long_range_reuse_count = 0;
  // Last production HydroStageCallback CFL record. This exposes workflow-level
  // physical identity diagnostics to integration tests without coupling callers
  // to HydroStageCallback implementation details.
  core::HydroCflDiagnostics final_hydro_cfl_diagnostics{};
  std::uint64_t final_hydro_imported_mpi_ghosts = 0;
  std::uint64_t final_hydro_remote_interface_faces = 0;
  std::uint64_t final_hydro_remote_stale_payloads = 0;
  std::size_t treepm_pm_grid = 0;
  std::size_t treepm_pm_grid_nx = 0;
  std::size_t treepm_pm_grid_ny = 0;
  std::size_t treepm_pm_grid_nz = 0;
  std::string treepm_pm_grid_shape;
  int treepm_update_cadence_steps = 1;
  // Recent stage/cadence details only (128 stages / 256 force decisions).
  // Lifetime counts remain in the profiler and refresh/reuse totals.
  std::vector<std::string> stage_sequence;
  std::vector<TreePmCadenceRecord> treepm_cadence_records;
  // Hierarchical coarse-block evidence (2026-10-26 stabilization; recent 512
  // executed blocks). Additive report fields; not persisted. Rejected
  // multirate attempts are rolled back and never appear here.
  struct HierarchicalBlockRecord {
    std::uint64_t step_index = 0U;
    double time_begin_code = 0.0;
    double time_end_code = 0.0;
    double coarse_interval_time_code = 0.0;
    double quantum_time_code = 0.0;
    // Global minimum per-particle criterion (gravity, encounter, growth
    // limiter) at the block's opening synchronization point.
    double global_min_particle_criterion_time_code = 0.0;
    std::string coarse_limiter;
    bool multirate = false;
    std::uint64_t pm_cycle_index = 0U;
    std::uint64_t pm_cycle_block_index = 0U;
    bool opens_pm_cycle = true;
    bool closes_pm_cycle = true;
    double pm_cycle_opening_time_code = 0.0;
    double pm_cycle_end_time_code = 0.0;
    unsigned pm_cycle_retries = 0U;
    // Global PM displacement bound (code time) evaluated at this block's
    // opening synchronization point; +infinity when no speed bounds it.
    // Multirate invariant: pm_cycle_end - pm_cycle_opening <= this bound for
    // every block of an accepted cycle.
    double pm_displacement_limit_time_code = 0.0;
    // Limiter that fixed this cycle's length T: target_blocks,
    // pm_displacement, output_or_end_boundary, rollback_retry (physical or
    // retry-selected) or rollback_memory (resource-induced single-block
    // fallback; see numerics.treepm_multirate_rollback_memory_policy).
    std::string pm_cycle_length_limiter;
  };
  std::vector<HierarchicalBlockRecord> hierarchical_block_records;
  std::uint64_t multirate_pm_cycles_completed = 0U;
  std::uint64_t multirate_pm_cycle_rollbacks = 0U;
  std::uint64_t multirate_forced_test_rollbacks = 0U;
  // Work of rejected (rolled-back) multirate attempts, excluded from the
  // trajectory diagnostics above (treepm_long_range_*_count, cadence records,
  // hierarchical_block_records). Profiler timings still include it (cost).
  std::uint64_t multirate_rejected_blocks = 0U;
  std::uint64_t multirate_rejected_pm_refreshes = 0U;
  std::uint64_t multirate_rejected_pm_reuses = 0U;
  std::uint64_t multirate_rejected_cadence_records = 0U;
  // Governed rollback-snapshot storage: high-water reservation of the segment
  // and the reservation still held when the segment returned (always 0).
  std::uint64_t multirate_rollback_snapshot_peak_bytes = 0U;
  std::uint64_t multirate_rollback_snapshot_retained_bytes = 0U;
  // Highest capacity actually retained by the rollback store (must never
  // exceed multirate_rollback_snapshot_peak_bytes) and the number of fresh
  // (re)allocations of the store (1 for a run whose extents never grow and
  // whose admissions are never refused: capacity is reused across cycles
  // and retries).
  std::uint64_t multirate_rollback_snapshot_peak_capacity_bytes = 0U;
  std::uint64_t multirate_rollback_snapshot_allocations = 0U;
  // Hierarchical scheduler-preparation caches (segment-local; additive,
  // not persisted): identity row list rebuilds/reuses and local maximum
  // mass scans/reuses. One rebuild and one scan per segment when the
  // particle rows and masses never change.
  std::uint64_t hierarchical_identity_row_rebuilds = 0U;
  std::uint64_t hierarchical_identity_row_reuses = 0U;
  std::uint64_t hierarchical_mass_max_scans = 0U;
  std::uint64_t hierarchical_mass_max_reuses = 0U;
  // Multi-block cycles opened as a single block because the rollback
  // snapshot was refused by the memory governor on some rank.
  std::uint64_t multirate_memory_limited_cycles = 0U;
  // Blocks executed beyond the requested step budget to close an open PM
  // cycle (soft stop); completed_steps includes them.
  std::uint64_t multirate_soft_stop_extension_blocks = 0U;
  std::uint64_t requested_steps = 0U;
  std::filesystem::path shared_run_directory;
  std::filesystem::path run_directory;
  std::filesystem::path normalized_config_snapshot_path;
  std::filesystem::path profiler_json_path;
  std::filesystem::path profiler_csv_path;
  std::filesystem::path operational_report_json_path;
  std::filesystem::path runtime_capability_report_path;
  std::filesystem::path ic_manifest_path;
  std::filesystem::path restart_path;
  std::filesystem::path snapshot_path;
  std::filesystem::path snapshot_set_path;
};

class ReferenceWorkflowRunner {
 public:
  explicit ReferenceWorkflowRunner(core::FrozenConfig frozen_config);

  [[nodiscard]] const core::FrozenConfig& frozenConfig() const noexcept;

  [[nodiscard]] ReferenceWorkflowReport run(
      const ReferenceWorkflowOptions& options = {}) const;

  // Testing/benchmark override for the output root only. The final run directory remains
  // <output_root_override>/<output.run_name>, so the typed config still governs run identity.
  [[nodiscard]] ReferenceWorkflowReport run(
      const std::filesystem::path& output_root_override,
      const ReferenceWorkflowOptions& options) const;

 private:
  [[nodiscard]] ReferenceWorkflowReport runImpl(
      const std::filesystem::path* output_root_override,
      const ReferenceWorkflowOptions& options) const;

  core::FrozenConfig m_frozen_config;
};

}  // namespace cosmosim::workflows
