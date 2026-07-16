#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/time_integration.hpp"

namespace cosmosim::core {
struct SimulationState;
}

namespace cosmosim::io {
struct RestartReadResult;
}

namespace cosmosim::workflows {

// Integration-layer ownership note:
// This runner intentionally lives outside core/ because it assembles analysis, physics,
// and I/O callbacks into a concrete end-to-end workflow.
struct ReferenceWorkflowOptions {
  std::uint64_t step_index = 0;
  double dt_time_code = 0.0;
  bool write_outputs = true;

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
};

struct ReferenceWorkflowReport {
  struct TreePmCadenceRecord {
    std::uint64_t step_index = 0;
    std::string stage_name;
    std::string pm_sync_surface;
    std::uint64_t gravity_kick_opportunity = 0;
    std::uint64_t field_version = 0;
    std::uint64_t last_refresh_opportunity = 0;
    std::uint64_t field_built_step_index = 0;
    double field_built_scale_factor = 1.0;
    std::uint64_t field_age_in_kick_opportunities = 0;
    std::uint64_t active_particles_kicked = 0;
    std::uint64_t inactive_particles_skipped = 0;
    bool refreshed_long_range_field = false;
  };

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
  std::vector<std::string> stage_sequence;
  std::vector<TreePmCadenceRecord> treepm_cadence_records;
  std::filesystem::path run_directory;
  std::filesystem::path normalized_config_snapshot_path;
  std::filesystem::path profiler_json_path;
  std::filesystem::path profiler_csv_path;
  std::filesystem::path operational_report_json_path;
  std::filesystem::path runtime_capability_report_path;
  std::filesystem::path ic_manifest_path;
  std::filesystem::path restart_path;
  std::filesystem::path snapshot_path;
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

namespace cosmosim::core {
// Transitional namespace compatibility aliases.
// TODO(architecture): remove these aliases after downstream callers migrate to cosmosim::workflows.
using ReferenceWorkflowOptions = workflows::ReferenceWorkflowOptions;
using ReferenceWorkflowReport = workflows::ReferenceWorkflowReport;
using ReferenceWorkflowRunner = workflows::ReferenceWorkflowRunner;
}  // namespace cosmosim::core
