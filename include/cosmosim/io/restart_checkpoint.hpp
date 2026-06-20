#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::io {

struct RestartSchema {
  // v19 persists the gas-cell scheduler independently of the particle scheduler.
  std::string name = "cosmosim_restart_v19";
  std::uint32_t version = 19;
};

[[nodiscard]] const RestartSchema& restartSchema();
[[nodiscard]] bool isRestartSchemaCompatible(std::uint32_t file_schema_version);
[[nodiscard]] const std::vector<std::string_view>& exactRestartCompletenessChecklist();

struct RestartWritePolicy {
  bool enable_fsync_finalize = false;
  std::string temporary_suffix = ".tmp";
};

struct RestartPersistentStateView {
  // Narrow restart-serialization root: persistent-only state ownership path.
  // TransientStepWorkspace/HydroScratch/PM/MPI/output scratch are intentionally
  // outside this view so restart traversal cannot accidentally serialize them.
  const core::SimulationState* simulation_state = nullptr;
};

struct OutputCadencePersistentState {
  bool output_enabled = false;
  bool write_restarts = false;
  bool snapshot_due = false;
  bool checkpoint_due = false;
  std::uint64_t last_completed_step_index = 0;
  std::uint64_t snapshot_interval_steps = 0;
  std::uint64_t next_snapshot_step_index = 0;
  std::string snapshot_stem;
  std::string restart_stem;
};

struct StochasticModulePersistentState {
  std::string module_name;
  std::uint32_t schema_version = 1;
  std::string rng_policy;
  std::uint64_t random_seed = 0;
  std::uint32_t rank_local_seed_offset = 0;
  std::uint64_t last_committed_step_index = 0;
  bool deterministic_from_serialized_inputs = true;
};

struct StochasticPersistentState {
  std::vector<StochasticModulePersistentState> modules;
};

struct RestartWritePayload {
  RestartPersistentStateView persistent_state;
  const core::IntegratorState* integrator_state = nullptr;
  // Particle scheduler; particles.time_bin is a derived compatibility mirror.
  const core::HierarchicalTimeBinScheduler* scheduler = nullptr;
  // Authoritative scheduler for gas-cell rows keyed by stable gas_cell_id in the
  // workflow-level persistence contract.  Required for v19 writes whenever gas
  // cells exist.
  const core::HierarchicalTimeBinScheduler* gas_cell_scheduler = nullptr;
  core::ProvenanceRecord provenance;
  std::string normalized_config_text;
  std::string normalized_config_hash_hex;
  parallel::DistributedRestartState distributed_gravity_state;
  OutputCadencePersistentState output_cadence_state;
  StochasticPersistentState stochastic_state;
};


struct RestartDiagnosticsSummary {
  std::string restart_schema_name;
  std::uint32_t restart_schema_version = 0;
  std::string current_boundary_kind;
  std::string last_completed_boundary_kind;
  bool restart_safe = false;
  std::uint64_t step_index = 0;
  std::uint64_t scheduler_current_tick = 0;
  std::uint32_t scheduler_max_bin = 0;
  std::uint64_t scheduler_element_count = 0;
  std::uint64_t scheduler_active_count = 0;
  std::uint64_t scheduler_pending_transition_count = 0;
  std::uint64_t gas_cell_scheduler_current_tick = 0;
  std::uint32_t gas_cell_scheduler_max_bin = 0;
  std::uint64_t gas_cell_scheduler_element_count = 0;
  std::uint64_t gas_cell_scheduler_active_count = 0;
  std::uint64_t gas_cell_scheduler_pending_transition_count = 0;
  std::uint64_t pm_cadence_steps = 0;
  std::uint64_t pm_gravity_kick_opportunity = 0;
  std::uint64_t pm_field_version = 0;
  std::uint64_t pm_last_refresh_opportunity = 0;
  std::uint64_t pm_last_refresh_step_index = 0;
  bool pm_refresh_commit_pending = false;
  bool pm_long_range_field_valid = false;
  bool output_enabled = false;
  bool output_snapshot_due = false;
  bool output_checkpoint_due = false;
  std::uint64_t output_last_completed_step_index = 0;
  std::uint64_t output_next_snapshot_step_index = 0;
  std::uint64_t stochastic_module_count = 0;
};

struct RestartReadResult {
  core::SimulationState state;
  core::IntegratorState integrator_state;
  core::TimeBinPersistentState scheduler_state;
  // v19+: state keyed by stable gas_cell_id.  For v14-v18 reads this is
  // reconstructed from the historic cell time_bin mirror for compatibility.
  core::TimeBinPersistentState gas_cell_scheduler_state;
  std::vector<std::uint64_t> gas_cell_scheduler_ids;
  core::ProvenanceRecord provenance;
  std::string normalized_config_text;
  std::string normalized_config_hash_hex;
  parallel::DistributedRestartState distributed_gravity_state;
  OutputCadencePersistentState output_cadence_state;
  StochasticPersistentState stochastic_state;
  RestartDiagnosticsSummary diagnostics;
  std::uint64_t payload_hash = 0;
  std::string payload_hash_hex;
};

[[nodiscard]] std::uint64_t restartPayloadIntegrityHash(const RestartWritePayload& payload);
[[nodiscard]] std::string restartPayloadIntegrityHashHex(const RestartWritePayload& payload);

void writeRestartCheckpointHdf5(
    const std::filesystem::path& output_path,
    const RestartWritePayload& payload,
    const RestartWritePolicy& policy = {});

[[nodiscard]] RestartReadResult readRestartCheckpointHdf5(const std::filesystem::path& input_path);

}  // namespace cosmosim::io
