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
  std::string name = "cosmosim_restart_v13";
  std::uint32_t version = 13;
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
  const core::HierarchicalTimeBinScheduler* scheduler = nullptr;
  core::ProvenanceRecord provenance;
  std::string normalized_config_text;
  std::string normalized_config_hash_hex;
  parallel::DistributedRestartState distributed_gravity_state;
  OutputCadencePersistentState output_cadence_state;
  StochasticPersistentState stochastic_state;
};

struct RestartReadResult {
  core::SimulationState state;
  core::IntegratorState integrator_state;
  core::TimeBinPersistentState scheduler_state;
  core::ProvenanceRecord provenance;
  std::string normalized_config_text;
  std::string normalized_config_hash_hex;
  parallel::DistributedRestartState distributed_gravity_state;
  OutputCadencePersistentState output_cadence_state;
  StochasticPersistentState stochastic_state;
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
