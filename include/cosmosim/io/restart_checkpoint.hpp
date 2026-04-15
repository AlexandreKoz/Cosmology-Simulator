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

namespace cosmosim::io {

struct RestartSchema {
  std::string name = "cosmosim_restart_v3";
  std::uint32_t version = 3;
};

[[nodiscard]] const RestartSchema& restartSchema();
[[nodiscard]] bool isRestartSchemaCompatible(std::uint32_t file_schema_version);
[[nodiscard]] const std::vector<std::string_view>& exactRestartCompletenessChecklist();

struct RestartWritePolicy {
  bool enable_fsync_finalize = false;
  std::string temporary_suffix = ".tmp";
};

struct RestartWritePayload {
  const core::SimulationState* state = nullptr;
  const core::IntegratorState* integrator_state = nullptr;
  const core::HierarchicalTimeBinScheduler* scheduler = nullptr;
  core::ProvenanceRecord provenance;
  std::string normalized_config_text;
  std::string normalized_config_hash_hex;
};

struct RestartReadResult {
  core::SimulationState state;
  core::IntegratorState integrator_state;
  core::TimeBinPersistentState scheduler_state;
  core::ProvenanceRecord provenance;
  std::string normalized_config_text;
  std::string normalized_config_hash_hex;
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
