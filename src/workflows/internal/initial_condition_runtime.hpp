#pragma once

#include <filesystem>

#include "cosmosim/core/config.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"
#include "cosmosim/workflows/runtime_services.hpp"

namespace cosmosim::workflows::internal {

struct InitialConditionStartupResult {
  core::SimulationState state;
  io::IcImportReport import_report;
  std::filesystem::path manifest_path;
  bool restoring_from_restart = false;
  bool already_partitioned = false;
};

// Owns the single supported IC/restart/test-override materialization path,
// including bounded distributed external IC ingestion and state finalization.
class InitialConditionRuntime {
 public:
  InitialConditionRuntime(
      const core::FrozenConfig& frozen_config,
      const RuntimeServices& services) noexcept;

  [[nodiscard]] InitialConditionStartupResult materialize(
      const ReferenceWorkflowOptions& options,
      const std::filesystem::path& run_directory) const;

 private:
  const core::FrozenConfig& m_frozen_config;
  const RuntimeServices& m_services;
};

}  // namespace cosmosim::workflows::internal
