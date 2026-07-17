#pragma once

#include <filesystem>

#include "cosmosim/core/config.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"

namespace cosmosim::workflows::internal {

struct InitialConditionStartupResult {
  core::SimulationState state;
  io::IcImportReport import_report;
  std::filesystem::path manifest_path;
  bool restoring_from_restart = false;
};

// Owns the single supported IC/restart/test-override materialization path and
// state finalization. Distributed IC materialization is deliberately absent.
class InitialConditionRuntime {
 public:
  explicit InitialConditionRuntime(const core::FrozenConfig& frozen_config) noexcept;

  [[nodiscard]] InitialConditionStartupResult materialize(
      const ReferenceWorkflowOptions& options,
      const std::filesystem::path& run_directory) const;

 private:
  const core::FrozenConfig& m_frozen_config;
};

}  // namespace cosmosim::workflows::internal
