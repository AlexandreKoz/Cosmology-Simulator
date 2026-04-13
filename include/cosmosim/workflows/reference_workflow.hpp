#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"

namespace cosmosim::workflows {

// Integration-layer ownership note:
// This runner intentionally lives outside core/ because it assembles analysis, physics,
// and I/O callbacks into a concrete end-to-end workflow.
struct ReferenceWorkflowOptions {
  std::uint64_t step_index = 0;
  double dt_time_code = 1.0e-3;
  bool write_outputs = true;
};

struct ReferenceWorkflowReport {
  bool config_compatible = false;
  bool schema_compatible = false;
  bool canonical_stage_order = false;
  bool restart_roundtrip_executed = false;
  bool snapshot_roundtrip_executed = false;
  bool restart_roundtrip_ok = false;
  bool snapshot_roundtrip_ok = false;
  std::vector<std::string> stage_sequence;
  std::filesystem::path profiler_json_path;
  std::filesystem::path profiler_csv_path;
  std::filesystem::path restart_path;
  std::filesystem::path snapshot_path;
};

class ReferenceWorkflowRunner {
 public:
  explicit ReferenceWorkflowRunner(core::FrozenConfig frozen_config);

  [[nodiscard]] const core::FrozenConfig& frozenConfig() const noexcept;

  [[nodiscard]] ReferenceWorkflowReport run(
      const std::filesystem::path& run_directory,
      const ReferenceWorkflowOptions& options = {}) const;

 private:
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
