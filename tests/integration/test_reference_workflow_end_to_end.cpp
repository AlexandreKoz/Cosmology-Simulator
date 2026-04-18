#include <cassert>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include "cosmosim/cosmosim.hpp"

int main() {
  std::stringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[mode]\n";
  stream << "mode = zoom_in\n";
  stream << "ic_file = generated\n";
  stream << "zoom_high_res_region = false\n\n";
  stream << "[numerics]\n";
  stream << "time_begin_code = 0.01\n";
  stream << "time_end_code = 0.0101\n";
  stream << "max_global_steps = 1\n";
  stream << "hierarchical_max_rung = 1\n\n";
  stream << "treepm_pm_grid = 9\n";
  stream << "treepm_asmth_cells = 1.75\n";
  stream << "treepm_rcut_cells = 6.0\n";
  stream << "treepm_assignment_scheme = cic\n";
  stream << "treepm_update_cadence_steps = 1\n\n";
  stream << "[output]\n";
  stream << "run_name = reference_integration_test\n";
  stream << "output_directory = integration_outputs\n";
  stream << "output_stem = snapshot\n";
  stream << "restart_stem = restart\n";

  const cosmosim::core::FrozenConfig frozen =
      cosmosim::core::loadFrozenConfigFromString(stream.str(), "test_reference_workflow_end_to_end");

  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  const std::filesystem::path output_dir =
      std::filesystem::temp_directory_path() / "cosmosim_reference_workflow_test";
  const cosmosim::workflows::ReferenceWorkflowReport report =
      runner.run(output_dir, cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = false});

  const std::filesystem::path expected_run_dir = output_dir / "reference_integration_test";

  assert(report.config_compatible);
  assert(report.schema_compatible);
  assert(report.canonical_stage_order);
  assert(report.stage_sequence.size() == 7);
  assert(report.stage_sequence.front() == "gravity_kick_pre");
  assert(report.stage_sequence.back() == "output_check");
  assert(report.completed_steps == 1);
  assert(report.treepm_pm_grid == 9);
  assert(report.run_directory == expected_run_dir);
  assert(report.normalized_config_snapshot_written);
  assert(std::filesystem::exists(report.normalized_config_snapshot_path));
  assert(std::filesystem::exists(report.profiler_json_path));
  assert(std::filesystem::exists(report.profiler_csv_path));
  assert(std::filesystem::exists(report.operational_report_json_path));

  std::ifstream op_in(report.operational_report_json_path);
  const std::string op_text((std::istreambuf_iterator<char>(op_in)), std::istreambuf_iterator<char>());
  assert(op_text.find("\"event_kind\": \"config.freeze\"") != std::string::npos);
  assert(op_text.find("\"provenance_config_hash_hex\"") != std::string::npos);
  assert(op_text.find("\"status\": \"ok\"") != std::string::npos);

  std::stringstream tsc_stream;
  tsc_stream << "schema_version = 1\n\n";
  tsc_stream << "[mode]\n";
  tsc_stream << "mode = zoom_in\n";
  tsc_stream << "ic_file = generated\n";
  tsc_stream << "zoom_high_res_region = false\n\n";
  tsc_stream << "[numerics]\n";
  tsc_stream << "time_begin_code = 0.01\n";
  tsc_stream << "time_end_code = 0.0101\n";
  tsc_stream << "max_global_steps = 1\n";
  tsc_stream << "hierarchical_max_rung = 1\n";
  tsc_stream << "treepm_pm_grid = 9\n";
  tsc_stream << "treepm_asmth_cells = 1.75\n";
  tsc_stream << "treepm_rcut_cells = 6.0\n";
  tsc_stream << "treepm_assignment_scheme = tsc\n";
  tsc_stream << "treepm_update_cadence_steps = 1\n\n";
  tsc_stream << "[output]\n";
  tsc_stream << "run_name = reference_integration_test_tsc\n";
  tsc_stream << "output_directory = integration_outputs\n";
  tsc_stream << "output_stem = snapshot\n";
  tsc_stream << "restart_stem = restart\n";

  const cosmosim::core::FrozenConfig tsc_frozen =
      cosmosim::core::loadFrozenConfigFromString(tsc_stream.str(), "test_reference_workflow_tsc");
  cosmosim::workflows::ReferenceWorkflowRunner tsc_runner(tsc_frozen);
  const cosmosim::workflows::ReferenceWorkflowReport tsc_report =
      tsc_runner.run(output_dir, cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = false});
  assert(tsc_report.completed_steps == 1);
  assert(tsc_report.treepm_pm_grid == 9);

  std::stringstream bad_stream;
  bad_stream << "schema_version = 1\n\n";
  bad_stream << "[mode]\n";
  bad_stream << "mode = zoom_in\n";
  bad_stream << "ic_file = generated\n";
  bad_stream << "zoom_high_res_region = false\n\n";
  bad_stream << "[numerics]\n";
  bad_stream << "time_begin_code = 0.01\n";
  bad_stream << "time_end_code = 0.0101\n";
  bad_stream << "max_global_steps = 1\n";
  bad_stream << "hierarchical_max_rung = 1\n";
  bad_stream << "treepm_pm_grid = 9\n";
  bad_stream << "treepm_asmth_cells = 1.75\n";
  bad_stream << "treepm_rcut_cells = 6.0\n";
  bad_stream << "treepm_assignment_scheme = cic\n";
  bad_stream << "treepm_update_cadence_steps = 2\n\n";
  bad_stream << "[output]\n";
  bad_stream << "run_name = reference_integration_test_bad_cadence\n";
  bad_stream << "output_directory = integration_outputs\n";
  bad_stream << "output_stem = snapshot\n";
  bad_stream << "restart_stem = restart\n";

  const cosmosim::core::FrozenConfig bad_frozen =
      cosmosim::core::loadFrozenConfigFromString(bad_stream.str(), "test_reference_workflow_bad_cadence");
  cosmosim::workflows::ReferenceWorkflowRunner bad_runner(bad_frozen);
  bool bad_cadence_threw = false;
  try {
    (void)bad_runner.run(output_dir, cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = false});
  } catch (const std::runtime_error&) {
    bad_cadence_threw = true;
  }
  assert(bad_cadence_threw);

  std::filesystem::remove_all(output_dir);
  return 0;
}
