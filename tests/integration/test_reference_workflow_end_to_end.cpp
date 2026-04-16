#include <cassert>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>

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

  std::filesystem::remove_all(output_dir);
  return 0;
}
