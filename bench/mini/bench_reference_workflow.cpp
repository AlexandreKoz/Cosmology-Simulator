#include <filesystem>
#include <sstream>

#include "bench/reporting/bench_report.hpp"
#include "cosmosim/cosmosim.hpp"

int main() {
  std::stringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[mode]\n";
  stream << "mode = zoom_in\n";
  stream << "ic_file = generated\n";

  const cosmosim::core::FrozenConfig frozen =
      cosmosim::core::loadFrozenConfigFromString(stream.str(), "bench_reference_workflow");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);

  const auto execution = cosmosim::bench::defaultExecutionConfig(1, 3);
  const std::filesystem::path run_directory = std::filesystem::temp_directory_path() / "cosmosim_bench_reference_workflow";

  for (std::size_t iter = 0; iter < execution.warmup_iterations; ++iter) {
    (void)runner.run(run_directory / ("warmup_" + std::to_string(iter)),
                     cosmosim::workflows::ReferenceWorkflowOptions{.step_index = iter, .write_outputs = false});
  }

  const auto begin = cosmosim::bench::BenchmarkClock::now();
  cosmosim::workflows::ReferenceWorkflowReport last_report;
  for (std::size_t iter = 0; iter < execution.measurement_iterations; ++iter) {
    last_report = runner.run(run_directory / ("measure_" + std::to_string(iter)),
                             cosmosim::workflows::ReferenceWorkflowOptions{.step_index = iter, .write_outputs = false});
  }
  const auto end = cosmosim::bench::BenchmarkClock::now();

  cosmosim::bench::BenchmarkReporter reporter("bench_reference_workflow");
  cosmosim::bench::addExecutionFields(reporter, execution);
  reporter.addField("measurement_ms", cosmosim::bench::BenchmarkClock::millisecondsBetween(begin, end));
  reporter.addField("stage_count", last_report.stage_sequence.size());
  reporter.addField("canonical_stage_order", last_report.canonical_stage_order ? 1 : 0);
  reporter.addField("config_compatible", last_report.config_compatible ? 1 : 0);
  reporter.addField("schema_compatible", last_report.schema_compatible ? 1 : 0);
  reporter.write();

  std::filesystem::remove_all(run_directory);
  return 0;
}
