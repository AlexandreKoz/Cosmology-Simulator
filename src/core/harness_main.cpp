#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

#include "cosmosim/cosmosim.hpp"

namespace {

int printUsage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " <config.param.txt>\n";
  return 2;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    return printUsage(argc > 0 ? argv[0] : "cosmosim_harness");
  }

  try {
    const std::filesystem::path config_path = argv[1];
    const cosmosim::core::FrozenConfig frozen = cosmosim::core::loadFrozenConfigFromFile(config_path, {});
    cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
    const cosmosim::workflows::ReferenceWorkflowReport report = runner.run();

    std::cout << cosmosim::core::projectName() << ' ' << cosmosim::core::versionString() << '\n';
    std::cout << "config=" << config_path.string() << '\n';
    std::cout << "run_directory=" << report.run_directory.string() << '\n';
    std::cout << "completed_steps=" << report.completed_steps << '\n';
    std::cout << "normalized_config=" << report.normalized_config_snapshot_path.string() << '\n';
    std::cout << "operational_report=" << report.operational_report_json_path.string() << '\n';
    if (!report.snapshot_path.empty()) {
      std::cout << "last_snapshot=" << report.snapshot_path.string() << '\n';
    }
    if (!report.restart_path.empty()) {
      std::cout << "last_restart=" << report.restart_path.string() << '\n';
    }
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "cosmosim runtime failed: " << ex.what() << '\n';
    return 1;
  }
}
