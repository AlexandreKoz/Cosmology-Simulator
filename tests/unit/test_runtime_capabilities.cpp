#include "cosmosim/workflows/runtime_capabilities.hpp"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>

namespace {

void testCapabilityTruth() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.numerics.hierarchical_max_rung = 0;
  const auto report = cosmosim::workflows::buildRuntimeCapabilityReport(config);
  assert(report.require("fixed_global_timestep").status ==
         cosmosim::workflows::RuntimeCapabilityStatus::kSupported);
  assert(report.require("production_hierarchical_local_timestep").status ==
         cosmosim::workflows::RuntimeCapabilityStatus::kUnsupported);
  assert(report.require("distributed_ic_import").status ==
         cosmosim::workflows::RuntimeCapabilityStatus::kUnsupported);
  assert(report.require("canonical_external_ic_import").status ==
         cosmosim::workflows::RuntimeCapabilityStatus::kProvisional);
  cosmosim::workflows::validateRequestedRuntimeCapabilities(config, report);

  config.numerics.hierarchical_max_rung = 1;
  bool threw = false;
  try {
    cosmosim::workflows::validateRequestedRuntimeCapabilities(config, report);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
}

void testJsonSerialization() {
  const auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  const auto report = cosmosim::workflows::buildRuntimeCapabilityReport(config);
  const std::string json =
      cosmosim::workflows::serializeRuntimeCapabilityReportJson(report);
  assert(json.find("\"schema_version\": 1") != std::string::npos);
  assert(json.find("\"name\": \"rank_remappable_restart\"") !=
         std::string::npos);
  assert(json.find("\"status\": \"unsupported\"") !=
         std::string::npos);

  const auto path = std::filesystem::current_path() /
      "chui_runtime_capabilities_test.json";
  cosmosim::workflows::writeRuntimeCapabilityReportJson(report, path);
  std::ifstream input(path, std::ios::binary);
  const std::string readback{
      std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
  assert(readback == json);
  std::filesystem::remove(path);
}

}  // namespace

int main() {
  testCapabilityTruth();
  testJsonSerialization();
  return 0;
}
