#include "cosmosim/workflows/runtime_capabilities.hpp"
#include "cosmosim/core/build_config.hpp"

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
  const auto& openmp = report.require("openmp_execution");
#if COSMOSIM_HAVE_OPENMP
  assert(openmp.compiled);
  assert(openmp.runtime_available);
#else
  assert(!openmp.compiled);
#endif
  assert(report.require("production_fof_halo_finder").status ==
         cosmosim::workflows::RuntimeCapabilityStatus::kSupported);
  assert(report.require("production_fft_power_spectrum").status ==
         cosmosim::workflows::RuntimeCapabilityStatus::kSupported);
  assert(report.require("bound_subhalo_finder").status ==
         cosmosim::workflows::RuntimeCapabilityStatus::kUnsupported);
  assert(report.require("fixed_global_timestep").status ==
         cosmosim::workflows::RuntimeCapabilityStatus::kSupported);
  assert(report.require("production_hierarchical_local_timestep").status ==
         cosmosim::workflows::RuntimeCapabilityStatus::kUnsupported);
#if COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI
  assert(report.require("distributed_ic_import").status ==
         cosmosim::workflows::RuntimeCapabilityStatus::kProvisional);
#else
  assert(report.require("distributed_ic_import").status ==
         cosmosim::workflows::RuntimeCapabilityStatus::kUnsupported);
#endif
  const auto canonical_status =
      report.require("canonical_external_ic_import").status;
  const auto multifile_status =
      report.require("multifile_external_ic_import").status;
  assert(canonical_status == multifile_status);
  assert(canonical_status !=
         cosmosim::workflows::RuntimeCapabilityStatus::kProvisional);
#if COSMOSIM_ENABLE_HDF5
  const auto& canonical_detail =
      report.require("canonical_external_ic_import").detail;
  assert(canonical_detail.find("manifest v4") != std::string::npos);
  assert(canonical_detail.find("source-chunk-to-canonical") !=
         std::string::npos);
#endif
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

void testPowerSpectrumCapabilitySchedulingTruth() {
  auto blocked = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  blocked.analysis.enable_diagnostics = true;
  blocked.analysis.diagnostics_execution_policy =
      cosmosim::core::AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthAndLightScience;
  const auto blocked_report = cosmosim::workflows::buildRuntimeCapabilityReport(blocked);
  const auto& blocked_fft = blocked_report.require("production_fft_power_spectrum");
  assert(blocked_fft.requested);
  assert(!blocked_fft.eligible);
  assert(!blocked_fft.scheduled);
  assert(!blocked_fft.active);
  assert(blocked_fft.implementation_maturity ==
         cosmosim::workflows::RuntimeImplementationMaturity::kProductionScalable);
  assert(blocked_fft.scientific_maturity ==
         cosmosim::workflows::RuntimeScientificMaturity::kProvisional);
  assert(blocked_fft.scalability ==
         cosmosim::workflows::RuntimeScalabilityClass::kScalableFft);

  auto allowed = blocked;
  allowed.analysis.diagnostics_execution_policy =
      cosmosim::core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional;
  const auto allowed_report = cosmosim::workflows::buildRuntimeCapabilityReport(allowed);
  const auto& allowed_fft = allowed_report.require("production_fft_power_spectrum");
  assert(allowed_fft.requested);
  assert(allowed_fft.eligible);
  assert(allowed_fft.scheduled);
  assert(allowed_fft.active);
}

void testJsonSerialization() {
  const auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  const auto report = cosmosim::workflows::buildRuntimeCapabilityReport(config);
  const std::string json =
      cosmosim::workflows::serializeRuntimeCapabilityReportJson(report);
  assert(json.find("\"schema_version\": 3") != std::string::npos);
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
  testPowerSpectrumCapabilitySchedulingTruth();
  testJsonSerialization();
  return 0;
}
