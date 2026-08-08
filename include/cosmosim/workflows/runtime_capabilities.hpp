#pragma once

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/config.hpp"

namespace cosmosim::workflows {

enum class RuntimeCapabilityStatus : unsigned char {
  kSupported = 0,
  kProvisional = 1,
  kUnsupported = 2,
};

enum class RuntimeImplementationMaturity : unsigned char {
  kNotApplicable = 0,
  kProductionScalable = 1,
  kReference = 2,
};

enum class RuntimeScientificMaturity : unsigned char {
  kNotApplicable = 0,
  kProvisional = 1,
  kValidated = 2,
};

enum class RuntimeScalabilityClass : unsigned char {
  kNotApplicable = 0,
  kLight = 1,
  kModerate = 2,
  kScalableFft = 3,
  kHeavyReference = 4,
};

struct RuntimeCapability {
  std::string name;
  RuntimeCapabilityStatus status = RuntimeCapabilityStatus::kUnsupported;
  bool requested = false;
  bool compiled = false;
  bool dependency_available = false;
  bool runtime_available = false;
  bool eligible = false;
  bool scheduled = false;
  // Compatibility summary for existing consumers: active means scheduled by the
  // current configuration/policy, not merely requested or compiled.
  bool active = false;
  RuntimeImplementationMaturity implementation_maturity =
      RuntimeImplementationMaturity::kNotApplicable;
  RuntimeScientificMaturity scientific_maturity =
      RuntimeScientificMaturity::kNotApplicable;
  RuntimeScalabilityClass scalability = RuntimeScalabilityClass::kNotApplicable;
  std::string detail;
};

struct RuntimeCapabilityReport {
  unsigned int schema_version = 3;
  std::vector<RuntimeCapability> capabilities;

  [[nodiscard]] const RuntimeCapability& require(std::string_view name) const;
};

[[nodiscard]] std::string_view runtimeCapabilityStatusName(
    RuntimeCapabilityStatus status) noexcept;
[[nodiscard]] std::string_view runtimeImplementationMaturityName(
    RuntimeImplementationMaturity maturity) noexcept;
[[nodiscard]] std::string_view runtimeScientificMaturityName(
    RuntimeScientificMaturity maturity) noexcept;
[[nodiscard]] std::string_view runtimeScalabilityClassName(
    RuntimeScalabilityClass scalability) noexcept;

[[nodiscard]] RuntimeCapabilityReport buildRuntimeCapabilityReport(
    const core::SimulationConfig& config);

// Reject configurations that explicitly request a capability whose production
// path is not represented by this build/runtime contract.
void validateRequestedRuntimeCapabilities(
    const core::SimulationConfig& config,
    const RuntimeCapabilityReport& report);

[[nodiscard]] std::string serializeRuntimeCapabilityReportJson(
    const RuntimeCapabilityReport& report);

void writeRuntimeCapabilityReportJson(
    const RuntimeCapabilityReport& report,
    const std::filesystem::path& output_path);

}  // namespace cosmosim::workflows
