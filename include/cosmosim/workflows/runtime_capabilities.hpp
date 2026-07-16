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

struct RuntimeCapability {
  std::string name;
  RuntimeCapabilityStatus status = RuntimeCapabilityStatus::kUnsupported;
  std::string detail;
};

struct RuntimeCapabilityReport {
  unsigned int schema_version = 1;
  std::vector<RuntimeCapability> capabilities;

  [[nodiscard]] const RuntimeCapability& require(std::string_view name) const;
};

[[nodiscard]] std::string_view runtimeCapabilityStatusName(
    RuntimeCapabilityStatus status) noexcept;

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
