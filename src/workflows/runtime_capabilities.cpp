#include "cosmosim/workflows/runtime_capabilities.hpp"

#include <fstream>
#include <stdexcept>

#include "cosmosim/core/build_config.hpp"

namespace cosmosim::workflows {
namespace {

[[nodiscard]] std::string escapeJson(std::string_view value) {
  std::string escaped;
  escaped.reserve(value.size());
  for (const char ch : value) {
    switch (ch) {
      case '\\':
        escaped += "\\\\";
        break;
      case '"':
        escaped += "\\\"";
        break;
      case '\n':
        escaped += "\\n";
        break;
      case '\r':
        escaped += "\\r";
        break;
      case '\t':
        escaped += "\\t";
        break;
      default:
        escaped.push_back(ch);
        break;
    }
  }
  return escaped;
}

}  // namespace

const RuntimeCapability& RuntimeCapabilityReport::require(
    std::string_view name) const {
  for (const RuntimeCapability& capability : capabilities) {
    if (capability.name == name) {
      return capability;
    }
  }
  throw std::out_of_range("runtime capability is absent from report: " +
                          std::string(name));
}

std::string_view runtimeCapabilityStatusName(
    RuntimeCapabilityStatus status) noexcept {
  switch (status) {
    case RuntimeCapabilityStatus::kSupported:
      return "supported";
    case RuntimeCapabilityStatus::kProvisional:
      return "provisional";
    case RuntimeCapabilityStatus::kUnsupported:
      return "unsupported";
  }
  return "unsupported";
}

RuntimeCapabilityReport buildRuntimeCapabilityReport(
    const core::SimulationConfig& config) {
  RuntimeCapabilityReport report;
  report.capabilities = {
      {"fixed_global_timestep", RuntimeCapabilityStatus::kSupported,
       "Production ReferenceWorkflow KDK with hierarchical_max_rung=0."},
      {"adaptive_global_timestep", RuntimeCapabilityStatus::kUnsupported,
       "Criteria currently constrain scheduler bins; they do not yet resize the global base interval."},
      {"production_hierarchical_local_timestep",
       RuntimeCapabilityStatus::kUnsupported,
       "Per-element drift/kick epochs are not yet complete in the production KDK path."},
      {"canonical_external_ic_import", RuntimeCapabilityStatus::kProvisional,
       "Versioned IcManifest conversion/validation and a fail-closed single-file bridge exist; canonical conversion tooling and typed manifest-file config are not complete."},
      {"distributed_ic_import", RuntimeCapabilityStatus::kUnsupported,
       "Every rank currently enters the same single-file import path."},
      {"rank_remappable_restart", RuntimeCapabilityStatus::kUnsupported,
       "Restart schema v20 supports same-world-size rank-local continuation only."},
      {"asynchronous_output", RuntimeCapabilityStatus::kUnsupported,
       "Snapshot and checkpoint writes are synchronous at restart-safe boundaries."},
      {"production_diagnostics", RuntimeCapabilityStatus::kSupported,
       "Run-health and validated lightweight diagnostics follow the typed execution policy."},
      {"provisional_diagnostics",
       config.analysis.diagnostics_execution_policy ==
               core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional
           ? RuntimeCapabilityStatus::kProvisional
           : RuntimeCapabilityStatus::kUnsupported,
       "Heavy reference diagnostics require explicit all_including_provisional opt-in."},
  };
  return report;
}

void validateRequestedRuntimeCapabilities(
    const core::SimulationConfig& config,
    const RuntimeCapabilityReport& report) {
  if (config.numerics.hierarchical_max_rung != 0 &&
      report.require("production_hierarchical_local_timestep").status !=
          RuntimeCapabilityStatus::kSupported) {
    throw std::invalid_argument(
        "runtime capability production_hierarchical_local_timestep is unsupported: " +
        report.require("production_hierarchical_local_timestep").detail);
  }
}

std::string serializeRuntimeCapabilityReportJson(
    const RuntimeCapabilityReport& report) {
  std::string json = "{\n  \"schema_version\": " +
      std::to_string(report.schema_version) + ",\n  \"capabilities\": [\n";
  for (std::size_t index = 0; index < report.capabilities.size(); ++index) {
    const RuntimeCapability& capability = report.capabilities[index];
    json += "    {\"name\": \"" + escapeJson(capability.name) +
        "\", \"status\": \"" +
        std::string(runtimeCapabilityStatusName(capability.status)) +
        "\", \"detail\": \"" + escapeJson(capability.detail) + "\"}";
    json += index + 1U == report.capabilities.size() ? "\n" : ",\n";
  }
  json += "  ]\n}\n";
  return json;
}

void writeRuntimeCapabilityReportJson(
    const RuntimeCapabilityReport& report,
    const std::filesystem::path& output_path) {
  std::ofstream output(output_path, std::ios::binary | std::ios::trunc);
  if (!output) {
    throw std::runtime_error(
        "failed to open runtime capability report: " + output_path.string());
  }
  output << serializeRuntimeCapabilityReportJson(report);
  if (!output) {
    throw std::runtime_error(
        "failed to write runtime capability report: " + output_path.string());
  }
}

}  // namespace cosmosim::workflows
