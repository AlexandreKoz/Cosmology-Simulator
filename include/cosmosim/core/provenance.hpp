#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace cosmosim::core {

// Canonical provenance payload written alongside run outputs.
struct ProvenanceRecord {
  std::string schema_version = "provenance_v1";
  std::string git_sha = "unknown";
  std::string compiler_id = "unknown";
  std::string compiler_version = "unknown";
  std::string build_preset = "manual";
  std::string enabled_features;
  std::string config_hash_hex;
  std::string timestamp_utc;
  std::string hardware_summary;
  int author_rank = 0;
};

[[nodiscard]] std::string collectCompilerId();
[[nodiscard]] std::string collectCompilerVersion();
[[nodiscard]] std::string collectHardwareSummary();
[[nodiscard]] std::string utcTimestampNowIso8601();

// Deterministic hashing for normalized configuration text.
[[nodiscard]] std::uint64_t stableConfigHash(const std::string& normalized_config_text);
[[nodiscard]] std::string stableConfigHashHex(const std::string& normalized_config_text);

[[nodiscard]] ProvenanceRecord makeProvenanceRecord(
    const std::string& config_hash_hex,
    const std::string& git_sha,
    int rank = 0);

void writeProvenanceRecord(
    const ProvenanceRecord& record,
    const std::filesystem::path& run_directory,
    const std::string& file_name = "provenance.meta.txt");

[[nodiscard]] ProvenanceRecord readProvenanceRecord(
    const std::filesystem::path& run_directory,
    const std::string& file_name = "provenance.meta.txt");

[[nodiscard]] std::string serializeProvenanceRecord(const ProvenanceRecord& record);
[[nodiscard]] ProvenanceRecord deserializeProvenanceRecord(std::string_view text);

}  // namespace cosmosim::core
