#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "io/internal/ic_reader_session.hpp"
#include "io/internal/ic_record_codec.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace cosmosim::io::file_set_internal {

inline constexpr std::size_t kParticleTypeCount = 6U;
inline constexpr std::uint32_t kInvalidIndex =
    std::numeric_limits<std::uint32_t>::max();
using ParticleRecord = internal::IcParticleRecord;

struct Inspection {
  IcManifest manifest;
  std::vector<IcSchemaSummary> schemas;
  IcImportCounters counters;
};

struct SourceFileInspection {
  std::filesystem::path path;
  IcSchemaSummary schema;
  std::uint64_t source_size_bytes = 0U;
  std::string source_sha256;
  std::string original_header_attributes;
  std::vector<IcFieldManifest> fields;
  std::vector<IcMissingFieldContract> missing_field_contracts;
  std::vector<std::string> defaulted_fields;
  std::vector<std::string> dropped_fields;
  std::vector<std::string> preserved_auxiliary_fields;
  std::vector<std::string> warnings;
  bool canonical_manifest_verified = false;
  std::string canonical_manifest_sha256;
  IcImportCounters counters;
};

[[nodiscard]] bool nearlyEqual(double lhs, double rhs);
[[nodiscard]] bool missingFieldContractsEqual(
    const std::vector<IcMissingFieldContract>& lhs,
    const std::vector<IcMissingFieldContract>& rhs);
[[nodiscard]] IcSpeciesPolicy mapConfiguredPolicy(
    core::InitialConditionSpeciesPolicy policy,
    std::size_t type_index);
[[nodiscard]] std::uint32_t speciesTag(IcSpeciesPolicy policy);

void checkedCounterAdd(
    std::uint64_t& destination,
    std::uint64_t value,
    std::string_view counter_name);

template <typename T>
[[nodiscard]] std::uint64_t vectorCapacityBytes(
    const std::vector<T>& values) {
  if (values.capacity() >
      std::numeric_limits<std::uint64_t>::max() / sizeof(T)) {
    throw std::overflow_error("IC vector capacity byte count overflow");
  }
  return static_cast<std::uint64_t>(values.capacity()) * sizeof(T);
}

#if COSMOSIM_ENABLE_HDF5
void readAttributeNonnegativeU32(
    hid_t group,
    const char* name,
    std::uint32_t& value);
[[nodiscard]] std::vector<std::filesystem::path> discoverFiles(
    const std::filesystem::path& requested,
    std::uint32_t count);
void validateCrossFileSchema(const IcManifest& manifest);
[[nodiscard]] SourceFileInspection inspectOneSourceFile(
    const std::filesystem::path& path,
    std::uint32_t file_index,
    IcDialect dialect,
    const std::array<IcSpeciesPolicy, kParticleTypeCount>& species_policy,
    const core::SimulationConfig& config,
    const IcImportOptions& options,
    bool has_authoritative_manifest);
[[nodiscard]] Inspection inspectFileSet(
    const std::filesystem::path& requested,
    const core::SimulationConfig& config,
    const IcImportOptions& options);
[[nodiscard]] std::vector<ParticleRecord> readRecordChunk(
    internal::IcReaderSession& session,
    const Inspection& inspection,
    std::size_t file_index,
    std::size_t type_index,
    std::size_t start,
    std::size_t count,
    const core::SimulationConfig& config,
    const IcImportOptions& options,
    IcImportCounters& counters);
void validateRecordScientificState(
    const ParticleRecord& record,
    IcSpeciesPolicy policy,
    std::uint64_t source_record_index);
void appendRecords(
    core::SimulationState& state,
    const std::vector<ParticleRecord>& records,
    std::uint32_t owner_rank);
void finalizeImportedState(
    core::SimulationState& state,
    const IcManifest& manifest,
    const core::SimulationConfig& config);
[[nodiscard]] double convertedBoxSizeCode(
    const IcManifest& manifest,
    const core::SimulationConfig& config);
void validateRuntimeCosmology(
    const IcManifest& manifest,
    const core::SimulationConfig& config);
void validateSerialCountsAndIds(
    const core::SimulationState& state,
    const IcManifest& manifest);
#endif

}  // namespace cosmosim::io::file_set_internal
