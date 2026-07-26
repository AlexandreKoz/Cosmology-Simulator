#include "cosmosim/io/ic_reader.hpp"
#include "io/internal/ic_file_set_common.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "io/internal/ic_byte_codec.hpp"
#include "io/internal/ic_canonical_bundle.hpp"
#include "io/internal/ic_conversion_catalog.hpp"
#include "io/internal/ic_distributed_audit.hpp"
#include "io/internal/ic_failure_protocol.hpp"
#include "io/internal/ic_hdf5_handle.hpp"
#include "io/internal/ic_mpi_collectives.hpp"
#include "io/internal/ic_reader_session.hpp"
#include "io/internal/ic_record_codec.hpp"
#include "io/internal/ic_stream_ingestion.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif
#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::io {

using namespace file_set_internal;
#if COSMOSIM_ENABLE_HDF5
using internal::Hdf5Handle;
using internal::IcReaderSession;
#endif
#if COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI
using distributed_audit_internal::exactDistributedChunkReconciliation;
using distributed_audit_internal::exactDistributedIdAudit;
using distributed_audit_internal::nestedByteCapacity;
using distributed_audit_internal::ownerForX;
using distributed_audit_internal::validateChunkCoverage;
using distributed_audit_internal::validateDistributedTotals;
using failure_protocol_internal::alltoallBytes;
using failure_protocol_internal::broadcastRootString;
using failure_protocol_internal::CollectivePhaseCounterScope;
using failure_protocol_internal::injectIcTestFault;
using failure_protocol_internal::mutateIcTestRoute;
using failure_protocol_internal::runCollectivePhase;
using failure_protocol_internal::runCollectivePhaseVoid;
using mpi_collective_internal::MpiCollectiveCallCounts;
using mpi_collective_internal::MpiCollectiveCounterScope;
using mpi_collective_internal::RoutingMpiCollectiveScope;
using mpi_collective_internal::mpiAllreduce;
using mpi_collective_internal::mpiGather;
using mpi_collective_internal::mpiGatherv;
#endif

namespace {
#if COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI

void populateSchemasFromManifest(Inspection& inspection) {
  inspection.schemas.clear();
  inspection.schemas.reserve(inspection.manifest.num_part_this_file.size());
  for (const auto& per_file_counts : inspection.manifest.num_part_this_file) {
    IcSchemaSummary schema;
    schema.count_by_type = per_file_counts;
    schema.total_count_by_type = inspection.manifest.num_part_total;
    schema.total_count_high_word = inspection.manifest.num_part_total_high_word;
    schema.mass_table = inspection.manifest.mass_table;
    schema.num_files_per_snapshot =
        inspection.manifest.num_files_per_snapshot;
    schema.box_size = inspection.manifest.box_size;
    schema.scale_factor = inspection.manifest.scale_factor;
    schema.redshift = inspection.manifest.redshift;
    schema.omega_matter = inspection.manifest.omega_matter;
    schema.omega_lambda = inspection.manifest.omega_lambda;
    schema.hubble_param = inspection.manifest.hubble_param;
    inspection.schemas.push_back(schema);
  }
}

[[nodiscard]] IcManifest makeTransferManifest(
    const SourceFileInspection& source,
    IcDialect dialect,
    const std::array<IcSpeciesPolicy, kParticleTypeCount>& species_policy) {
  IcManifest manifest;
  manifest.dialect = dialect;
  manifest.dialect_version = "1";
  manifest.converter_version = "chui_distributed_inspector_fragment_v3";
  manifest.source_files = {source.path};
  manifest.source_provenance_ids = {"sha256:" + source.source_sha256};
  manifest.source_file_sizes_bytes = {source.source_size_bytes};
  manifest.source_sha256 = {source.source_sha256};
  manifest.original_header_attributes = {source.original_header_attributes};
  manifest.num_files_per_snapshot = 1U;
  manifest.num_part_this_file = {source.schema.count_by_type};
  manifest.num_part_total = source.schema.count_by_type;
  for (std::size_t type = 0; type < kParticleTypeCount; ++type) {
    manifest.num_part_total_high_word[type] =
        manifest.num_part_total[type] >> 32U;
  }
  manifest.mass_table = source.schema.mass_table;
  manifest.box_size = source.schema.box_size;
  manifest.scale_factor = source.schema.scale_factor;
  manifest.redshift = source.schema.redshift;
  manifest.omega_matter = source.schema.omega_matter;
  manifest.omega_lambda = source.schema.omega_lambda;
  manifest.hubble_param = source.schema.hubble_param;
  manifest.species_policy = species_policy;
  manifest.fields = source.fields;
  for (IcFieldManifest& field : manifest.fields) {
    field.source_file_index = 0U;
    if (field.disposition == IcFieldDisposition::kConverted) {
      manifest.converted_fields.push_back(field.dataset_path);
      if (std::find(
              manifest.conversion_equations.begin(),
              manifest.conversion_equations.end(), field.conversion_equation) ==
          manifest.conversion_equations.end()) {
        manifest.conversion_equations.push_back(field.conversion_equation);
      }
    }
  }
  manifest.missing_field_contracts = source.missing_field_contracts;
  for (IcMissingFieldContract& contract : manifest.missing_field_contracts) {
    contract.source_file_index = 0U;
  }
  manifest.defaulted_fields = source.defaulted_fields;
  manifest.dropped_fields = source.dropped_fields;
  manifest.preserved_auxiliary_fields = source.preserved_auxiliary_fields;
  manifest.warnings = source.warnings;
  manifest.canonical_source_manifest_verified =
      source.canonical_manifest_verified;
  manifest.canonical_source_manifest_sha256 =
      source.canonical_manifest_sha256;
  validateIcManifest(manifest);
  return manifest;
}

void appendSchemaSummary(
    std::vector<std::uint8_t>& output,
    const IcSchemaSummary& schema) {
  for (std::uint64_t value : schema.count_by_type) {
    internal::appendLe64(output, value);
  }
  for (std::uint64_t value : schema.total_count_by_type) {
    internal::appendLe64(output, value);
  }
  for (std::uint32_t value : schema.total_count_high_word) {
    internal::appendLe32(output, value);
  }
  for (double value : schema.mass_table) {
    internal::appendDouble(output, value);
  }
  internal::appendLe32(output, schema.num_files_per_snapshot);
  internal::appendDouble(output, schema.box_size);
  internal::appendDouble(output, schema.scale_factor);
  internal::appendDouble(output, schema.redshift);
  internal::appendDouble(output, schema.omega_matter);
  internal::appendDouble(output, schema.omega_lambda);
  internal::appendDouble(output, schema.hubble_param);
}

[[nodiscard]] IcSchemaSummary readSchemaSummary(
    std::span<const std::uint8_t> input,
    std::size_t& offset) {
  IcSchemaSummary schema;
  for (std::uint64_t& value : schema.count_by_type) {
    value = internal::readLe64(input, offset);
  }
  for (std::uint64_t& value : schema.total_count_by_type) {
    value = internal::readLe64(input, offset);
  }
  for (std::uint32_t& value : schema.total_count_high_word) {
    value = internal::readLe32(input, offset);
  }
  for (double& value : schema.mass_table) {
    value = internal::readDouble(input, offset);
  }
  schema.num_files_per_snapshot = internal::readLe32(input, offset);
  schema.box_size = internal::readDouble(input, offset);
  schema.scale_factor = internal::readDouble(input, offset);
  schema.redshift = internal::readDouble(input, offset);
  schema.omega_matter = internal::readDouble(input, offset);
  schema.omega_lambda = internal::readDouble(input, offset);
  schema.hubble_param = internal::readDouble(input, offset);
  return schema;
}

void appendLengthPrefixedString(
    std::vector<std::uint8_t>& output,
    std::string_view value) {
  internal::appendLe64(output, value.size());
  output.insert(output.end(), value.begin(), value.end());
}

[[nodiscard]] std::string readLengthPrefixedString(
    std::span<const std::uint8_t> input,
    std::size_t& offset) {
  const std::uint64_t length = internal::readLe64(input, offset);
  if (length > input.size() - offset) {
    throw std::runtime_error("truncated distributed IC metadata string");
  }
  std::string value(
      reinterpret_cast<const char*>(input.data() + offset),
      static_cast<std::size_t>(length));
  offset += static_cast<std::size_t>(length);
  return value;
}

[[nodiscard]] std::vector<std::uint8_t> gatherRootBytes(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::uint8_t>& local_bytes) {
  const int local_count = runCollectivePhase<int>(
      mpi_context, "IC metadata gather local-count validation", [&]() {
        if (local_bytes.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max())) {
          throw std::overflow_error(
              "IC metadata fragment exceeds MPI int count");
        }
        return static_cast<int>(local_bytes.size());
      });
  std::vector<int> counts = runCollectivePhase<std::vector<int>>(
      mpi_context, "IC metadata gather count allocation", [&]() {
        if (!mpi_context.isRoot()) {
          return std::vector<int>{};
        }
        return std::vector<int>(
            static_cast<std::size_t>(mpi_context.worldSize()));
      });
  mpiGather(
      &local_count, 1, MPI_INT,
      mpi_context.isRoot() ? counts.data() : nullptr, 1, MPI_INT, 0,
      MPI_COMM_WORLD);

  struct RootGatherLayout {
    std::vector<int> displacements;
    std::vector<std::uint8_t> bytes;
  };
  RootGatherLayout layout = runCollectivePhase<RootGatherLayout>(
      mpi_context, "IC metadata gather root-layout preparation", [&]() {
        RootGatherLayout prepared;
        if (!mpi_context.isRoot()) {
          return prepared;
        }
        prepared.displacements.resize(counts.size());
        std::uint64_t total = 0U;
        for (std::size_t rank = 0; rank < counts.size(); ++rank) {
          if (counts[rank] < 0 ||
              total >
                  static_cast<std::uint64_t>(std::numeric_limits<int>::max()) -
                      static_cast<std::uint64_t>(counts[rank])) {
            throw std::overflow_error(
                "IC metadata gather exceeds MPI int displacement range");
          }
          prepared.displacements[rank] = static_cast<int>(total);
          total += static_cast<std::uint64_t>(counts[rank]);
        }
        prepared.bytes.resize(static_cast<std::size_t>(total));
        return prepared;
      });

  mpiGatherv(
      local_bytes.empty() ? nullptr : local_bytes.data(), local_count, MPI_BYTE,
      mpi_context.isRoot() && !layout.bytes.empty() ? layout.bytes.data() : nullptr,
      mpi_context.isRoot() ? counts.data() : nullptr,
      mpi_context.isRoot() ? layout.displacements.data() : nullptr,
      MPI_BYTE, 0, MPI_COMM_WORLD);
  return std::move(layout.bytes);
}

[[nodiscard]] std::string encodeDiscoveredPaths(
    const std::vector<std::filesystem::path>& paths) {
  std::vector<std::uint8_t> bytes;
  internal::appendLe64(bytes, paths.size());
  for (const auto& path : paths) {
    appendLengthPrefixedString(bytes, path.lexically_normal().generic_string());
  }
  return std::string(
      reinterpret_cast<const char*>(bytes.data()), bytes.size());
}

[[nodiscard]] std::vector<std::filesystem::path> decodeDiscoveredPaths(
    std::string_view encoded) {
  const auto* begin = reinterpret_cast<const std::uint8_t*>(encoded.data());
  std::span<const std::uint8_t> bytes(begin, encoded.size());
  std::size_t offset = 0U;
  const std::uint64_t count = internal::readLe64(bytes, offset);
  if (count > std::numeric_limits<std::size_t>::max()) {
    throw std::overflow_error("discovered IC file count exceeds size_t");
  }
  std::vector<std::filesystem::path> paths;
  paths.reserve(static_cast<std::size_t>(count));
  for (std::uint64_t index = 0U; index < count; ++index) {
    paths.emplace_back(readLengthPrefixedString(bytes, offset));
  }
  if (offset != bytes.size()) {
    throw std::runtime_error("distributed IC path metadata has trailing bytes");
  }
  return paths;
}

void appendManifestLists(IcManifest& destination, IcManifest&& source) {
  destination.fields.insert(
      destination.fields.end(),
      std::make_move_iterator(source.fields.begin()),
      std::make_move_iterator(source.fields.end()));
  destination.missing_field_contracts.insert(
      destination.missing_field_contracts.end(),
      std::make_move_iterator(source.missing_field_contracts.begin()),
      std::make_move_iterator(source.missing_field_contracts.end()));
  destination.defaulted_fields.insert(
      destination.defaulted_fields.end(),
      std::make_move_iterator(source.defaulted_fields.begin()),
      std::make_move_iterator(source.defaulted_fields.end()));
  destination.dropped_fields.insert(
      destination.dropped_fields.end(),
      std::make_move_iterator(source.dropped_fields.begin()),
      std::make_move_iterator(source.dropped_fields.end()));
  destination.preserved_auxiliary_fields.insert(
      destination.preserved_auxiliary_fields.end(),
      std::make_move_iterator(source.preserved_auxiliary_fields.begin()),
      std::make_move_iterator(source.preserved_auxiliary_fields.end()));
  destination.warnings.insert(
      destination.warnings.end(),
      std::make_move_iterator(source.warnings.begin()),
      std::make_move_iterator(source.warnings.end()));
}

[[nodiscard]] Inspection inspectFileSetDistributed(
    const std::filesystem::path& requested,
    const core::SimulationConfig& config,
    const parallel::MpiContext& mpi_context,
    const IcImportOptions& options) {
  const bool has_authoritative_manifest = options.manifest != nullptr;
  const IcDialect dialect = has_authoritative_manifest
      ? options.manifest->dialect
      : (config.mode.ic_convention ==
                 core::InitialConditionConvention::kChuiCanonicalV1
             ? IcDialect::kChuiCanonicalV1
             : IcDialect::kGadgetArepoBridgeV1);
  std::array<IcSpeciesPolicy, kParticleTypeCount> species_policy{
      IcSpeciesPolicy::kGas,
      IcSpeciesPolicy::kDarkMatter,
      mapConfiguredPolicy(config.mode.ic_part_type2_policy, 2U),
      mapConfiguredPolicy(config.mode.ic_part_type3_policy, 3U),
      IcSpeciesPolicy::kStar,
      IcSpeciesPolicy::kBlackHole};
  if (options.manifest != nullptr) {
    species_policy = options.manifest->species_policy;
  }

  std::string root_paths = runCollectivePhase<std::string>(
      mpi_context, "IC distributed file discovery", [&]() {
        if (!mpi_context.isRoot()) {
          return std::string{};
        }
        const std::filesystem::path first_source =
            has_authoritative_manifest &&
                !options.manifest->source_files.empty()
            ? options.manifest->source_files.front()
            : requested;
        Hdf5Handle file(
            H5Fopen(first_source.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        if (!file.valid()) {
          throw std::runtime_error(
              "failed to open first IC source for distributed discovery: " +
              first_source.string());
        }
        Hdf5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
        if (!header.valid()) {
          throw std::runtime_error("IC source is missing /Header");
        }
        std::uint32_t file_count = 0U;
        readAttributeNonnegativeU32(
            header.get(), "NumFilesPerSnapshot", file_count);
        std::vector<std::filesystem::path> paths;
        if (has_authoritative_manifest &&
            !options.manifest->source_files.empty()) {
          paths = options.manifest->source_files;
        } else {
          paths = discoverFiles(requested, file_count);
        }
        if (paths.size() != file_count) {
          throw std::runtime_error(
              "distributed IC discovery disagrees with "
              "NumFilesPerSnapshot");
        }
        return encodeDiscoveredPaths(paths);
      });
  root_paths = broadcastRootString(mpi_context, std::move(root_paths));
  const std::vector<std::filesystem::path> paths =
      runCollectivePhase<std::vector<std::filesystem::path>>(
          mpi_context, "IC distributed path metadata decoding", [&]() {
            return decodeDiscoveredPaths(root_paths);
          });

  struct LocalInspectionCollection {
    std::vector<std::pair<std::uint32_t, SourceFileInspection>> sources;
    IcImportCounters counters;
  };
  LocalInspectionCollection local =
      runCollectivePhase<LocalInspectionCollection>(
          mpi_context, "IC assigned-file inspection and hashing", [&]() {
            LocalInspectionCollection collection;
            for (std::size_t file_index = 0; file_index < paths.size();
                 ++file_index) {
              if (static_cast<int>(file_index %
                                   static_cast<std::size_t>(
                                       mpi_context.worldSize())) !=
                  mpi_context.worldRank()) {
                continue;
              }
              SourceFileInspection source = inspectOneSourceFile(
                  paths[file_index], static_cast<std::uint32_t>(file_index),
                  dialect, species_policy, config,
                  options, has_authoritative_manifest);
              checkedCounterAdd(
                  collection.counters.metadata_bytes_read,
                  source.counters.metadata_bytes_read,
                  "metadata_bytes_read");
              checkedCounterAdd(
                  collection.counters.hash_bytes_read,
                  source.counters.hash_bytes_read, "hash_bytes_read");
              checkedCounterAdd(
                  collection.counters.full_file_hash_pass_count,
                  source.counters.full_file_hash_pass_count,
                  "full_file_hash_pass_count");
              collection.sources.emplace_back(
                  static_cast<std::uint32_t>(file_index), std::move(source));
            }
            collection.counters.files_assigned = collection.sources.size();
            collection.counters.bytes_read =
                collection.counters.metadata_bytes_read +
                collection.counters.hash_bytes_read;
            return collection;
          });

  std::vector<std::uint8_t> local_metadata =
      runCollectivePhase<std::vector<std::uint8_t>>(
          mpi_context, "IC manifest-fragment serialization", [&]() {
            std::vector<std::uint8_t> bytes;
            internal::appendLe64(bytes, local.sources.size());
            for (const auto& [file_index, source] : local.sources) {
              internal::appendLe32(bytes, file_index);
              appendSchemaSummary(bytes, source.schema);
              const std::string json = serializeIcManifestJson(
                  makeTransferManifest(source, dialect, species_policy));
              appendLengthPrefixedString(bytes, json);
            }
            return bytes;
          });
  runCollectivePhaseVoid(
      mpi_context, "IC manifest-fragment accounting", [&]() {
        checkedCounterAdd(
            local.counters.manifest_metadata_bytes_communicated,
            local_metadata.size(),
            "manifest_metadata_bytes_communicated");
      });
  const std::vector<std::uint8_t> gathered =
      gatherRootBytes(mpi_context, local_metadata);
  runCollectivePhaseVoid(
      mpi_context, "IC gathered-manifest accounting", [&]() {
        if (mpi_context.isRoot()) {
          checkedCounterAdd(
              local.counters.manifest_metadata_bytes_communicated,
              gathered.size(), "manifest_metadata_bytes_communicated");
        }
      });

  Inspection root_inspection = runCollectivePhase<Inspection>(
      mpi_context, "IC distributed manifest assembly", [&]() {
        Inspection assembled;
        assembled.counters = local.counters;
        if (!mpi_context.isRoot()) {
          return assembled;
        }
        std::vector<std::optional<std::pair<IcSchemaSummary, IcManifest>>>
            fragments(paths.size());
        std::size_t offset = 0U;
        for (int rank = 0; rank < mpi_context.worldSize(); ++rank) {
          if (offset >= gathered.size()) {
            throw std::runtime_error(
                "distributed IC metadata is missing a rank fragment");
          }
          const std::uint64_t entry_count = internal::readLe64(gathered, offset);
          for (std::uint64_t entry = 0U; entry < entry_count; ++entry) {
            const std::uint32_t file_index = internal::readLe32(gathered, offset);
            if (file_index >= fragments.size() || fragments[file_index]) {
              throw std::runtime_error(
                  "distributed IC file fragment is duplicated or out of range");
            }
            IcSchemaSummary schema = readSchemaSummary(gathered, offset);
            IcManifest fragment = deserializeIcManifestJson(
                readLengthPrefixedString(gathered, offset));
            fragments[file_index].emplace(
                std::move(schema), std::move(fragment));
          }
        }
        if (offset != gathered.size()) {
          throw std::runtime_error(
              "distributed IC metadata has trailing bytes");
        }
        if (std::any_of(
                fragments.begin(), fragments.end(),
                [](const auto& fragment) { return !fragment.has_value(); })) {
          throw std::runtime_error(
              "distributed IC inspection did not cover every source file");
        }

        const IcSchemaSummary& baseline = fragments.front()->first;
        IcManifest& manifest = assembled.manifest;
        manifest.dialect = dialect;
        manifest.dialect_version = "1";
        manifest.converter_version = "chui_distributed_inspector_v3";
        manifest.num_files_per_snapshot =
            static_cast<std::uint32_t>(paths.size());
        manifest.source_files = paths;
        manifest.species_policy = species_policy;
        manifest.num_part_total = baseline.total_count_by_type;
        manifest.num_part_total_high_word = baseline.total_count_high_word;
        manifest.mass_table = baseline.mass_table;
        manifest.box_size = baseline.box_size;
        manifest.scale_factor = baseline.scale_factor;
        manifest.redshift = baseline.redshift;
        manifest.omega_matter = baseline.omega_matter;
        manifest.omega_lambda = baseline.omega_lambda;
        manifest.hubble_param = baseline.hubble_param;

        std::array<std::uint64_t, kParticleTypeCount> summed{};
        for (std::size_t file_index = 0; file_index < fragments.size();
             ++file_index) {
          IcSchemaSummary schema = fragments[file_index]->first;
          IcManifest fragment = std::move(fragments[file_index]->second);
          if (fragment.canonical_source_manifest_verified) {
            if (manifest.canonical_source_manifest_verified &&
                manifest.canonical_source_manifest_sha256 !=
                    fragment.canonical_source_manifest_sha256) {
              throw std::runtime_error(
                  "distributed canonical IC fragments disagree on manifest binding");
            }
            manifest.canonical_source_manifest_verified = true;
            manifest.canonical_source_manifest_sha256 =
                fragment.canonical_source_manifest_sha256;
          }
          if (schema.num_files_per_snapshot != paths.size() ||
              schema.total_count_by_type != baseline.total_count_by_type ||
              schema.total_count_high_word != baseline.total_count_high_word ||
              schema.mass_table != baseline.mass_table ||
              !nearlyEqual(schema.box_size, baseline.box_size) ||
              !nearlyEqual(schema.scale_factor, baseline.scale_factor) ||
              !nearlyEqual(schema.redshift, baseline.redshift) ||
              !nearlyEqual(schema.omega_matter, baseline.omega_matter) ||
              !nearlyEqual(schema.omega_lambda, baseline.omega_lambda) ||
              !nearlyEqual(schema.hubble_param, baseline.hubble_param)) {
            throw std::runtime_error(
                "distributed IC source files disagree on scientific header "
                "or totals");
          }
          for (std::size_t type = 0; type < kParticleTypeCount; ++type) {
            if (summed[type] >
                std::numeric_limits<std::uint64_t>::max() -
                    schema.count_by_type[type]) {
              throw std::overflow_error(
                  "distributed IC particle-count sum overflow");
            }
            summed[type] += schema.count_by_type[type];
          }
          assembled.schemas.push_back(schema);
          manifest.num_part_this_file.push_back(schema.count_by_type);
          manifest.source_file_sizes_bytes.push_back(
              fragment.source_file_sizes_bytes.front());
          manifest.source_sha256.push_back(fragment.source_sha256.front());
          manifest.source_provenance_ids.push_back(
              fragment.source_provenance_ids.front());
          manifest.original_header_attributes.push_back(
              fragment.original_header_attributes.front());
          for (IcFieldManifest& field : fragment.fields) {
            field.source_file_index = static_cast<std::uint32_t>(file_index);
          }
          for (IcMissingFieldContract& contract :
               fragment.missing_field_contracts) {
            contract.source_file_index =
                static_cast<std::uint32_t>(file_index);
          }
          appendManifestLists(manifest, std::move(fragment));
        }
        if (summed != manifest.num_part_total) {
          throw std::runtime_error(
              "distributed IC per-file counts do not equal declared totals");
        }
        validateCrossFileSchema(manifest);
        for (const IcFieldManifest& field : manifest.fields) {
          if (field.disposition == IcFieldDisposition::kConverted) {
            manifest.converted_fields.push_back(field.dataset_path);
            if (std::find(
                    manifest.conversion_equations.begin(),
                    manifest.conversion_equations.end(),
                    field.conversion_equation) ==
                manifest.conversion_equations.end()) {
              manifest.conversion_equations.push_back(
                  field.conversion_equation);
            }
          }
        }
        if (has_authoritative_manifest) {
          const IcManifest& supplied = *options.manifest;
          if (supplied.source_files != manifest.source_files ||
              supplied.source_sha256 != manifest.source_sha256 ||
              supplied.source_file_sizes_bytes !=
                  manifest.source_file_sizes_bytes ||
              supplied.num_part_this_file != manifest.num_part_this_file ||
              supplied.num_part_total != manifest.num_part_total ||
              supplied.species_policy != manifest.species_policy ||
              !missingFieldContractsEqual(
                  supplied.missing_field_contracts,
                  manifest.missing_field_contracts) ||
              supplied.fields.size() != manifest.fields.size()) {
            throw std::runtime_error(
                "supplied distributed IC manifest does not match inspected "
                "source provenance or schema");
          }
          for (const IcFieldManifest& expected : supplied.fields) {
            const auto observed = std::find_if(
                manifest.fields.begin(), manifest.fields.end(),
                [&](const IcFieldManifest& actual) {
                  return actual.source_file_index ==
                             expected.source_file_index &&
                      actual.dataset_path == expected.dataset_path &&
                      actual.selected_alias == expected.selected_alias;
                });
            if (observed == manifest.fields.end() ||
                observed->scalar_type != expected.scalar_type ||
                observed->scalar_class != expected.scalar_class ||
                observed->byte_width != expected.byte_width ||
                observed->is_signed != expected.is_signed ||
                observed->byte_order != expected.byte_order ||
                observed->rank != expected.rank ||
                observed->dimensions != expected.dimensions ||
                observed->disposition != expected.disposition) {
              throw std::runtime_error(
                  "supplied distributed IC manifest field does not match " +
                  expected.dataset_path);
            }
          }
          manifest = supplied;
        }
        validateIcManifest(manifest);
        return assembled;
      });

  std::string manifest_json = runCollectivePhase<std::string>(
      mpi_context, "IC root manifest serialization", [&]() {
        return mpi_context.isRoot()
            ? serializeIcManifestJson(root_inspection.manifest)
            : std::string{};
      });
  manifest_json = broadcastRootString(mpi_context, std::move(manifest_json));
  runCollectivePhaseVoid(
      mpi_context, "IC broadcast-manifest accounting", [&]() {
        checkedCounterAdd(
            local.counters.manifest_metadata_bytes_communicated,
            manifest_json.size(), "manifest_metadata_bytes_communicated");
      });
  Inspection inspection = runCollectivePhase<Inspection>(
      mpi_context, "IC distributed manifest decode", [&]() {
        Inspection decoded;
        decoded.counters = local.counters;
        decoded.manifest = deserializeIcManifestJson(manifest_json);
        populateSchemasFromManifest(decoded);
        return decoded;
      });
  return inspection;
}


#endif  // COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI

}  // namespace

IcReadResult readDistributedGadgetArepoHdf5Ic(
    const std::filesystem::path& ic_path,
    const core::SimulationConfig& config,
    const parallel::MpiContext& mpi_context,
    const IcImportOptions& options) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(ic_path);static_cast<void>(config);static_cast<void>(mpi_context);static_cast<void>(options);throw std::runtime_error("COSMOSIM_ENABLE_HDF5=OFF: distributed IC reader unavailable");
#elif !COSMOSIM_ENABLE_MPI
  if(mpi_context.worldSize()!=1)throw std::runtime_error("distributed IC reader requires an MPI-enabled build for world_size > 1");return readGadgetArepoHdf5Ic(ic_path,config,options);
#else
  if (!mpi_context.isEnabled() || mpi_context.worldSize() == 1) {
    return readGadgetArepoHdf5Ic(ic_path, config, options);
  }
  std::uint64_t collective_phase_count = 0U;
  CollectivePhaseCounterScope collective_counter_scope(collective_phase_count);
  MpiCollectiveCallCounts mpi_collective_counts;
  MpiCollectiveCounterScope mpi_collective_counter_scope(mpi_collective_counts);
  runCollectivePhaseVoid(
      mpi_context, "IC distributed configuration validation", [&]() {
        if (options.chunk_particle_count == 0U) {
          throw std::invalid_argument("chunk_particle_count must be positive");
        }
        if (config.mode.ic_staging_particle_count == 0U) {
          throw std::invalid_argument(
              "distributed IC staging particle count must be positive");
        }
        if (options.manifest != nullptr) {
          validateIcManifest(*options.manifest);
        }
      });

  Inspection inspection =
      inspectFileSetDistributed(ic_path, config, mpi_context, options);
  if (options.validate_runtime_cosmology) {
    runCollectivePhaseVoid(
        mpi_context, "IC distributed runtime-cosmology validation", [&]() {
          validateRuntimeCosmology(inspection.manifest, config);
        });
  }

  const std::string local_manifest_json =
      runCollectivePhase<std::string>(
          mpi_context, "IC local manifest serialization", [&]() {
            return serializeIcManifestJson(inspection.manifest);
          });
  const std::string local_manifest_digest =
      runCollectivePhase<std::string>(
          mpi_context, "IC manifest digest", [&]() {
            injectIcTestFault(mpi_context, "manifest_digest");
            return icSha256Hex(std::string_view(local_manifest_json));
          });
  std::string root_manifest_digest =
      mpi_context.isRoot() ? local_manifest_digest : std::string{};
  root_manifest_digest =
      broadcastRootString(mpi_context, std::move(root_manifest_digest));
  const int local_digest_mismatch =
      local_manifest_digest == root_manifest_digest ? 0 : 1;
  int any_digest_mismatch = 0;
  mpiAllreduce(
      &local_digest_mismatch, &any_digest_mismatch, 1, MPI_INT, MPI_MAX,
      MPI_COMM_WORLD);
  if (any_digest_mismatch != 0) {
    throw std::runtime_error(
        "IC manifest SHA-256 digest is inconsistent across ranks");
  }

  IcReadResult result;
  result.report.counters = inspection.counters;
  result.report.manifest = inspection.manifest;
  result.report.already_partitioned = true;
  result.report.schema.count_by_type =
      inspection.manifest.num_part_this_file.front();
  result.report.schema.total_count_by_type =
      inspection.manifest.num_part_total;
  result.report.schema.total_count_high_word =
      inspection.manifest.num_part_total_high_word;
  result.report.schema.mass_table = inspection.manifest.mass_table;
  result.report.schema.num_files_per_snapshot =
      inspection.manifest.num_files_per_snapshot;
  result.report.schema.box_size = inspection.manifest.box_size;
  result.report.schema.scale_factor = inspection.manifest.scale_factor;
  result.report.schema.redshift = inspection.manifest.redshift;
  result.report.schema.omega_matter = inspection.manifest.omega_matter;
  result.report.schema.omega_lambda = inspection.manifest.omega_lambda;
  result.report.schema.hubble_param = inspection.manifest.hubble_param;
  result.report.defaulted_fields = inspection.manifest.defaulted_fields;
  for (const auto& value : result.report.defaulted_fields) {
    const auto equals = value.find('=');
    result.report.missing_optional_fields.push_back(
        value.substr(0, equals));
  }
  result.report.unsupported_fields = inspection.manifest.dropped_fields;
  result.report.manifest_verified =
      inspection.manifest.canonical_source_manifest_verified;
  result.report.verified_manifest_sha256 =
      inspection.manifest.canonical_source_manifest_sha256;
  result.report.provenance_authority = options.manifest != nullptr
      ? "supplied_manifest_v1"
      : (inspection.manifest.canonical_source_manifest_verified
            ? "canonical_embedded_manifest_v2"
            : "runtime_config_and_inspected_source");

  struct RoutingConfiguration {
    double box_size = 0.0;
    std::size_t chunk_particle_count = 0U;
    std::size_t batch_particle_count = 0U;
  };
  const RoutingConfiguration routing =
      runCollectivePhase<RoutingConfiguration>(
          mpi_context, "IC distributed routing configuration", [&]() {
            RoutingConfiguration values;
            values.box_size = convertedBoxSizeCode(
                inspection.manifest, config);
            values.chunk_particle_count = std::min(
                options.chunk_particle_count,
                config.mode.ic_staging_particle_count);
            values.batch_particle_count =
                config.mode.ic_staging_particle_count;
            if (values.chunk_particle_count == 0U ||
                values.batch_particle_count == 0U) {
              throw std::invalid_argument(
                  "distributed IC staging particle count must be positive");
            }
            return values;
          });

  std::set<std::size_t> assigned_files;
  std::array<double, 5> local_source_mass{};
  std::unordered_map<std::size_t, IcReaderSession> reader_sessions;
  for (std::size_t file_index = 0;
       file_index < inspection.manifest.source_files.size(); ++file_index) {
    const int reader_rank = static_cast<int>(
        file_index % static_cast<std::size_t>(mpi_context.worldSize()));
    const std::uint64_t file_record_count = std::accumulate(
        inspection.manifest.num_part_this_file[file_index].begin(),
        inspection.manifest.num_part_this_file[file_index].end(),
        std::uint64_t{0});
    for (std::size_t type_index = 0; type_index < kParticleTypeCount;
         ++type_index) {
      const std::size_t total = static_cast<std::size_t>(
          inspection.manifest.num_part_this_file[file_index][type_index]);
      for (std::size_t batch_start = 0U; batch_start < total;
           batch_start += routing.batch_particle_count) {
        const std::size_t batch_count = std::min(
            routing.batch_particle_count, total - batch_start);
        RoutingMpiCollectiveScope routing_mpi_collective_scope;
        const std::uint64_t batch_collective_begin = collective_phase_count;

        std::vector<ParticleRecord> records =
            runCollectivePhase<std::vector<ParticleRecord>>(
                mpi_context, "IC batch read and scientific conversion", [&]() {
                  if (mpi_context.worldRank() != reader_rank) {
                    return std::vector<ParticleRecord>{};
                  }
                  auto session = reader_sessions.find(file_index);
                  if (session == reader_sessions.end()) {
                    injectIcTestFault(
                        mpi_context, "reader_session_creation");
                    session = reader_sessions.try_emplace(
                        file_index,
                        inspection.manifest.source_files[file_index],
                        inspection.manifest.source_file_sizes_bytes[file_index],
                        inspection.manifest.source_sha256[file_index],
                        result.report.counters).first;
                  }
                  std::vector<ParticleRecord> local_records;
                  local_records.reserve(batch_count);
                  for (std::size_t chunk_start = batch_start;
                       chunk_start < batch_start + batch_count;
                       chunk_start += routing.chunk_particle_count) {
                    const std::size_t chunk_count = std::min(
                        routing.chunk_particle_count,
                        batch_start + batch_count - chunk_start);
                    std::vector<ParticleRecord> chunk = readRecordChunk(
                        session->second, inspection, file_index, type_index,
                        chunk_start, chunk_count, config, options,
                        result.report.counters);
                    local_records.insert(
                        local_records.end(),
                        std::make_move_iterator(chunk.begin()),
                        std::make_move_iterator(chunk.end()));
                    checkedCounterAdd(
                        result.report.counters.chunks_assigned, 1U,
                        "chunks_assigned");
                  }
                  assigned_files.insert(file_index);
                  checkedCounterAdd(
                      result.report.counters.routing_batch_count, 1U,
                      "routing_batch_count");
                  checkedCounterAdd(
                      result.report.counters.reader_batches_assigned, 1U,
                      "reader_batches_assigned");
                  checkedCounterAdd(
                      result.report.counters.reader_records_assigned,
                      local_records.size(), "reader_records_assigned");
                  for (const ParticleRecord& record : local_records) {
                    if (record.species >= local_source_mass.size()) {
                      throw std::runtime_error(
                          "source IC record has invalid species tag");
                    }
                    local_source_mass[record.species] += record.mass;
                  }
                  return local_records;
                });
        runCollectivePhaseVoid(
            mpi_context, "IC batch coverage validation", [&]() {
              validateChunkCoverage(
                  mpi_context, file_index, type_index, batch_start,
                  batch_count, reader_rank);
            });

        std::vector<std::vector<std::uint8_t>> per_rank =
            runCollectivePhase<std::vector<std::vector<std::uint8_t>>>(
                mpi_context, "IC batched owner calculation and serialization",
                [&]() {
                  injectIcTestFault(mpi_context, "owner_serialization");
                  std::vector<std::vector<std::uint8_t>> buckets(
                      static_cast<std::size_t>(mpi_context.worldSize()));
                  if (mpi_context.worldRank() == reader_rank) {
                    for (const ParticleRecord& record : records) {
                      const int owner = ownerForX(
                          record.x, routing.box_size,
                          mpi_context.worldSize());
                      internal::serializeIcRecord(
                          record, buckets[static_cast<std::size_t>(owner)]);
                    }
                    checkedCounterAdd(
                        result.report.counters.records_routed,
                        records.size(), "records_routed");
                    mutateIcTestRoute(mpi_context, buckets);
                  }
                  return buckets;
                });
        const std::uint64_t serialized_capacity =
            runCollectivePhase<std::uint64_t>(
                mpi_context, "IC serialized-capacity accounting", [&]() {
                  return nestedByteCapacity(per_rank);
                });
        runCollectivePhaseVoid(
            mpi_context, "IC serialization accounting", [&]() {
              injectIcTestFault(mpi_context, "serialization_accounting");
              std::uint64_t serialized_size = 0U;
              for (const auto& bucket : per_rank) {
                checkedCounterAdd(
                    serialized_size, bucket.size(),
                    "serialized payload size");
              }
              checkedCounterAdd(
                  result.report.counters.bytes_serialized, serialized_size,
                  "bytes_serialized");
            });

        runCollectivePhaseVoid(
            mpi_context, "IC send-layout fault-injection boundary", [&]() {
              injectIcTestFault(mpi_context, "send_layout");
            });
        std::uint64_t sent = 0U;
        std::uint64_t received_bytes = 0U;
        std::uint64_t exchange_peak = 0U;
        const std::vector<std::uint8_t> inbound_bytes = alltoallBytes(
            mpi_context, per_rank, sent, received_bytes, &exchange_peak);
        runCollectivePhaseVoid(
            mpi_context, "IC post-exchange accounting", [&]() {
              injectIcTestFault(mpi_context, "post_exchange_accounting");
              if (mpi_context.isRoot()) {
                checkedCounterAdd(
                    result.report.counters.main_exchange_count, 1U,
                    "main_exchange_count");
              }
              checkedCounterAdd(
                  result.report.counters.bytes_sent, sent, "bytes_sent");
              checkedCounterAdd(
                  result.report.counters.bytes_received, received_bytes,
                  "bytes_received");
            });

        std::vector<ParticleRecord> inbound =
            runCollectivePhase<std::vector<ParticleRecord>>(
                mpi_context, "IC wire validation and deserialization", [&]() {
                  injectIcTestFault(mpi_context, "payload_validation");
                  if (inbound_bytes.size() % internal::kIcWireRecordBytes != 0U) {
                    throw std::runtime_error(
                        "distributed IC wire payload has invalid length");
                  }
                  std::vector<ParticleRecord> decoded;
                  decoded.reserve(inbound_bytes.size() / internal::kIcWireRecordBytes);
                  std::size_t offset = 0U;
                  while (offset < inbound_bytes.size()) {
                    injectIcTestFault(mpi_context, "deserialization");
                    ParticleRecord record =
                        internal::deserializeIcRecord(inbound_bytes, offset);
                    if (record.species >= 5U) {
                      throw std::runtime_error(
                          "distributed IC wire record has invalid species");
                    }
                    IcSpeciesPolicy policy = IcSpeciesPolicy::kDarkMatter;
                    if (record.species == static_cast<std::uint32_t>(
                                              core::ParticleSpecies::kGas)) {
                      policy = IcSpeciesPolicy::kGas;
                    } else if (record.species == static_cast<std::uint32_t>(
                                                     core::ParticleSpecies::kStar)) {
                      policy = IcSpeciesPolicy::kStar;
                    } else if (record.species == static_cast<std::uint32_t>(
                                                     core::ParticleSpecies::kBlackHole)) {
                      policy = IcSpeciesPolicy::kBlackHole;
                    } else if (record.species == static_cast<std::uint32_t>(
                                                     core::ParticleSpecies::kTracer)) {
                      policy = IcSpeciesPolicy::kTracer;
                    }
                    validateRecordScientificState(
                        record, policy, routing.box_size);
                    if (ownerForX(
                            record.x, routing.box_size,
                            mpi_context.worldSize()) !=
                        mpi_context.worldRank()) {
                      throw std::runtime_error(
                          "distributed IC record arrived at the wrong owner");
                    }
                    decoded.push_back(record);
                  }
                  return decoded;
                });

        const std::uint64_t reconciliation_peak =
            exactDistributedChunkReconciliation(
                mpi_context, records, inbound, result.report.counters);
        runCollectivePhaseVoid(
            mpi_context, "IC authoritative-state append", [&]() {
              injectIcTestFault(mpi_context, "sidecar_append");
              appendRecords(
                  result.state, inbound,
                  static_cast<std::uint32_t>(mpi_context.worldRank()));
            });

        const std::uint64_t routing_peak =
            runCollectivePhase<std::uint64_t>(
                mpi_context, "IC routing staging accounting", [&]() {
                  std::uint64_t peak = 0U;
                  checkedCounterAdd(
                      peak, vectorCapacityBytes(records),
                      "routing staging bytes");
                  checkedCounterAdd(
                      peak, serialized_capacity, "routing staging bytes");
                  checkedCounterAdd(
                      peak, exchange_peak, "routing staging bytes");
                  checkedCounterAdd(
                      peak, inbound_bytes.capacity(),
                      "routing staging bytes");
                  checkedCounterAdd(
                      peak, vectorCapacityBytes(inbound),
                      "routing staging bytes");
                  checkedCounterAdd(
                      peak, reconciliation_peak, "routing staging bytes");
                  return peak;
                });
        result.report.counters.peak_staging_bytes = std::max(
            result.report.counters.peak_staging_bytes, routing_peak);
        runCollectivePhaseVoid(
            mpi_context, "IC routing collective accounting", [&]() {
              const std::uint64_t batch_collective_end =
                  collective_phase_count;
              if (batch_collective_end < batch_collective_begin) {
                throw std::overflow_error(
                    "IC routing collective counter wrapped");
              }
              if (mpi_context.isRoot()) {
                checkedCounterAdd(
                    result.report.counters.routing_collective_phase_count,
                    batch_collective_end - batch_collective_begin,
                    "routing_collective_phase_count");
              }
            });
      }
    }
    runCollectivePhaseVoid(
        mpi_context, "IC reader-session completion validation", [&]() {
          injectIcTestFault(mpi_context, "source_identity_revalidation");
          const auto session = reader_sessions.find(file_index);
          if (file_record_count == 0U) {
            return;
          }
          if (mpi_context.worldRank() == reader_rank) {
            if (session == reader_sessions.end()) {
              throw std::runtime_error(
                  "assigned IC reader session was not constructed");
            }
            session->second.revalidateSourceIdentity(
                result.report.counters);
          }
        });
    reader_sessions.erase(file_index);
  }

  result.report.counters.files_assigned = assigned_files.size();
  const std::array<std::uint64_t, 2> reader_record_bounds =
      runCollectivePhase<std::array<std::uint64_t, 2>>(
          mpi_context, "IC reader-work imbalance accounting", [&]() {
            const std::uint64_t local_records =
                result.report.counters.reader_records_assigned;
            std::uint64_t minimum = 0U;
            std::uint64_t maximum = 0U;
            mpiAllreduce(
                &local_records, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                MPI_COMM_WORLD);
            mpiAllreduce(
                &local_records, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                MPI_COMM_WORLD);
            return std::array<std::uint64_t, 2>{minimum, maximum};
          });
  if (mpi_context.isRoot()) {
    result.report.counters.reader_record_imbalance =
        reader_record_bounds[1] - reader_record_bounds[0];
  }
  runCollectivePhaseVoid(
      mpi_context, "IC final authoritative-state construction", [&]() {
        injectIcTestFault(mpi_context, "final_state");
        finalizeImportedState(result.state, inspection.manifest, config);
      });
  result.report.counters.peak_staging_bytes = std::max(
      result.report.counters.peak_staging_bytes,
      exactDistributedIdAudit(
          mpi_context, result.state.particle_sidecar.particle_id,
          config.mode.ic_staging_particle_count,
          result.report.counters));
  validateDistributedTotals(
      mpi_context, result.state, inspection.manifest, local_source_mass);

  result.report.counters.final_local_particle_count =
      result.state.particles.size();
  result.report.counters.final_local_gas_cell_count = result.state.cells.size();
  result.report.counters.final_local_star_count =
      result.state.star_particles.size();
  result.report.counters.final_local_black_hole_count =
      result.state.black_holes.size();
  result.report.counters.final_local_tracer_count = result.state.tracers.size();
  result.report.counters.bytes_read =
      result.report.counters.metadata_bytes_read +
      result.report.counters.hash_bytes_read +
      result.report.counters.payload_bytes_read;
  const std::uint64_t local_routed_records =
      result.report.counters.records_routed;
  std::uint64_t global_routed_records = 0U;
  mpiAllreduce(
      &local_routed_records, &global_routed_records, 1, MPI_UINT64_T, MPI_SUM,
      MPI_COMM_WORLD);
  result.report.counters.logical_consensus_phase_count =
      collective_phase_count;
  result.report.counters.collective_phase_count = collective_phase_count;
  result.report.counters.routing_logical_consensus_phase_count =
      result.report.counters.routing_collective_phase_count;
  result.report.counters.mpi_collective_call_count =
      mpi_collective_counts.total;
  result.report.counters.routing_mpi_collective_call_count =
      mpi_collective_counts.routing_total;
  result.report.counters.nonrouting_mpi_collective_call_count =
      mpi_collective_counts.total >= mpi_collective_counts.routing_total
      ? mpi_collective_counts.total - mpi_collective_counts.routing_total
      : 0U;
  result.report.counters.mpi_allreduce_call_count =
      mpi_collective_counts.allreduce;
  result.report.counters.mpi_bcast_call_count = mpi_collective_counts.bcast;
  result.report.counters.mpi_gather_call_count = mpi_collective_counts.gather;
  result.report.counters.mpi_gatherv_call_count = mpi_collective_counts.gatherv;
  result.report.counters.mpi_alltoall_call_count = mpi_collective_counts.alltoall;
  result.report.counters.mpi_alltoallv_call_count =
      mpi_collective_counts.alltoallv;
  if (global_routed_records > 0U) {
    result.report.counters.collectives_per_million_records =
        static_cast<double>(mpi_collective_counts.total) * 1.0e6 /
        static_cast<double>(global_routed_records);
  }
  return result;
#endif
}

}  // namespace cosmosim::io
