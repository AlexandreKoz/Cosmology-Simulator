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
#include "io/internal/ic_hdf5_handle.hpp"
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
using internal::IcReaderSession;
#endif

IcReadResult readGadgetArepoHdf5Ic(
    const std::filesystem::path& ic_path,
    const core::SimulationConfig& config,
    const IcImportOptions& options) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(ic_path);
  static_cast<void>(config);
  static_cast<void>(options);
  throw std::runtime_error(
      "COSMOSIM_ENABLE_HDF5=OFF: GADGET/AREPO HDF5 IC reader unavailable in this build");
#else
  if (options.chunk_particle_count == 0U) {
    throw std::invalid_argument("chunk_particle_count must be positive");
  }
  if (options.manifest != nullptr) {
    validateIcManifest(*options.manifest);
  }
  Inspection inspection = inspectFileSet(ic_path, config, options);
  if (options.validate_runtime_cosmology) {
    validateRuntimeCosmology(inspection.manifest, config);
  }

  IcReadResult result;
  result.report.counters = inspection.counters;
  result.report.manifest = inspection.manifest;
  result.report.schema = inspection.schemas.front();
  result.report.defaulted_fields = inspection.manifest.defaulted_fields;
  for (const auto& value : result.report.defaulted_fields) {
    const auto separator = value.find('=');
    result.report.missing_optional_fields.push_back(
        value.substr(0U, separator));
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

  for (std::size_t file_index = 0U;
       file_index < inspection.manifest.source_files.size(); ++file_index) {
    ++result.report.counters.files_assigned;
    IcReaderSession session(
        inspection.manifest.source_files[file_index],
        inspection.manifest.source_file_sizes_bytes[file_index],
        inspection.manifest.source_sha256[file_index],
        result.report.counters);
    for (std::size_t type_index = 0U; type_index < 6U; ++type_index) {
      const std::size_t total = static_cast<std::size_t>(
          inspection.manifest.num_part_this_file[file_index][type_index]);
      for (std::size_t start = 0U; start < total;
           start += options.chunk_particle_count) {
        const std::size_t count =
            std::min(options.chunk_particle_count, total - start);
        auto records = readRecordChunk(
            session, inspection, file_index, type_index, start, count, config, options,
            result.report.counters);
        appendRecords(result.state, records, 0U);
        ++result.report.counters.chunks_assigned;
        result.report.counters.records_routed += records.size();
      }
    }
    session.revalidateSourceIdentity(result.report.counters);
  }

  validateSerialCountsAndIds(result.state, inspection.manifest);
  finalizeImportedState(result.state, inspection.manifest, config);
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
  return result;
#endif
}


IcImportReport internal::streamGadgetArepoHdf5Ic(
    const std::filesystem::path& ic_path,
    const core::SimulationConfig& config,
    const IcImportOptions& options,
    const IcManifestReadyCallback& on_manifest_ready,
    const IcRecordBatchCallback& on_record_batch) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(ic_path);
  static_cast<void>(config);
  static_cast<void>(options);
  static_cast<void>(on_manifest_ready);
  static_cast<void>(on_record_batch);
  throw std::runtime_error(
      "COSMOSIM_ENABLE_HDF5=OFF: streaming GADGET/AREPO IC ingestion unavailable");
#else
  if (options.chunk_particle_count == 0U) {
    throw std::invalid_argument("chunk_particle_count must be positive");
  }
  if (!on_manifest_ready || !on_record_batch) {
    throw std::invalid_argument(
        "streaming IC ingestion requires manifest and record callbacks");
  }
  if (options.manifest != nullptr) {
    validateIcManifest(*options.manifest);
  }
  Inspection inspection = inspectFileSet(ic_path, config, options);
  if (options.validate_runtime_cosmology) {
    validateRuntimeCosmology(inspection.manifest, config);
  }

  IcImportReport report;
  report.counters = inspection.counters;
  report.manifest = inspection.manifest;
  report.schema = inspection.schemas.front();
  report.defaulted_fields = inspection.manifest.defaulted_fields;
  for (const auto& value : report.defaulted_fields) {
    const auto separator = value.find('=');
    report.missing_optional_fields.push_back(value.substr(0U, separator));
  }
  report.unsupported_fields = inspection.manifest.dropped_fields;
  report.manifest_verified =
      inspection.manifest.canonical_source_manifest_verified;
  report.verified_manifest_sha256 =
      inspection.manifest.canonical_source_manifest_sha256;
  report.provenance_authority = options.manifest != nullptr
      ? "supplied_manifest_v1"
      : (inspection.manifest.canonical_source_manifest_verified
            ? "canonical_embedded_manifest_v2"
            : "runtime_config_and_inspected_source");
  on_manifest_ready(inspection.manifest);

  std::array<std::uint64_t, 6> observed_by_source_type{};
  for (std::size_t file_index = 0U;
       file_index < inspection.manifest.source_files.size(); ++file_index) {
    ++report.counters.files_assigned;
    IcReaderSession session(
        inspection.manifest.source_files[file_index],
        inspection.manifest.source_file_sizes_bytes[file_index],
        inspection.manifest.source_sha256[file_index],
        report.counters);
    for (std::size_t type_index = 0U; type_index < 6U; ++type_index) {
      const std::size_t total = static_cast<std::size_t>(
          inspection.manifest.num_part_this_file[file_index][type_index]);
      for (std::size_t start = 0U; start < total;
           start += options.chunk_particle_count) {
        const std::size_t count =
            std::min(options.chunk_particle_count, total - start);
        auto records = readRecordChunk(
            session, inspection, file_index, type_index, start, count, config,
            options, report.counters);
        on_record_batch(records);
        checkedCounterAdd(
            observed_by_source_type[type_index], records.size(),
            "streamed_source_type_count");
        ++report.counters.chunks_assigned;
        checkedCounterAdd(
            report.counters.records_routed, records.size(),
            "records_routed");
      }
    }
    session.revalidateSourceIdentity(report.counters);
  }
  for (std::size_t type_index = 0U; type_index < 6U; ++type_index) {
    if (observed_by_source_type[type_index] !=
        inspection.manifest.num_part_total[type_index]) {
      throw std::runtime_error(
          "streaming IC ingestion source coverage disagrees with the manifest");
    }
  }
  report.counters.final_local_particle_count =
      std::accumulate(
          observed_by_source_type.begin(), observed_by_source_type.end(),
          std::uint64_t{0});
  report.counters.bytes_read =
      report.counters.metadata_bytes_read +
      report.counters.hash_bytes_read +
      report.counters.payload_bytes_read;
  return report;
#endif
}


}  // namespace cosmosim::io
