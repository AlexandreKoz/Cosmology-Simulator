#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::parallel {
class MpiContext;
}

namespace cosmosim::io {

enum class IcParticleKind : std::uint8_t {
  kGas = 0,
  kDarkMatter = 1,
  kDisk = 2,
  kBulge = 3,
  kStar = 4,
  kBoundary = 5,
};

enum class IcDialect : std::uint8_t {
  kChuiCanonicalV1 = 0,
  kGadgetArepoBridgeV1 = 1,
};

enum class IcCoordinateFrame : std::uint8_t {
  kComoving = 0,
  kPhysical = 1,
};

enum class IcVelocityConvention : std::uint8_t {
  kNotVelocity = 0,
  kPhysicalPeculiar = 1,
  kSqrtAScaledPeculiar = 2,
  kComovingCoordinateRate = 3,
};

enum class IcFieldSemantics : std::uint8_t {
  kCoordinate = 0,
  kVelocity = 1,
  kExtensive = 2,
  kIntensive = 3,
  kSpecific = 4,
  kIdentifier = 5,
};

enum class IcFieldDisposition : std::uint8_t {
  kConverted = 0,
  kPreserved = 1,
  kDefaulted = 2,
  kDropped = 3,
  kRejected = 4,
};

enum class IcSpeciesPolicy : std::uint8_t {
  kReject = 0,
  kGas = 1,
  kDarkMatter = 2,
  kCollisionlessFamily2AsDarkMatter = 3,
  kCollisionlessFamily3AsDarkMatter = 4,
  kStar = 5,
  kBlackHole = 6,
  kTracer = 7,
};

enum class IcScalarClass : std::uint8_t {
  kFloatingPoint = 0,
  kInteger = 1,
};

enum class IcByteOrder : std::uint8_t {
  kNotApplicable = 0,
  kLittleEndian = 1,
  kBigEndian = 2,
  kNative = 3,
};

struct IcFieldManifest {
  std::uint32_t source_file_index = 0;
  std::string dataset_path;
  std::string selected_alias;
  std::string scalar_type;
  IcScalarClass scalar_class = IcScalarClass::kFloatingPoint;
  std::uint8_t byte_width = 0;
  bool is_signed = false;
  IcByteOrder byte_order = IcByteOrder::kNotApplicable;
  std::uint8_t rank = 0;
  std::vector<std::uint64_t> dimensions;
  std::uint64_t record_count = 0;
  // Explicit SI multiplier for one stored unit before h/a scaling.
  double base_unit_to_si = 1.0;
  double hubble_exponent = 0.0;
  double scale_factor_exponent = 0.0;
  IcCoordinateFrame coordinate_frame = IcCoordinateFrame::kComoving;
  IcVelocityConvention velocity_convention =
      IcVelocityConvention::kNotVelocity;
  IcFieldSemantics semantics = IcFieldSemantics::kIntensive;
  IcFieldDisposition disposition = IcFieldDisposition::kConverted;
  std::string source_unit;
  std::string target_unit;
  std::string conversion_equation;
};

struct IcManifest {
  std::string schema_name = "chui_ic_audit_manifest";
  std::uint32_t schema_version = 2;
  std::string converter_version = "chui_ic_converter_v1";
  IcDialect dialect = IcDialect::kGadgetArepoBridgeV1;
  std::string dialect_version = "1";
  std::vector<std::filesystem::path> source_files;
  std::vector<std::string> source_provenance_ids;
  std::vector<std::uint64_t> source_file_sizes_bytes;
  std::vector<std::string> source_sha256;
  std::vector<std::string> original_header_attributes;
  std::vector<std::array<std::uint64_t, 6>> num_part_this_file;
  std::array<std::uint64_t, 6> num_part_total{};
  std::array<std::uint32_t, 6> num_part_total_high_word{};
  std::uint32_t num_files_per_snapshot = 1;
  std::array<double, 6> mass_table{};
  double box_size = 0.0;
  double scale_factor = 1.0;
  double redshift = 0.0;
  double omega_matter = 0.0;
  double omega_lambda = 0.0;
  double hubble_param = 0.0;
  std::vector<IcFieldManifest> fields;
  std::array<IcSpeciesPolicy, 6> species_policy{
      IcSpeciesPolicy::kGas,
      IcSpeciesPolicy::kDarkMatter,
      IcSpeciesPolicy::kReject,
      IcSpeciesPolicy::kReject,
      IcSpeciesPolicy::kStar,
      IcSpeciesPolicy::kBlackHole};
  std::vector<std::string> defaulted_fields;
  std::vector<std::string> converted_fields;
  std::vector<std::string> dropped_fields;
  std::vector<std::string> rejected_fields;
  std::vector<std::string> preserved_auxiliary_fields;
  std::vector<std::string> conversion_equations;
  std::vector<std::string> warnings;
};

void validateIcManifest(const IcManifest& manifest);
[[nodiscard]] double icStoredToSiMultiplier(
    const IcFieldManifest& field,
    double hubble_param,
    double scale_factor);
[[nodiscard]] double icVelocityConventionMultiplier(
    IcVelocityConvention convention,
    double scale_factor);
[[nodiscard]] std::string serializeIcManifestJson(const IcManifest& manifest);
[[nodiscard]] IcManifest deserializeIcManifestJson(std::string_view json_text);
[[nodiscard]] IcManifest readIcManifestJson(
    const std::filesystem::path& input_path);
void writeIcManifestJson(
    const IcManifest& manifest,
    const std::filesystem::path& output_path);

// A compact per-file schema summary for transparent import auditing.
struct IcSchemaSummary {
  std::array<std::uint64_t, 6> count_by_type{};
  std::array<std::uint64_t, 6> total_count_by_type{};
  std::array<std::uint32_t, 6> total_count_high_word{};
  std::array<double, 6> mass_table{};
  std::uint32_t num_files_per_snapshot = 1;
  double box_size = 0.0;
  double scale_factor = 1.0;
  double redshift = 0.0;
  double omega_matter = 0.0;
  double omega_lambda = 0.0;
  double hubble_param = 0.0;
  bool velocities_are_peculiar = true;
};

// Build the explicit recognized direct-import contract. This helper does not
// inspect the filesystem; callers must populate the schema from audited header
// metadata and validate the returned manifest before use.
[[nodiscard]] IcManifest makeGadgetArepoBridgeV1Manifest(
    const std::filesystem::path& ic_path,
    const IcSchemaSummary& schema);

struct IcImportOptions {
  bool require_velocities = true;
  bool require_particle_ids = true;
  bool allow_mass_table_fallback = true;
  std::size_t chunk_particle_count = 1u << 16;
  // Null selects the recognized, versioned gadget_arepo_bridge_v1 contract.
  // A non-null manifest is caller-owned for the synchronous read.
  const IcManifest* manifest = nullptr;
};

struct IcImportCounters {
  std::uint64_t files_assigned = 0;
  std::uint64_t chunks_assigned = 0;
  std::uint64_t bytes_read = 0;
  std::uint64_t records_read = 0;
  std::uint64_t records_converted = 0;
  std::uint64_t records_routed = 0;
  std::uint64_t bytes_sent = 0;
  std::uint64_t bytes_received = 0;
  std::uint64_t peak_staging_bytes = 0;
  std::uint64_t final_local_particle_count = 0;
  std::uint64_t final_local_gas_cell_count = 0;
  std::uint64_t final_local_star_count = 0;
  std::uint64_t final_local_black_hole_count = 0;
  std::uint64_t final_local_tracer_count = 0;
};

struct IcImportReport {
  IcSchemaSummary schema;
  std::optional<IcManifest> manifest;
  std::vector<std::string> present_aliases;
  std::vector<std::string> defaulted_fields;
  std::vector<std::string> missing_optional_fields;
  std::vector<std::string> unsupported_fields;
  IcImportCounters counters;
  bool already_partitioned = false;
};

struct IcReadResult {
  core::SimulationState state;
  IcImportReport report;
};

// Read common GADGET/AREPO-style HDF5 initial conditions into SimulationState.
// This API is species-aware and writes directly into SoA buffers in chunked batches.
[[nodiscard]] IcReadResult readGadgetArepoHdf5Ic(
    const std::filesystem::path& ic_path,
    const core::SimulationConfig& config,
    const IcImportOptions& options = {});

// MPI-aware bounded reader. Source chunks are assigned once, converted once,
// and routed directly to deterministic initial owners. The returned state is
// already partitioned and must not pass through replicated-state compaction.
[[nodiscard]] IcReadResult readDistributedGadgetArepoHdf5Ic(
    const std::filesystem::path& ic_path,
    const core::SimulationConfig& config,
    const parallel::MpiContext& mpi_context,
    const IcImportOptions& options = {});

// Internal helper path for generated isolated test problems without external files.
[[nodiscard]] IcReadResult buildGeneratedIsolatedIc(
    const core::SimulationConfig& config,
    std::size_t dark_matter_particle_count,
    std::size_t gas_particle_count,
    std::uint64_t particle_id_seed = 1u);

// A small deterministic converter for compatibility tooling and regression tests.
[[nodiscard]] IcReadResult convertGeneratedIsolatedIcToState(
    const core::SimulationConfig& config,
    std::size_t particle_count_per_axis);

}  // namespace cosmosim::io
