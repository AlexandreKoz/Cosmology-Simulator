#include "cosmosim/io/ic_reader.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/units.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace cosmosim::io {
namespace {

constexpr std::uint32_t k_species_dark_matter =
    static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
constexpr std::uint32_t k_species_gas = static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
constexpr std::uint32_t k_species_star = static_cast<std::uint32_t>(core::ParticleSpecies::kStar);
constexpr std::uint32_t k_species_black_hole =
    static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole);

[[nodiscard]] std::uint32_t mapTypeIndexToSpeciesTag(
    std::size_t type_index,
    const IcManifest& manifest) {
  switch (manifest.species_policy.at(type_index)) {
    case IcSpeciesPolicy::kGas:
      return k_species_gas;
    case IcSpeciesPolicy::kDarkMatter:
    case IcSpeciesPolicy::kCollisionlessFamily2AsDarkMatter:
    case IcSpeciesPolicy::kCollisionlessFamily3AsDarkMatter:
      return k_species_dark_matter;
    case IcSpeciesPolicy::kStar:
      return k_species_star;
    case IcSpeciesPolicy::kBlackHole:
      return k_species_black_hole;
    case IcSpeciesPolicy::kTracer:
      return static_cast<std::uint32_t>(core::ParticleSpecies::kTracer);
    case IcSpeciesPolicy::kReject:
      break;
  }
  throw std::runtime_error(
      "IC particle type " + std::to_string(type_index) +
      " has no accepted species mapping policy");
}

[[nodiscard]] bool isFinite(double value) { return std::isfinite(value) != 0; }

[[nodiscard]] std::string fileFnv1aProvenanceId(
    const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error(
        "failed to hash IC source file: " + path.string());
  }
  std::uint64_t hash = 1469598103934665603ULL;
  std::array<char, 1U << 15U> buffer{};
  while (input) {
    input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const std::streamsize count = input.gcount();
    for (std::streamsize i = 0; i < count; ++i) {
      hash ^= static_cast<unsigned char>(buffer[static_cast<std::size_t>(i)]);
      hash *= 1099511628211ULL;
    }
  }
  if (!input.eof()) {
    throw std::runtime_error(
        "failed while hashing IC source file: " + path.string());
  }
  std::ostringstream encoded;
  encoded << "fnv1a64:" << std::hex << std::setw(16) << std::setfill('0')
          << hash;
  return encoded.str();
}

void validateFiniteField(const std::vector<double>& values, std::string_view field_name) {
  for (double value : values) {
    if (!isFinite(value)) {
      throw std::runtime_error("non-finite value in field: " + std::string(field_name));
    }
  }
}

void applyLengthFrameConversion(
    std::vector<double>& coordinate_values,
    double scale_factor,
    std::string_view source_frame,
    std::string_view target_frame) {
  if (source_frame == target_frame) {
    return;
  }
  if (source_frame == "physical" && target_frame == "comoving") {
    for (double& value : coordinate_values) {
      value = core::physicalToComovingLength(value, scale_factor);
    }
    return;
  }
  if (source_frame == "comoving" && target_frame == "physical") {
    for (double& value : coordinate_values) {
      value = core::comovingToPhysicalLength(value, scale_factor);
    }
    return;
  }
  throw std::runtime_error(
      "unsupported coordinate frame conversion from " + std::string(source_frame) + " to " +
      std::string(target_frame));
}

void fillSpeciesLedger(core::SimulationState& state) {
  state.species.count_by_species = {};
  for (std::uint32_t species_tag : state.particle_sidecar.species_tag) {
    if (species_tag < state.species.count_by_species.size()) {
      state.species.count_by_species[species_tag] += 1;
    }
  }
}

[[nodiscard]] IcManifest buildRecognizedSingleFileManifest(
    const std::filesystem::path& ic_path,
    const IcSchemaSummary& schema,
    const core::UnitSystem& source_units) {
  IcManifest manifest;
  manifest.dialect = IcDialect::kGadgetArepoBridgeV1;
  manifest.dialect_version = "1";
  manifest.source_files = {ic_path.lexically_normal()};
  manifest.num_part_this_file = {schema.count_by_type};
  manifest.num_part_total = schema.total_count_by_type;
  manifest.num_part_total_high_word = schema.total_count_high_word;
  manifest.num_files_per_snapshot = schema.num_files_per_snapshot;
  manifest.mass_table = schema.mass_table;
  manifest.box_size = schema.box_size;
  manifest.scale_factor = schema.scale_factor;
  manifest.redshift = schema.redshift;
  manifest.omega_matter = schema.omega_matter;
  manifest.omega_lambda = schema.omega_lambda;
  manifest.hubble_param = schema.hubble_param;

  const auto add_field = [&](std::string path,
                             std::string scalar_type,
                             std::vector<std::uint64_t> dimensions,
                             double base_unit_to_si,
                             IcFieldSemantics semantics,
                             IcVelocityConvention velocity_convention =
                                 IcVelocityConvention::kNotVelocity) {
    const std::uint64_t records = dimensions.front();
    manifest.fields.push_back(IcFieldManifest{
        .dataset_path = std::move(path),
        .scalar_type = std::move(scalar_type),
        .rank = static_cast<std::uint8_t>(dimensions.size()),
        .dimensions = std::move(dimensions),
        .record_count = records,
        .base_unit_to_si = base_unit_to_si,
        .hubble_exponent = 0.0,
        .scale_factor_exponent = 0.0,
        .coordinate_frame = IcCoordinateFrame::kComoving,
        .velocity_convention = velocity_convention,
        .semantics = semantics,
        .disposition = IcFieldDisposition::kConverted});
  };
  add_field("/Header/MassTable", "float64", {6U},
            source_units.mass_si_per_code,
            IcFieldSemantics::kExtensive);
  add_field("/Header/BoxSize", "float64", {1U},
            source_units.length_si_per_code,
            IcFieldSemantics::kCoordinate);
  for (std::size_t type = 0; type < schema.count_by_type.size(); ++type) {
    const std::uint64_t count = schema.count_by_type[type];
    if (count == 0U) {
      continue;
    }
    const std::string prefix = "/PartType" + std::to_string(type) + "/";
    add_field(prefix + "Coordinates", "float64", {count, 3U},
              source_units.length_si_per_code,
              IcFieldSemantics::kCoordinate);
    add_field(prefix + "Velocities", "float64", {count, 3U},
              source_units.velocity_si_per_code,
              IcFieldSemantics::kVelocity,
              IcVelocityConvention::kPhysicalPeculiar);
    add_field(prefix + "ParticleIDs", "uint64", {count}, 1.0,
              IcFieldSemantics::kIdentifier);
    if (schema.mass_table[type] == 0.0) {
      add_field(prefix + "Masses", "float64", {count},
                source_units.mass_si_per_code,
                IcFieldSemantics::kExtensive);
    }
    if (type == 0U) {
      add_field(prefix + "InternalEnergy", "float64", {count},
                source_units.velocity_si_per_code *
                    source_units.velocity_si_per_code,
                IcFieldSemantics::kSpecific);
      add_field(prefix + "Density", "float64", {count},
                source_units.mass_si_per_code /
                    std::pow(source_units.length_si_per_code, 3.0),
                IcFieldSemantics::kIntensive);
    } else if (type == 5U) {
      add_field(prefix + "BH_Mass", "float64", {count},
                source_units.mass_si_per_code,
                IcFieldSemantics::kExtensive);
      add_field(prefix + "BH_Mdot", "float64", {count},
                source_units.mass_si_per_code /
                    source_units.timeSiPerCode(),
                IcFieldSemantics::kIntensive);
    }
  }
  return manifest;
}

[[nodiscard]] const IcFieldManifest& requireFieldConvention(
    const IcManifest& manifest,
    std::string_view canonical_path,
    IcFieldSemantics semantics) {
  for (const IcFieldManifest& field : manifest.fields) {
    if (field.dataset_path == canonical_path) {
      if (field.semantics != semantics) {
        throw std::runtime_error(
            "IC manifest field semantics mismatch for " +
            std::string(canonical_path));
      }
      return field;
    }
  }
  throw std::runtime_error(
      "IC manifest lacks required conversion entry for " +
      std::string(canonical_path));
}

void applyManifestUnitConversion(
    std::vector<double>& values,
    const IcFieldManifest& field,
    const IcManifest& manifest,
    double target_si_per_code) {
  if (!(target_si_per_code > 0.0)) {
    throw std::invalid_argument("target SI-per-code conversion must be positive");
  }
  double multiplier = icStoredToSiMultiplier(
      field, manifest.hubble_param, manifest.scale_factor) /
      target_si_per_code;
  if (field.semantics == IcFieldSemantics::kVelocity) {
    multiplier *= icVelocityConventionMultiplier(
        field.velocity_convention, manifest.scale_factor);
  }
  for (double& value : values) {
    value *= multiplier;
  }
}

[[nodiscard]] const std::vector<std::string>& gasFieldAliases(std::string_view canonical_key) {
  static const std::vector<std::string> k_internal_energy_aliases = {
      "InternalEnergy", "U", "Internal_Energy"};
  static const std::vector<std::string> k_density_aliases = {"Density", "Rho"};
  static const std::vector<std::string> k_metallicity_aliases = {"Metallicity", "GFM_Metallicity"};
  static const std::vector<std::string> k_smoothing_length_aliases = {
      "SmoothingLength", "Hsml", "Smoothing_Length"};
  if (canonical_key == "InternalEnergy") {
    return k_internal_energy_aliases;
  }
  if (canonical_key == "Density") {
    return k_density_aliases;
  }
  if (canonical_key == "Metallicity") {
    return k_metallicity_aliases;
  }
  if (canonical_key == "SmoothingLength") {
    return k_smoothing_length_aliases;
  }
  throw std::runtime_error("unknown gas field alias key: " + std::string(canonical_key));
}

#if COSMOSIM_ENABLE_HDF5

class Hdf5Handle {
 public:
  explicit Hdf5Handle(hid_t handle = -1) : m_handle(handle) {}
  Hdf5Handle(const Hdf5Handle&) = delete;
  Hdf5Handle& operator=(const Hdf5Handle&) = delete;
  Hdf5Handle(Hdf5Handle&& other) noexcept : m_handle(other.m_handle) { other.m_handle = -1; }
  Hdf5Handle& operator=(Hdf5Handle&& other) noexcept {
    if (this != &other) {
      close();
      m_handle = other.m_handle;
      other.m_handle = -1;
    }
    return *this;
  }
  ~Hdf5Handle() { close(); }

  [[nodiscard]] hid_t get() const { return m_handle; }
  [[nodiscard]] bool valid() const { return m_handle >= 0; }

 private:
  void close() {
    if (m_handle >= 0) {
      const H5I_type_t type = H5Iget_type(m_handle);
      if (type == H5I_FILE) {
        H5Fclose(m_handle);
      } else if (type == H5I_GROUP) {
        H5Gclose(m_handle);
      } else if (type == H5I_DATASET) {
        H5Dclose(m_handle);
      } else if (type == H5I_DATASPACE) {
        H5Sclose(m_handle);
      } else if (type == H5I_ATTR) {
        H5Aclose(m_handle);
      } else if (type == H5I_DATATYPE) {
        H5Tclose(m_handle);
      }
      m_handle = -1;
    }
  }

  hid_t m_handle = -1;
};

[[nodiscard]] bool hdf5PathExists(hid_t parent, const std::string& path) {
  return H5Lexists(parent, path.c_str(), H5P_DEFAULT) > 0;
}

[[nodiscard]] std::string pickAlias(
    hid_t group,
    const std::vector<std::string>& aliases,
    IcImportReport& report,
    const std::string& fallback_key,
    bool required) {
  for (const std::string& alias : aliases) {
    if (hdf5PathExists(group, alias)) {
      report.present_aliases.push_back(fallback_key + "=" + alias);
      return alias;
    }
  }
  if (required) {
    throw std::runtime_error("required dataset missing for key: " + fallback_key);
  }
  report.missing_optional_fields.push_back(fallback_key);
  return {};
}

void readHeaderArrays(hid_t header_group, IcSchemaSummary& schema) {
  Hdf5Handle attr_count(H5Aopen(header_group, "NumPart_ThisFile", H5P_DEFAULT));
  if (!attr_count.valid()) {
    throw std::runtime_error("Header/NumPart_ThisFile attribute is required");
  }
  std::array<std::uint32_t, 6> raw_count{};
  if (H5Aread(attr_count.get(), H5T_NATIVE_UINT32, raw_count.data()) < 0) {
    throw std::runtime_error("failed to read NumPart_ThisFile");
  }
  for (std::size_t i = 0; i < 6; ++i) {
    schema.count_by_type[i] = raw_count[i];
  }

  Hdf5Handle attr_total(H5Aopen(header_group, "NumPart_Total", H5P_DEFAULT));
  if (!attr_total.valid()) {
    throw std::runtime_error("Header/NumPart_Total attribute is required");
  }
  std::array<std::uint32_t, 6> raw_total{};
  if (H5Aread(attr_total.get(), H5T_NATIVE_UINT32, raw_total.data()) < 0) {
    throw std::runtime_error("failed to read NumPart_Total");
  }

  Hdf5Handle attr_high_word(
      H5Aopen(header_group, "NumPart_Total_HighWord", H5P_DEFAULT));
  if (attr_high_word.valid() &&
      H5Aread(attr_high_word.get(), H5T_NATIVE_UINT32,
              schema.total_count_high_word.data()) < 0) {
    throw std::runtime_error("failed to read NumPart_Total_HighWord");
  }
  for (std::size_t i = 0; i < 6; ++i) {
    schema.total_count_by_type[i] =
        static_cast<std::uint64_t>(raw_total[i]) |
        (static_cast<std::uint64_t>(schema.total_count_high_word[i]) << 32U);
  }

  Hdf5Handle attr_mass(H5Aopen(header_group, "MassTable", H5P_DEFAULT));
  if (!attr_mass.valid()) {
    throw std::runtime_error("Header/MassTable attribute is required");
  }
  if (H5Aread(attr_mass.get(), H5T_NATIVE_DOUBLE, schema.mass_table.data()) < 0) {
    throw std::runtime_error("failed to read MassTable");
  }

  Hdf5Handle attr_time(H5Aopen(header_group, "Time", H5P_DEFAULT));
  if (!attr_time.valid() ||
      H5Aread(attr_time.get(), H5T_NATIVE_DOUBLE, &schema.scale_factor) < 0) {
    throw std::runtime_error("failed to read Header/Time");
  }

  const auto read_required_double = [&](const char* name, double& value) {
    Hdf5Handle attribute(H5Aopen(header_group, name, H5P_DEFAULT));
    if (!attribute.valid() ||
        H5Aread(attribute.get(), H5T_NATIVE_DOUBLE, &value) < 0) {
      throw std::runtime_error(
          std::string("failed to read required Header/") + name);
    }
  };
  read_required_double("Redshift", schema.redshift);
  read_required_double("BoxSize", schema.box_size);
  read_required_double("Omega0", schema.omega_matter);
  read_required_double("OmegaLambda", schema.omega_lambda);
  read_required_double("HubbleParam", schema.hubble_param);

  Hdf5Handle attr_num_files(
      H5Aopen(header_group, "NumFilesPerSnapshot", H5P_DEFAULT));
  if (!attr_num_files.valid() ||
      H5Aread(attr_num_files.get(), H5T_NATIVE_UINT32,
              &schema.num_files_per_snapshot) < 0) {
    throw std::runtime_error(
        "failed to read required Header/NumFilesPerSnapshot");
  }
}

void readDatasetChunk1d(
    hid_t group,
    const std::string& dataset_name,
    std::size_t start,
    std::size_t count,
    std::vector<double>& out) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("failed to open dataset: " + dataset_name);
  }
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  if (!file_space.valid()) {
    throw std::runtime_error("failed to get dataspace for: " + dataset_name);
  }
  if (H5Sget_simple_extent_ndims(file_space.get()) != 1) {
    throw std::runtime_error("dataset must have rank 1: " + dataset_name);
  }
  hsize_t file_dims[1]{};
  if (H5Sget_simple_extent_dims(file_space.get(), file_dims, nullptr) < 0 ||
      start + count > file_dims[0]) {
    throw std::runtime_error("dataset extent/count mismatch: " + dataset_name);
  }
  Hdf5Handle datatype(H5Dget_type(dataset.get()));
  if (!datatype.valid() || H5Tget_class(datatype.get()) != H5T_FLOAT) {
    throw std::runtime_error("dataset must use a floating-point type: " + dataset_name);
  }

  hsize_t file_offset[1] = {static_cast<hsize_t>(start)};
  hsize_t file_count[1] = {static_cast<hsize_t>(count)};
  if (H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr) <
      0) {
    throw std::runtime_error("failed to select hyperslab for: " + dataset_name);
  }

  hsize_t mem_dims[1] = {static_cast<hsize_t>(count)};
  Hdf5Handle mem_space(H5Screate_simple(1, mem_dims, nullptr));
  out.resize(count);
  if (H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, mem_space.get(), file_space.get(), H5P_DEFAULT,
              out.data()) < 0) {
    throw std::runtime_error("failed to read dataset chunk: " + dataset_name);
  }
}

void readDatasetChunk2d(
    hid_t group,
    const std::string& dataset_name,
    std::size_t start,
    std::size_t count,
    std::vector<double>& out) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("failed to open dataset: " + dataset_name);
  }
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  if (!file_space.valid()) {
    throw std::runtime_error("failed to get dataspace for: " + dataset_name);
  }
  if (H5Sget_simple_extent_ndims(file_space.get()) != 2) {
    throw std::runtime_error("vector dataset must have rank 2: " + dataset_name);
  }
  hsize_t file_dims[2]{};
  if (H5Sget_simple_extent_dims(file_space.get(), file_dims, nullptr) < 0 ||
      start + count > file_dims[0] || file_dims[1] != 3U) {
    throw std::runtime_error(
        "vector dataset must have dimensions [record_count,3]: " +
        dataset_name);
  }
  Hdf5Handle datatype(H5Dget_type(dataset.get()));
  if (!datatype.valid() || H5Tget_class(datatype.get()) != H5T_FLOAT) {
    throw std::runtime_error(
        "vector dataset must use a floating-point type: " + dataset_name);
  }

  hsize_t file_offset[2] = {static_cast<hsize_t>(start), 0};
  hsize_t file_count[2] = {static_cast<hsize_t>(count), 3};
  if (H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr) <
      0) {
    throw std::runtime_error("failed to select hyperslab for: " + dataset_name);
  }

  hsize_t mem_dims[2] = {static_cast<hsize_t>(count), 3};
  Hdf5Handle mem_space(H5Screate_simple(2, mem_dims, nullptr));
  out.resize(count * 3);
  if (H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, mem_space.get(), file_space.get(), H5P_DEFAULT,
              out.data()) < 0) {
    throw std::runtime_error("failed to read dataset chunk: " + dataset_name);
  }
}

void readDatasetChunkIds(
    hid_t group,
    const std::string& dataset_name,
    std::size_t start,
    std::size_t count,
    std::vector<std::uint64_t>& out) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("failed to open dataset: " + dataset_name);
  }
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  if (!file_space.valid() || H5Sget_simple_extent_ndims(file_space.get()) != 1) {
    throw std::runtime_error("particle ID dataset must have rank 1");
  }
  hsize_t file_dims[1]{};
  if (H5Sget_simple_extent_dims(file_space.get(), file_dims, nullptr) < 0 ||
      start + count > file_dims[0]) {
    throw std::runtime_error("particle ID dataset extent/count mismatch");
  }
  Hdf5Handle datatype(H5Dget_type(dataset.get()));
  if (!datatype.valid() || H5Tget_class(datatype.get()) != H5T_INTEGER) {
    throw std::runtime_error("particle ID dataset must use an integer type");
  }
  hsize_t file_offset[1] = {static_cast<hsize_t>(start)};
  hsize_t file_count[1] = {static_cast<hsize_t>(count)};
  if (H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr) <
      0) {
    throw std::runtime_error("failed to select hyperslab for IDs");
  }

  hsize_t mem_dims[1] = {static_cast<hsize_t>(count)};
  Hdf5Handle mem_space(H5Screate_simple(1, mem_dims, nullptr));
  out.resize(count);
  if (H5Dread(dataset.get(), H5T_NATIVE_UINT64, mem_space.get(), file_space.get(), H5P_DEFAULT,
              out.data()) < 0) {
    throw std::runtime_error("failed to read particle IDs");
  }
}

#endif

}  // namespace

IcManifest makeGadgetArepoBridgeV1Manifest(
    const std::filesystem::path& ic_path,
    const IcSchemaSummary& schema) {
  const core::UnitSystem source_units =
      core::makeUnitSystem("kpc", "msun", "km_s");
  return buildRecognizedSingleFileManifest(ic_path, schema, source_units);
}

IcReadResult buildGeneratedIsolatedIc(
    const core::SimulationConfig& config,
    std::size_t dark_matter_particle_count,
    std::size_t gas_particle_count,
    std::uint64_t particle_id_seed) {
  IcReadResult result;
  const std::size_t total_particle_count = dark_matter_particle_count + gas_particle_count;
  result.state.resizeParticles(total_particle_count);
  result.state.resizeCells(gas_particle_count);
  result.state.metadata.scale_factor = 1.0;
  result.state.metadata.run_name = config.output.run_name;
  result.state.metadata.snapshot_stem = config.output.output_stem;
  result.state.metadata.restart_stem = config.output.restart_stem;

  const core::UnitSystem target_units =
      core::makeUnitSystem(config.units.length_unit, config.units.mass_unit, config.units.velocity_unit);

  std::size_t global_index = 0;
  for (std::size_t i = 0; i < gas_particle_count; ++i, ++global_index) {
    const double x = static_cast<double>(i) * 0.01;
    result.state.particle_sidecar.particle_id[global_index] = particle_id_seed + global_index;
    result.state.particle_sidecar.species_tag[global_index] = k_species_gas;
    result.state.particle_sidecar.owning_rank[global_index] = 0;
    result.state.particles.position_x_comoving[global_index] = x;
    result.state.particles.position_y_comoving[global_index] = 0.5;
    result.state.particles.position_z_comoving[global_index] = 0.5;
    result.state.particles.velocity_x_peculiar[global_index] = 10.0;
    result.state.particles.velocity_y_peculiar[global_index] = 0.0;
    result.state.particles.velocity_z_peculiar[global_index] = 0.0;
    result.state.particles.mass_code[global_index] = target_units.massSiToCode(1.0e34);
    result.state.particles.time_bin[global_index] = 0;

    result.state.cells.center_x_comoving[i] = x;
    result.state.cells.center_y_comoving[i] = 0.5;
    result.state.cells.center_z_comoving[i] = 0.5;
    result.state.cells.mass_code[i] = result.state.particles.mass_code[global_index];
    result.state.cells.time_bin[i] = 0;
    result.state.cells.patch_index[i] = 0;
    result.state.gas_cells.density_code[i] = 1.0;
    result.state.gas_cells.pressure_code[i] = 1.0;
    result.state.gas_cells.internal_energy_code[i] = 1.5;
    result.state.gas_cells.temperature_code[i] = 1.0e4;
    result.state.gas_cells.sound_speed_code[i] = 1.0;
  }

  for (std::size_t i = 0; i < dark_matter_particle_count; ++i, ++global_index) {
    result.state.particle_sidecar.particle_id[global_index] = particle_id_seed + global_index;
    result.state.particle_sidecar.species_tag[global_index] = k_species_dark_matter;
    result.state.particle_sidecar.owning_rank[global_index] = 0;
    result.state.particles.position_x_comoving[global_index] = static_cast<double>(i) * 0.01;
    result.state.particles.position_y_comoving[global_index] = 0.0;
    result.state.particles.position_z_comoving[global_index] = 0.0;
    result.state.particles.velocity_x_peculiar[global_index] = 0.0;
    result.state.particles.velocity_y_peculiar[global_index] = 0.0;
    result.state.particles.velocity_z_peculiar[global_index] = 0.0;
    result.state.particles.mass_code[global_index] = target_units.massSiToCode(1.0e36);
    result.state.particles.time_bin[global_index] = 0;
  }

  fillSpeciesLedger(result.state);
  result.state.rebuildSpeciesIndex();
  if (gas_particle_count > 0) {
    result.state.refreshGasCellIdentityFromParticleOrder();
  }
  return result;
}

IcReadResult convertGeneratedIsolatedIcToState(
    const core::SimulationConfig& config,
    std::size_t particle_count_per_axis) {
  const std::size_t dm_count = particle_count_per_axis * particle_count_per_axis;
  const std::size_t gas_count = particle_count_per_axis;
  IcReadResult result = buildGeneratedIsolatedIc(config, dm_count, gas_count);
  result.report.defaulted_fields.push_back("gas/internal_energy_code=uniform_default");
  return result;
}

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
  if (options.chunk_particle_count == 0) {
    throw std::runtime_error("chunk_particle_count must be positive");
  }

  IcReadResult result;
  Hdf5Handle file(H5Fopen(ic_path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed to open IC file: " + ic_path.string());
  }

  Hdf5Handle header_group(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
  if (!header_group.valid()) {
    throw std::runtime_error("IC file missing /Header group");
  }
  readHeaderArrays(header_group.get(), result.report.schema);
  if (!(result.report.schema.scale_factor > 0.0)) {
    throw std::runtime_error("Header/Time (scale factor) must be positive");
  }

  // The recognized bridge is a versioned convention, not an inference from a
  // filename: kpc, Msun, km/s, comoving coordinates, physical peculiar
  // velocities, and zero h/a exponents. Callers can replace every convention
  // through an explicit validated IcManifest.
  const core::UnitSystem source_units = core::makeUnitSystem("kpc", "msun", "km_s");
  const core::UnitSystem target_units =
      core::makeUnitSystem(config.units.length_unit, config.units.mass_unit, config.units.velocity_unit);
  IcManifest manifest = options.manifest != nullptr
      ? *options.manifest
      : makeGadgetArepoBridgeV1Manifest(ic_path, result.report.schema);
  if (manifest.source_provenance_ids.empty()) {
    manifest.source_provenance_ids.reserve(manifest.source_files.size());
    for (const auto& source_file : manifest.source_files) {
      manifest.source_provenance_ids.push_back(
          fileFnv1aProvenanceId(source_file));
    }
  }
  validateIcManifest(manifest);
  if (manifest.source_files.size() != 1U ||
      manifest.num_files_per_snapshot != 1U) {
    throw std::runtime_error(
        "multifile IC manifests require the distributed IC service; convert to canonical single-file CHUI IC or use a supported distributed build");
  }
  if (manifest.source_files.front().lexically_normal() !=
      ic_path.lexically_normal()) {
    throw std::runtime_error(
        "IcManifest source file does not match the requested IC path");
  }
  if (manifest.num_part_this_file.front() != result.report.schema.count_by_type ||
      manifest.num_part_total != result.report.schema.total_count_by_type ||
      manifest.mass_table != result.report.schema.mass_table ||
      manifest.num_part_total_high_word !=
          result.report.schema.total_count_high_word) {
    throw std::runtime_error(
        "IcManifest particle counts/high words/MassTable do not match the HDF5 header");
  }
  const auto nearly_equal = [](double lhs, double rhs) {
    return std::abs(lhs - rhs) <=
        1.0e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
  };
  const core::UnitSystem megaparsec_units =
      core::makeUnitSystem("mpc", "msun", "km_s");
  const IcFieldManifest& box_convention = requireFieldConvention(
      manifest, "/Header/BoxSize", IcFieldSemantics::kCoordinate);
  double box_size_mpc = result.report.schema.box_size *
      icStoredToSiMultiplier(
          box_convention, manifest.hubble_param, manifest.scale_factor) /
      megaparsec_units.length_si_per_code;
  if (box_convention.coordinate_frame == IcCoordinateFrame::kPhysical) {
    box_size_mpc = core::physicalToComovingLength(
        box_size_mpc, manifest.scale_factor);
  }
  if (!nearly_equal(manifest.scale_factor, config.numerics.a_begin) ||
      !nearly_equal(manifest.omega_matter, config.cosmology.omega_matter) ||
      !nearly_equal(manifest.omega_lambda, config.cosmology.omega_lambda) ||
      !nearly_equal(manifest.hubble_param, config.cosmology.hubble_param) ||
      !nearly_equal(box_size_mpc, config.cosmology.box_size_x_mpc_comoving) ||
      !nearly_equal(box_size_mpc, config.cosmology.box_size_y_mpc_comoving) ||
      !nearly_equal(box_size_mpc, config.cosmology.box_size_z_mpc_comoving)) {
    throw std::runtime_error(
        "IC manifest cosmology/BoxSize/start epoch does not match the frozen runtime configuration");
  }
  result.report.manifest = manifest;

  std::size_t total_count = 0;
  for (std::uint64_t count : result.report.schema.count_by_type) {
    total_count += static_cast<std::size_t>(count);
  }
  result.state.resizeParticles(total_count);
  const std::size_t gas_cell_count = static_cast<std::size_t>(result.report.schema.count_by_type[0]);
  result.state.resizeCells(gas_cell_count);
  const std::size_t black_hole_count =
      static_cast<std::size_t>(result.report.schema.count_by_type[5]);
  result.state.black_holes.resize(black_hole_count);
  result.state.metadata.run_name = config.output.run_name;
  result.state.metadata.scale_factor = result.report.schema.scale_factor;

  std::size_t global_offset = 0;
  for (std::size_t type_index = 0; type_index < result.report.schema.count_by_type.size(); ++type_index) {
    const std::size_t local_count = static_cast<std::size_t>(result.report.schema.count_by_type[type_index]);
    if (local_count == 0) {
      continue;
    }

    const std::string group_name = "/PartType" + std::to_string(type_index);
    Hdf5Handle group(H5Gopen2(file.get(), group_name.c_str(), H5P_DEFAULT));
    if (!group.valid()) {
      throw std::runtime_error("missing particle group: " + group_name);
    }

    const std::string coordinates_name = pickAlias(
        group.get(), {"Coordinates", "Position", "POS"}, result.report,
        group_name + "/Coordinates", true);
    const std::string velocities_name = pickAlias(
        group.get(), {"Velocities", "Velocity", "VEL"}, result.report,
        group_name + "/Velocities", options.require_velocities);
    const std::string ids_name = pickAlias(
        group.get(), {"ParticleIDs", "ParticleID", "ID"}, result.report,
        group_name + "/ParticleIDs", options.require_particle_ids);
    const std::string masses_name = pickAlias(
        group.get(), {"Masses", "Mass"}, result.report, group_name + "/Masses", false);
    std::string internal_energy_name;
    std::string density_name;
    std::string metallicity_name;
    std::string smoothing_length_name;
    std::string black_hole_mass_name;
    std::string black_hole_mdot_name;
    if (type_index == 0) {
      internal_energy_name = pickAlias(
          group.get(),
          gasFieldAliases("InternalEnergy"),
          result.report,
          group_name + "/InternalEnergy",
          false);
      density_name = pickAlias(
          group.get(),
          gasFieldAliases("Density"),
          result.report,
          group_name + "/Density",
          false);
      metallicity_name = pickAlias(
          group.get(),
          gasFieldAliases("Metallicity"),
          result.report,
          group_name + "/Metallicity",
          false);
      smoothing_length_name = pickAlias(
          group.get(),
          gasFieldAliases("SmoothingLength"),
          result.report,
          group_name + "/SmoothingLength",
          false);
    } else if (type_index == 5U) {
      black_hole_mass_name = pickAlias(
          group.get(), {"BH_Mass"}, result.report,
          group_name + "/BH_Mass", true);
      black_hole_mdot_name = pickAlias(
          group.get(), {"BH_Mdot"}, result.report,
          group_name + "/BH_Mdot", false);
    }

    std::vector<double> coordinates_chunk;
    std::vector<double> velocity_chunk;
    std::vector<double> mass_chunk;
    std::vector<double> internal_energy_chunk;
    std::vector<double> density_chunk;
    std::vector<double> black_hole_mass_chunk;
    std::vector<double> black_hole_mdot_chunk;
    std::vector<std::uint64_t> ids_chunk;

    for (std::size_t local_start = 0; local_start < local_count; local_start += options.chunk_particle_count) {
      const std::size_t chunk_count = std::min(options.chunk_particle_count, local_count - local_start);
      readDatasetChunk2d(group.get(), coordinates_name, local_start, chunk_count, coordinates_chunk);
      validateFiniteField(coordinates_chunk, group_name + "/" + coordinates_name);
      const IcFieldManifest& coordinate_convention = requireFieldConvention(
          manifest, group_name + "/Coordinates", IcFieldSemantics::kCoordinate);
      applyLengthFrameConversion(
          coordinates_chunk,
          result.report.schema.scale_factor,
          coordinate_convention.coordinate_frame == IcCoordinateFrame::kPhysical
              ? "physical"
              : "comoving",
          "comoving");
      applyManifestUnitConversion(
          coordinates_chunk,
          coordinate_convention,
          manifest,
          target_units.length_si_per_code);

      if (!velocities_name.empty()) {
        readDatasetChunk2d(group.get(), velocities_name, local_start, chunk_count, velocity_chunk);
        validateFiniteField(velocity_chunk, group_name + "/" + velocities_name);
        applyManifestUnitConversion(
            velocity_chunk,
            requireFieldConvention(
                manifest,
                group_name + "/Velocities",
                IcFieldSemantics::kVelocity),
            manifest,
            target_units.velocity_si_per_code);
      } else {
        velocity_chunk.assign(chunk_count * 3, 0.0);
        result.report.defaulted_fields.push_back(group_name + "/Velocities=zero");
      }

      if (!ids_name.empty()) {
        readDatasetChunkIds(group.get(), ids_name, local_start, chunk_count, ids_chunk);
      } else {
        ids_chunk.resize(chunk_count);
        for (std::size_t i = 0; i < chunk_count; ++i) {
          ids_chunk[i] = static_cast<std::uint64_t>(global_offset + local_start + i + 1);
        }
        result.report.defaulted_fields.push_back(group_name + "/ParticleIDs=generated");
      }

      if (!masses_name.empty()) {
        readDatasetChunk1d(group.get(), masses_name, local_start, chunk_count, mass_chunk);
        validateFiniteField(mass_chunk, group_name + "/" + masses_name);
        applyManifestUnitConversion(
            mass_chunk,
            requireFieldConvention(
                manifest,
                group_name + "/Masses",
                IcFieldSemantics::kExtensive),
            manifest,
            target_units.mass_si_per_code);
      } else {
        const double constant_mass = result.report.schema.mass_table[type_index];
        if (!(constant_mass > 0.0) || !options.allow_mass_table_fallback) {
          throw std::runtime_error(
              "missing Masses dataset and invalid/non-enabled MassTable fallback in " + group_name);
        }
        mass_chunk.assign(chunk_count, constant_mass);
        applyManifestUnitConversion(
            mass_chunk,
            requireFieldConvention(
                manifest, "/Header/MassTable",
                IcFieldSemantics::kExtensive),
            manifest,
            target_units.mass_si_per_code);
        result.report.defaulted_fields.push_back(group_name + "/Masses=MassTable");
      }

      if (type_index == 0) {
        if (!internal_energy_name.empty()) {
          readDatasetChunk1d(group.get(), internal_energy_name, local_start, chunk_count, internal_energy_chunk);
          validateFiniteField(internal_energy_chunk, group_name + "/" + internal_energy_name);
          applyManifestUnitConversion(
              internal_energy_chunk,
              requireFieldConvention(
                  manifest,
                  group_name + "/InternalEnergy",
                  IcFieldSemantics::kSpecific),
              manifest,
              target_units.velocity_si_per_code *
                  target_units.velocity_si_per_code);
        } else {
          internal_energy_chunk.assign(chunk_count, 0.0);
          result.report.defaulted_fields.push_back(group_name + "/InternalEnergy=zero");
        }

        if (!density_name.empty()) {
          readDatasetChunk1d(group.get(), density_name, local_start, chunk_count, density_chunk);
          validateFiniteField(density_chunk, group_name + "/" + density_name);
          applyManifestUnitConversion(
              density_chunk,
              requireFieldConvention(
                  manifest,
                  group_name + "/Density",
                  IcFieldSemantics::kIntensive),
              manifest,
              target_units.mass_si_per_code /
                  std::pow(target_units.length_si_per_code, 3.0));
        } else {
          density_chunk.assign(chunk_count, 0.0);
          result.report.defaulted_fields.push_back(group_name + "/Density=zero");
        }
      } else if (type_index == 5U) {
        readDatasetChunk1d(
            group.get(), black_hole_mass_name, local_start, chunk_count,
            black_hole_mass_chunk);
        validateFiniteField(
            black_hole_mass_chunk, group_name + "/" + black_hole_mass_name);
        applyManifestUnitConversion(
            black_hole_mass_chunk,
            requireFieldConvention(
                manifest, group_name + "/BH_Mass",
                IcFieldSemantics::kExtensive),
            manifest,
            target_units.mass_si_per_code);
        if (!black_hole_mdot_name.empty()) {
          readDatasetChunk1d(
              group.get(), black_hole_mdot_name, local_start, chunk_count,
              black_hole_mdot_chunk);
          validateFiniteField(
              black_hole_mdot_chunk,
              group_name + "/" + black_hole_mdot_name);
          applyManifestUnitConversion(
              black_hole_mdot_chunk,
              requireFieldConvention(
                  manifest, group_name + "/BH_Mdot",
                  IcFieldSemantics::kIntensive),
              manifest,
              target_units.mass_si_per_code / target_units.timeSiPerCode());
        } else {
          black_hole_mdot_chunk.assign(chunk_count, 0.0);
          result.report.defaulted_fields.push_back(
              group_name + "/BH_Mdot=zero");
        }
      }

      for (std::size_t local_i = 0; local_i < chunk_count; ++local_i) {
        const std::size_t global_i = global_offset + local_start + local_i;
        const std::size_t coord_offset = local_i * 3;
        result.state.particles.position_x_comoving[global_i] = coordinates_chunk[coord_offset + 0];
        result.state.particles.position_y_comoving[global_i] = coordinates_chunk[coord_offset + 1];
        result.state.particles.position_z_comoving[global_i] = coordinates_chunk[coord_offset + 2];

        result.state.particles.velocity_x_peculiar[global_i] = velocity_chunk[coord_offset + 0];
        result.state.particles.velocity_y_peculiar[global_i] = velocity_chunk[coord_offset + 1];
        result.state.particles.velocity_z_peculiar[global_i] = velocity_chunk[coord_offset + 2];

        result.state.particles.mass_code[global_i] = mass_chunk[local_i];
        result.state.particles.time_bin[global_i] = 0;

        result.state.particle_sidecar.particle_id[global_i] = ids_chunk[local_i];
        if (!(mass_chunk[local_i] > 0.0) ||
            !std::isfinite(mass_chunk[local_i])) {
          throw std::runtime_error(
              "IC particle masses must be finite and positive");
        }
        result.state.particle_sidecar.species_tag[global_i] =
            mapTypeIndexToSpeciesTag(type_index, manifest);
        result.state.particle_sidecar.owning_rank[global_i] = 0;

        if (type_index == 0) {
          const std::size_t cell_i = local_start + local_i;
          result.state.cells.center_x_comoving[cell_i] = coordinates_chunk[coord_offset + 0];
          result.state.cells.center_y_comoving[cell_i] = coordinates_chunk[coord_offset + 1];
          result.state.cells.center_z_comoving[cell_i] = coordinates_chunk[coord_offset + 2];
          result.state.cells.mass_code[cell_i] = mass_chunk[local_i];
          result.state.cells.time_bin[cell_i] = 0;
          result.state.cells.patch_index[cell_i] = 0;

          result.state.gas_cells.density_code[cell_i] = density_chunk[local_i];
          result.state.gas_cells.internal_energy_code[cell_i] = internal_energy_chunk[local_i];
          result.state.gas_cells.pressure_code[cell_i] = 0.0;
          result.state.gas_cells.temperature_code[cell_i] = 0.0;
          result.state.gas_cells.sound_speed_code[cell_i] = 0.0;
        } else if (type_index == 5U) {
          const std::size_t black_hole_row = local_start + local_i;
          if (!(black_hole_mass_chunk[local_i] > 0.0) ||
              !std::isfinite(black_hole_mass_chunk[local_i]) ||
              black_hole_mdot_chunk[local_i] < 0.0 ||
              !std::isfinite(black_hole_mdot_chunk[local_i])) {
            throw std::runtime_error(
                "PartType5 BH_Mass/BH_Mdot must be finite with positive mass and non-negative accretion rate");
          }
          result.state.black_holes.particle_index[black_hole_row] =
              static_cast<std::uint32_t>(global_i);
          result.state.black_holes.host_cell_index[black_hole_row] =
              std::numeric_limits<std::uint32_t>::max();
          result.state.black_holes.subgrid_mass_code[black_hole_row] =
              black_hole_mass_chunk[local_i];
          result.state.black_holes.accretion_rate_code[black_hole_row] =
              black_hole_mdot_chunk[local_i];
        }
      }
    }

    if (type_index == 0) {
      if (!metallicity_name.empty()) {
        result.report.unsupported_fields.push_back(
            "PartType0/Metallicity present but unmapped: SimulationState has no gas metallicity lane");
      }
      if (!smoothing_length_name.empty()) {
        result.report.unsupported_fields.push_back(
            "PartType0/SmoothingLength present but unmapped: SimulationState has no gas smoothing-length lane");
      }
    }

    global_offset += local_count;
  }

  fillSpeciesLedger(result.state);
  {
    std::vector<std::uint64_t> sorted_ids(
        result.state.particle_sidecar.particle_id.begin(),
        result.state.particle_sidecar.particle_id.end());
    std::sort(sorted_ids.begin(), sorted_ids.end());
    if (std::adjacent_find(sorted_ids.begin(), sorted_ids.end()) !=
        sorted_ids.end()) {
      throw std::runtime_error("IC import contains duplicate particle IDs");
    }
  }
  result.state.rebuildSpeciesIndex();
  if (gas_cell_count > 0) {
    result.state.refreshGasCellIdentityFromParticleOrder();
  }
  if (!result.state.validateOwnershipInvariants()) {
    throw std::runtime_error("IC import produced invalid ownership invariants");
  }

  return result;
#endif
}

}  // namespace cosmosim::io
