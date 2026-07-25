#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include <hdf5.h>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/io/ic_reader.hpp"

namespace {

class Hdf5Handle {
 public:
  explicit Hdf5Handle(hid_t handle = -1) : m_handle(handle) {}
  Hdf5Handle(const Hdf5Handle&) = delete;
  Hdf5Handle& operator=(const Hdf5Handle&) = delete;
  Hdf5Handle(Hdf5Handle&& other) noexcept : m_handle(other.m_handle) {
    other.m_handle = -1;
  }
  Hdf5Handle& operator=(Hdf5Handle&& other) noexcept {
    if (this != &other) {
      close();
      m_handle = other.m_handle;
      other.m_handle = -1;
    }
    return *this;
  }
  ~Hdf5Handle() { close(); }
  [[nodiscard]] hid_t get() const noexcept { return m_handle; }
  [[nodiscard]] bool valid() const noexcept { return m_handle >= 0; }

 private:
  void close() {
    if (m_handle < 0) return;
    switch (H5Iget_type(m_handle)) {
      case H5I_FILE: H5Fclose(m_handle); break;
      case H5I_GROUP: H5Gclose(m_handle); break;
      case H5I_DATASET: H5Dclose(m_handle); break;
      case H5I_DATASPACE: H5Sclose(m_handle); break;
      case H5I_ATTR: H5Aclose(m_handle); break;
      case H5I_DATATYPE: H5Tclose(m_handle); break;
      default: break;
    }
    m_handle = -1;
  }
  hid_t m_handle = -1;
};

struct Arguments {
  std::filesystem::path input;
  std::filesystem::path source_manifest;
  std::filesystem::path output;
  std::filesystem::path manifest_output;
  std::size_t chunk_particles = 65536U;
  double source_length_si = 0.0;
  double source_mass_si = 0.0;
  double source_velocity_si = 0.0;
  std::string coordinate_frame;
  std::string velocity_convention;
  double length_h = std::numeric_limits<double>::quiet_NaN();
  double length_a = std::numeric_limits<double>::quiet_NaN();
  double mass_h = std::numeric_limits<double>::quiet_NaN();
  double mass_a = std::numeric_limits<double>::quiet_NaN();
  double velocity_h = std::numeric_limits<double>::quiet_NaN();
  double velocity_a = std::numeric_limits<double>::quiet_NaN();
  std::string part_type2_policy = "reject";
  std::string part_type3_policy = "reject";
  std::string target_length_unit = "kpc";
  std::string target_mass_unit = "msun";
  std::string target_velocity_unit = "km_s";
};

[[nodiscard]] double parseDouble(std::string_view text, std::string_view name) {
  std::string owned(text);
  std::size_t consumed = 0;
  const double value = std::stod(owned, &consumed);
  if (consumed != owned.size() || !std::isfinite(value)) {
    throw std::invalid_argument(std::string(name) + " must be finite");
  }
  return value;
}

[[nodiscard]] std::size_t parseSize(
    std::string_view text, std::string_view name) {
  std::uint64_t value = 0;
  const auto [end, error] = std::from_chars(
      text.data(), text.data() + text.size(), value);
  if (error != std::errc{} || end != text.data() + text.size() || value == 0U ||
      value > std::numeric_limits<std::size_t>::max()) {
    throw std::invalid_argument(std::string(name) + " must be a positive size");
  }
  return static_cast<std::size_t>(value);
}

[[nodiscard]] std::string usage() {
  return
      "Usage: cosmosim_convert_ic --input FILE --output FILE --manifest FILE "
      "--source-convention gadget_arepo_bridge_v1 "
      "--source-length-unit-to-si X --source-mass-unit-to-si X "
      "--source-velocity-unit-to-si X --coordinate-frame MODE "
      "--velocity-convention MODE --length-h-exponent X "
      "--length-a-exponent X --mass-h-exponent X --mass-a-exponent X "
      "--velocity-h-exponent X --velocity-a-exponent X [options]\n"
      "   or: cosmosim_convert_ic --source-manifest FILE --output FILE "
      "--manifest FILE [options]\n";
}

[[nodiscard]] Arguments parseArguments(int argc, char** argv) {
  Arguments arguments;
  bool source_convention_selected = false;
  for (int i = 1; i < argc; ++i) {
    const std::string key = argv[i];
    const auto value = [&]() -> std::string {
      if (i + 1 >= argc) {
        throw std::invalid_argument("missing value for " + key);
      }
      return argv[++i];
    };
    if (key == "--help" || key == "-h") {
      std::cout << usage();
      std::exit(0);
    } else if (key == "--input") arguments.input = value();
    else if (key == "--source-manifest") arguments.source_manifest = value();
    else if (key == "--output") arguments.output = value();
    else if (key == "--manifest") arguments.manifest_output = value();
    else if (key == "--source-convention") {
      const std::string selected = value();
      if (selected != "gadget_arepo_bridge_v1") {
        throw std::invalid_argument("unsupported source convention: " + selected);
      }
      source_convention_selected = true;
    } else if (key == "--chunk-particles") {
      arguments.chunk_particles = parseSize(value(), key);
    } else if (key == "--source-length-unit-to-si") {
      arguments.source_length_si = parseDouble(value(), key);
    } else if (key == "--source-mass-unit-to-si") {
      arguments.source_mass_si = parseDouble(value(), key);
    } else if (key == "--source-velocity-unit-to-si") {
      arguments.source_velocity_si = parseDouble(value(), key);
    } else if (key == "--coordinate-frame") arguments.coordinate_frame = value();
    else if (key == "--velocity-convention") arguments.velocity_convention = value();
    else if (key == "--length-h-exponent") arguments.length_h = parseDouble(value(), key);
    else if (key == "--length-a-exponent") arguments.length_a = parseDouble(value(), key);
    else if (key == "--mass-h-exponent") arguments.mass_h = parseDouble(value(), key);
    else if (key == "--mass-a-exponent") arguments.mass_a = parseDouble(value(), key);
    else if (key == "--velocity-h-exponent") arguments.velocity_h = parseDouble(value(), key);
    else if (key == "--velocity-a-exponent") arguments.velocity_a = parseDouble(value(), key);
    else if (key == "--part-type2-policy") arguments.part_type2_policy = value();
    else if (key == "--part-type3-policy") arguments.part_type3_policy = value();
    else if (key == "--target-length-unit") arguments.target_length_unit = value();
    else if (key == "--target-mass-unit") arguments.target_mass_unit = value();
    else if (key == "--target-velocity-unit") arguments.target_velocity_unit = value();
    else if (key == "--target-length-unit-to-si" ||
             key == "--target-mass-unit-to-si" ||
             key == "--target-velocity-unit-to-si") {
      // Retained only as a compatibility check. CHUÍ target units remain typed.
      static_cast<void>(parseDouble(value(), key));
    } else {
      throw std::invalid_argument("unknown argument: " + key);
    }
  }
  if (arguments.output.empty() || arguments.manifest_output.empty()) {
    throw std::invalid_argument("--output and --manifest are required");
  }
  const bool manifest_mode = !arguments.source_manifest.empty();
  if (manifest_mode == source_convention_selected) {
    throw std::invalid_argument(
        "select exactly one of --source-manifest or --source-convention");
  }
  if (!manifest_mode) {
    if (arguments.input.empty()) {
      throw std::invalid_argument("--input is required for direct conversion");
    }
    const bool complete = arguments.source_length_si > 0.0 &&
        arguments.source_mass_si > 0.0 && arguments.source_velocity_si > 0.0 &&
        !arguments.coordinate_frame.empty() &&
        !arguments.velocity_convention.empty() &&
        std::isfinite(arguments.length_h) && std::isfinite(arguments.length_a) &&
        std::isfinite(arguments.mass_h) && std::isfinite(arguments.mass_a) &&
        std::isfinite(arguments.velocity_h) &&
        std::isfinite(arguments.velocity_a);
    if (!complete) {
      throw std::invalid_argument(
          "direct conversion requires the complete explicit source convention");
    }
  }
  return arguments;
}

[[nodiscard]] std::string configText(const Arguments& arguments, bool manifest_mode) {
  std::string text;
  text += "[mode]\nmode = zoom_in\n";
  text += "ic_file = input.hdf5\n";
  text += "ic_chunk_particle_count = " +
      std::to_string(arguments.chunk_particles) + "\n";
  text += "ic_part_type2_policy = " + arguments.part_type2_policy + "\n";
  text += "ic_part_type3_policy = " + arguments.part_type3_policy + "\n";
  if (manifest_mode) {
    text += "ic_convention = manifest_v1\n";
    text += "ic_manifest_file = supplied.audit.json\n";
  } else {
    text += "ic_convention = gadget_arepo_bridge_v1\n";
    text += "ic_bridge_source_length_unit_to_si = " +
        std::to_string(arguments.source_length_si) + "\n";
    text += "ic_bridge_source_mass_unit_to_si = " +
        std::to_string(arguments.source_mass_si) + "\n";
    text += "ic_bridge_source_velocity_unit_to_si = " +
        std::to_string(arguments.source_velocity_si) + "\n";
    text += "ic_bridge_coordinate_frame = " + arguments.coordinate_frame + "\n";
    text += "ic_bridge_velocity_convention = " +
        arguments.velocity_convention + "\n";
    text += "ic_bridge_length_hubble_exponent = " +
        std::to_string(arguments.length_h) + "\n";
    text += "ic_bridge_length_scale_factor_exponent = " +
        std::to_string(arguments.length_a) + "\n";
    text += "ic_bridge_mass_hubble_exponent = " +
        std::to_string(arguments.mass_h) + "\n";
    text += "ic_bridge_mass_scale_factor_exponent = " +
        std::to_string(arguments.mass_a) + "\n";
    text += "ic_bridge_velocity_hubble_exponent = " +
        std::to_string(arguments.velocity_h) + "\n";
    text += "ic_bridge_velocity_scale_factor_exponent = " +
        std::to_string(arguments.velocity_a) + "\n";
  }
  text += "[units]\nlength_unit = " + arguments.target_length_unit + "\n";
  text += "mass_unit = " + arguments.target_mass_unit + "\n";
  text += "velocity_unit = " + arguments.target_velocity_unit + "\n";
  return text;
}

void writeAttributeU32(hid_t group, const char* name, std::uint32_t value) {
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attribute(H5Acreate2(
      group, name, H5T_STD_U32LE, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attribute.valid() ||
      H5Awrite(attribute.get(), H5T_NATIVE_UINT32, &value) < 0) {
    throw std::runtime_error(std::string("failed to write Header/") + name);
  }
}

void writeAttributeF64(hid_t group, const char* name, double value) {
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attribute(H5Acreate2(
      group, name, H5T_IEEE_F64LE, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attribute.valid() ||
      H5Awrite(attribute.get(), H5T_NATIVE_DOUBLE, &value) < 0) {
    throw std::runtime_error(std::string("failed to write Header/") + name);
  }
}

template <typename T>
void writeAttributeArray6(
    hid_t group, const char* name, hid_t file_type, hid_t memory_type,
    const std::array<T, 6>& values) {
  hsize_t dimensions[1]{6U};
  Hdf5Handle space(H5Screate_simple(1, dimensions, nullptr));
  Hdf5Handle attribute(H5Acreate2(
      group, name, file_type, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attribute.valid() ||
      H5Awrite(attribute.get(), memory_type, values.data()) < 0) {
    throw std::runtime_error(std::string("failed to write Header/") + name);
  }
}

void writeAttributeString(
    hid_t group, const char* name, const std::string& value) {
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle type(H5Tcopy(H5T_C_S1));
  if (!type.valid() || H5Tset_size(type.get(), value.size() + 1U) < 0 ||
      H5Tset_strpad(type.get(), H5T_STR_NULLTERM) < 0) {
    throw std::runtime_error("failed to construct canonical string type");
  }
  Hdf5Handle attribute(H5Acreate2(
      group, name, type.get(), space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attribute.valid() ||
      H5Awrite(attribute.get(), type.get(), value.c_str()) < 0) {
    throw std::runtime_error(std::string("failed to write Header/") + name);
  }
}

template <typename T>
void writeDataset1d(
    hid_t group, const char* name, hid_t file_type, hid_t memory_type,
    const std::vector<T>& values) {
  hsize_t dimensions[1]{static_cast<hsize_t>(values.size())};
  Hdf5Handle space(H5Screate_simple(1, dimensions, nullptr));
  Hdf5Handle dataset(H5Dcreate2(
      group, name, file_type, space.get(), H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT));
  if (!dataset.valid() ||
      (!values.empty() && H5Dwrite(
          dataset.get(), memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          values.data()) < 0)) {
    throw std::runtime_error(std::string("failed to write dataset ") + name);
  }
}

void writeDatasetVec3(
    hid_t group, const char* name, const std::vector<double>& values) {
  if (values.size() % 3U != 0U) {
    throw std::logic_error("canonical vector dataset is not N x 3");
  }
  hsize_t dimensions[2]{
      static_cast<hsize_t>(values.size() / 3U), 3U};
  Hdf5Handle space(H5Screate_simple(2, dimensions, nullptr));
  Hdf5Handle dataset(H5Dcreate2(
      group, name, H5T_IEEE_F64LE, space.get(), H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT));
  if (!dataset.valid() ||
      (!values.empty() && H5Dwrite(
          dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          values.data()) < 0)) {
    throw std::runtime_error(std::string("failed to write dataset ") + name);
  }
}

[[nodiscard]] std::size_t canonicalType(std::uint32_t species) {
  using cosmosim::core::ParticleSpecies;
  if (species == static_cast<std::uint32_t>(ParticleSpecies::kGas)) return 0U;
  if (species == static_cast<std::uint32_t>(ParticleSpecies::kDarkMatter)) return 1U;
  if (species == static_cast<std::uint32_t>(ParticleSpecies::kStar)) return 4U;
  if (species == static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) return 5U;
  throw std::runtime_error(
      "canonical CHUÍ v1 conversion does not yet define a tracer output family");
}

void writeCanonicalFile(
    const std::filesystem::path& output,
    const cosmosim::io::IcReadResult& result,
    const cosmosim::core::SimulationConfig& config,
    const std::string& manifest_digest) {
  const std::filesystem::path partial = output.string() + ".part";
  std::filesystem::remove(partial);
  Hdf5Handle file(H5Fcreate(
      partial.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  if (!file.valid()) throw std::runtime_error("failed to create canonical output");
  Hdf5Handle header(H5Gcreate2(
      file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  if (!header.valid()) throw std::runtime_error("failed to create /Header");

  std::array<std::uint64_t, 6> counts64{};
  std::array<std::vector<std::size_t>, 6> indices;
  for (std::size_t particle = 0; particle < result.state.particles.size();
       ++particle) {
    const std::size_t type =
        canonicalType(result.state.particle_sidecar.species_tag[particle]);
    indices[type].push_back(particle);
    ++counts64[type];
  }
  std::array<std::uint32_t, 6> low{};
  std::array<std::uint32_t, 6> high{};
  for (std::size_t type = 0; type < 6; ++type) {
    low[type] = static_cast<std::uint32_t>(counts64[type] & 0xffffffffULL);
    high[type] = static_cast<std::uint32_t>(counts64[type] >> 32U);
  }
  writeAttributeArray6(
      header.get(), "NumPart_ThisFile", H5T_STD_U32LE,
      H5T_NATIVE_UINT32, low);
  writeAttributeArray6(
      header.get(), "NumPart_Total", H5T_STD_U32LE,
      H5T_NATIVE_UINT32, low);
  writeAttributeArray6(
      header.get(), "NumPart_Total_HighWord", H5T_STD_U32LE,
      H5T_NATIVE_UINT32, high);
  writeAttributeArray6(
      header.get(), "MassTable", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
      std::array<double, 6>{});
  const auto& manifest = *result.report.manifest;
  const auto target_units = cosmosim::core::makeUnitSystem(
      config.units.length_unit, config.units.mass_unit,
      config.units.velocity_unit);
  const auto box_field = std::find_if(
      manifest.fields.begin(), manifest.fields.end(), [](const auto& field) {
        return field.source_file_index == 0U &&
            field.dataset_path == "/Header/BoxSize";
      });
  if (box_field == manifest.fields.end()) {
    throw std::runtime_error("source manifest lacks /Header/BoxSize");
  }
  const double box_code = manifest.box_size *
      cosmosim::io::icFieldConversionMultiplier(
          *box_field, manifest, target_units);
  writeAttributeF64(header.get(), "Time", manifest.scale_factor);
  writeAttributeF64(header.get(), "Redshift", manifest.redshift);
  writeAttributeF64(header.get(), "BoxSize", box_code);
  writeAttributeF64(header.get(), "Omega0", manifest.omega_matter);
  writeAttributeF64(header.get(), "OmegaLambda", manifest.omega_lambda);
  writeAttributeF64(header.get(), "HubbleParam", manifest.hubble_param);
  writeAttributeU32(header.get(), "NumFilesPerSnapshot", 1U);
  writeAttributeString(header.get(), "ChuiIcSchemaName", "chui_canonical_v1");
  writeAttributeU32(header.get(), "ChuiIcSchemaVersion", 1U);
  writeAttributeF64(
      header.get(), "ChuiLengthUnitToSI", target_units.length_si_per_code);
  writeAttributeF64(
      header.get(), "ChuiMassUnitToSI", target_units.mass_si_per_code);
  writeAttributeF64(
      header.get(), "ChuiVelocityUnitToSI", target_units.velocity_si_per_code);
  writeAttributeString(header.get(), "ChuiCoordinateFrame", "comoving");
  writeAttributeString(
      header.get(), "ChuiVelocityConvention", "physical_peculiar");
  writeAttributeString(
      header.get(), "ConversionManifestSha256", manifest_digest);

  std::unordered_map<std::uint32_t, std::size_t> star_rows;
  for (std::size_t row = 0; row < result.state.star_particles.size(); ++row) {
    star_rows.emplace(result.state.star_particles.particle_index[row], row);
  }
  std::unordered_map<std::uint32_t, std::size_t> bh_rows;
  for (std::size_t row = 0; row < result.state.black_holes.size(); ++row) {
    bh_rows.emplace(result.state.black_holes.particle_index[row], row);
  }
  std::unordered_map<std::uint64_t, std::size_t> gas_rows;
  for (std::size_t row = 0; row < result.state.gas_cells.size(); ++row) {
    gas_rows.emplace(result.state.gas_cells.parent_particle_id[row], row);
  }

  for (std::size_t type : {0U, 1U, 4U, 5U}) {
    if (indices[type].empty()) continue;
    const std::string group_name = "/PartType" + std::to_string(type);
    Hdf5Handle group(H5Gcreate2(
        file.get(), group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT,
        H5P_DEFAULT));
    std::vector<double> coordinates;
    std::vector<double> velocities;
    std::vector<double> masses;
    std::vector<std::uint64_t> ids;
    coordinates.reserve(indices[type].size() * 3U);
    velocities.reserve(indices[type].size() * 3U);
    masses.reserve(indices[type].size());
    ids.reserve(indices[type].size());
    for (const std::size_t particle : indices[type]) {
      coordinates.insert(coordinates.end(), {
          result.state.particles.position_x_comoving[particle],
          result.state.particles.position_y_comoving[particle],
          result.state.particles.position_z_comoving[particle]});
      velocities.insert(velocities.end(), {
          result.state.particles.velocity_x_peculiar[particle],
          result.state.particles.velocity_y_peculiar[particle],
          result.state.particles.velocity_z_peculiar[particle]});
      masses.push_back(result.state.particles.mass_code[particle]);
      ids.push_back(result.state.particle_sidecar.particle_id[particle]);
    }
    writeDatasetVec3(group.get(), "Coordinates", coordinates);
    writeDatasetVec3(group.get(), "Velocities", velocities);
    writeDataset1d(
        group.get(), "Masses", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
        masses);
    writeDataset1d(
        group.get(), "ParticleIDs", H5T_STD_U64LE, H5T_NATIVE_UINT64,
        ids);

    if (type == 0U) {
      std::vector<double> internal_energy;
      std::vector<double> density;
      for (const std::size_t particle : indices[type]) {
        const auto row = gas_rows.find(
            result.state.particle_sidecar.particle_id[particle]);
        if (row == gas_rows.end()) {
          throw std::runtime_error("gas particle lacks canonical sidecar row");
        }
        internal_energy.push_back(
            result.state.gas_cells.internal_energy_code[row->second]);
        density.push_back(result.state.gas_cells.density_code[row->second]);
      }
      writeDataset1d(
          group.get(), "InternalEnergy", H5T_IEEE_F64LE,
          H5T_NATIVE_DOUBLE, internal_energy);
      writeDataset1d(
          group.get(), "Density", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
          density);
    } else if (type == 4U) {
      std::vector<double> formation;
      std::vector<double> initial_mass;
      std::vector<double> metallicity;
      for (const std::size_t particle : indices[type]) {
        const auto row = star_rows.find(static_cast<std::uint32_t>(particle));
        if (row == star_rows.end()) {
          throw std::runtime_error("star particle lacks canonical sidecar row");
        }
        formation.push_back(
            result.state.star_particles.formation_scale_factor[row->second]);
        initial_mass.push_back(
            result.state.star_particles.birth_mass_code[row->second]);
        metallicity.push_back(
            result.state.star_particles.metallicity_mass_fraction[row->second]);
      }
      writeDataset1d(
          group.get(), "StellarFormationTime", H5T_IEEE_F64LE,
          H5T_NATIVE_DOUBLE, formation);
      writeDataset1d(
          group.get(), "InitialMass", H5T_IEEE_F64LE,
          H5T_NATIVE_DOUBLE, initial_mass);
      writeDataset1d(
          group.get(), "Metallicity", H5T_IEEE_F64LE,
          H5T_NATIVE_DOUBLE, metallicity);
    } else if (type == 5U) {
      std::vector<double> subgrid_mass;
      std::vector<double> accretion_rate;
      for (const std::size_t particle : indices[type]) {
        const auto row = bh_rows.find(static_cast<std::uint32_t>(particle));
        if (row == bh_rows.end()) {
          throw std::runtime_error(
              "black-hole particle lacks canonical sidecar row");
        }
        subgrid_mass.push_back(
            result.state.black_holes.subgrid_mass_code[row->second]);
        accretion_rate.push_back(
            result.state.black_holes.accretion_rate_code[row->second]);
      }
      writeDataset1d(
          group.get(), "BH_Mass", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
          subgrid_mass);
      writeDataset1d(
          group.get(), "BH_Mdot", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
          accretion_rate);
    }
  }
  file = Hdf5Handle{};
  std::filesystem::remove(output);
  std::filesystem::rename(partial, output);
}

int run(const Arguments& arguments) {
  const bool manifest_mode = !arguments.source_manifest.empty();
  cosmosim::io::IcManifest supplied_manifest;
  const cosmosim::io::IcManifest* manifest_pointer = nullptr;
  std::filesystem::path source = arguments.input;
  if (manifest_mode) {
    supplied_manifest =
        cosmosim::io::readIcManifestJson(arguments.source_manifest);
    if (supplied_manifest.source_files.empty()) {
      throw std::runtime_error("source manifest contains no source files");
    }
    source = supplied_manifest.source_files.front();
    manifest_pointer = &supplied_manifest;
  }
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      configText(arguments, manifest_mode), "cosmosim_convert_ic");
  cosmosim::io::IcImportOptions options;
  options.chunk_particle_count = arguments.chunk_particles;
  options.manifest = manifest_pointer;
  options.validate_runtime_cosmology = false;
  cosmosim::io::IcReadResult result =
      cosmosim::io::readGadgetArepoHdf5Ic(source, frozen.config, options);
  if (!result.report.manifest.has_value()) {
    throw std::logic_error("IC converter did not receive an audit manifest");
  }
  result.report.manifest->converter_version = "chui_ic_converter_v3";
  if (manifest_mode) {
    result.report.manifest->source_manifest_file =
        std::filesystem::absolute(arguments.source_manifest).lexically_normal().string();
    result.report.manifest->source_manifest_sha256 =
        cosmosim::io::icSha256FileHex(arguments.source_manifest);
  }
  const std::string manifest_json =
      cosmosim::io::serializeIcManifestJson(*result.report.manifest);
  const std::string digest = cosmosim::io::icSha256Hex(manifest_json);
  writeCanonicalFile(
      arguments.output, result, frozen.config, digest);

  const std::filesystem::path manifest_partial =
      arguments.manifest_output.string() + ".part";
  std::filesystem::remove(manifest_partial);
  cosmosim::io::writeIcManifestJson(
      *result.report.manifest, manifest_partial);
  std::filesystem::remove(arguments.manifest_output);
  std::filesystem::rename(manifest_partial, arguments.manifest_output);
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    return run(parseArguments(argc, argv));
  } catch (const std::exception& error) {
    std::cerr << "cosmosim_convert_ic: " << error.what() << '\n';
    return 1;
  }
}
