#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <limits>
#include <memory>
#include <optional>
#include <queue>
#include <span>
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
#include "io/internal/ic_stream_ingestion.hpp"
#include "io/internal/ic_canonical_limits.hpp"

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
  std::string gas_internal_energy_policy = "reject";
  std::string gas_density_policy = "reject";
  std::string star_formation_time_policy = "reject";
  std::string star_initial_mass_policy = "reject";
  std::string star_metallicity_policy = "reject";
  std::string bh_mdot_policy = "reject";
  double gas_internal_energy_value_code = 0.0;
  double gas_density_value_code = 0.0;
  double star_formation_time_value = 0.0;
  double star_initial_mass_value_code = 0.0;
  double star_metallicity_value = 0.0;
  double bh_mdot_value_code = 0.0;
  bool policy_argument_supplied = false;
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
      "--manifest FILE [options]\n"
      "Missing-field options: --gas-internal-energy-policy POLICY "
      "--gas-density-policy POLICY --star-formation-time-policy POLICY "
      "--star-initial-mass-policy POLICY --star-metallicity-policy POLICY "
      "--bh-mdot-policy POLICY, with matching --*-value options for "
      "use_config_value.\n";
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
    else if (key == "--part-type2-policy") {
      arguments.part_type2_policy = value();
      arguments.policy_argument_supplied = true;
    } else if (key == "--part-type3-policy") {
      arguments.part_type3_policy = value();
      arguments.policy_argument_supplied = true;
    }
    else if (key == "--gas-internal-energy-policy") {
      arguments.gas_internal_energy_policy = value();
      arguments.policy_argument_supplied = true;
    } else if (key == "--gas-density-policy") {
      arguments.gas_density_policy = value();
      arguments.policy_argument_supplied = true;
    } else if (key == "--star-formation-time-policy") {
      arguments.star_formation_time_policy = value();
      arguments.policy_argument_supplied = true;
    } else if (key == "--star-initial-mass-policy") {
      arguments.star_initial_mass_policy = value();
      arguments.policy_argument_supplied = true;
    } else if (key == "--star-metallicity-policy") {
      arguments.star_metallicity_policy = value();
      arguments.policy_argument_supplied = true;
    } else if (key == "--bh-mdot-policy") {
      arguments.bh_mdot_policy = value();
      arguments.policy_argument_supplied = true;
    } else if (key == "--gas-internal-energy-value-code") {
      arguments.gas_internal_energy_value_code = parseDouble(value(), key);
      arguments.policy_argument_supplied = true;
    } else if (key == "--gas-density-value-code") {
      arguments.gas_density_value_code = parseDouble(value(), key);
      arguments.policy_argument_supplied = true;
    } else if (key == "--star-formation-time-value") {
      arguments.star_formation_time_value = parseDouble(value(), key);
      arguments.policy_argument_supplied = true;
    } else if (key == "--star-initial-mass-value-code") {
      arguments.star_initial_mass_value_code = parseDouble(value(), key);
      arguments.policy_argument_supplied = true;
    } else if (key == "--star-metallicity-value") {
      arguments.star_metallicity_value = parseDouble(value(), key);
      arguments.policy_argument_supplied = true;
    } else if (key == "--bh-mdot-value-code") {
      arguments.bh_mdot_value_code = parseDouble(value(), key);
      arguments.policy_argument_supplied = true;
    } else if (key == "--target-length-unit") arguments.target_length_unit = value();
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
  if (manifest_mode && arguments.policy_argument_supplied) {
    throw std::invalid_argument(
        "manifest mode derives species and missing-field policies from the "
        "manifest; policy/value CLI overrides are not allowed");
  }
  if (arguments.part_type2_policy == "tracer" ||
      arguments.part_type3_policy == "tracer") {
    throw std::invalid_argument(
        "canonical CHUÍ v1 has no tracer output family; tracer policy is unsupported");
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

[[nodiscard]] std::string missingFieldPolicyName(
    cosmosim::io::IcMissingFieldPolicy policy) {
  using Policy = cosmosim::io::IcMissingFieldPolicy;
  switch (policy) {
    case Policy::kReject: return "reject";
    case Policy::kReconstruct: return "reconstruct";
    case Policy::kUseConfigValue: return "use_config_value";
    case Policy::kDialectDefinedDefault: return "dialect_defined_default";
    case Policy::kPreserveUnavailable: return "preserve_unavailable";
  }
  throw std::invalid_argument("unknown manifest missing-field policy");
}

void applyManifestMissingFieldPolicies(
    Arguments& arguments,
    const cosmosim::io::IcManifest& manifest) {
  struct SelectedPolicy {
    bool assigned = false;
    cosmosim::io::IcMissingFieldPolicy policy =
        cosmosim::io::IcMissingFieldPolicy::kReject;
    double configured_value_code = 0.0;
  };
  std::array<SelectedPolicy, 6> selected{};
  const auto index_for_path = [](std::string_view path) -> std::size_t {
    if (path.ends_with("/InternalEnergy")) return 0U;
    if (path.ends_with("/Density")) return 1U;
    if (path.ends_with("/StellarFormationTime")) return 2U;
    if (path.ends_with("/InitialMass")) return 3U;
    if (path.ends_with("/Metallicity")) return 4U;
    if (path.ends_with("/BH_Mdot")) return 5U;
    throw std::invalid_argument(
        "manifest contains an unsupported missing-field contract: " +
        std::string(path));
  };
  for (const auto& contract : manifest.missing_field_contracts) {
    SelectedPolicy& destination = selected[index_for_path(contract.field_path)];
    if (destination.assigned &&
        (destination.policy != contract.policy ||
         destination.configured_value_code != contract.configured_value_code)) {
      throw std::invalid_argument(
          "manifest contains conflicting missing-field policies for " +
          contract.field_path);
    }
    destination.assigned = true;
    destination.policy = contract.policy;
    destination.configured_value_code = contract.configured_value_code;
  }
  const auto policy_name = [&](std::size_t index) {
    return selected[index].assigned
        ? missingFieldPolicyName(selected[index].policy)
        : std::string("reject");
  };
  arguments.gas_internal_energy_policy = policy_name(0U);
  arguments.gas_density_policy = policy_name(1U);
  arguments.star_formation_time_policy = policy_name(2U);
  arguments.star_initial_mass_policy = policy_name(3U);
  arguments.star_metallicity_policy = policy_name(4U);
  arguments.bh_mdot_policy = policy_name(5U);
  arguments.gas_internal_energy_value_code =
      selected[0].configured_value_code;
  arguments.gas_density_value_code = selected[1].configured_value_code;
  arguments.star_formation_time_value = selected[2].configured_value_code;
  arguments.star_initial_mass_value_code = selected[3].configured_value_code;
  arguments.star_metallicity_value = selected[4].configured_value_code;
  arguments.bh_mdot_value_code = selected[5].configured_value_code;
}

[[nodiscard]] std::string configText(const Arguments& arguments, bool manifest_mode) {
  std::string text;
  text += "[mode]\nmode = zoom_in\n";
  if (!manifest_mode) text += "ic_file = input.hdf5\n";
  text += "ic_chunk_particle_count = " +
      std::to_string(arguments.chunk_particles) + "\n";
  text += "ic_part_type2_policy = " + arguments.part_type2_policy + "\n";
  text += "ic_part_type3_policy = " + arguments.part_type3_policy + "\n";
  text += "ic_gas_internal_energy_policy = " +
      arguments.gas_internal_energy_policy + "\n";
  text += "ic_gas_density_policy = " + arguments.gas_density_policy + "\n";
  text += "ic_star_formation_time_policy = " +
      arguments.star_formation_time_policy + "\n";
  text += "ic_star_initial_mass_policy = " +
      arguments.star_initial_mass_policy + "\n";
  text += "ic_star_metallicity_policy = " +
      arguments.star_metallicity_policy + "\n";
  text += "ic_bh_mdot_policy = " + arguments.bh_mdot_policy + "\n";
  text += "ic_gas_internal_energy_value_code = " +
      std::to_string(arguments.gas_internal_energy_value_code) + "\n";
  text += "ic_gas_density_value_code = " +
      std::to_string(arguments.gas_density_value_code) + "\n";
  text += "ic_star_formation_time_value = " +
      std::to_string(arguments.star_formation_time_value) + "\n";
  text += "ic_star_initial_mass_value_code = " +
      std::to_string(arguments.star_initial_mass_value_code) + "\n";
  text += "ic_star_metallicity_value = " +
      std::to_string(arguments.star_metallicity_value) + "\n";
  text += "ic_bh_mdot_value_code = " +
      std::to_string(arguments.bh_mdot_value_code) + "\n";
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

void writeAttributeU64(hid_t group, const char* name, std::uint64_t value) {
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attribute(H5Acreate2(
      group, name, H5T_STD_U64LE, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attribute.valid() ||
      H5Awrite(attribute.get(), H5T_NATIVE_UINT64, &value) < 0) {
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

void writeByteDataset(
    hid_t group, const char* name, const std::string& value) {
  hsize_t dimensions[1]{static_cast<hsize_t>(value.size())};
  Hdf5Handle space(H5Screate_simple(1, dimensions, nullptr));
  Hdf5Handle dataset(H5Dcreate2(
      group, name, H5T_STD_U8LE, space.get(), H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT));
  if (!dataset.valid() ||
      (!value.empty() && H5Dwrite(
          dataset.get(), H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          value.data()) < 0)) {
    throw std::runtime_error(std::string("failed to write dataset ") + name);
  }
}

void writeTextExact(
    const std::filesystem::path& path, const std::string& contents) {
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output) {
    throw std::runtime_error("failed to create " + path.string());
  }
  output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
  output.flush();
  if (!output) {
    throw std::runtime_error("failed to write " + path.string());
  }
}

[[nodiscard]] bool injectFault(std::string_view phase) {
  const char* value = std::getenv("COSMOSIM_CONVERTER_TEST_FAULT");
  return value != nullptr && phase == value;
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
      "canonical CHUÍ v1 conversion does not define a tracer output family");
}

[[nodiscard]] std::size_t canonicalTypeForPolicy(
    cosmosim::io::IcSpeciesPolicy policy) {
  using cosmosim::io::IcSpeciesPolicy;
  switch (policy) {
    case IcSpeciesPolicy::kGas:
      return 0U;
    case IcSpeciesPolicy::kDarkMatter:
    case IcSpeciesPolicy::kCollisionlessFamily2AsDarkMatter:
    case IcSpeciesPolicy::kCollisionlessFamily3AsDarkMatter:
      return 1U;
    case IcSpeciesPolicy::kStar:
      return 4U;
    case IcSpeciesPolicy::kBlackHole:
      return 5U;
    case IcSpeciesPolicy::kTracer:
      throw std::runtime_error(
          "canonical CHUÍ v1 conversion does not define a tracer output family");
    case IcSpeciesPolicy::kReject:
      break;
  }
  throw std::logic_error("populated source type has reject species policy");
}

[[nodiscard]] Hdf5Handle createDataset1d(
    hid_t group, const char* name, hid_t file_type, std::uint64_t count) {
  const hsize_t dimensions[1]{static_cast<hsize_t>(count)};
  Hdf5Handle space(H5Screate_simple(1, dimensions, nullptr));
  Hdf5Handle dataset(H5Dcreate2(
      group, name, file_type, space.get(), H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error(std::string("failed to create dataset ") + name);
  }
  return dataset;
}

[[nodiscard]] Hdf5Handle createDatasetVec3(
    hid_t group, const char* name, std::uint64_t count) {
  const hsize_t dimensions[2]{static_cast<hsize_t>(count), 3U};
  Hdf5Handle space(H5Screate_simple(2, dimensions, nullptr));
  Hdf5Handle dataset(H5Dcreate2(
      group, name, H5T_IEEE_F64LE, space.get(), H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error(std::string("failed to create dataset ") + name);
  }
  return dataset;
}

void writeDataset1dSlice(
    hid_t dataset,
    hid_t memory_type,
    std::uint64_t row_offset,
    std::span<const double> values) {
  if (values.empty()) return;
  Hdf5Handle file_space(H5Dget_space(dataset));
  const hsize_t offset[1]{static_cast<hsize_t>(row_offset)};
  const hsize_t extent[1]{static_cast<hsize_t>(values.size())};
  if (!file_space.valid() ||
      H5Sselect_hyperslab(
          file_space.get(), H5S_SELECT_SET, offset, nullptr, extent, nullptr) < 0) {
    throw std::runtime_error("failed to select canonical scalar hyperslab");
  }
  Hdf5Handle memory_space(H5Screate_simple(1, extent, nullptr));
  if (!memory_space.valid() ||
      H5Dwrite(
          dataset, memory_type, memory_space.get(), file_space.get(),
          H5P_DEFAULT, values.data()) < 0) {
    throw std::runtime_error("failed to append canonical scalar dataset");
  }
}

void writeDatasetU64Slice(
    hid_t dataset,
    std::uint64_t row_offset,
    std::span<const std::uint64_t> values) {
  if (values.empty()) return;
  Hdf5Handle file_space(H5Dget_space(dataset));
  const hsize_t offset[1]{static_cast<hsize_t>(row_offset)};
  const hsize_t extent[1]{static_cast<hsize_t>(values.size())};
  if (!file_space.valid() ||
      H5Sselect_hyperslab(
          file_space.get(), H5S_SELECT_SET, offset, nullptr, extent, nullptr) < 0) {
    throw std::runtime_error("failed to select canonical ID hyperslab");
  }
  Hdf5Handle memory_space(H5Screate_simple(1, extent, nullptr));
  if (!memory_space.valid() ||
      H5Dwrite(
          dataset, H5T_NATIVE_UINT64, memory_space.get(), file_space.get(),
          H5P_DEFAULT, values.data()) < 0) {
    throw std::runtime_error("failed to append canonical ID dataset");
  }
}

void writeDatasetVec3Slice(
    hid_t dataset,
    std::uint64_t row_offset,
    std::span<const double> values) {
  if (values.empty()) return;
  if (values.size() % 3U != 0U) {
    throw std::logic_error("canonical vector batch is not N x 3");
  }
  Hdf5Handle file_space(H5Dget_space(dataset));
  const hsize_t offset[2]{static_cast<hsize_t>(row_offset), 0U};
  const hsize_t extent[2]{static_cast<hsize_t>(values.size() / 3U), 3U};
  if (!file_space.valid() ||
      H5Sselect_hyperslab(
          file_space.get(), H5S_SELECT_SET, offset, nullptr, extent, nullptr) < 0) {
    throw std::runtime_error("failed to select canonical vector hyperslab");
  }
  Hdf5Handle memory_space(H5Screate_simple(2, extent, nullptr));
  if (!memory_space.valid() ||
      H5Dwrite(
          dataset, H5T_NATIVE_DOUBLE, memory_space.get(), file_space.get(),
          H5P_DEFAULT, values.data()) < 0) {
    throw std::runtime_error("failed to append canonical vector dataset");
  }
}

class ExternalIdAudit {
 public:
  explicit ExternalIdAudit(std::filesystem::path directory)
      : m_directory(std::move(directory)) {
    std::error_code error;
    std::filesystem::remove_all(m_directory, error);
    if (!std::filesystem::create_directories(m_directory)) {
      throw std::runtime_error(
          "failed to create bounded duplicate-ID audit directory");
    }
  }

  ExternalIdAudit(const ExternalIdAudit&) = delete;
  ExternalIdAudit& operator=(const ExternalIdAudit&) = delete;

  ~ExternalIdAudit() {
    std::error_code error;
    std::filesystem::remove_all(m_directory, error);
  }

  void addBatch(
      std::span<const cosmosim::io::internal::IcParticleRecord> records) {
    if (records.empty()) return;
    std::vector<std::uint64_t> ids;
    ids.reserve(records.size());
    for (const auto& record : records) ids.push_back(record.id);
    std::sort(ids.begin(), ids.end());
    if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) {
      throw std::runtime_error(
          "duplicate particle ID detected within a conversion batch");
    }
    const auto path = m_directory /
        ("run_0_" + std::to_string(m_next_run++) + ".bin");
    writeRun(path, ids);
    m_runs.push_back(path);
  }

  void finish() {
    std::size_t generation = 1U;
    while (m_runs.size() > 1U) {
      std::vector<std::filesystem::path> next;
      for (std::size_t begin = 0U; begin < m_runs.size();
           begin += kMergeFanIn) {
        const std::size_t end =
            std::min(begin + kMergeFanIn, m_runs.size());
        const auto output = m_directory /
            ("run_" + std::to_string(generation) + "_" +
             std::to_string(next.size()) + ".bin");
        mergeRuns(
            std::span<const std::filesystem::path>(
                m_runs.data() + begin, end - begin),
            output);
        for (std::size_t i = begin; i < end; ++i) {
          std::filesystem::remove(m_runs[i]);
        }
        next.push_back(output);
      }
      m_runs = std::move(next);
      ++generation;
    }
  }

 private:
  static constexpr std::size_t kMergeFanIn = 32U;

  struct RunReader {
    std::ifstream input;
    std::uint64_t current = 0U;
  };

  static void writeRun(
      const std::filesystem::path& path,
      std::span<const std::uint64_t> ids) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("failed to create ID audit run");
    output.write(
        reinterpret_cast<const char*>(ids.data()),
        static_cast<std::streamsize>(ids.size_bytes()));
    output.flush();
    if (!output) throw std::runtime_error("failed to write ID audit run");
  }

  static bool readNext(RunReader& reader) {
    reader.input.read(
        reinterpret_cast<char*>(&reader.current),
        static_cast<std::streamsize>(sizeof(reader.current)));
    if (reader.input.gcount() == 0 && reader.input.eof()) return false;
    if (reader.input.gcount() !=
        static_cast<std::streamsize>(sizeof(reader.current))) {
      throw std::runtime_error("truncated duplicate-ID audit run");
    }
    return true;
  }

  static void mergeRuns(
      std::span<const std::filesystem::path> inputs,
      const std::filesystem::path& output_path) {
    std::vector<RunReader> readers(inputs.size());
    using QueueEntry = std::pair<std::uint64_t, std::size_t>;
    std::priority_queue<
        QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
    for (std::size_t i = 0U; i < inputs.size(); ++i) {
      readers[i].input.open(inputs[i], std::ios::binary);
      if (!readers[i].input) {
        throw std::runtime_error("failed to open duplicate-ID audit run");
      }
      if (readNext(readers[i])) queue.emplace(readers[i].current, i);
    }
    std::ofstream output(output_path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("failed to create merged ID audit run");
    std::optional<std::uint64_t> previous;
    while (!queue.empty()) {
      const auto [id, reader_index] = queue.top();
      queue.pop();
      if (previous.has_value() && *previous == id) {
        throw std::runtime_error(
            "duplicate particle ID detected by bounded external audit");
      }
      output.write(
          reinterpret_cast<const char*>(&id),
          static_cast<std::streamsize>(sizeof(id)));
      if (!output) throw std::runtime_error("failed to merge ID audit runs");
      previous = id;
      if (readNext(readers[reader_index])) {
        queue.emplace(readers[reader_index].current, reader_index);
      }
    }
  }

  std::filesystem::path m_directory;
  std::vector<std::filesystem::path> m_runs;
  std::size_t m_next_run = 0U;
};

class CanonicalStreamingWriter {
 public:
  CanonicalStreamingWriter(
      const std::filesystem::path& output_partial,
      const cosmosim::io::IcManifest& manifest,
      const cosmosim::core::SimulationConfig& config,
      const std::string& manifest_json,
      const std::string& manifest_digest,
      const std::string& manifest_filename,
      const std::string& marker_filename,
      std::size_t chunk_capacity)
      : m_expected_counts(canonicalCounts(manifest)),
        m_chunk_capacity(chunk_capacity),
        m_file(H5Fcreate(
            output_partial.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT)),
        m_header(H5Gcreate2(
            m_file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) {
    if (!m_file.valid() || !m_header.valid()) {
      throw std::runtime_error("failed to create canonical streaming output");
    }
    writeHeader(
        manifest, config, manifest_json, manifest_digest, manifest_filename,
        marker_filename);
    createParticleDatasets();
  }

  void append(
      std::span<const cosmosim::io::internal::IcParticleRecord> records) {
    m_peak_batch_records = std::max(m_peak_batch_records, records.size());
    std::array<std::vector<const cosmosim::io::internal::IcParticleRecord*>, 6>
        by_type;
    for (const auto& record : records) {
      by_type[canonicalType(record.species)].push_back(&record);
    }
    for (std::size_t type : {0U, 1U, 4U, 5U}) {
      if (by_type[type].empty()) continue;
      appendType(type, by_type[type]);
    }
  }

  void finish(std::uint64_t reader_peak_staging_bytes) {
    for (std::size_t type : {0U, 1U, 4U, 5U}) {
      if (m_offsets[type] != m_expected_counts[type]) {
        throw std::runtime_error(
            "canonical streaming writer count disagrees with manifest");
      }
    }
    writeAttributeU64(
        m_header.get(), "ChuiConverterPeakBatchRecords",
        static_cast<std::uint64_t>(m_peak_batch_records));
    writeAttributeU64(
        m_header.get(), "ChuiConverterChunkParticleCapacity",
        static_cast<std::uint64_t>(m_chunk_capacity));
    writeAttributeU64(
        m_header.get(), "ChuiConverterPeakReaderStagingBytes",
        reader_peak_staging_bytes);
    writeAttributeU32(
        m_header.get(), "ChuiConverterFullStateMaterialized", 0U);
    writeAttributeString(
        m_header.get(), "ChuiConverterFlow",
        "source_chunk_to_validate_convert_append");
    m_header = Hdf5Handle{};
    m_file = Hdf5Handle{};
  }

 private:
  struct DatasetSet {
    Hdf5Handle group;
    Hdf5Handle coordinates;
    Hdf5Handle velocities;
    Hdf5Handle masses;
    Hdf5Handle ids;
    Hdf5Handle first_sidecar;
    Hdf5Handle second_sidecar;
    Hdf5Handle third_sidecar;
  };

  static std::array<std::uint64_t, 6> canonicalCounts(
      const cosmosim::io::IcManifest& manifest) {
    std::array<std::uint64_t, 6> counts{};
    for (std::size_t source_type = 0U; source_type < 6U; ++source_type) {
      const std::uint64_t count = manifest.num_part_total[source_type];
      if (count == 0U) continue;
      const std::size_t target =
          canonicalTypeForPolicy(manifest.species_policy[source_type]);
      if (counts[target] > std::numeric_limits<std::uint64_t>::max() - count) {
        throw std::overflow_error("canonical particle-count overflow");
      }
      counts[target] += count;
    }
    cosmosim::io::internal::validateCanonicalSingleFileCounts(counts);
    return counts;
  }

  void writeHeader(
      const cosmosim::io::IcManifest& manifest,
      const cosmosim::core::SimulationConfig& config,
      const std::string& manifest_json,
      const std::string& manifest_digest,
      const std::string& manifest_filename,
      const std::string& marker_filename) {
    std::array<std::uint32_t, 6> low{};
    std::array<std::uint32_t, 6> high{};
    for (std::size_t type = 0U; type < 6U; ++type) {
      low[type] = static_cast<std::uint32_t>(
          m_expected_counts[type] & 0xffffffffULL);
      high[type] = static_cast<std::uint32_t>(
          m_expected_counts[type] >> 32U);
    }
    writeAttributeArray6(
        m_header.get(), "NumPart_ThisFile", H5T_STD_U32LE,
        H5T_NATIVE_UINT32, low);
    writeAttributeArray6(
        m_header.get(), "NumPart_Total", H5T_STD_U32LE,
        H5T_NATIVE_UINT32, low);
    writeAttributeArray6(
        m_header.get(), "NumPart_Total_HighWord", H5T_STD_U32LE,
        H5T_NATIVE_UINT32, high);
    writeAttributeArray6(
        m_header.get(), "MassTable", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
        std::array<double, 6>{});

    const auto target_units = cosmosim::core::makeUnitSystem(
        config.units.length_unit, config.units.mass_unit,
        config.units.velocity_unit);
    const auto box_field = std::find_if(
        manifest.fields.begin(), manifest.fields.end(),
        [](const cosmosim::io::IcFieldManifest& field) {
          return field.dataset_path == "/Header/BoxSize" &&
              field.disposition == cosmosim::io::IcFieldDisposition::kConverted;
        });
    if (box_field == manifest.fields.end()) {
      throw std::runtime_error(
          "canonical streaming conversion requires Header/BoxSize contract");
    }
    const double box_code = manifest.box_size *
        cosmosim::io::icFieldConversionMultiplier(
            *box_field, manifest, target_units);
    writeAttributeF64(m_header.get(), "Time", manifest.scale_factor);
    writeAttributeF64(m_header.get(), "Redshift", manifest.redshift);
    writeAttributeF64(m_header.get(), "BoxSize", box_code);
    writeAttributeF64(m_header.get(), "Omega0", manifest.omega_matter);
    writeAttributeF64(m_header.get(), "OmegaLambda", manifest.omega_lambda);
    writeAttributeF64(m_header.get(), "HubbleParam", manifest.hubble_param);
    writeAttributeU32(m_header.get(), "NumFilesPerSnapshot", 1U);
    writeAttributeString(m_header.get(), "ChuiIcSchemaName", "chui_canonical_v1");
    writeAttributeU32(m_header.get(), "ChuiIcSchemaVersion", 2U);
    writeAttributeF64(
        m_header.get(), "ChuiLengthUnitToSI", target_units.length_si_per_code);
    writeAttributeF64(
        m_header.get(), "ChuiMassUnitToSI", target_units.mass_si_per_code);
    writeAttributeF64(
        m_header.get(), "ChuiVelocityUnitToSI", target_units.velocity_si_per_code);
    writeAttributeString(m_header.get(), "ChuiCoordinateFrame", "comoving");
    writeAttributeString(
        m_header.get(), "ChuiVelocityConvention", "physical_peculiar");
    writeAttributeString(
        m_header.get(), "ConversionManifestSha256", manifest_digest);
    writeAttributeString(
        m_header.get(), "ConversionManifestSidecar", manifest_filename);
    writeAttributeString(
        m_header.get(), "ConversionBundleMarker", marker_filename);
    Hdf5Handle provenance(H5Gcreate2(
        m_file.get(), "/Provenance", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    if (!provenance.valid()) {
      throw std::runtime_error("failed to create /Provenance");
    }
    writeByteDataset(
        provenance.get(), "ConversionManifestJson", manifest_json);
  }

  void createParticleDatasets() {
    for (std::size_t type : {0U, 1U, 4U, 5U}) {
      if (m_expected_counts[type] == 0U) continue;
      auto datasets = std::make_unique<DatasetSet>();
      const std::string group_name = "/PartType" + std::to_string(type);
      datasets->group = Hdf5Handle(H5Gcreate2(
          m_file.get(), group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT,
          H5P_DEFAULT));
      if (!datasets->group.valid()) {
        throw std::runtime_error("failed to create " + group_name);
      }
      datasets->coordinates = createDatasetVec3(
          datasets->group.get(), "Coordinates", m_expected_counts[type]);
      datasets->velocities = createDatasetVec3(
          datasets->group.get(), "Velocities", m_expected_counts[type]);
      datasets->masses = createDataset1d(
          datasets->group.get(), "Masses", H5T_IEEE_F64LE,
          m_expected_counts[type]);
      datasets->ids = createDataset1d(
          datasets->group.get(), "ParticleIDs", H5T_STD_U64LE,
          m_expected_counts[type]);
      if (type == 0U) {
        datasets->first_sidecar = createDataset1d(
            datasets->group.get(), "InternalEnergy", H5T_IEEE_F64LE,
            m_expected_counts[type]);
        datasets->second_sidecar = createDataset1d(
            datasets->group.get(), "Density", H5T_IEEE_F64LE,
            m_expected_counts[type]);
      } else if (type == 4U) {
        datasets->first_sidecar = createDataset1d(
            datasets->group.get(), "StellarFormationTime", H5T_IEEE_F64LE,
            m_expected_counts[type]);
        datasets->second_sidecar = createDataset1d(
            datasets->group.get(), "InitialMass", H5T_IEEE_F64LE,
            m_expected_counts[type]);
        datasets->third_sidecar = createDataset1d(
            datasets->group.get(), "Metallicity", H5T_IEEE_F64LE,
            m_expected_counts[type]);
      } else if (type == 5U) {
        datasets->first_sidecar = createDataset1d(
            datasets->group.get(), "BH_Mass", H5T_IEEE_F64LE,
            m_expected_counts[type]);
        datasets->second_sidecar = createDataset1d(
            datasets->group.get(), "BH_Mdot", H5T_IEEE_F64LE,
            m_expected_counts[type]);
      }
      m_datasets[type] = std::move(datasets);
    }
  }

  void appendType(
      std::size_t type,
      const std::vector<const cosmosim::io::internal::IcParticleRecord*>&
          records) {
    DatasetSet& datasets = *m_datasets[type];
    const std::uint64_t offset = m_offsets[type];
    if (records.size() > m_expected_counts[type] - offset) {
      throw std::runtime_error("canonical streaming append exceeds declared count");
    }
    std::vector<double> coordinates;
    std::vector<double> velocities;
    std::vector<double> masses;
    std::vector<std::uint64_t> ids;
    std::vector<double> first;
    std::vector<double> second;
    std::vector<double> third;
    coordinates.reserve(records.size() * 3U);
    velocities.reserve(records.size() * 3U);
    masses.reserve(records.size());
    ids.reserve(records.size());
    first.reserve(records.size());
    second.reserve(records.size());
    third.reserve(records.size());
    for (const auto* record : records) {
      coordinates.insert(
          coordinates.end(), {record->x, record->y, record->z});
      velocities.insert(
          velocities.end(), {record->vx, record->vy, record->vz});
      masses.push_back(record->mass);
      ids.push_back(record->id);
      if (type == 0U) {
        first.push_back(record->gas_internal_energy);
        second.push_back(record->gas_density);
      } else if (type == 4U) {
        first.push_back(record->star_formation);
        second.push_back(record->star_birth_mass);
        third.push_back(record->star_metallicity);
      } else if (type == 5U) {
        first.push_back(record->bh_mass);
        second.push_back(record->bh_mdot);
      }
    }
    writeDatasetVec3Slice(datasets.coordinates.get(), offset, coordinates);
    writeDatasetVec3Slice(datasets.velocities.get(), offset, velocities);
    writeDataset1dSlice(
        datasets.masses.get(), H5T_NATIVE_DOUBLE, offset, masses);
    writeDatasetU64Slice(datasets.ids.get(), offset, ids);
    if (!first.empty()) {
      writeDataset1dSlice(
          datasets.first_sidecar.get(), H5T_NATIVE_DOUBLE, offset, first);
    }
    if (!second.empty()) {
      writeDataset1dSlice(
          datasets.second_sidecar.get(), H5T_NATIVE_DOUBLE, offset, second);
    }
    if (!third.empty()) {
      writeDataset1dSlice(
          datasets.third_sidecar.get(), H5T_NATIVE_DOUBLE, offset, third);
    }
    m_offsets[type] += records.size();
  }

  std::array<std::uint64_t, 6> m_expected_counts{};
  std::array<std::uint64_t, 6> m_offsets{};
  std::array<std::unique_ptr<DatasetSet>, 6> m_datasets;
  std::size_t m_chunk_capacity = 0U;
  std::size_t m_peak_batch_records = 0U;
  Hdf5Handle m_file;
  Hdf5Handle m_header;
};

int run(const Arguments& arguments) {
  const bool manifest_mode = !arguments.source_manifest.empty();
  Arguments effective_arguments = arguments;
  cosmosim::io::IcManifest supplied_manifest;
  const cosmosim::io::IcManifest* manifest_pointer = nullptr;
  std::filesystem::path source = arguments.input;
  if (manifest_mode) {
    supplied_manifest =
        cosmosim::io::readIcManifestJson(arguments.source_manifest);
    if (supplied_manifest.source_files.empty()) {
      throw std::runtime_error("source manifest contains no source files");
    }
    if (supplied_manifest.species_policy[2] ==
            cosmosim::io::IcSpeciesPolicy::kTracer ||
        supplied_manifest.species_policy[3] ==
            cosmosim::io::IcSpeciesPolicy::kTracer) {
      throw std::invalid_argument(
          "canonical CHUÍ v1 has no tracer output family; tracer policy in "
          "the supplied manifest is unsupported");
    }
    applyManifestMissingFieldPolicies(effective_arguments, supplied_manifest);
    const std::filesystem::path manifest_directory =
        std::filesystem::absolute(arguments.source_manifest)
            .lexically_normal()
            .parent_path();
    for (auto& source_file : supplied_manifest.source_files) {
      std::filesystem::path source_path(source_file);
      if (source_path.is_relative()) {
        source_path = manifest_directory / source_path;
      }
      source_file = std::filesystem::absolute(source_path).lexically_normal().string();
    }
    source = supplied_manifest.source_files.front();
    manifest_pointer = &supplied_manifest;
  }
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      configText(effective_arguments, manifest_mode), "cosmosim_convert_ic");
  cosmosim::io::IcImportOptions options;
  options.chunk_particle_count = arguments.chunk_particles;
  options.manifest = manifest_pointer;
  options.validate_runtime_cosmology = false;

  const std::filesystem::path output_absolute =
      std::filesystem::absolute(arguments.output).lexically_normal();
  const std::filesystem::path manifest_absolute =
      std::filesystem::absolute(arguments.manifest_output).lexically_normal();
  if (output_absolute.parent_path() != manifest_absolute.parent_path()) {
    throw std::invalid_argument(
        "canonical HDF5 and audit manifest must share one directory");
  }
  const std::filesystem::path marker_absolute =
      output_absolute.string() + ".complete";
  const std::filesystem::path output_partial =
      output_absolute.string() + ".part";
  const std::filesystem::path manifest_partial =
      manifest_absolute.string() + ".part";
  const std::filesystem::path marker_partial =
      marker_absolute.string() + ".part";
  for (const auto& path : {output_partial, manifest_partial, marker_partial}) {
    std::error_code error;
    std::filesystem::remove(path, error);
  }
  for (const auto& path : {output_absolute, manifest_absolute, marker_absolute}) {
    if (std::filesystem::exists(path)) {
      throw std::runtime_error(
          "refusing to overwrite existing canonical bundle member: " +
          path.string());
    }
  }

  bool manifest_finalized = false;
  bool output_finalized = false;
  bool marker_finalized = false;
  const auto cleanup = [&]() noexcept {
    for (const auto& path : {output_partial, manifest_partial, marker_partial}) {
      std::error_code error;
      std::filesystem::remove(path, error);
    }
    if (marker_finalized) {
      std::error_code error;
      std::filesystem::remove(marker_absolute, error);
    }
    if (output_finalized) {
      std::error_code error;
      std::filesystem::remove(output_absolute, error);
    }
    if (manifest_finalized) {
      std::error_code error;
      std::filesystem::remove(manifest_absolute, error);
    }
  };

  std::optional<cosmosim::io::IcManifest> final_manifest;
  std::string manifest_json;
  std::string digest;
  std::string hdf5_bound_digest;
  std::unique_ptr<CanonicalStreamingWriter> writer;
  ExternalIdAudit id_audit(
      output_absolute.parent_path() /
      ("." + output_absolute.filename().string() + ".id_audit"));

  try {
    const cosmosim::io::IcImportReport report =
        cosmosim::io::internal::streamGadgetArepoHdf5Ic(
            source, frozen.config, options,
            [&](const cosmosim::io::IcManifest& inspected_manifest) {
              final_manifest = inspected_manifest;
              final_manifest->converter_version = "chui_ic_converter_v4";
              if (manifest_mode) {
                final_manifest->source_manifest_file =
                    std::filesystem::absolute(arguments.source_manifest)
                        .lexically_normal()
                        .string();
                final_manifest->source_manifest_sha256 =
                    cosmosim::io::icSha256FileHex(arguments.source_manifest);
              }
              manifest_json =
                  cosmosim::io::serializeIcManifestJson(*final_manifest);
              digest = cosmosim::io::icSha256Hex(manifest_json);
              if (injectFault("manifest_write")) {
                throw std::runtime_error(
                    "injected audit manifest write failure");
              }
              writeTextExact(manifest_partial, manifest_json);
              hdf5_bound_digest = digest;
              if (injectFault("hash_mismatch")) {
                hdf5_bound_digest.assign(64U, '0');
              }
              writer = std::make_unique<CanonicalStreamingWriter>(
                  output_partial, *final_manifest, frozen.config,
                  manifest_json, hdf5_bound_digest,
                  manifest_absolute.filename().string(),
                  marker_absolute.filename().string(),
                  arguments.chunk_particles);
            },
            [&](std::span<const cosmosim::io::internal::IcParticleRecord>
                    records) {
              if (!writer) {
                throw std::logic_error(
                    "streaming converter received records before the manifest");
              }
              id_audit.addBatch(records);
              writer->append(records);
            });
    if (!writer || !final_manifest.has_value()) {
      throw std::logic_error("streaming converter did not initialize output");
    }
    id_audit.finish();
    writer->finish(report.counters.peak_staging_bytes);
    writer.reset();
    if (injectFault("hdf5_write")) {
      throw std::runtime_error("injected canonical HDF5 write failure");
    }
    if (cosmosim::io::icSha256FileHex(manifest_partial) != digest) {
      throw std::runtime_error(
          "canonical bundle manifest SHA-256 failed pre-commit verification");
    }
    if (hdf5_bound_digest != digest) {
      throw std::runtime_error(
          "canonical HDF5 manifest binding disagrees with the audit manifest");
    }
    const std::string marker_text =
        "chui_ic_bundle_v1\nsha256=" + digest +
        "\ncanonical=" + output_absolute.filename().string() +
        "\nmanifest=" + manifest_absolute.filename().string() + "\n";
    writeTextExact(marker_partial, marker_text);

    if (injectFault("rename_manifest")) {
      throw std::runtime_error("injected manifest finalization failure");
    }
    std::filesystem::rename(manifest_partial, manifest_absolute);
    manifest_finalized = true;
    if (injectFault("rename_hdf5")) {
      throw std::runtime_error("injected HDF5 finalization failure");
    }
    std::filesystem::rename(output_partial, output_absolute);
    output_finalized = true;
    if (injectFault("rename_marker")) {
      throw std::runtime_error(
          "injected bundle-marker finalization failure");
    }
    std::filesystem::rename(marker_partial, marker_absolute);
    marker_finalized = true;
  } catch (...) {
    writer.reset();
    cleanup();
    throw;
  }
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
