#include "cosmosim/core/config.hpp"

#include <algorithm>
#include <charconv>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_mode.hpp"

namespace cosmosim::core {
namespace {

struct ParsedEntry {
  std::string value;
  int line_number = 0;
};

[[nodiscard]] std::string trim(const std::string& input) {
  const auto begin = std::find_if_not(input.begin(), input.end(), [](unsigned char c) {
    return std::isspace(c) != 0;
  });
  const auto end = std::find_if_not(input.rbegin(), input.rend(), [](unsigned char c) {
    return std::isspace(c) != 0;
  }).base();

  if (begin >= end) {
    return {};
  }

  return std::string(begin, end);
}

[[nodiscard]] std::string toLower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return value;
}

[[nodiscard]] std::string removeInlineComment(const std::string& line) {
  std::size_t comment_pos = std::string::npos;
  for (const std::string marker : {"#", ";", "//"}) {
    const std::size_t pos = line.find(marker);
    if (pos != std::string::npos && (comment_pos == std::string::npos || pos < comment_pos)) {
      comment_pos = pos;
    }
  }

  if (comment_pos == std::string::npos) {
    return line;
  }
  return line.substr(0, comment_pos);
}

[[nodiscard]] std::map<std::string, ParsedEntry> parseEntries(const std::string& text) {
  std::istringstream stream(text);
  std::string raw_line;
  int line_number = 0;
  std::string current_section;
  std::map<std::string, ParsedEntry> entries;

  while (std::getline(stream, raw_line)) {
    ++line_number;
    const std::string line = trim(removeInlineComment(raw_line));
    if (line.empty()) {
      continue;
    }

    if (line.front() == '[' && line.back() == ']') {
      current_section = toLower(trim(line.substr(1, line.size() - 2)));
      if (current_section.empty()) {
        throw ConfigError("line " + std::to_string(line_number) + ": empty section name");
      }
      continue;
    }

    std::string key;
    std::string value;
    const std::size_t equal_pos = line.find('=');
    if (equal_pos != std::string::npos) {
      key = trim(line.substr(0, equal_pos));
      value = trim(line.substr(equal_pos + 1));
    } else {
      std::istringstream token_stream(line);
      token_stream >> key;
      std::getline(token_stream, value);
      value = trim(value);
    }

    key = toLower(trim(key));
    if (key.empty() || value.empty()) {
      throw ConfigError("line " + std::to_string(line_number) + ": expected key and value");
    }

    if (!current_section.empty() && key.find('.') == std::string::npos) {
      key = current_section + "." + key;
    }

    auto [it, inserted] = entries.emplace(key, ParsedEntry{value, line_number});
    if (!inserted) {
      throw ConfigError("line " + std::to_string(line_number) + ": duplicate key '" + key +
                        "' (first at line " + std::to_string(it->second.line_number) + ")");
    }
  }

  return entries;
}

[[nodiscard]] bool parseBool(const std::string& value, const std::string& key) {
  const std::string lower = toLower(trim(value));
  if (lower == "1" || lower == "true" || lower == "on" || lower == "yes") {
    return true;
  }
  if (lower == "0" || lower == "false" || lower == "off" || lower == "no") {
    return false;
  }
  throw ConfigError("key '" + key + "': invalid boolean value '" + value + "'");
}

template <typename T>
[[nodiscard]] T parseNumber(const std::string& value, const std::string& key) {
  const std::string trimmed = trim(value);
  T number{};
  const char* begin = trimmed.data();
  const char* end = trimmed.data() + trimmed.size();
  const auto [ptr, ec] = std::from_chars(begin, end, number);
  if (ec != std::errc{} || ptr != end) {
    throw ConfigError("key '" + key + "': invalid numeric value '" + value + "'");
  }
  return number;
}

[[nodiscard]] double parseFloating(const std::string& value, const std::string& key) {
  const std::string trimmed = trim(value);
  char* consumed = nullptr;
  const double number = std::strtod(trimmed.c_str(), &consumed);
  if (consumed != trimmed.c_str() + static_cast<std::ptrdiff_t>(trimmed.size())) {
    throw ConfigError("key '" + key + "': invalid floating-point value '" + value + "'");
  }
  return number;
}

[[nodiscard]] std::pair<double, std::string> splitMagnitudeUnit(const std::string& value) {
  std::istringstream stream(value);
  double magnitude = 0.0;
  stream >> magnitude;
  if (stream.fail()) {
    throw ConfigError("invalid value '" + value + "': expected numeric magnitude");
  }

  std::string unit;
  stream >> unit;
  unit = toLower(unit);
  return {magnitude, unit};
}

[[nodiscard]] double lengthToMpc(double value, const std::string& unit) {
  const std::string normalized = toLower(unit);
  if (normalized.empty() || normalized == "mpc") {
    return value;
  }
  if (normalized == "kpc") {
    return value * 1.0e-3;
  }
  throw ConfigError("unsupported length unit '" + unit + "' (supported: mpc, kpc)");
}

[[nodiscard]] double parseLengthMpc(
    const std::string& value,
    const std::string& default_unit,
    const std::string& key) {
  const auto [magnitude, provided_unit] = splitMagnitudeUnit(value);
  const std::string unit = provided_unit.empty() ? default_unit : provided_unit;
  try {
    return lengthToMpc(magnitude, unit);
  } catch (const ConfigError&) {
    throw ConfigError("key '" + key + "': unsupported unit in value '" + value + "'");
  }
}

[[nodiscard]] double parseLengthKpc(
    const std::string& value,
    const std::string& default_unit,
    const std::string& key) {
  const double in_mpc = parseLengthMpc(value, default_unit, key);
  return in_mpc * 1.0e3;
}

[[nodiscard]] std::string requireString(
    std::map<std::string, ParsedEntry>& entries,
    std::set<std::string>& consumed,
    const std::string& key,
    const std::string& default_value) {
  const auto it = entries.find(key);
  if (it == entries.end()) {
    return default_value;
  }
  consumed.insert(key);
  return trim(it->second.value);
}

[[nodiscard]] std::string sanitizeStem(const std::string& value, const std::string& key) {
  if (value.empty()) {
    throw ConfigError("key '" + key + "': value cannot be empty");
  }
  for (const char c : value) {
    const bool valid = std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-';
    if (!valid) {
      throw ConfigError(
          "key '" + key + "': only [a-zA-Z0-9_-] are allowed for stable output naming");
    }
  }
  return value;
}

[[nodiscard]] std::string modeToLowerString(SimulationMode mode) {
  switch (mode) {
    case SimulationMode::kCosmoCube:
      return "cosmo_cube";
    case SimulationMode::kZoomIn:
      return "zoom_in";
    case SimulationMode::kIsolatedGalaxy:
      return "isolated_galaxy";
    case SimulationMode::kIsolatedCluster:
      return "isolated_cluster";
  }
  return "unknown";
}

[[nodiscard]] SimulationMode parseMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "cosmo_cube") {
    return SimulationMode::kCosmoCube;
  }
  if (lower == "zoom_in") {
    return SimulationMode::kZoomIn;
  }
  if (lower == "isolated_galaxy") {
    return SimulationMode::kIsolatedGalaxy;
  }
  if (lower == "isolated_cluster") {
    return SimulationMode::kIsolatedCluster;
  }
  throw ConfigError("key 'mode.mode': invalid mode '" + value + "'");
}

[[nodiscard]] GravitySolver parseGravitySolver(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "treepm") {
    return GravitySolver::kTreePm;
  }
  throw ConfigError("key 'numerics.gravity_solver': invalid value '" + value + "'");
}

[[nodiscard]] HydroSolver parseHydroSolver(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "godunov_fv") {
    return HydroSolver::kGodunovFv;
  }
  throw ConfigError("key 'numerics.hydro_solver': invalid value '" + value + "'");
}

[[nodiscard]] CoordinateFrame parseCoordinateFrame(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "comoving") {
    return CoordinateFrame::kComoving;
  }
  if (lower == "physical") {
    return CoordinateFrame::kPhysical;
  }
  throw ConfigError("key 'units.coordinate_frame': invalid value '" + value + "'");
}

[[nodiscard]] ModeHydroBoundary parseModeHydroBoundary(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "auto") {
    return ModeHydroBoundary::kAuto;
  }
  if (lower == "periodic") {
    return ModeHydroBoundary::kPeriodic;
  }
  if (lower == "open") {
    return ModeHydroBoundary::kOpen;
  }
  if (lower == "reflective") {
    return ModeHydroBoundary::kReflective;
  }
  throw ConfigError("key 'mode.hydro_boundary': invalid value '" + value + "'");
}

[[nodiscard]] ModeGravityBoundary parseModeGravityBoundary(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "auto") {
    return ModeGravityBoundary::kAuto;
  }
  if (lower == "periodic") {
    return ModeGravityBoundary::kPeriodic;
  }
  if (lower == "isolated_monopole") {
    return ModeGravityBoundary::kIsolatedMonopole;
  }
  throw ConfigError("key 'mode.gravity_boundary': invalid value '" + value + "'");
}

[[nodiscard]] FeedbackMode parseFeedbackMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "thermal") {
    return FeedbackMode::kThermal;
  }
  if (lower == "kinetic") {
    return FeedbackMode::kKinetic;
  }
  if (lower == "momentum") {
    return FeedbackMode::kMomentum;
  }
  if (lower == "thermal_kinetic_momentum") {
    return FeedbackMode::kThermalKineticMomentum;
  }
  throw ConfigError(
      "physics.fb_mode must be one of: thermal, kinetic, momentum, thermal_kinetic_momentum");
}

[[nodiscard]] FeedbackVariant parseFeedbackVariant(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "none") {
    return FeedbackVariant::kNone;
  }
  if (lower == "delayed_cooling") {
    return FeedbackVariant::kDelayedCooling;
  }
  if (lower == "stochastic") {
    return FeedbackVariant::kStochastic;
  }
  throw ConfigError("physics.fb_variant must be one of: none, delayed_cooling, stochastic");
}

[[nodiscard]] UvBackgroundModel parseUvBackgroundModel(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "none") {
    return UvBackgroundModel::kNone;
  }
  if (lower == "hm12") {
    return UvBackgroundModel::kHm12;
  }
  if (lower == "fg20") {
    return UvBackgroundModel::kFg20;
  }
  throw ConfigError("key 'physics.uv_background_model': invalid value '" + value + "'");
}

[[nodiscard]] SelfShieldingModel parseSelfShieldingModel(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "none") {
    return SelfShieldingModel::kNone;
  }
  if (lower == "rahmati13_like") {
    return SelfShieldingModel::kRahmati13Like;
  }
  throw ConfigError("key 'physics.self_shielding_model': invalid value '" + value + "'");
}

[[nodiscard]] CoolingModel parseCoolingModel(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "primordial") {
    return CoolingModel::kPrimordial;
  }
  if (lower == "primordial_metal_line") {
    return CoolingModel::kPrimordialMetalLine;
  }
  throw ConfigError("key 'physics.cooling_model': invalid value '" + value + "'");
}

struct ConfigKeySpec {
  const char* key;
  const char* default_value;
};

[[nodiscard]] const std::vector<ConfigKeySpec>& configKeyRegistry() {
  static const std::vector<ConfigKeySpec> keys = {
      {"schema_version", "1"},
      {"units.length_unit", "mpc"},
      {"units.mass_unit", "msun"},
      {"units.velocity_unit", "km_s"},
      {"units.coordinate_frame", "comoving"},
      {"mode.mode", "zoom_in"},
      {"mode.ic_file", "ics.hdf5"},
      {"mode.zoom_high_res_region", "false"},
      {"mode.zoom_region_file", ""},
      {"mode.hydro_boundary", "auto"},
      {"mode.gravity_boundary", "auto"},
      {"cosmology.omega_matter", "0.315"},
      {"cosmology.omega_lambda", "0.685"},
      {"cosmology.omega_baryon", "0.049"},
      {"cosmology.hubble_param", "0.674"},
      {"cosmology.sigma8", "0.811"},
      {"cosmology.scalar_index_ns", "0.965"},
      {"cosmology.box_size", "50.0"},
      {"numerics.time_begin_code", "0.0"},
      {"numerics.time_end_code", "1.0"},
      {"numerics.max_global_steps", "1024"},
      {"numerics.hierarchical_max_rung", "12"},
      {"numerics.amr_max_level", "10"},
      {"numerics.gravity_softening", "1.0 kpc"},
      {"numerics.gravity_solver", "treepm"},
      {"numerics.hydro_solver", "godunov_fv"},
      {"physics.enable_cooling", "true"},
      {"physics.enable_star_formation", "true"},
      {"physics.enable_feedback", "true"},
      {"physics.enable_stellar_evolution", "true"},
      {"physics.reionization_model", "hm12"},
      {"physics.uv_background_model", "hm12"},
      {"physics.self_shielding_model", "none"},
      {"physics.cooling_model", "primordial"},
      {"physics.metal_line_table_path", ""},
      {"physics.temperature_floor_k", "100.0"},
      {"physics.sf_density_threshold_code", "10.0"},
      {"physics.sf_temperature_threshold_k", "1.0e4"},
      {"physics.sf_min_converging_flow_rate_code", "0.0"},
      {"physics.sf_epsilon_ff", "0.01"},
      {"physics.sf_min_star_particle_mass_code", "0.1"},
      {"physics.sf_stochastic_spawning", "true"},
      {"physics.sf_random_seed", "123456789"},
      {"physics.fb_mode", "thermal_kinetic_momentum"},
      {"physics.fb_variant", "none"},
      {"physics.fb_use_returned_mass_budget", "true"},
      {"physics.fb_epsilon_thermal", "0.6"},
      {"physics.fb_epsilon_kinetic", "0.3"},
      {"physics.fb_epsilon_momentum", "0.1"},
      {"physics.fb_sn_energy_erg_per_mass_code", "1.0e49"},
      {"physics.fb_momentum_code_per_mass_code", "3.0e3"},
      {"physics.fb_neighbor_count", "8"},
      {"physics.fb_delayed_cooling_time_code", "0.0"},
      {"physics.fb_stochastic_event_probability", "0.25"},
      {"physics.fb_random_seed", "42424242"},
      {"physics.stellar_evolution_table_path", ""},
      {"physics.stellar_evolution_hubble_time_years", "1.44e10"},
      {"physics.enable_black_hole_agn", "false"},
      {"physics.bh_seed_halo_mass_threshold_code", "1.0e3"},
      {"physics.bh_seed_mass_code", "1.0"},
      {"physics.bh_seed_max_per_cell", "1"},
      {"physics.bh_alpha_bondi", "1.0"},
      {"physics.bh_use_eddington_cap", "true"},
      {"physics.bh_epsilon_r", "0.1"},
      {"physics.bh_epsilon_f", "0.05"},
      {"physics.bh_feedback_coupling_efficiency", "1.0"},
      {"physics.bh_duty_cycle_active_edd_ratio_threshold", "0.01"},
      {"physics.bh_proton_mass_si", "1.67262192369e-27"},
      {"physics.bh_thomson_cross_section_si", "6.6524587321e-29"},
      {"physics.bh_newton_g_si", "6.67430e-11"},
      {"physics.bh_speed_of_light_si", "2.99792458e8"},
      {"physics.enable_tracers", "false"},
      {"physics.tracer_track_mass", "true"},
      {"physics.tracer_min_host_mass_code", "0.0"},
      {"output.run_name", "cosmosim_run"},
      {"output.output_directory", "outputs"},
      {"output.output_stem", "snapshot"},
      {"output.restart_stem", "restart"},
      {"output.snapshot_interval_steps", "64"},
      {"output.write_restarts", "true"},
      {"parallel.mpi_ranks_expected", "1"},
      {"parallel.omp_threads", "1"},
      {"parallel.gpu_devices", "0"},
      {"parallel.deterministic_reduction", "true"},
      {"analysis.enable_diagnostics", "true"},
      {"analysis.enable_halo_workflow", "false"},
      {"analysis.halo_on_the_fly", "false"},
      {"analysis.run_health_interval_steps", "1"},
      {"analysis.science_light_interval_steps", "8"},
      {"analysis.science_heavy_interval_steps", "64"},
      {"analysis.retention_bundle_count", "8"},
      {"analysis.power_spectrum_mesh_n", "16"},
      {"analysis.power_spectrum_bin_count", "12"},
      {"analysis.sf_history_bin_count", "16"},
      {"analysis.quicklook_grid_n", "32"},
      {"analysis.diagnostics_stem", "diagnostics"},
      {"analysis.halo_catalog_stem", "halo_catalog"},
      {"analysis.merger_tree_stem", "merger_tree_plan"},
      {"analysis.halo_fof_linking_length_factor", "0.2"},
      {"analysis.halo_fof_min_group_size", "16"},
      {"analysis.halo_include_gas", "true"},
      {"analysis.halo_include_stars", "true"},
      {"analysis.halo_include_black_holes", "true"},
      {"compatibility.allow_unknown_keys", "false"},
  };
  return keys;
}

[[nodiscard]] const std::map<std::string, std::string>& deprecatedAliasRegistry() {
  static const std::map<std::string, std::string> aliases = {
      {"omega0", "cosmology.omega_matter"},
      {"omegalambda", "cosmology.omega_lambda"},
      {"hubbleparam", "cosmology.hubble_param"},
      {"timemax", "numerics.time_end_code"},
      {"mode", "mode.mode"},
      {"run_name", "output.run_name"},
  };
  return aliases;
}

[[nodiscard]] std::string defaultFor(const std::string& key) {
  for (const auto& entry : configKeyRegistry()) {
    if (entry.key == key) {
      return entry.default_value;
    }
  }
  throw ConfigError("internal config key registry missing default for key '" + key + "'");
}

void validateConfig(const SimulationConfig& config) {
  if (config.schema_version != 1) {
    throw ConfigError("schema_version must be 1 for this build");
  }
  if (config.numerics.time_end_code <= config.numerics.time_begin_code) {
    throw ConfigError("numerics.time_end_code must be greater than numerics.time_begin_code");
  }
  if (config.numerics.max_global_steps <= 0) {
    throw ConfigError("numerics.max_global_steps must be > 0");
  }
  if (config.output.snapshot_interval_steps <= 0) {
    throw ConfigError("output.snapshot_interval_steps must be > 0");
  }
  if (config.parallel.mpi_ranks_expected <= 0 || config.parallel.omp_threads <= 0) {
    throw ConfigError("parallel settings require positive mpi_ranks_expected and omp_threads");
  }
  if (config.analysis.run_health_interval_steps <= 0 ||
      config.analysis.science_light_interval_steps <= 0 ||
      config.analysis.science_heavy_interval_steps <= 0) {
    throw ConfigError("analysis cadence intervals must be > 0");
  }
  if (config.analysis.retention_bundle_count < 1) {
    throw ConfigError("analysis.retention_bundle_count must be >= 1");
  }
  if (config.analysis.power_spectrum_mesh_n < 4 ||
      config.analysis.power_spectrum_bin_count < 1 ||
      config.analysis.sf_history_bin_count < 1 ||
      config.analysis.quicklook_grid_n < 4) {
    throw ConfigError("analysis mesh/bin settings must be within conservative minimum bounds");
  }
  if (config.analysis.halo_fof_linking_length_factor <= 0.0 ||
      config.analysis.halo_fof_linking_length_factor > 1.0) {
    throw ConfigError("analysis.halo_fof_linking_length_factor must be in (0, 1]");
  }
  if (config.analysis.halo_fof_min_group_size < 2) {
    throw ConfigError("analysis.halo_fof_min_group_size must be >= 2");
  }
  if (config.cosmology.omega_matter <= 0.0 || config.cosmology.omega_lambda < 0.0) {
    throw ConfigError("cosmology requires omega_matter > 0 and omega_lambda >= 0");
  }
  if (config.physics.temperature_floor_k <= 0.0) {
    throw ConfigError("physics.temperature_floor_k must be > 0");
  }
  if (config.physics.sf_density_threshold_code <= 0.0 ||
      config.physics.sf_temperature_threshold_k <= 0.0 ||
      config.physics.sf_min_star_particle_mass_code <= 0.0) {
    throw ConfigError("star formation thresholds and min star mass must be > 0");
  }
  if (config.physics.sf_epsilon_ff < 0.0 || config.physics.sf_epsilon_ff > 1.0) {
    throw ConfigError("physics.sf_epsilon_ff must be in [0, 1]");
  }
  if (config.physics.fb_epsilon_thermal < 0.0 || config.physics.fb_epsilon_kinetic < 0.0 ||
      config.physics.fb_epsilon_momentum < 0.0) {
    throw ConfigError("physics feedback efficiencies must be >= 0");
  }
  if (config.physics.fb_sn_energy_erg_per_mass_code <= 0.0 ||
      config.physics.fb_momentum_code_per_mass_code < 0.0) {
    throw ConfigError("physics feedback budget scales must be physically non-negative");
  }
  if (config.physics.fb_neighbor_count == 0) {
    throw ConfigError("physics.fb_neighbor_count must be > 0");
  }
  if (config.physics.fb_stochastic_event_probability <= 0.0 ||
      config.physics.fb_stochastic_event_probability > 1.0) {
    throw ConfigError("physics.fb_stochastic_event_probability must be in (0, 1]");
  }
  if (config.physics.stellar_evolution_hubble_time_years <= 0.0) {
    throw ConfigError("physics.stellar_evolution_hubble_time_years must be > 0");
  }
  if (config.physics.bh_seed_halo_mass_threshold_code <= 0.0 || config.physics.bh_seed_mass_code <= 0.0 ||
      config.physics.bh_seed_max_per_cell == 0) {
    throw ConfigError("physics BH seeding parameters must be > 0");
  }
  if (config.physics.bh_alpha_bondi <= 0.0) {
    throw ConfigError("physics.bh_alpha_bondi must be > 0");
  }
  if (config.physics.bh_epsilon_r <= 0.0 || config.physics.bh_epsilon_r > 1.0 ||
      config.physics.bh_epsilon_f < 0.0 || config.physics.bh_epsilon_f > 1.0 ||
      config.physics.bh_feedback_coupling_efficiency < 0.0 ||
      config.physics.bh_feedback_coupling_efficiency > 1.0) {
    throw ConfigError("physics BH efficiencies must be within conservative physical bounds");
  }
  if (config.physics.bh_duty_cycle_active_edd_ratio_threshold < 0.0 ||
      config.physics.bh_duty_cycle_active_edd_ratio_threshold > 1.0) {
    throw ConfigError("physics.bh_duty_cycle_active_edd_ratio_threshold must be in [0, 1]");
  }
  if (config.physics.bh_proton_mass_si <= 0.0 || config.physics.bh_thomson_cross_section_si <= 0.0 ||
      config.physics.bh_newton_g_si <= 0.0 || config.physics.bh_speed_of_light_si <= 0.0) {
    throw ConfigError("physics BH constants must be > 0");
  }
  if (config.physics.tracer_min_host_mass_code < 0.0) {
    throw ConfigError("physics.tracer_min_host_mass_code must be >= 0");
  }
  const ModePolicy policy = buildModePolicy(config.mode);
  validateModePolicy(config, policy);
}


[[nodiscard]] std::string buildNormalizedText(const FrozenConfig& frozen) {
  std::ostringstream stream;
  stream << "schema_version = " << frozen.config.schema_version << '\n';
  stream << "provenance.config_hash = " << frozen.provenance.config_hash_hex << '\n';
  stream << "\n[units]\n";
  stream << "length_unit = " << frozen.config.units.length_unit << '\n';
  stream << "mass_unit = " << frozen.config.units.mass_unit << '\n';
  stream << "velocity_unit = " << frozen.config.units.velocity_unit << '\n';
  stream << "coordinate_frame = " << coordinateFrameToString(frozen.config.units.coordinate_frame) << '\n';
  stream << "\n[mode]\n";
  stream << "mode = " << modeToLowerString(frozen.config.mode.mode) << '\n';
  stream << "ic_file = " << frozen.config.mode.ic_file << '\n';
  stream << "zoom_high_res_region = " << (frozen.config.mode.zoom_high_res_region ? "true" : "false")
         << '\n';
  stream << "zoom_region_file = " << frozen.config.mode.zoom_region_file << '\n';
  stream << "hydro_boundary = " << modeHydroBoundaryToString(frozen.config.mode.hydro_boundary) << '\n';
  stream << "gravity_boundary = " << modeGravityBoundaryToString(frozen.config.mode.gravity_boundary) << '\n';
  stream << "\n[cosmology]\n";
  stream << "omega_matter = " << frozen.config.cosmology.omega_matter << '\n';
  stream << "omega_lambda = " << frozen.config.cosmology.omega_lambda << '\n';
  stream << "omega_baryon = " << frozen.config.cosmology.omega_baryon << '\n';
  stream << "hubble_param = " << frozen.config.cosmology.hubble_param << '\n';
  stream << "sigma8 = " << frozen.config.cosmology.sigma8 << '\n';
  stream << "scalar_index_ns = " << frozen.config.cosmology.scalar_index_ns << '\n';
  stream << "box_size_mpc_comoving = " << frozen.config.cosmology.box_size_mpc_comoving << '\n';
  stream << "\n[numerics]\n";
  stream << "time_begin_code = " << frozen.config.numerics.time_begin_code << '\n';
  stream << "time_end_code = " << frozen.config.numerics.time_end_code << '\n';
  stream << "max_global_steps = " << frozen.config.numerics.max_global_steps << '\n';
  stream << "hierarchical_max_rung = " << frozen.config.numerics.hierarchical_max_rung << '\n';
  stream << "amr_max_level = " << frozen.config.numerics.amr_max_level << '\n';
  stream << "gravity_softening_kpc_comoving = "
         << frozen.config.numerics.gravity_softening_kpc_comoving << '\n';
  stream << "gravity_solver = " << gravitySolverToString(frozen.config.numerics.gravity_solver) << '\n';
  stream << "hydro_solver = " << hydroSolverToString(frozen.config.numerics.hydro_solver) << '\n';
  stream << "\n[physics]\n";
  stream << "enable_cooling = " << (frozen.config.physics.enable_cooling ? "true" : "false") << '\n';
  stream << "enable_star_formation = "
         << (frozen.config.physics.enable_star_formation ? "true" : "false") << '\n';
  stream << "enable_feedback = " << (frozen.config.physics.enable_feedback ? "true" : "false") << '\n';
  stream << "enable_stellar_evolution = "
         << (frozen.config.physics.enable_stellar_evolution ? "true" : "false") << '\n';
  stream << "reionization_model = " << frozen.config.physics.reionization_model << '\n';
  stream << "uv_background_model = " << uvBackgroundModelToString(frozen.config.physics.uv_background_model) << '\n';
  stream << "self_shielding_model = " << selfShieldingModelToString(frozen.config.physics.self_shielding_model) << '\n';
  stream << "cooling_model = " << coolingModelToString(frozen.config.physics.cooling_model) << '\n';
  stream << "metal_line_table_path = " << frozen.config.physics.metal_line_table_path << '\n';
  stream << "temperature_floor_k = " << frozen.config.physics.temperature_floor_k << '\n';
  stream << "sf_density_threshold_code = " << frozen.config.physics.sf_density_threshold_code << '\n';
  stream << "sf_temperature_threshold_k = " << frozen.config.physics.sf_temperature_threshold_k << '\n';
  stream << "sf_min_converging_flow_rate_code = " << frozen.config.physics.sf_min_converging_flow_rate_code
         << '\n';
  stream << "sf_epsilon_ff = " << frozen.config.physics.sf_epsilon_ff << '\n';
  stream << "sf_min_star_particle_mass_code = "
         << frozen.config.physics.sf_min_star_particle_mass_code << '\n';
  stream << "sf_stochastic_spawning = "
         << (frozen.config.physics.sf_stochastic_spawning ? "true" : "false") << '\n';
  stream << "sf_random_seed = " << frozen.config.physics.sf_random_seed << '\n';
  stream << "fb_mode = " << feedbackModeToString(frozen.config.physics.fb_mode) << '\n';
  stream << "fb_variant = " << feedbackVariantToString(frozen.config.physics.fb_variant) << '\n';
  stream << "fb_use_returned_mass_budget = "
         << (frozen.config.physics.fb_use_returned_mass_budget ? "true" : "false") << '\n';
  stream << "fb_epsilon_thermal = " << frozen.config.physics.fb_epsilon_thermal << '\n';
  stream << "fb_epsilon_kinetic = " << frozen.config.physics.fb_epsilon_kinetic << '\n';
  stream << "fb_epsilon_momentum = " << frozen.config.physics.fb_epsilon_momentum << '\n';
  stream << "fb_sn_energy_erg_per_mass_code = " << frozen.config.physics.fb_sn_energy_erg_per_mass_code << '\n';
  stream << "fb_momentum_code_per_mass_code = " << frozen.config.physics.fb_momentum_code_per_mass_code << '\n';
  stream << "fb_neighbor_count = " << frozen.config.physics.fb_neighbor_count << '\n';
  stream << "fb_delayed_cooling_time_code = " << frozen.config.physics.fb_delayed_cooling_time_code << '\n';
  stream << "fb_stochastic_event_probability = " << frozen.config.physics.fb_stochastic_event_probability << '\n';
  stream << "fb_random_seed = " << frozen.config.physics.fb_random_seed << '\n';
  stream << "stellar_evolution_table_path = " << frozen.config.physics.stellar_evolution_table_path << '\n';
  stream << "stellar_evolution_hubble_time_years = "
         << frozen.config.physics.stellar_evolution_hubble_time_years << '\n';
  stream << "enable_black_hole_agn = " << (frozen.config.physics.enable_black_hole_agn ? "true" : "false")
         << '\n';
  stream << "bh_seed_halo_mass_threshold_code = " << frozen.config.physics.bh_seed_halo_mass_threshold_code
         << '\n';
  stream << "bh_seed_mass_code = " << frozen.config.physics.bh_seed_mass_code << '\n';
  stream << "bh_seed_max_per_cell = " << frozen.config.physics.bh_seed_max_per_cell << '\n';
  stream << "bh_alpha_bondi = " << frozen.config.physics.bh_alpha_bondi << '\n';
  stream << "bh_use_eddington_cap = " << (frozen.config.physics.bh_use_eddington_cap ? "true" : "false")
         << '\n';
  stream << "bh_epsilon_r = " << frozen.config.physics.bh_epsilon_r << '\n';
  stream << "bh_epsilon_f = " << frozen.config.physics.bh_epsilon_f << '\n';
  stream << "bh_feedback_coupling_efficiency = "
         << frozen.config.physics.bh_feedback_coupling_efficiency << '\n';
  stream << "bh_duty_cycle_active_edd_ratio_threshold = "
         << frozen.config.physics.bh_duty_cycle_active_edd_ratio_threshold << '\n';
  stream << "bh_proton_mass_si = " << frozen.config.physics.bh_proton_mass_si << '\n';
  stream << "bh_thomson_cross_section_si = " << frozen.config.physics.bh_thomson_cross_section_si << '\n';
  stream << "bh_newton_g_si = " << frozen.config.physics.bh_newton_g_si << '\n';
  stream << "bh_speed_of_light_si = " << frozen.config.physics.bh_speed_of_light_si << '\n';
  stream << "enable_tracers = " << (frozen.config.physics.enable_tracers ? "true" : "false") << '\n';
  stream << "tracer_track_mass = " << (frozen.config.physics.tracer_track_mass ? "true" : "false") << '\n';
  stream << "tracer_min_host_mass_code = " << frozen.config.physics.tracer_min_host_mass_code << '\n';
  stream << "\n[output]\n";
  stream << "run_name = " << frozen.config.output.run_name << '\n';
  stream << "output_directory = " << frozen.config.output.output_directory << '\n';
  stream << "output_stem = " << frozen.config.output.output_stem << '\n';
  stream << "restart_stem = " << frozen.config.output.restart_stem << '\n';
  stream << "snapshot_interval_steps = " << frozen.config.output.snapshot_interval_steps << '\n';
  stream << "write_restarts = " << (frozen.config.output.write_restarts ? "true" : "false") << '\n';
  stream << "\n[parallel]\n";
  stream << "mpi_ranks_expected = " << frozen.config.parallel.mpi_ranks_expected << '\n';
  stream << "omp_threads = " << frozen.config.parallel.omp_threads << '\n';
  stream << "gpu_devices = " << frozen.config.parallel.gpu_devices << '\n';
  stream << "deterministic_reduction = "
         << (frozen.config.parallel.deterministic_reduction ? "true" : "false") << '\n';
  stream << "\n[analysis]\n";
  stream << "enable_diagnostics = " << (frozen.config.analysis.enable_diagnostics ? "true" : "false")
         << '\n';
  stream << "enable_halo_workflow = " << (frozen.config.analysis.enable_halo_workflow ? "true" : "false")
         << '\n';
  stream << "halo_on_the_fly = " << (frozen.config.analysis.halo_on_the_fly ? "true" : "false") << '\n';
  stream << "run_health_interval_steps = " << frozen.config.analysis.run_health_interval_steps << '\n';
  stream << "science_light_interval_steps = " << frozen.config.analysis.science_light_interval_steps << '\n';
  stream << "science_heavy_interval_steps = " << frozen.config.analysis.science_heavy_interval_steps << '\n';
  stream << "retention_bundle_count = " << frozen.config.analysis.retention_bundle_count << '\n';
  stream << "power_spectrum_mesh_n = " << frozen.config.analysis.power_spectrum_mesh_n << '\n';
  stream << "power_spectrum_bin_count = " << frozen.config.analysis.power_spectrum_bin_count << '\n';
  stream << "sf_history_bin_count = " << frozen.config.analysis.sf_history_bin_count << '\n';
  stream << "quicklook_grid_n = " << frozen.config.analysis.quicklook_grid_n << '\n';
  stream << "diagnostics_stem = " << frozen.config.analysis.diagnostics_stem << '\n';
  stream << "halo_catalog_stem = " << frozen.config.analysis.halo_catalog_stem << '\n';
  stream << "merger_tree_stem = " << frozen.config.analysis.merger_tree_stem << '\n';
  stream << "halo_fof_linking_length_factor = " << frozen.config.analysis.halo_fof_linking_length_factor << '\n';
  stream << "halo_fof_min_group_size = " << frozen.config.analysis.halo_fof_min_group_size << '\n';
  stream << "halo_include_gas = " << (frozen.config.analysis.halo_include_gas ? "true" : "false") << '\n';
  stream << "halo_include_stars = " << (frozen.config.analysis.halo_include_stars ? "true" : "false") << '\n';
  stream << "halo_include_black_holes = " << (frozen.config.analysis.halo_include_black_holes ? "true" : "false")
         << '\n';
  stream << "\n[compatibility]\n";
  stream << "allow_unknown_keys = "
         << (frozen.config.compatibility.allow_unknown_keys ? "true" : "false") << '\n';
  return stream.str();
}

[[nodiscard]] FrozenConfig normalizeValidateFreeze(
    const std::map<std::string, ParsedEntry>& parsed_entries,
    const std::string& source_name,
    const ParseOptions& options) {
  std::map<std::string, ParsedEntry> entries = parsed_entries;
  std::set<std::string> consumed;
  FrozenConfig frozen;
  frozen.provenance.source_name = source_name;

  for (const auto& [legacy_key, canonical_key] : deprecatedAliasRegistry()) {
    const auto it = entries.find(legacy_key);
    if (it == entries.end()) {
      continue;
    }
    if (entries.contains(canonical_key)) {
      throw ConfigError("deprecated key '" + legacy_key + "' cannot be combined with '" +
                        canonical_key + "'");
    }
    entries.emplace(canonical_key, it->second);
    consumed.insert(legacy_key);
    frozen.provenance.deprecation_warnings.push_back(
        "deprecated key '" + legacy_key + "' mapped to '" + canonical_key + "'");
  }

  frozen.config.schema_version = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "schema_version", "1"), "schema_version"));

  frozen.config.units.length_unit =
      toLower(requireString(entries, consumed, "units.length_unit", frozen.config.units.length_unit));
  frozen.config.units.mass_unit =
      toLower(requireString(entries, consumed, "units.mass_unit", frozen.config.units.mass_unit));
  frozen.config.units.velocity_unit = toLower(
      requireString(entries, consumed, "units.velocity_unit", frozen.config.units.velocity_unit));
  frozen.config.units.coordinate_frame = parseCoordinateFrame(
      requireString(entries, consumed, "units.coordinate_frame", defaultFor("units.coordinate_frame")));

  frozen.config.mode.mode = parseMode(requireString(entries, consumed, "mode.mode", defaultFor("mode.mode")));
  frozen.config.mode.ic_file = requireString(entries, consumed, "mode.ic_file", defaultFor("mode.ic_file"));
  frozen.config.mode.zoom_high_res_region = parseBool(
      requireString(entries, consumed, "mode.zoom_high_res_region", defaultFor("mode.zoom_high_res_region")),
      "mode.zoom_high_res_region");
  frozen.config.mode.zoom_region_file =
      requireString(entries, consumed, "mode.zoom_region_file", defaultFor("mode.zoom_region_file"));
  frozen.config.mode.hydro_boundary = parseModeHydroBoundary(
      requireString(entries, consumed, "mode.hydro_boundary", defaultFor("mode.hydro_boundary")));
  frozen.config.mode.gravity_boundary = parseModeGravityBoundary(
      requireString(entries, consumed, "mode.gravity_boundary", defaultFor("mode.gravity_boundary")));

  frozen.config.cosmology.omega_matter = parseFloating(
      requireString(entries, consumed, "cosmology.omega_matter", "0.315"),
      "cosmology.omega_matter");
  frozen.config.cosmology.omega_lambda = parseFloating(
      requireString(entries, consumed, "cosmology.omega_lambda", "0.685"),
      "cosmology.omega_lambda");
  frozen.config.cosmology.omega_baryon = parseFloating(
      requireString(entries, consumed, "cosmology.omega_baryon", "0.049"),
      "cosmology.omega_baryon");
  frozen.config.cosmology.hubble_param = parseFloating(
      requireString(entries, consumed, "cosmology.hubble_param", "0.674"),
      "cosmology.hubble_param");
  frozen.config.cosmology.sigma8 =
      parseFloating(requireString(entries, consumed, "cosmology.sigma8", "0.811"), "cosmology.sigma8");
  frozen.config.cosmology.scalar_index_ns = parseFloating(
      requireString(entries, consumed, "cosmology.scalar_index_ns", "0.965"),
      "cosmology.scalar_index_ns");
  frozen.config.cosmology.box_size_mpc_comoving = parseLengthMpc(
      requireString(entries, consumed, "cosmology.box_size", "50.0"),
      frozen.config.units.length_unit,
      "cosmology.box_size");

  frozen.config.numerics.time_begin_code = parseFloating(
      requireString(entries, consumed, "numerics.time_begin_code", "0.0"),
      "numerics.time_begin_code");
  frozen.config.numerics.time_end_code = parseFloating(
      requireString(entries, consumed, "numerics.time_end_code", "1.0"),
      "numerics.time_end_code");
  frozen.config.numerics.max_global_steps = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "numerics.max_global_steps", "1024"),
      "numerics.max_global_steps"));
  frozen.config.numerics.hierarchical_max_rung = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "numerics.hierarchical_max_rung", "12"),
      "numerics.hierarchical_max_rung"));
  frozen.config.numerics.amr_max_level = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "numerics.amr_max_level", "10"), "numerics.amr_max_level"));
  frozen.config.numerics.gravity_softening_kpc_comoving = parseLengthKpc(
      requireString(entries, consumed, "numerics.gravity_softening", "1.0 kpc"),
      frozen.config.units.length_unit,
      "numerics.gravity_softening");
  frozen.config.numerics.gravity_solver = parseGravitySolver(
      requireString(entries, consumed, "numerics.gravity_solver", defaultFor("numerics.gravity_solver")));
  frozen.config.numerics.hydro_solver = parseHydroSolver(
      requireString(entries, consumed, "numerics.hydro_solver", defaultFor("numerics.hydro_solver")));

  frozen.config.physics.enable_cooling = parseBool(
      requireString(entries, consumed, "physics.enable_cooling", "true"), "physics.enable_cooling");
  frozen.config.physics.enable_star_formation =
      parseBool(requireString(entries, consumed, "physics.enable_star_formation", "true"),
                "physics.enable_star_formation");
  frozen.config.physics.enable_feedback = parseBool(
      requireString(entries, consumed, "physics.enable_feedback", "true"), "physics.enable_feedback");
  frozen.config.physics.enable_stellar_evolution = parseBool(
      requireString(entries, consumed, "physics.enable_stellar_evolution", "true"),
      "physics.enable_stellar_evolution");
  frozen.config.physics.reionization_model =
      requireString(entries, consumed, "physics.reionization_model", "hm12");
  frozen.config.physics.uv_background_model = parseUvBackgroundModel(
      requireString(entries, consumed, "physics.uv_background_model", defaultFor("physics.uv_background_model")));
  frozen.config.physics.self_shielding_model = parseSelfShieldingModel(
      requireString(entries, consumed, "physics.self_shielding_model", defaultFor("physics.self_shielding_model")));
  frozen.config.physics.cooling_model = parseCoolingModel(
      requireString(entries, consumed, "physics.cooling_model", defaultFor("physics.cooling_model")));
  frozen.config.physics.metal_line_table_path =
      requireString(entries, consumed, "physics.metal_line_table_path", "");
  frozen.config.physics.temperature_floor_k = parseFloating(
      requireString(entries, consumed, "physics.temperature_floor_k", "100.0"),
      "physics.temperature_floor_k");
  frozen.config.physics.sf_density_threshold_code = parseFloating(
      requireString(entries, consumed, "physics.sf_density_threshold_code", "10.0"),
      "physics.sf_density_threshold_code");
  frozen.config.physics.sf_temperature_threshold_k = parseFloating(
      requireString(entries, consumed, "physics.sf_temperature_threshold_k", "1.0e4"),
      "physics.sf_temperature_threshold_k");
  frozen.config.physics.sf_min_converging_flow_rate_code = parseFloating(
      requireString(entries, consumed, "physics.sf_min_converging_flow_rate_code", "0.0"),
      "physics.sf_min_converging_flow_rate_code");
  frozen.config.physics.sf_epsilon_ff = parseFloating(
      requireString(entries, consumed, "physics.sf_epsilon_ff", "0.01"),
      "physics.sf_epsilon_ff");
  frozen.config.physics.sf_min_star_particle_mass_code = parseFloating(
      requireString(entries, consumed, "physics.sf_min_star_particle_mass_code", "0.1"),
      "physics.sf_min_star_particle_mass_code");
  frozen.config.physics.sf_stochastic_spawning = parseBool(
      requireString(entries, consumed, "physics.sf_stochastic_spawning", "true"),
      "physics.sf_stochastic_spawning");
  frozen.config.physics.sf_random_seed = static_cast<std::uint64_t>(parseNumber<unsigned long long>(
      requireString(entries, consumed, "physics.sf_random_seed", "123456789"),
      "physics.sf_random_seed"));
  frozen.config.physics.fb_mode = parseFeedbackMode(
      requireString(entries, consumed, "physics.fb_mode", defaultFor("physics.fb_mode")));
  frozen.config.physics.fb_variant = parseFeedbackVariant(
      requireString(entries, consumed, "physics.fb_variant", defaultFor("physics.fb_variant")));
  frozen.config.physics.fb_use_returned_mass_budget = parseBool(
      requireString(entries, consumed, "physics.fb_use_returned_mass_budget", "true"),
      "physics.fb_use_returned_mass_budget");
  frozen.config.physics.fb_epsilon_thermal = parseFloating(
      requireString(entries, consumed, "physics.fb_epsilon_thermal", "0.6"),
      "physics.fb_epsilon_thermal");
  frozen.config.physics.fb_epsilon_kinetic = parseFloating(
      requireString(entries, consumed, "physics.fb_epsilon_kinetic", "0.3"),
      "physics.fb_epsilon_kinetic");
  frozen.config.physics.fb_epsilon_momentum = parseFloating(
      requireString(entries, consumed, "physics.fb_epsilon_momentum", "0.1"),
      "physics.fb_epsilon_momentum");
  frozen.config.physics.fb_sn_energy_erg_per_mass_code = parseFloating(
      requireString(entries, consumed, "physics.fb_sn_energy_erg_per_mass_code", "1.0e49"),
      "physics.fb_sn_energy_erg_per_mass_code");
  frozen.config.physics.fb_momentum_code_per_mass_code = parseFloating(
      requireString(entries, consumed, "physics.fb_momentum_code_per_mass_code", "3.0e3"),
      "physics.fb_momentum_code_per_mass_code");
  frozen.config.physics.fb_neighbor_count = static_cast<std::uint32_t>(parseNumber<unsigned>(
      requireString(entries, consumed, "physics.fb_neighbor_count", "8"),
      "physics.fb_neighbor_count"));
  frozen.config.physics.fb_delayed_cooling_time_code = parseFloating(
      requireString(entries, consumed, "physics.fb_delayed_cooling_time_code", "0.0"),
      "physics.fb_delayed_cooling_time_code");
  frozen.config.physics.fb_stochastic_event_probability = parseFloating(
      requireString(entries, consumed, "physics.fb_stochastic_event_probability", "0.25"),
      "physics.fb_stochastic_event_probability");
  frozen.config.physics.fb_random_seed = static_cast<std::uint64_t>(parseNumber<unsigned long long>(
      requireString(entries, consumed, "physics.fb_random_seed", "42424242"),
      "physics.fb_random_seed"));
  frozen.config.physics.stellar_evolution_table_path =
      requireString(entries, consumed, "physics.stellar_evolution_table_path", "");
  frozen.config.physics.stellar_evolution_hubble_time_years = parseFloating(
      requireString(entries, consumed, "physics.stellar_evolution_hubble_time_years", "1.44e10"),
      "physics.stellar_evolution_hubble_time_years");
  frozen.config.physics.enable_black_hole_agn = parseBool(
      requireString(entries, consumed, "physics.enable_black_hole_agn", "false"),
      "physics.enable_black_hole_agn");
  frozen.config.physics.bh_seed_halo_mass_threshold_code = parseFloating(
      requireString(entries, consumed, "physics.bh_seed_halo_mass_threshold_code", "1.0e3"),
      "physics.bh_seed_halo_mass_threshold_code");
  frozen.config.physics.bh_seed_mass_code = parseFloating(
      requireString(entries, consumed, "physics.bh_seed_mass_code", "1.0"),
      "physics.bh_seed_mass_code");
  frozen.config.physics.bh_seed_max_per_cell = static_cast<std::uint32_t>(parseNumber<unsigned>(
      requireString(entries, consumed, "physics.bh_seed_max_per_cell", "1"),
      "physics.bh_seed_max_per_cell"));
  frozen.config.physics.bh_alpha_bondi = parseFloating(
      requireString(entries, consumed, "physics.bh_alpha_bondi", "1.0"),
      "physics.bh_alpha_bondi");
  frozen.config.physics.bh_use_eddington_cap = parseBool(
      requireString(entries, consumed, "physics.bh_use_eddington_cap", "true"),
      "physics.bh_use_eddington_cap");
  frozen.config.physics.bh_epsilon_r = parseFloating(
      requireString(entries, consumed, "physics.bh_epsilon_r", "0.1"),
      "physics.bh_epsilon_r");
  frozen.config.physics.bh_epsilon_f = parseFloating(
      requireString(entries, consumed, "physics.bh_epsilon_f", "0.05"),
      "physics.bh_epsilon_f");
  frozen.config.physics.bh_feedback_coupling_efficiency = parseFloating(
      requireString(entries, consumed, "physics.bh_feedback_coupling_efficiency", "1.0"),
      "physics.bh_feedback_coupling_efficiency");
  frozen.config.physics.bh_duty_cycle_active_edd_ratio_threshold = parseFloating(
      requireString(entries, consumed, "physics.bh_duty_cycle_active_edd_ratio_threshold", "0.01"),
      "physics.bh_duty_cycle_active_edd_ratio_threshold");
  frozen.config.physics.bh_proton_mass_si = parseFloating(
      requireString(entries, consumed, "physics.bh_proton_mass_si", "1.67262192369e-27"),
      "physics.bh_proton_mass_si");
  frozen.config.physics.bh_thomson_cross_section_si = parseFloating(
      requireString(entries, consumed, "physics.bh_thomson_cross_section_si", "6.6524587321e-29"),
      "physics.bh_thomson_cross_section_si");
  frozen.config.physics.bh_newton_g_si = parseFloating(
      requireString(entries, consumed, "physics.bh_newton_g_si", "6.67430e-11"),
      "physics.bh_newton_g_si");
  frozen.config.physics.bh_speed_of_light_si = parseFloating(
      requireString(entries, consumed, "physics.bh_speed_of_light_si", "2.99792458e8"),
      "physics.bh_speed_of_light_si");
  frozen.config.physics.enable_tracers = parseBool(
      requireString(entries, consumed, "physics.enable_tracers", "false"),
      "physics.enable_tracers");
  frozen.config.physics.tracer_track_mass = parseBool(
      requireString(entries, consumed, "physics.tracer_track_mass", "true"),
      "physics.tracer_track_mass");
  frozen.config.physics.tracer_min_host_mass_code = parseFloating(
      requireString(entries, consumed, "physics.tracer_min_host_mass_code", "0.0"),
      "physics.tracer_min_host_mass_code");

  frozen.config.output.run_name =
      requireString(entries, consumed, "output.run_name", frozen.config.output.run_name);
  frozen.config.output.output_directory =
      requireString(entries, consumed, "output.output_directory", frozen.config.output.output_directory);
  frozen.config.output.output_stem = sanitizeStem(
      requireString(entries, consumed, "output.output_stem", frozen.config.output.output_stem),
      "output.output_stem");
  frozen.config.output.restart_stem = sanitizeStem(
      requireString(entries, consumed, "output.restart_stem", frozen.config.output.restart_stem),
      "output.restart_stem");
  frozen.config.output.snapshot_interval_steps = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "output.snapshot_interval_steps", "64"),
      "output.snapshot_interval_steps"));
  frozen.config.output.write_restarts = parseBool(
      requireString(entries, consumed, "output.write_restarts", "true"), "output.write_restarts");

  frozen.config.parallel.mpi_ranks_expected = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "parallel.mpi_ranks_expected", "1"),
      "parallel.mpi_ranks_expected"));
  frozen.config.parallel.omp_threads = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "parallel.omp_threads", "1"), "parallel.omp_threads"));
  frozen.config.parallel.gpu_devices = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "parallel.gpu_devices", "0"), "parallel.gpu_devices"));
  frozen.config.parallel.deterministic_reduction =
      parseBool(requireString(entries, consumed, "parallel.deterministic_reduction", "true"),
                "parallel.deterministic_reduction");

  frozen.config.analysis.enable_diagnostics = parseBool(
      requireString(entries, consumed, "analysis.enable_diagnostics", "true"),
      "analysis.enable_diagnostics");
  frozen.config.analysis.enable_halo_workflow = parseBool(
      requireString(entries, consumed, "analysis.enable_halo_workflow", "false"),
      "analysis.enable_halo_workflow");
  frozen.config.analysis.halo_on_the_fly = parseBool(
      requireString(entries, consumed, "analysis.halo_on_the_fly", "false"),
      "analysis.halo_on_the_fly");
  frozen.config.analysis.run_health_interval_steps = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "analysis.run_health_interval_steps", "1"),
      "analysis.run_health_interval_steps"));
  frozen.config.analysis.science_light_interval_steps = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "analysis.science_light_interval_steps", "8"),
      "analysis.science_light_interval_steps"));
  frozen.config.analysis.science_heavy_interval_steps = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "analysis.science_heavy_interval_steps", "64"),
      "analysis.science_heavy_interval_steps"));
  frozen.config.analysis.retention_bundle_count = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "analysis.retention_bundle_count", "8"),
      "analysis.retention_bundle_count"));
  frozen.config.analysis.power_spectrum_mesh_n = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "analysis.power_spectrum_mesh_n", "16"),
      "analysis.power_spectrum_mesh_n"));
  frozen.config.analysis.power_spectrum_bin_count = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "analysis.power_spectrum_bin_count", "12"),
      "analysis.power_spectrum_bin_count"));
  frozen.config.analysis.sf_history_bin_count = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "analysis.sf_history_bin_count", "16"),
      "analysis.sf_history_bin_count"));
  frozen.config.analysis.quicklook_grid_n = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "analysis.quicklook_grid_n", "32"),
      "analysis.quicklook_grid_n"));
  frozen.config.analysis.diagnostics_stem = sanitizeStem(
      requireString(entries, consumed, "analysis.diagnostics_stem", frozen.config.analysis.diagnostics_stem),
      "analysis.diagnostics_stem");
  frozen.config.analysis.halo_catalog_stem = sanitizeStem(
      requireString(entries, consumed, "analysis.halo_catalog_stem", frozen.config.analysis.halo_catalog_stem),
      "analysis.halo_catalog_stem");
  frozen.config.analysis.merger_tree_stem = sanitizeStem(
      requireString(entries, consumed, "analysis.merger_tree_stem", frozen.config.analysis.merger_tree_stem),
      "analysis.merger_tree_stem");
  frozen.config.analysis.halo_fof_linking_length_factor = parseFloating(
      requireString(
          entries,
          consumed,
          "analysis.halo_fof_linking_length_factor",
          "0.2"),
      "analysis.halo_fof_linking_length_factor");
  frozen.config.analysis.halo_fof_min_group_size = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "analysis.halo_fof_min_group_size", "16"),
      "analysis.halo_fof_min_group_size"));
  frozen.config.analysis.halo_include_gas = parseBool(
      requireString(entries, consumed, "analysis.halo_include_gas", "true"),
      "analysis.halo_include_gas");
  frozen.config.analysis.halo_include_stars = parseBool(
      requireString(entries, consumed, "analysis.halo_include_stars", "true"),
      "analysis.halo_include_stars");
  frozen.config.analysis.halo_include_black_holes = parseBool(
      requireString(entries, consumed, "analysis.halo_include_black_holes", "true"),
      "analysis.halo_include_black_holes");

  const bool compatible_by_file = parseBool(
      requireString(entries, consumed, "compatibility.allow_unknown_keys", "false"),
      "compatibility.allow_unknown_keys");
  frozen.config.compatibility.allow_unknown_keys = options.allow_unknown_keys || compatible_by_file;

  validateConfig(frozen.config);

  std::vector<std::string> unknown;
  std::set<std::string> known_keys;
  for (const auto& entry : configKeyRegistry()) {
    known_keys.insert(entry.key);
  }
  for (const auto& [legacy_key, _] : deprecatedAliasRegistry()) {
    known_keys.insert(legacy_key);
  }

  for (const auto& [key, _] : entries) {
    if (consumed.contains(key)) {
      continue;
    }
    if (!known_keys.contains(key)) {
      unknown.push_back(key);
      continue;
    }
    throw ConfigError("internal config parsing error: known key '" + key + "' was not consumed");
  }

  if (!unknown.empty() && !frozen.config.compatibility.allow_unknown_keys) {
    std::ostringstream stream;
    stream << "unknown parameter keys:";
    for (const std::string& key : unknown) {
      stream << ' ' << key;
    }
    stream << " (set compatibility.allow_unknown_keys=true to bypass)";
    throw ConfigError(stream.str());
  }

  frozen.provenance.config_hash = 0;
  frozen.provenance.config_hash_hex = "0000000000000000";
  frozen.normalized_text = buildNormalizedText(frozen);

  frozen.provenance.config_hash = stableConfigHash(frozen.normalized_text);
  frozen.provenance.config_hash_hex = stableConfigHashHex(frozen.normalized_text);
  frozen.normalized_text = buildNormalizedText(frozen);

  return frozen;
}

}  // namespace

ConfigError::ConfigError(const std::string& message) : std::runtime_error(message) {}

FrozenConfig loadFrozenConfigFromFile(const std::filesystem::path& path, const ParseOptions& options) {
  std::ifstream input(path);
  if (!input) {
    throw ConfigError("failed to open config file: " + path.string());
  }
  std::ostringstream contents;
  contents << input.rdbuf();
  return loadFrozenConfigFromString(contents.str(), path.string(), options);
}

FrozenConfig loadFrozenConfigFromString(
    const std::string& config_text,
    const std::string& source_name,
    const ParseOptions& options) {
  const auto entries = parseEntries(config_text);
  return normalizeValidateFreeze(entries, source_name, options);
}

void writeNormalizedConfigSnapshot(
    const FrozenConfig& frozen_config,
    const std::filesystem::path& run_directory) {
  std::filesystem::create_directories(run_directory);
  const std::filesystem::path snapshot_path = run_directory / "normalized_config.param.txt";
  std::ofstream output(snapshot_path);
  if (!output) {
    throw ConfigError("failed to write normalized config snapshot: " + snapshot_path.string());
  }
  output << "# normalized_config generated by cosmosim\n";
  output << "# source = " << frozen_config.provenance.source_name << '\n';
  output << frozen_config.normalized_text;
}

std::string modeToString(SimulationMode mode) {
  return modeToLowerString(mode);
}

std::string gravitySolverToString(GravitySolver solver) {
  switch (solver) {
    case GravitySolver::kTreePm:
      return "treepm";
  }
  throw ConfigError("unhandled GravitySolver enum value during serialization");
}

std::string hydroSolverToString(HydroSolver solver) {
  switch (solver) {
    case HydroSolver::kGodunovFv:
      return "godunov_fv";
  }
  throw ConfigError("unhandled HydroSolver enum value during serialization");
}

std::string coordinateFrameToString(CoordinateFrame frame) {
  switch (frame) {
    case CoordinateFrame::kComoving:
      return "comoving";
    case CoordinateFrame::kPhysical:
      return "physical";
  }
  throw ConfigError("unhandled CoordinateFrame enum value during serialization");
}

std::string modeHydroBoundaryToString(ModeHydroBoundary boundary) {
  switch (boundary) {
    case ModeHydroBoundary::kAuto:
      return "auto";
    case ModeHydroBoundary::kPeriodic:
      return "periodic";
    case ModeHydroBoundary::kOpen:
      return "open";
    case ModeHydroBoundary::kReflective:
      return "reflective";
  }
  throw ConfigError("unhandled ModeHydroBoundary enum value during serialization");
}

std::string modeGravityBoundaryToString(ModeGravityBoundary boundary) {
  switch (boundary) {
    case ModeGravityBoundary::kAuto:
      return "auto";
    case ModeGravityBoundary::kPeriodic:
      return "periodic";
    case ModeGravityBoundary::kIsolatedMonopole:
      return "isolated_monopole";
  }
  throw ConfigError("unhandled ModeGravityBoundary enum value during serialization");
}

std::string feedbackModeToString(FeedbackMode mode) {
  switch (mode) {
    case FeedbackMode::kThermal:
      return "thermal";
    case FeedbackMode::kKinetic:
      return "kinetic";
    case FeedbackMode::kMomentum:
      return "momentum";
    case FeedbackMode::kThermalKineticMomentum:
      return "thermal_kinetic_momentum";
  }
  throw ConfigError("unhandled FeedbackMode enum value during serialization");
}

std::string feedbackVariantToString(FeedbackVariant variant) {
  switch (variant) {
    case FeedbackVariant::kNone:
      return "none";
    case FeedbackVariant::kDelayedCooling:
      return "delayed_cooling";
    case FeedbackVariant::kStochastic:
      return "stochastic";
  }
  throw ConfigError("unhandled FeedbackVariant enum value during serialization");
}


std::string uvBackgroundModelToString(UvBackgroundModel model) {
  switch (model) {
    case UvBackgroundModel::kNone:
      return "none";
    case UvBackgroundModel::kHm12:
      return "hm12";
    case UvBackgroundModel::kFg20:
      return "fg20";
  }
  throw ConfigError("unhandled UvBackgroundModel enum value during serialization");
}

std::string selfShieldingModelToString(SelfShieldingModel model) {
  switch (model) {
    case SelfShieldingModel::kNone:
      return "none";
    case SelfShieldingModel::kRahmati13Like:
      return "rahmati13_like";
  }
  throw ConfigError("unhandled SelfShieldingModel enum value during serialization");
}

std::string coolingModelToString(CoolingModel model) {
  switch (model) {
    case CoolingModel::kPrimordial:
      return "primordial";
    case CoolingModel::kPrimordialMetalLine:
      return "primordial_metal_line";
  }
  throw ConfigError("unhandled CoolingModel enum value during serialization");
}

}  // namespace cosmosim::core
