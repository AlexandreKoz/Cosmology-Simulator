#include "cosmosim/core/config.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_mode.hpp"

namespace cosmosim::core {
namespace {

struct ParsedEntry {
  std::string value;
  int line_number = 0;
};

struct SectionRequirement {
  std::string name;
  std::vector<std::string> required_fields;
};

[[nodiscard]] const std::map<std::string, std::string>& deprecatedAliasRegistry();

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
  bool in_single_quotes = false;
  bool in_double_quotes = false;

  const auto is_whitespace_boundary = [&line](std::size_t pos) {
    return pos == 0 || std::isspace(static_cast<unsigned char>(line[pos - 1])) != 0;
  };

  for (std::size_t pos = 0; pos < line.size(); ++pos) {
    const char current = line[pos];
    if (current == '"' && !in_single_quotes) {
      in_double_quotes = !in_double_quotes;
      continue;
    }
    if (current == '\'' && !in_double_quotes) {
      in_single_quotes = !in_single_quotes;
      continue;
    }
    if (in_single_quotes || in_double_quotes) {
      continue;
    }

    if ((current == '#' || current == ';') && is_whitespace_boundary(pos)) {
      return line.substr(0, pos);
    }
    if (current == '/' && pos + 1 < line.size() && line[pos + 1] == '/' &&
        is_whitespace_boundary(pos)) {
      return line.substr(0, pos);
    }
  }

  return line;
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
    const bool allow_empty_value = (equal_pos != std::string::npos);
    if (key.empty() || (value.empty() && !allow_empty_value)) {
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

[[nodiscard]] std::set<std::string> inputSections(const std::map<std::string, ParsedEntry>& entries) {
  std::set<std::string> sections;
  for (const auto& [key, _] : entries) {
    const std::size_t dot = key.find('.');
    if (dot == std::string::npos) {
      continue;
    }
    sections.insert(key.substr(0, dot));
  }
  return sections;
}

[[nodiscard]] const std::vector<SectionRequirement>& executableSchemaRequirements() {
  static const std::vector<SectionRequirement> requirements = {
      {"mode", {"mode.mode"}},
  };
  return requirements;
}

void validateRequiredSchema(const std::map<std::string, ParsedEntry>& entries) {
  const std::set<std::string> sections = inputSections(entries);
  for (const SectionRequirement& requirement : executableSchemaRequirements()) {
    const bool has_legacy_mode_root = (requirement.name == "mode" && entries.contains("mode"));
    if (!sections.contains(requirement.name) && !has_legacy_mode_root) {
      throw ConfigError("missing required section '" + requirement.name + "'");
    }
    for (const std::string& field : requirement.required_fields) {
      const bool has_legacy_mode_field = (field == "mode.mode" && entries.contains("mode"));
      if (!entries.contains(field) && !has_legacy_mode_field) {
        throw ConfigError("missing required field '" + field + "'");
      }
    }
  }
}

void applyDeprecatedAliases(
    std::map<std::string, ParsedEntry>& entries,
    std::set<std::string>& consumed,
    FrozenConfig& frozen) {
  for (const auto& [legacy_key, canonical_key] : deprecatedAliasRegistry()) {
    const auto legacy_it = entries.find(legacy_key);
    if (legacy_it == entries.end()) {
      continue;
    }
    const auto canonical_it = entries.find(canonical_key);
    if (canonical_it != entries.end()) {
      if (trim(canonical_it->second.value) != trim(legacy_it->second.value)) {
        throw ConfigError(
            "alias conflict at path '" + canonical_key + "': keys '" + legacy_key + "' and '" +
            canonical_key + "' provide different values");
      }
      throw ConfigError("deprecated key '" + legacy_key + "' cannot be combined with '" +
                        canonical_key + "'");
    }
    entries.emplace(canonical_key, legacy_it->second);
    consumed.insert(legacy_key);
    const std::string note =
        "deprecated key '" + legacy_key + "' mapped to '" + canonical_key + "'";
    frozen.provenance.deprecation_warnings.push_back(note);
    frozen.user_config.alias_resolution_notes.push_back(note);
  }
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
  if (trimmed.empty()) {
    throw ConfigError("key '" + key + "': invalid finite floating-point value '" + value + "'");
  }
  char* consumed = nullptr;
  const double number = std::strtod(trimmed.c_str(), &consumed);
  if (consumed != trimmed.c_str() + static_cast<std::ptrdiff_t>(trimmed.size()) ||
      !std::isfinite(number)) {
    throw ConfigError("key '" + key + "': invalid finite floating-point value '" + value + "'");
  }
  return number;
}

void requireFinite(double value, const std::string& key) {
  if (!std::isfinite(value)) {
    throw ConfigError("key '" + key + "': value must be finite");
  }
}

void requireAllFinite(std::initializer_list<std::pair<double, const char*>> values) {
  for (const auto& [value, key] : values) {
    requireFinite(value, key);
  }
}

[[nodiscard]] std::pair<double, std::string> splitMagnitudeUnit(const std::string& value) {
  std::istringstream stream(value);
  double magnitude = 0.0;
  stream >> magnitude;
  if (stream.fail() || !std::isfinite(magnitude)) {
    throw ConfigError("invalid value '" + value + "': expected finite numeric magnitude");
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
  throw ConfigError("unhandled SimulationMode enum value during serialization");
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

[[nodiscard]] InitialConditionConvention parseInitialConditionConvention(
    const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "generated") {
    return InitialConditionConvention::kGenerated;
  }
  if (lower == "chui_canonical_v1") {
    return InitialConditionConvention::kChuiCanonicalV1;
  }
  if (lower == "gadget_arepo_bridge_v1") {
    return InitialConditionConvention::kGadgetArepoBridgeV1;
  }
  if (lower == "manifest_v1") {
    return InitialConditionConvention::kManifestV1;
  }
  throw ConfigError(
      "key 'mode.ic_convention': invalid value '" + value +
      "' (supported: generated, chui_canonical_v1, gadget_arepo_bridge_v1, manifest_v1)");
}

[[nodiscard]] InitialConditionSpeciesPolicy parseInitialConditionSpeciesPolicy(
    const std::string& value,
    const std::string& key) {
  const std::string lower = toLower(trim(value));
  if (lower == "reject") {
    return InitialConditionSpeciesPolicy::kReject;
  }
  if (lower == "dark_matter") {
    return InitialConditionSpeciesPolicy::kDarkMatter;
  }
  if (lower == "star") {
    return InitialConditionSpeciesPolicy::kStar;
  }
  if (lower == "black_hole") {
    return InitialConditionSpeciesPolicy::kBlackHole;
  }
  if (lower == "tracer") {
    return InitialConditionSpeciesPolicy::kTracer;
  }
  throw ConfigError(
      "key '" + key + "': invalid value '" + value +
      "' (supported: reject, dark_matter, star, black_hole, tracer)");
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

[[nodiscard]] TreePmAssignmentScheme parseTreePmAssignmentScheme(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "cic") {
    return TreePmAssignmentScheme::kCic;
  }
  if (lower == "tsc") {
    return TreePmAssignmentScheme::kTsc;
  }
  throw ConfigError(
      "key 'numerics.treepm_assignment_scheme': invalid value '" + value +
      "' (supported: cic, tsc)");
}

[[nodiscard]] TreePmOpeningCriterion parseTreePmOpeningCriterion(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "geometric") {
    return TreePmOpeningCriterion::kGeometric;
  }
  if (lower == "com_distance") {
    return TreePmOpeningCriterion::kComDistance;
  }
  if (lower == "relative_force_error") {
    return TreePmOpeningCriterion::kRelativeForceError;
  }
  throw ConfigError(
      "key 'numerics.treepm_tree_opening_criterion': invalid value '" + value +
      "' (supported: geometric, com_distance, relative_force_error)");
}

[[nodiscard]] std::string treePmOpeningCriterionToString(TreePmOpeningCriterion criterion) {
  switch (criterion) {
    case TreePmOpeningCriterion::kGeometric:
      return "geometric";
    case TreePmOpeningCriterion::kComDistance:
      return "com_distance";
    case TreePmOpeningCriterion::kRelativeForceError:
      return "relative_force_error";
  }
  throw ConfigError("unhandled TreePmOpeningCriterion enum value during serialization");
}

[[nodiscard]] PmDecompositionMode parsePmDecompositionMode(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "slab") {
    return PmDecompositionMode::kSlab;
  }
  if (lower == "pencil") {
    return PmDecompositionMode::kPencil;
  }
  throw ConfigError(
      "key 'numerics.treepm_pm_decomposition_mode': invalid value '" + value +
      "' (supported: slab, pencil)");
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

[[nodiscard]] ZoomLongRangeStrategy parseZoomLongRangeStrategy(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "disabled") {
    return ZoomLongRangeStrategy::kDisabled;
  }
  if (lower == "global_coarse_plus_focused_highres_correction") {
    return ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection;
  }
  throw ConfigError("key 'mode.zoom_long_range_strategy': invalid value '" + value + "'");
}

[[nodiscard]] std::string zoomLongRangeStrategyToString(ZoomLongRangeStrategy strategy) {
  switch (strategy) {
    case ZoomLongRangeStrategy::kDisabled:
      return "disabled";
    case ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection:
      return "global_coarse_plus_focused_highres_correction";
  }
  throw ConfigError("unhandled ZoomLongRangeStrategy enum value during serialization");
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

[[nodiscard]] IntegratorTimeVariable parseIntegratorTimeVariable(const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "scale_factor" || lower == "a") {
    return IntegratorTimeVariable::kScaleFactor;
  }
  if (lower == "ln_a" || lower == "log_scale_factor") {
    return IntegratorTimeVariable::kLogScaleFactor;
  }
  if (lower == "t_code" || lower == "code_time") {
    return IntegratorTimeVariable::kCodeTime;
  }
  if (lower == "t_phys" || lower == "physical_time") {
    return IntegratorTimeVariable::kPhysicalTime;
  }
  throw ConfigError(
      "key 'numerics.integrator_time_variable': invalid value '" + value +
      "' (supported: scale_factor, ln_a, t_code, t_phys)");
}

[[nodiscard]] double redshiftToScaleFactor(double z, const std::string& key) {
  if (!std::isfinite(z) || z <= -1.0) {
    throw ConfigError("key '" + key + "': redshift must be finite and > -1");
  }
  return 1.0 / (1.0 + z);
}

[[nodiscard]] double scaleFactorToRedshift(double a, const std::string& key) {
  if (!std::isfinite(a) || a <= 0.0) {
    throw ConfigError("key '" + key + "': scale factor must be finite and > 0");
  }
  return (1.0 / a) - 1.0;
}

void validateScaleRedshiftPair(double a, double z, const std::string& a_key, const std::string& z_key) {
  const double expected_a = redshiftToScaleFactor(z, z_key);
  if (!std::isfinite(a) || a <= 0.0) {
    throw ConfigError("key '" + a_key + "': scale factor must be finite and > 0");
  }
  const double tolerance = 16.0 * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::abs(a), std::abs(expected_a)});
  if (std::abs(expected_a - a) > tolerance) {
    throw ConfigError("keys '" + a_key + "' and '" + z_key +
                      "' are inconsistent: expected " + a_key + " = 1/(1+" + z_key + ")");
  }
}

void normalizeScaleRedshiftPair(
    std::map<std::string, ParsedEntry>& entries,
    std::set<std::string>& consumed,
    const std::string& a_key,
    const std::string& z_key,
    double default_a,
    double default_z,
    double& out_a,
    double& out_z) {
  const bool has_a = entries.contains(a_key);
  const bool has_z = entries.contains(z_key);

  if (has_a && has_z) {
    out_a = parseFloating(requireString(entries, consumed, a_key, std::to_string(default_a)), a_key);
    out_z = parseFloating(requireString(entries, consumed, z_key, std::to_string(default_z)), z_key);
    validateScaleRedshiftPair(out_a, out_z, a_key, z_key);
    return;
  }

  if (has_z) {
    out_z = parseFloating(requireString(entries, consumed, z_key, std::to_string(default_z)), z_key);
    out_a = redshiftToScaleFactor(out_z, z_key);
    return;
  }

  if (has_a) {
    out_a = parseFloating(requireString(entries, consumed, a_key, std::to_string(default_a)), a_key);
    out_z = scaleFactorToRedshift(out_a, a_key);
    return;
  }

  out_a = default_a;
  out_z = default_z;
  validateScaleRedshiftPair(out_a, out_z, a_key, z_key);
}

[[nodiscard]] AnalysisConfig::DiagnosticsExecutionPolicy parseDiagnosticsExecutionPolicy(
    const std::string& value) {
  const std::string lower = toLower(trim(value));
  if (lower == "run_health_only") {
    return AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthOnly;
  }
  if (lower == "run_health_and_light_science") {
    return AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthAndLightScience;
  }
  if (lower == "all_including_provisional") {
    return AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional;
  }
  throw ConfigError(
      "analysis.diagnostics_execution_policy must be one of: run_health_only, "
      "run_health_and_light_science, all_including_provisional");
}

[[nodiscard]] std::string diagnosticsExecutionPolicyToString(
    AnalysisConfig::DiagnosticsExecutionPolicy policy) {
  switch (policy) {
    case AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthOnly:
      return "run_health_only";
    case AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthAndLightScience:
      return "run_health_and_light_science";
    case AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional:
      return "all_including_provisional";
  }
  throw ConfigError("unhandled DiagnosticsExecutionPolicy enum value during serialization");
}

[[nodiscard]] std::string pmDecompositionModeToString(PmDecompositionMode mode) {
  switch (mode) {
    case PmDecompositionMode::kSlab:
      return "slab";
    case PmDecompositionMode::kPencil:
      return "pencil";
  }
  throw ConfigError("unhandled PmDecompositionMode enum value during serialization");
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
      {"mode.ic_file", "generated"},
      {"mode.ic_convention", "generated"},
      {"mode.ic_manifest_file", ""},
      {"mode.ic_chunk_particle_count", "65536"},
      {"mode.ic_staging_particle_count", "65536"},
      {"mode.ic_part_type2_policy", "reject"},
      {"mode.ic_part_type3_policy", "reject"},
      {"mode.zoom_high_res_region", "false"},
      {"mode.zoom_region_file", ""},
      {"mode.zoom_long_range_strategy", "disabled"},
      {"mode.zoom_region_center_x", "0.5"},
      {"mode.zoom_region_center_y", "0.5"},
      {"mode.zoom_region_center_z", "0.5"},
      {"mode.zoom_region_radius", "0.0"},
      {"mode.zoom_focused_pm_grid_nx", "0"},
      {"mode.zoom_focused_pm_grid_ny", "0"},
      {"mode.zoom_focused_pm_grid_nz", "0"},
      {"mode.zoom_contamination_radius", "0.0"},
      {"mode.hydro_boundary", "auto"},
      {"mode.gravity_boundary", "auto"},
      {"cosmology.omega_matter", "0.315"},
      {"cosmology.omega_lambda", "0.685"},
      {"cosmology.omega_baryon", "0.049"},
      {"cosmology.hubble_param", "0.674"},
      {"cosmology.sigma8", "0.811"},
      {"cosmology.scalar_index_ns", "0.965"},
      {"cosmology.box_size_x", "50.0"},
      {"cosmology.box_size_y", "50.0"},
      {"cosmology.box_size_z", "50.0"},
      {"cosmology.box_size", "50.0"},
      {"numerics.a_begin", "1.0"},
      {"numerics.a_end", "1.0"},
      {"numerics.z_begin", "0.0"},
      {"numerics.z_end", "0.0"},
      {"numerics.t_code_begin", "0.0"},
      {"numerics.t_code_end", "1.0"},
      {"numerics.t_phys_begin", "0.0"},
      {"numerics.t_phys_end", "0.0"},
      {"numerics.integrator_time_variable", "scale_factor"},
      {"numerics.cosmology_max_delta_ln_a", "1.0e-2"},
      {"numerics.cosmology_max_hubble_time_fraction", "1.0e-2"},
      {"numerics.source_max_fractional_change", "0.1"},
      {"numerics.max_global_steps", "1024"},
      {"numerics.hierarchical_max_rung", "0"},
      {"numerics.amr_max_level", "10"},
      {"numerics.gravity_softening", "1.0 kpc"},
      {"numerics.gravity_softening_gas", ""},
      {"numerics.gravity_softening_dark_matter", ""},
      {"numerics.gravity_softening_star", ""},
      {"numerics.gravity_softening_black_hole", ""},
      {"numerics.gravity_softening_tracer", ""},
      {"numerics.gravity_solver", "treepm"},
      {"numerics.hydro_solver", "godunov_fv"},
      {"numerics.treepm_pm_grid_nx", "16"},
      {"numerics.treepm_pm_grid_ny", "16"},
      {"numerics.treepm_pm_grid_nz", "16"},
      {"numerics.treepm_pm_grid", "16"},
      {"numerics.treepm_asmth_cells", "1.25"},
      {"numerics.treepm_rcut_cells", "6.25"},
      {"numerics.treepm_tree_opening_criterion", "com_distance"},
      {"numerics.treepm_tree_opening_theta", "0.7"},
      {"numerics.treepm_tree_relative_force_tolerance", "0.005"},
      {"numerics.treepm_tree_relative_force_acceleration_floor", "1.0e-30"},
      {"numerics.treepm_assignment_scheme", "tsc"},
      {"numerics.treepm_enable_window_deconvolution", "true"},
      {"numerics.treepm_update_cadence_steps", "1"},
      {"numerics.treepm_pm_decomposition_mode", "slab"},
      {"numerics.treepm_tree_exchange_batch_bytes", "4194304"},
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
      {"output.snapshot_interval_time_code", "0.0"},
      {"output.write_restarts", "true"},
      {"parallel.mpi_ranks_expected", "1"},
      {"parallel.omp_threads", "1"},
      {"parallel.gpu_devices", "0"},
      {"parallel.deterministic_reduction", "true"},
      {"parallel.decomposition_particle_count_weight", "1.0"},
      {"parallel.decomposition_gas_cell_weight", "1.5"},
      {"parallel.decomposition_tree_interaction_weight", "1.0"},
      {"parallel.decomposition_pm_mesh_weight", "0.25"},
      {"parallel.decomposition_amr_patch_weight", "1.0"},
      {"parallel.decomposition_active_fraction_weight", "2.0"},
      {"parallel.decomposition_memory_pressure_weight", "9.5367431640625e-7"},
      {"parallel.decomposition_gpu_occupancy_weight", "0.0"},
      {"parallel.decomposition_generic_work_weight", "0.5"},
      {"parallel.decomposition_runtime_rebalance_enabled", "true"},
      {"parallel.decomposition_debug_exact_ownership_audit", "false"},
      {"parallel.decomposition_rebalance_imbalance_trigger", "1.25"},
      {"parallel.decomposition_rebalance_memory_trigger", "1.50"},
      {"parallel.decomposition_rebalance_max_migrated_load_fraction", "0.25"},
      {"parallel.decomposition_measured_tree_pair_weight", "1.0"},
      {"parallel.decomposition_measured_pm_cell_weight", "1.0"},
      {"parallel.decomposition_measured_amr_cell_weight", "1.0"},
      {"parallel.decomposition_measured_hydro_face_weight", "1.0"},
      {"parallel.decomposition_measured_wall_ms_weight", "1.0"},
      {"parallel.isolated_pm_root_workspace_limit_bytes", "268435456"},
      {"parallel.zoom_high_res_allgather_limit_bytes", "268435456"},
      {"analysis.enable_diagnostics", "true"},
      {"analysis.enable_halo_workflow", "false"},
      {"analysis.halo_on_the_fly", "false"},
      {"analysis.diagnostics_execution_policy", "run_health_and_light_science"},
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
      {"timemax", "numerics.t_code_end"},
      {"mode", "mode.mode"},
      {"run_name", "output.run_name"},
      {"numerics.initial_scale_factor", "numerics.a_begin"},
      {"numerics.initial_redshift", "numerics.z_begin"},
      {"numerics.time_begin_code", "numerics.t_code_begin"},
      {"numerics.time_end_code", "numerics.t_code_end"},
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

  requireAllFinite({
      {config.cosmology.omega_matter, "cosmology.omega_matter"},
      {config.cosmology.omega_lambda, "cosmology.omega_lambda"},
      {config.cosmology.omega_baryon, "cosmology.omega_baryon"},
      {config.cosmology.hubble_param, "cosmology.hubble_param"},
      {config.cosmology.sigma8, "cosmology.sigma8"},
      {config.cosmology.scalar_index_ns, "cosmology.scalar_index_ns"},
      {config.cosmology.box_size_x_mpc_comoving, "cosmology.box_size_x"},
      {config.cosmology.box_size_y_mpc_comoving, "cosmology.box_size_y"},
      {config.cosmology.box_size_z_mpc_comoving, "cosmology.box_size_z"},
      {config.cosmology.box_size_mpc_comoving, "cosmology.box_size"},
      {config.numerics.a_begin, "numerics.a_begin"},
      {config.numerics.a_end, "numerics.a_end"},
      {config.numerics.z_begin, "numerics.z_begin"},
      {config.numerics.z_end, "numerics.z_end"},
      {config.numerics.t_code_begin, "numerics.t_code_begin"},
      {config.numerics.t_code_end, "numerics.t_code_end"},
      {config.numerics.t_phys_begin, "numerics.t_phys_begin"},
      {config.numerics.t_phys_end, "numerics.t_phys_end"},
      {config.numerics.cosmology_max_delta_ln_a, "numerics.cosmology_max_delta_ln_a"},
      {config.numerics.cosmology_max_hubble_time_fraction, "numerics.cosmology_max_hubble_time_fraction"},
      {config.numerics.source_max_fractional_change, "numerics.source_max_fractional_change"},
      {config.numerics.gravity_softening_kpc_comoving, "numerics.gravity_softening"},
      {config.numerics.gravity_softening_gas_kpc_comoving, "numerics.gravity_softening_gas"},
      {config.numerics.gravity_softening_dark_matter_kpc_comoving, "numerics.gravity_softening_dark_matter"},
      {config.numerics.gravity_softening_star_kpc_comoving, "numerics.gravity_softening_star"},
      {config.numerics.gravity_softening_black_hole_kpc_comoving, "numerics.gravity_softening_black_hole"},
      {config.numerics.gravity_softening_tracer_kpc_comoving, "numerics.gravity_softening_tracer"},
      {config.numerics.treepm_asmth_cells, "numerics.treepm_asmth_cells"},
      {config.numerics.treepm_rcut_cells, "numerics.treepm_rcut_cells"},
      {config.numerics.treepm_tree_opening_theta, "numerics.treepm_tree_opening_theta"},
      {config.numerics.treepm_tree_relative_force_tolerance,
       "numerics.treepm_tree_relative_force_tolerance"},
      {config.numerics.treepm_tree_relative_force_acceleration_floor,
       "numerics.treepm_tree_relative_force_acceleration_floor"},
      {config.mode.zoom_region_center_x_mpc_comoving, "mode.zoom_region_center_x"},
      {config.mode.zoom_region_center_y_mpc_comoving, "mode.zoom_region_center_y"},
      {config.mode.zoom_region_center_z_mpc_comoving, "mode.zoom_region_center_z"},
      {config.mode.zoom_region_radius_mpc_comoving, "mode.zoom_region_radius"},
      {config.mode.zoom_contamination_radius_mpc_comoving, "mode.zoom_contamination_radius"},
      {config.physics.temperature_floor_k, "physics.temperature_floor_k"},
      {config.physics.sf_density_threshold_code, "physics.sf_density_threshold_code"},
      {config.physics.sf_temperature_threshold_k, "physics.sf_temperature_threshold_k"},
      {config.physics.sf_min_converging_flow_rate_code, "physics.sf_min_converging_flow_rate_code"},
      {config.physics.sf_epsilon_ff, "physics.sf_epsilon_ff"},
      {config.physics.sf_min_star_particle_mass_code, "physics.sf_min_star_particle_mass_code"},
      {config.physics.fb_epsilon_thermal, "physics.fb_epsilon_thermal"},
      {config.physics.fb_epsilon_kinetic, "physics.fb_epsilon_kinetic"},
      {config.physics.fb_epsilon_momentum, "physics.fb_epsilon_momentum"},
      {config.physics.fb_sn_energy_erg_per_mass_code, "physics.fb_sn_energy_erg_per_mass_code"},
      {config.physics.fb_momentum_code_per_mass_code, "physics.fb_momentum_code_per_mass_code"},
      {config.physics.fb_delayed_cooling_time_code, "physics.fb_delayed_cooling_time_code"},
      {config.physics.fb_stochastic_event_probability, "physics.fb_stochastic_event_probability"},
      {config.physics.stellar_evolution_hubble_time_years, "physics.stellar_evolution_hubble_time_years"},
      {config.physics.bh_seed_halo_mass_threshold_code, "physics.bh_seed_halo_mass_threshold_code"},
      {config.physics.bh_seed_mass_code, "physics.bh_seed_mass_code"},
      {config.physics.bh_alpha_bondi, "physics.bh_alpha_bondi"},
      {config.physics.bh_epsilon_r, "physics.bh_epsilon_r"},
      {config.physics.bh_epsilon_f, "physics.bh_epsilon_f"},
      {config.physics.bh_feedback_coupling_efficiency, "physics.bh_feedback_coupling_efficiency"},
      {config.physics.bh_duty_cycle_active_edd_ratio_threshold, "physics.bh_duty_cycle_active_edd_ratio_threshold"},
      {config.physics.bh_proton_mass_si, "physics.bh_proton_mass_si"},
      {config.physics.bh_thomson_cross_section_si, "physics.bh_thomson_cross_section_si"},
      {config.physics.bh_newton_g_si, "physics.bh_newton_g_si"},
      {config.physics.bh_speed_of_light_si, "physics.bh_speed_of_light_si"},
      {config.physics.tracer_min_host_mass_code, "physics.tracer_min_host_mass_code"},
      {config.analysis.halo_fof_linking_length_factor, "analysis.halo_fof_linking_length_factor"},
  });

  if (config.numerics.t_code_end <= config.numerics.t_code_begin) {
    throw ConfigError("numerics.t_code_end must be greater than numerics.t_code_begin");
  }
  if (config.numerics.t_phys_end < config.numerics.t_phys_begin) {
    throw ConfigError("numerics.t_phys_end must be greater than or equal to numerics.t_phys_begin");
  }
  validateScaleRedshiftPair(
      config.numerics.a_begin,
      config.numerics.z_begin,
      "numerics.a_begin",
      "numerics.z_begin");
  validateScaleRedshiftPair(
      config.numerics.a_end,
      config.numerics.z_end,
      "numerics.a_end",
      "numerics.z_end");
  if (config.numerics.a_end < config.numerics.a_begin) {
    throw ConfigError(
        "numerics.a_end must be >= numerics.a_begin for forward cosmological integration; "
        "reverse integration requires an explicit future config mode");
  }
  if (config.numerics.z_end > config.numerics.z_begin) {
    throw ConfigError(
        "numerics.z_end must be <= numerics.z_begin for forward cosmological integration; "
        "reverse integration requires an explicit future config mode");
  }
  if (config.numerics.cosmology_max_delta_ln_a <= 0.0 ||
      config.numerics.cosmology_max_hubble_time_fraction <= 0.0 ||
      config.numerics.source_max_fractional_change <= 0.0 ||
      config.numerics.source_max_fractional_change > 1.0) {
    throw ConfigError("numerics cosmology/source timestep limits must be positive; source_max_fractional_change must be <= 1");
  }
  if (config.numerics.max_global_steps <= 0) {
    throw ConfigError("numerics.max_global_steps must be > 0");
  }
  if (config.numerics.hierarchical_max_rung != 0) {
    throw ConfigError(
        "numerics.hierarchical_max_rung must be 0: production ReferenceWorkflow "
        "does not yet carry per-element kick/drift epochs for mixed-rung KDK integration");
  }
  if (config.output.snapshot_interval_steps < 0) {
    throw ConfigError("output.snapshot_interval_steps must be >= 0");
  }
  if (!std::isfinite(config.output.snapshot_interval_time_code) ||
      config.output.snapshot_interval_time_code < 0.0) {
    throw ConfigError("output.snapshot_interval_time_code must be finite and >= 0");
  }
  if (config.output.snapshot_interval_steps == 0 &&
      config.output.snapshot_interval_time_code == 0.0) {
    throw ConfigError(
        "output requires snapshot_interval_steps > 0 or snapshot_interval_time_code > 0");
  }
  if (config.parallel.mpi_ranks_expected <= 0 || config.parallel.omp_threads <= 0) {
    throw ConfigError("parallel settings require positive mpi_ranks_expected and omp_threads");
  }
  if (config.parallel.gpu_devices < 0) {
    throw ConfigError("parallel.gpu_devices must be >= 0");
  }
  if (config.parallel.decomposition_particle_count_weight < 0.0 ||
      config.parallel.decomposition_gas_cell_weight < 0.0 ||
      config.parallel.decomposition_tree_interaction_weight < 0.0 ||
      config.parallel.decomposition_pm_mesh_weight < 0.0 ||
      config.parallel.decomposition_amr_patch_weight < 0.0 ||
      config.parallel.decomposition_active_fraction_weight < 0.0 ||
      config.parallel.decomposition_memory_pressure_weight < 0.0 ||
      config.parallel.decomposition_gpu_occupancy_weight < 0.0 ||
      config.parallel.decomposition_generic_work_weight < 0.0 ||
      config.parallel.decomposition_rebalance_imbalance_trigger < 1.0 ||
      config.parallel.decomposition_rebalance_memory_trigger < 1.0 ||
      config.parallel.decomposition_rebalance_max_migrated_load_fraction < 0.0 ||
      config.parallel.decomposition_rebalance_max_migrated_load_fraction > 1.0 ||
      config.parallel.decomposition_measured_tree_pair_weight < 0.0 ||
      config.parallel.decomposition_measured_pm_cell_weight < 0.0 ||
      config.parallel.decomposition_measured_amr_cell_weight < 0.0 ||
      config.parallel.decomposition_measured_hydro_face_weight < 0.0 ||
      config.parallel.decomposition_measured_wall_ms_weight < 0.0) {
    throw ConfigError("parallel decomposition work/rebalance weights must be non-negative; rebalance triggers must be >= 1");
  }
  if (config.parallel.isolated_pm_root_workspace_limit_bytes == 0U ||
      config.parallel.zoom_high_res_allgather_limit_bytes == 0U) {
    throw ConfigError("parallel PM gather guard byte limits must be > 0");
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
  if (config.cosmology.omega_matter <= 0.0 || config.cosmology.omega_lambda < 0.0 ||
      config.cosmology.omega_baryon < 0.0 || config.cosmology.hubble_param <= 0.0 ||
      config.cosmology.sigma8 <= 0.0 || config.cosmology.scalar_index_ns <= 0.0) {
    throw ConfigError(
        "cosmology requires omega_matter > 0, omega_lambda >= 0, omega_baryon >= 0, "
        "hubble_param > 0, sigma8 > 0, and scalar_index_ns > 0");
  }
  if (config.cosmology.omega_baryon > config.cosmology.omega_matter) {
    throw ConfigError("cosmology.omega_baryon must be <= cosmology.omega_matter");
  }
  auto requireOptionalPositiveSoftening = [](double value, const char* key) {
    if (value == 0.0) {
      throw ConfigError(std::string(key) + " must be > 0 when specified");
    }
  };
  requireOptionalPositiveSoftening(config.numerics.gravity_softening_gas_kpc_comoving, "numerics.gravity_softening_gas");
  requireOptionalPositiveSoftening(config.numerics.gravity_softening_dark_matter_kpc_comoving, "numerics.gravity_softening_dark_matter");
  requireOptionalPositiveSoftening(config.numerics.gravity_softening_star_kpc_comoving, "numerics.gravity_softening_star");
  requireOptionalPositiveSoftening(config.numerics.gravity_softening_black_hole_kpc_comoving, "numerics.gravity_softening_black_hole");
  requireOptionalPositiveSoftening(config.numerics.gravity_softening_tracer_kpc_comoving, "numerics.gravity_softening_tracer");
  if (config.physics.enable_feedback && !config.physics.enable_star_formation) {
    throw ConfigError(
        "physics.enable_feedback requires physics.enable_star_formation=true or an explicit future external-source policy");
  }
  if (config.physics.enable_black_hole_agn && !config.physics.enable_feedback) {
    throw ConfigError("physics.enable_black_hole_agn requires physics.enable_feedback=true");
  }
  if (config.physics.enable_cooling &&
      config.physics.cooling_model == CoolingModel::kPrimordialMetalLine &&
      trim(config.physics.metal_line_table_path).empty()) {
    throw ConfigError(
        "physics.cooling_model=primordial_metal_line requires physics.metal_line_table_path");
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
  if (config.cosmology.box_size_x_mpc_comoving <= 0.0 ||
      config.cosmology.box_size_y_mpc_comoving <= 0.0 ||
      config.cosmology.box_size_z_mpc_comoving <= 0.0) {
    throw ConfigError("cosmology.box_size_{x,y,z} must all be > 0");
  }
  if (config.numerics.treepm_pm_grid_nx <= 0 ||
      config.numerics.treepm_pm_grid_ny <= 0 ||
      config.numerics.treepm_pm_grid_nz <= 0) {
    throw ConfigError("numerics.treepm_pm_grid_n{xyz} must all be > 0");
  }
  if (!std::isfinite(config.numerics.treepm_asmth_cells) ||
      config.numerics.treepm_asmth_cells <= 0.0) {
    throw ConfigError("numerics.treepm_asmth_cells must be finite and > 0");
  }
  if (!std::isfinite(config.numerics.treepm_rcut_cells) ||
      config.numerics.treepm_rcut_cells <= 0.0) {
    throw ConfigError("numerics.treepm_rcut_cells must be finite and > 0");
  }
  if (!std::isfinite(config.numerics.treepm_tree_opening_theta) ||
      config.numerics.treepm_tree_opening_theta <= 0.0) {
    throw ConfigError("numerics.treepm_tree_opening_theta must be finite and > 0");
  }
  if (!std::isfinite(config.numerics.treepm_tree_relative_force_tolerance) ||
      config.numerics.treepm_tree_relative_force_tolerance <= 0.0) {
    throw ConfigError(
        "numerics.treepm_tree_relative_force_tolerance must be finite and > 0");
  }
  if (!std::isfinite(config.numerics.treepm_tree_relative_force_acceleration_floor) ||
      config.numerics.treepm_tree_relative_force_acceleration_floor <= 0.0) {
    throw ConfigError(
        "numerics.treepm_tree_relative_force_acceleration_floor must be finite and > 0");
  }
  if (config.numerics.treepm_update_cadence_steps != 1) {
    throw ConfigError(
        "numerics.treepm_update_cadence_steps must be exactly 1 in production: PM mesh reuse "
        "does not yet validate source, scale-factor, or decomposition epochs");
  }
  if (config.numerics.treepm_tree_exchange_batch_bytes == 0) {
    throw ConfigError("numerics.treepm_tree_exchange_batch_bytes must be > 0");
  }
  if (config.mode.ic_chunk_particle_count == 0U) {
    throw ConfigError("mode.ic_chunk_particle_count must be > 0");
  }
  if (config.mode.ic_staging_particle_count == 0U) {
    throw ConfigError("mode.ic_staging_particle_count must be > 0");
  }
  if (config.mode.ic_convention == InitialConditionConvention::kGenerated) {
    if (config.mode.ic_file != "generated") {
      throw ConfigError(
          "mode.ic_convention=generated requires mode.ic_file=generated; external HDF5 input must select an explicit versioned convention");
    }
    if (!config.mode.ic_manifest_file.empty()) {
      throw ConfigError(
          "mode.ic_manifest_file must be empty when mode.ic_convention=generated");
    }
  } else {
    if (config.mode.ic_file.empty() || config.mode.ic_file == "generated") {
      throw ConfigError(
          "external initial-condition conventions require a non-generated mode.ic_file path");
    }
    if (config.mode.ic_convention == InitialConditionConvention::kManifestV1 &&
        config.mode.ic_manifest_file.empty()) {
      throw ConfigError(
          "mode.ic_convention=manifest_v1 requires mode.ic_manifest_file");
    }
  }
  if (config.mode.zoom_region_radius_mpc_comoving < 0.0) {
    throw ConfigError("mode.zoom_region_radius must be >= 0");
  }
  if (config.mode.zoom_contamination_radius_mpc_comoving < 0.0) {
    throw ConfigError("mode.zoom_contamination_radius must be >= 0");
  }
  if (config.mode.zoom_high_res_region &&
      config.mode.zoom_region_radius_mpc_comoving <= 0.0) {
    throw ConfigError("mode.zoom_region_radius must be > 0 when mode.zoom_high_res_region is true");
  }
  if (config.mode.zoom_long_range_strategy ==
      ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection) {
    if (!config.mode.zoom_high_res_region) {
      throw ConfigError(
          "mode.zoom_long_range_strategy requires mode.zoom_high_res_region=true");
    }
    if (config.mode.zoom_focused_pm_grid_nx <= 0 ||
        config.mode.zoom_focused_pm_grid_ny <= 0 ||
        config.mode.zoom_focused_pm_grid_nz <= 0) {
      throw ConfigError(
          "mode.zoom_focused_pm_grid_n{xyz} must be > 0 for focused zoom PM correction");
    }
  }
  const ModePolicy policy = buildModePolicy(config.mode);
  validateModePolicy(config, policy);
  if (policy.gravity_boundary == GravityBoundaryModel::kPeriodicPoisson) {
    const double mesh_spacing_x = config.cosmology.box_size_x_mpc_comoving /
        static_cast<double>(config.numerics.treepm_pm_grid_nx);
    const double mesh_spacing_y = config.cosmology.box_size_y_mpc_comoving /
        static_cast<double>(config.numerics.treepm_pm_grid_ny);
    const double mesh_spacing_z = config.cosmology.box_size_z_mpc_comoving /
        static_cast<double>(config.numerics.treepm_pm_grid_nz);
    const double representative_mesh_spacing =
        std::cbrt(mesh_spacing_x * mesh_spacing_y * mesh_spacing_z);
    const double cutoff_radius =
        config.numerics.treepm_rcut_cells * representative_mesh_spacing;
    const double shortest_box_axis = std::min({
        config.cosmology.box_size_x_mpc_comoving,
        config.cosmology.box_size_y_mpc_comoving,
        config.cosmology.box_size_z_mpc_comoving,
    });
    const double half_shortest_axis = 0.5 * shortest_box_axis;
    if (cutoff_radius >= std::nextafter(half_shortest_axis, 0.0)) {
      throw ConfigError(
          "periodic TreePM requires treepm_rcut_cells times the representative mesh spacing "
          "to be < half the shortest box axis; cutoff_radius=" +
          std::to_string(cutoff_radius) + ", half_shortest_axis=" +
          std::to_string(half_shortest_axis) + ", pm_grid=" +
          std::to_string(config.numerics.treepm_pm_grid_nx) + "x" +
          std::to_string(config.numerics.treepm_pm_grid_ny) + "x" +
          std::to_string(config.numerics.treepm_pm_grid_nz) +
          "; increase the PM grid or reduce treepm_rcut_cells");
    }
  }
  if (policy.gravity_boundary ==
          GravityBoundaryModel::kIsolatedMonopoleDirichlet &&
      config.numerics.treepm_enable_window_deconvolution) {
    throw ConfigError(
        "numerics.treepm_enable_window_deconvolution must be false for "
        "isolated/open gravity boundaries because no periodic assignment "
        "window is defined for that solver path");
  }
}


[[nodiscard]] std::string buildNormalizedText(const FrozenConfig& frozen) {
  std::ostringstream stream;
  stream << "schema_version = " << frozen.config.schema_version << '\n';
  stream << "\n[units]\n";
  stream << "length_unit = " << frozen.config.units.length_unit << '\n';
  stream << "mass_unit = " << frozen.config.units.mass_unit << '\n';
  stream << "velocity_unit = " << frozen.config.units.velocity_unit << '\n';
  stream << "coordinate_frame = " << coordinateFrameToString(frozen.config.units.coordinate_frame) << '\n';
  stream << "\n[mode]\n";
  stream << "mode = " << modeToLowerString(frozen.config.mode.mode) << '\n';
  stream << "ic_file = " << frozen.config.mode.ic_file << '\n';
  stream << "ic_convention = "
         << initialConditionConventionToString(frozen.config.mode.ic_convention)
         << '\n';
  stream << "ic_manifest_file = " << frozen.config.mode.ic_manifest_file << '\n';
  stream << "ic_chunk_particle_count = "
         << frozen.config.mode.ic_chunk_particle_count << '\n';
  stream << "ic_staging_particle_count = "
         << frozen.config.mode.ic_staging_particle_count << '\n';
  stream << "ic_part_type2_policy = "
         << initialConditionSpeciesPolicyToString(
                frozen.config.mode.ic_part_type2_policy)
         << '\n';
  stream << "ic_part_type3_policy = "
         << initialConditionSpeciesPolicyToString(
                frozen.config.mode.ic_part_type3_policy)
         << '\n';
  stream << "zoom_high_res_region = " << (frozen.config.mode.zoom_high_res_region ? "true" : "false")
         << '\n';
  stream << "zoom_region_file = " << frozen.config.mode.zoom_region_file << '\n';
  stream << "zoom_long_range_strategy = " << zoomLongRangeStrategyToString(frozen.config.mode.zoom_long_range_strategy) << '\n';
  stream << "zoom_region_center_x = " << frozen.config.mode.zoom_region_center_x_mpc_comoving << '\n';
  stream << "zoom_region_center_y = " << frozen.config.mode.zoom_region_center_y_mpc_comoving << '\n';
  stream << "zoom_region_center_z = " << frozen.config.mode.zoom_region_center_z_mpc_comoving << '\n';
  stream << "zoom_region_radius = " << frozen.config.mode.zoom_region_radius_mpc_comoving << '\n';
  stream << "zoom_focused_pm_grid_nx = " << frozen.config.mode.zoom_focused_pm_grid_nx << '\n';
  stream << "zoom_focused_pm_grid_ny = " << frozen.config.mode.zoom_focused_pm_grid_ny << '\n';
  stream << "zoom_focused_pm_grid_nz = " << frozen.config.mode.zoom_focused_pm_grid_nz << '\n';
  stream << "zoom_contamination_radius = " << frozen.config.mode.zoom_contamination_radius_mpc_comoving << '\n';
  stream << "hydro_boundary = " << modeHydroBoundaryToString(frozen.config.mode.hydro_boundary) << '\n';
  stream << "gravity_boundary = " << modeGravityBoundaryToString(frozen.config.mode.gravity_boundary) << '\n';
  stream << "\n[cosmology]\n";
  stream << "omega_matter = " << frozen.config.cosmology.omega_matter << '\n';
  stream << "omega_lambda = " << frozen.config.cosmology.omega_lambda << '\n';
  stream << "omega_baryon = " << frozen.config.cosmology.omega_baryon << '\n';
  stream << "hubble_param = " << frozen.config.cosmology.hubble_param << '\n';
  stream << "sigma8 = " << frozen.config.cosmology.sigma8 << '\n';
  stream << "scalar_index_ns = " << frozen.config.cosmology.scalar_index_ns << '\n';
  stream << "box_size_x = " << frozen.config.cosmology.box_size_x_mpc_comoving << " mpc\n";
  stream << "box_size_y = " << frozen.config.cosmology.box_size_y_mpc_comoving << " mpc\n";
  stream << "box_size_z = " << frozen.config.cosmology.box_size_z_mpc_comoving << " mpc\n";
  stream << "\n[numerics]\n";
  stream << "a_begin = " << frozen.config.numerics.a_begin << '\n';
  stream << "a_end = " << frozen.config.numerics.a_end << '\n';
  stream << "z_begin = " << frozen.config.numerics.z_begin << '\n';
  stream << "z_end = " << frozen.config.numerics.z_end << '\n';
  stream << "t_code_begin = " << frozen.config.numerics.t_code_begin << '\n';
  stream << "t_code_end = " << frozen.config.numerics.t_code_end << '\n';
  stream << "t_phys_begin = " << frozen.config.numerics.t_phys_begin << '\n';
  stream << "t_phys_end = " << frozen.config.numerics.t_phys_end << '\n';
  stream << "integrator_time_variable = "
         << integratorTimeVariableToString(frozen.config.numerics.integrator_time_variable) << '\n';
  stream << "cosmology_max_delta_ln_a = " << frozen.config.numerics.cosmology_max_delta_ln_a << '\n';
  stream << "cosmology_max_hubble_time_fraction = " << frozen.config.numerics.cosmology_max_hubble_time_fraction << '\n';
  stream << "source_max_fractional_change = " << frozen.config.numerics.source_max_fractional_change << '\n';
  stream << "max_global_steps = " << frozen.config.numerics.max_global_steps << '\n';
  stream << "hierarchical_max_rung = " << frozen.config.numerics.hierarchical_max_rung << '\n';
  stream << "amr_max_level = " << frozen.config.numerics.amr_max_level << '\n';
  stream << "gravity_softening = " << frozen.config.numerics.gravity_softening_kpc_comoving << " kpc\n";
  auto emitOptionalSoftening = [&](std::string_view key, double value_kpc_comoving) {
    stream << key << " = ";
    if (value_kpc_comoving > 0.0) {
      stream << value_kpc_comoving << " kpc";
    }
    stream << '\n';
  };
  emitOptionalSoftening("gravity_softening_gas", frozen.config.numerics.gravity_softening_gas_kpc_comoving);
  emitOptionalSoftening("gravity_softening_dark_matter", frozen.config.numerics.gravity_softening_dark_matter_kpc_comoving);
  emitOptionalSoftening("gravity_softening_star", frozen.config.numerics.gravity_softening_star_kpc_comoving);
  emitOptionalSoftening("gravity_softening_black_hole", frozen.config.numerics.gravity_softening_black_hole_kpc_comoving);
  emitOptionalSoftening("gravity_softening_tracer", frozen.config.numerics.gravity_softening_tracer_kpc_comoving);
  stream << "gravity_solver = " << gravitySolverToString(frozen.config.numerics.gravity_solver) << '\n';
  stream << "hydro_solver = " << hydroSolverToString(frozen.config.numerics.hydro_solver) << '\n';
  stream << "treepm_pm_grid_nx = " << frozen.config.numerics.treepm_pm_grid_nx << '\n';
  stream << "treepm_pm_grid_ny = " << frozen.config.numerics.treepm_pm_grid_ny << '\n';
  stream << "treepm_pm_grid_nz = " << frozen.config.numerics.treepm_pm_grid_nz << '\n';
  stream << "treepm_asmth_cells = " << frozen.config.numerics.treepm_asmth_cells << '\n';
  stream << "treepm_rcut_cells = " << frozen.config.numerics.treepm_rcut_cells << '\n';
  stream << "treepm_tree_opening_criterion = "
         << treePmOpeningCriterionToString(frozen.config.numerics.treepm_tree_opening_criterion)
         << '\n';
  stream << "treepm_tree_opening_theta = "
         << frozen.config.numerics.treepm_tree_opening_theta << '\n';
  stream << "treepm_tree_relative_force_tolerance = "
         << frozen.config.numerics.treepm_tree_relative_force_tolerance << '\n';
  stream << "treepm_tree_relative_force_acceleration_floor = "
         << frozen.config.numerics.treepm_tree_relative_force_acceleration_floor << '\n';
  stream << "treepm_assignment_scheme = "
         << (frozen.config.numerics.treepm_assignment_scheme == TreePmAssignmentScheme::kCic ? "cic" : "tsc")
         << '\n';
  stream << "treepm_enable_window_deconvolution = "
         << (frozen.config.numerics.treepm_enable_window_deconvolution ? "true" : "false") << '\n';
  stream << "treepm_update_cadence_steps = " << frozen.config.numerics.treepm_update_cadence_steps << '\n';
  stream << "treepm_pm_decomposition_mode = "
         << pmDecompositionModeToString(frozen.config.numerics.treepm_pm_decomposition_mode) << '\n';
  stream << "treepm_tree_exchange_batch_bytes = "
         << frozen.config.numerics.treepm_tree_exchange_batch_bytes << '\n';
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
  stream << "snapshot_interval_time_code = "
         << frozen.config.output.snapshot_interval_time_code << '\n';
  stream << "write_restarts = " << (frozen.config.output.write_restarts ? "true" : "false") << '\n';
  stream << "\n[parallel]\n";
  stream << "mpi_ranks_expected = " << frozen.config.parallel.mpi_ranks_expected << '\n';
  stream << "omp_threads = " << frozen.config.parallel.omp_threads << '\n';
  stream << "gpu_devices = " << frozen.config.parallel.gpu_devices << '\n';
  stream << "deterministic_reduction = "
         << (frozen.config.parallel.deterministic_reduction ? "true" : "false") << '\n';
  stream << "decomposition_particle_count_weight = "
         << frozen.config.parallel.decomposition_particle_count_weight << '\n';
  stream << "decomposition_gas_cell_weight = "
         << frozen.config.parallel.decomposition_gas_cell_weight << '\n';
  stream << "decomposition_tree_interaction_weight = "
         << frozen.config.parallel.decomposition_tree_interaction_weight << '\n';
  stream << "decomposition_pm_mesh_weight = "
         << frozen.config.parallel.decomposition_pm_mesh_weight << '\n';
  stream << "decomposition_amr_patch_weight = "
         << frozen.config.parallel.decomposition_amr_patch_weight << '\n';
  stream << "decomposition_active_fraction_weight = "
         << frozen.config.parallel.decomposition_active_fraction_weight << '\n';
  stream << "decomposition_memory_pressure_weight = "
         << frozen.config.parallel.decomposition_memory_pressure_weight << '\n';
  stream << "decomposition_gpu_occupancy_weight = "
         << frozen.config.parallel.decomposition_gpu_occupancy_weight << '\n';
  stream << "decomposition_generic_work_weight = "
         << frozen.config.parallel.decomposition_generic_work_weight << '\n';
  stream << "decomposition_runtime_rebalance_enabled = "
         << (frozen.config.parallel.decomposition_runtime_rebalance_enabled ? "true" : "false") << '\n';
  stream << "decomposition_debug_exact_ownership_audit = "
         << (frozen.config.parallel.decomposition_debug_exact_ownership_audit ? "true" : "false") << '\n';
  stream << "decomposition_rebalance_imbalance_trigger = "
         << frozen.config.parallel.decomposition_rebalance_imbalance_trigger << '\n';
  stream << "decomposition_rebalance_memory_trigger = "
         << frozen.config.parallel.decomposition_rebalance_memory_trigger << '\n';
  stream << "decomposition_rebalance_max_migrated_load_fraction = "
         << frozen.config.parallel.decomposition_rebalance_max_migrated_load_fraction << '\n';
  stream << "decomposition_measured_tree_pair_weight = "
         << frozen.config.parallel.decomposition_measured_tree_pair_weight << '\n';
  stream << "decomposition_measured_pm_cell_weight = "
         << frozen.config.parallel.decomposition_measured_pm_cell_weight << '\n';
  stream << "decomposition_measured_amr_cell_weight = "
         << frozen.config.parallel.decomposition_measured_amr_cell_weight << '\n';
  stream << "decomposition_measured_hydro_face_weight = "
         << frozen.config.parallel.decomposition_measured_hydro_face_weight << '\n';
  stream << "decomposition_measured_wall_ms_weight = "
         << frozen.config.parallel.decomposition_measured_wall_ms_weight << '\n';
  stream << "isolated_pm_root_workspace_limit_bytes = "
         << frozen.config.parallel.isolated_pm_root_workspace_limit_bytes << '\n';
  stream << "zoom_high_res_allgather_limit_bytes = "
         << frozen.config.parallel.zoom_high_res_allgather_limit_bytes << '\n';
  stream << "\n[analysis]\n";
  stream << "enable_diagnostics = " << (frozen.config.analysis.enable_diagnostics ? "true" : "false")
         << '\n';
  stream << "enable_halo_workflow = " << (frozen.config.analysis.enable_halo_workflow ? "true" : "false")
         << '\n';
  stream << "halo_on_the_fly = " << (frozen.config.analysis.halo_on_the_fly ? "true" : "false") << '\n';
  stream << "diagnostics_execution_policy = "
         << diagnosticsExecutionPolicyToString(frozen.config.analysis.diagnostics_execution_policy) << '\n';
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
    const std::string& raw_text,
    const std::string& source_name,
    const ParseOptions& options) {
  std::map<std::string, ParsedEntry> entries = parsed_entries;
  std::set<std::string> consumed;
  FrozenConfig frozen;
  frozen.raw_text = raw_text;
  frozen.user_config.source_name = source_name;
  frozen.provenance.source_name = source_name;
  validateRequiredSchema(entries);
  applyDeprecatedAliases(entries, consumed, frozen);

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
  const bool has_explicit_ic_convention = entries.contains("mode.ic_convention");
  frozen.config.mode.ic_convention = parseInitialConditionConvention(
      requireString(
          entries, consumed, "mode.ic_convention",
          defaultFor("mode.ic_convention")));
  frozen.config.mode.ic_manifest_file = requireString(
      entries, consumed, "mode.ic_manifest_file",
      defaultFor("mode.ic_manifest_file"));
  frozen.config.mode.ic_chunk_particle_count = parseNumber<std::uint64_t>(
      requireString(
          entries, consumed, "mode.ic_chunk_particle_count",
          defaultFor("mode.ic_chunk_particle_count")),
      "mode.ic_chunk_particle_count");
  frozen.config.mode.ic_staging_particle_count = parseNumber<std::uint64_t>(
      requireString(
          entries, consumed, "mode.ic_staging_particle_count",
          defaultFor("mode.ic_staging_particle_count")),
      "mode.ic_staging_particle_count");
  frozen.config.mode.ic_part_type2_policy =
      parseInitialConditionSpeciesPolicy(
          requireString(
              entries, consumed, "mode.ic_part_type2_policy",
              defaultFor("mode.ic_part_type2_policy")),
          "mode.ic_part_type2_policy");
  frozen.config.mode.ic_part_type3_policy =
      parseInitialConditionSpeciesPolicy(
          requireString(
              entries, consumed, "mode.ic_part_type3_policy",
              defaultFor("mode.ic_part_type3_policy")),
          "mode.ic_part_type3_policy");
  if (!has_explicit_ic_convention && frozen.config.mode.ic_file != "generated") {
    throw ConfigError(
        "mode.ic_convention is required for external initial conditions; the runtime will not guess unit, frame, scale-factor, velocity, or species conventions from mode.ic_file");
  }
  frozen.config.mode.zoom_high_res_region = parseBool(
      requireString(entries, consumed, "mode.zoom_high_res_region", defaultFor("mode.zoom_high_res_region")),
      "mode.zoom_high_res_region");
  frozen.config.mode.zoom_region_file =
      requireString(entries, consumed, "mode.zoom_region_file", defaultFor("mode.zoom_region_file"));
  frozen.config.mode.zoom_long_range_strategy = parseZoomLongRangeStrategy(requireString(
      entries,
      consumed,
      "mode.zoom_long_range_strategy",
      defaultFor("mode.zoom_long_range_strategy")));
  frozen.config.mode.zoom_region_center_x_mpc_comoving = parseFloating(
      requireString(entries, consumed, "mode.zoom_region_center_x", defaultFor("mode.zoom_region_center_x")),
      "mode.zoom_region_center_x");
  frozen.config.mode.zoom_region_center_y_mpc_comoving = parseFloating(
      requireString(entries, consumed, "mode.zoom_region_center_y", defaultFor("mode.zoom_region_center_y")),
      "mode.zoom_region_center_y");
  frozen.config.mode.zoom_region_center_z_mpc_comoving = parseFloating(
      requireString(entries, consumed, "mode.zoom_region_center_z", defaultFor("mode.zoom_region_center_z")),
      "mode.zoom_region_center_z");
  frozen.config.mode.zoom_region_radius_mpc_comoving = parseFloating(
      requireString(entries, consumed, "mode.zoom_region_radius", defaultFor("mode.zoom_region_radius")),
      "mode.zoom_region_radius");
  frozen.config.mode.zoom_focused_pm_grid_nx = parseNumber<int>(
      requireString(entries, consumed, "mode.zoom_focused_pm_grid_nx", defaultFor("mode.zoom_focused_pm_grid_nx")),
      "mode.zoom_focused_pm_grid_nx");
  frozen.config.mode.zoom_focused_pm_grid_ny = parseNumber<int>(
      requireString(entries, consumed, "mode.zoom_focused_pm_grid_ny", defaultFor("mode.zoom_focused_pm_grid_ny")),
      "mode.zoom_focused_pm_grid_ny");
  frozen.config.mode.zoom_focused_pm_grid_nz = parseNumber<int>(
      requireString(entries, consumed, "mode.zoom_focused_pm_grid_nz", defaultFor("mode.zoom_focused_pm_grid_nz")),
      "mode.zoom_focused_pm_grid_nz");
  frozen.config.mode.zoom_contamination_radius_mpc_comoving = parseFloating(
      requireString(
          entries,
          consumed,
          "mode.zoom_contamination_radius",
          defaultFor("mode.zoom_contamination_radius")),
      "mode.zoom_contamination_radius");
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
  const bool has_box_size_x = entries.contains("cosmology.box_size_x");
  const bool has_box_size_y = entries.contains("cosmology.box_size_y");
  const bool has_box_size_z = entries.contains("cosmology.box_size_z");
  const bool has_box_size_scalar = entries.contains("cosmology.box_size");
  if ((has_box_size_x || has_box_size_y || has_box_size_z) &&
      !(has_box_size_x && has_box_size_y && has_box_size_z)) {
    throw ConfigError(
        "cosmology.box_size_x, cosmology.box_size_y, and cosmology.box_size_z must be specified together");
  }
  if (has_box_size_x) {
    frozen.config.cosmology.box_size_x_mpc_comoving = parseLengthMpc(
        requireString(entries, consumed, "cosmology.box_size_x", defaultFor("cosmology.box_size_x")),
        frozen.config.units.length_unit,
        "cosmology.box_size_x");
    frozen.config.cosmology.box_size_y_mpc_comoving = parseLengthMpc(
        requireString(entries, consumed, "cosmology.box_size_y", defaultFor("cosmology.box_size_y")),
        frozen.config.units.length_unit,
        "cosmology.box_size_y");
    frozen.config.cosmology.box_size_z_mpc_comoving = parseLengthMpc(
        requireString(entries, consumed, "cosmology.box_size_z", defaultFor("cosmology.box_size_z")),
        frozen.config.units.length_unit,
        "cosmology.box_size_z");
    if (has_box_size_scalar) {
      consumed.insert("cosmology.box_size");
    }
  } else {
    const double box_size_scalar_mpc = parseLengthMpc(
        requireString(entries, consumed, "cosmology.box_size", defaultFor("cosmology.box_size")),
        frozen.config.units.length_unit,
        "cosmology.box_size");
    frozen.config.cosmology.box_size_x_mpc_comoving = box_size_scalar_mpc;
    frozen.config.cosmology.box_size_y_mpc_comoving = box_size_scalar_mpc;
    frozen.config.cosmology.box_size_z_mpc_comoving = box_size_scalar_mpc;
  }
  frozen.config.cosmology.box_size_mpc_comoving = frozen.config.cosmology.box_size_x_mpc_comoving;

  normalizeScaleRedshiftPair(
      entries,
      consumed,
      "numerics.a_begin",
      "numerics.z_begin",
      parseFloating(defaultFor("numerics.a_begin"), "numerics.a_begin"),
      parseFloating(defaultFor("numerics.z_begin"), "numerics.z_begin"),
      frozen.config.numerics.a_begin,
      frozen.config.numerics.z_begin);
  normalizeScaleRedshiftPair(
      entries,
      consumed,
      "numerics.a_end",
      "numerics.z_end",
      parseFloating(defaultFor("numerics.a_end"), "numerics.a_end"),
      parseFloating(defaultFor("numerics.z_end"), "numerics.z_end"),
      frozen.config.numerics.a_end,
      frozen.config.numerics.z_end);
  frozen.config.numerics.t_code_begin = parseFloating(
      requireString(entries, consumed, "numerics.t_code_begin", defaultFor("numerics.t_code_begin")),
      "numerics.t_code_begin");
  frozen.config.numerics.t_code_end = parseFloating(
      requireString(entries, consumed, "numerics.t_code_end", defaultFor("numerics.t_code_end")),
      "numerics.t_code_end");
  frozen.config.numerics.t_phys_begin = parseFloating(
      requireString(entries, consumed, "numerics.t_phys_begin", defaultFor("numerics.t_phys_begin")),
      "numerics.t_phys_begin");
  frozen.config.numerics.t_phys_end = parseFloating(
      requireString(entries, consumed, "numerics.t_phys_end", defaultFor("numerics.t_phys_end")),
      "numerics.t_phys_end");
  frozen.config.numerics.integrator_time_variable = parseIntegratorTimeVariable(requireString(
      entries,
      consumed,
      "numerics.integrator_time_variable",
      defaultFor("numerics.integrator_time_variable")));
  frozen.config.numerics.cosmology_max_delta_ln_a = parseFloating(
      requireString(entries, consumed, "numerics.cosmology_max_delta_ln_a", defaultFor("numerics.cosmology_max_delta_ln_a")),
      "numerics.cosmology_max_delta_ln_a");
  frozen.config.numerics.cosmology_max_hubble_time_fraction = parseFloating(
      requireString(entries, consumed, "numerics.cosmology_max_hubble_time_fraction", defaultFor("numerics.cosmology_max_hubble_time_fraction")),
      "numerics.cosmology_max_hubble_time_fraction");
  frozen.config.numerics.source_max_fractional_change = parseFloating(
      requireString(entries, consumed, "numerics.source_max_fractional_change", defaultFor("numerics.source_max_fractional_change")),
      "numerics.source_max_fractional_change");
  frozen.config.numerics.max_global_steps = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "numerics.max_global_steps", "1024"),
      "numerics.max_global_steps"));
  frozen.config.numerics.hierarchical_max_rung = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "numerics.hierarchical_max_rung", "0"),
      "numerics.hierarchical_max_rung"));
  frozen.config.numerics.amr_max_level = static_cast<int>(parseNumber<long>(
      requireString(entries, consumed, "numerics.amr_max_level", "10"), "numerics.amr_max_level"));
  frozen.config.numerics.gravity_softening_kpc_comoving = parseLengthKpc(
      requireString(entries, consumed, "numerics.gravity_softening", "1.0 kpc"),
      frozen.config.units.length_unit,
      "numerics.gravity_softening");
  auto parseOptionalSofteningKpc = [&](const char* key) -> double {
    const std::string value = requireString(entries, consumed, key, defaultFor(key));
    if (trim(value).empty()) {
      return -1.0;
    }
    return parseLengthKpc(value, frozen.config.units.length_unit, key);
  };
  frozen.config.numerics.gravity_softening_gas_kpc_comoving = parseOptionalSofteningKpc("numerics.gravity_softening_gas");
  frozen.config.numerics.gravity_softening_dark_matter_kpc_comoving = parseOptionalSofteningKpc("numerics.gravity_softening_dark_matter");
  frozen.config.numerics.gravity_softening_star_kpc_comoving = parseOptionalSofteningKpc("numerics.gravity_softening_star");
  frozen.config.numerics.gravity_softening_black_hole_kpc_comoving = parseOptionalSofteningKpc("numerics.gravity_softening_black_hole");
  frozen.config.numerics.gravity_softening_tracer_kpc_comoving = parseOptionalSofteningKpc("numerics.gravity_softening_tracer");
  frozen.config.numerics.gravity_solver = parseGravitySolver(
      requireString(entries, consumed, "numerics.gravity_solver", defaultFor("numerics.gravity_solver")));
  frozen.config.numerics.hydro_solver = parseHydroSolver(
      requireString(entries, consumed, "numerics.hydro_solver", defaultFor("numerics.hydro_solver")));
  const bool has_pm_grid_nx = entries.contains("numerics.treepm_pm_grid_nx");
  const bool has_pm_grid_ny = entries.contains("numerics.treepm_pm_grid_ny");
  const bool has_pm_grid_nz = entries.contains("numerics.treepm_pm_grid_nz");
  const bool has_pm_grid_scalar = entries.contains("numerics.treepm_pm_grid");
  if ((has_pm_grid_nx || has_pm_grid_ny || has_pm_grid_nz) &&
      !(has_pm_grid_nx && has_pm_grid_ny && has_pm_grid_nz)) {
    throw ConfigError(
        "numerics.treepm_pm_grid_nx, numerics.treepm_pm_grid_ny, and numerics.treepm_pm_grid_nz must be specified together");
  }
  if (has_pm_grid_nx) {
    frozen.config.numerics.treepm_pm_grid_nx = static_cast<int>(parseNumber<long>(
        requireString(entries, consumed, "numerics.treepm_pm_grid_nx", defaultFor("numerics.treepm_pm_grid_nx")),
        "numerics.treepm_pm_grid_nx"));
    frozen.config.numerics.treepm_pm_grid_ny = static_cast<int>(parseNumber<long>(
        requireString(entries, consumed, "numerics.treepm_pm_grid_ny", defaultFor("numerics.treepm_pm_grid_ny")),
        "numerics.treepm_pm_grid_ny"));
    frozen.config.numerics.treepm_pm_grid_nz = static_cast<int>(parseNumber<long>(
        requireString(entries, consumed, "numerics.treepm_pm_grid_nz", defaultFor("numerics.treepm_pm_grid_nz")),
        "numerics.treepm_pm_grid_nz"));
    if (has_pm_grid_scalar) {
      consumed.insert("numerics.treepm_pm_grid");
    }
  } else {
    const int pm_grid_scalar = static_cast<int>(parseNumber<long>(
        requireString(entries, consumed, "numerics.treepm_pm_grid", defaultFor("numerics.treepm_pm_grid")),
        "numerics.treepm_pm_grid"));
    frozen.config.numerics.treepm_pm_grid_nx = pm_grid_scalar;
    frozen.config.numerics.treepm_pm_grid_ny = pm_grid_scalar;
    frozen.config.numerics.treepm_pm_grid_nz = pm_grid_scalar;
  }
  frozen.config.numerics.treepm_pm_grid = frozen.config.numerics.treepm_pm_grid_nx;
  frozen.config.numerics.treepm_asmth_cells = parseFloating(
      requireString(entries, consumed, "numerics.treepm_asmth_cells", defaultFor("numerics.treepm_asmth_cells")),
      "numerics.treepm_asmth_cells");
  frozen.config.numerics.treepm_rcut_cells = parseFloating(
      requireString(entries, consumed, "numerics.treepm_rcut_cells", defaultFor("numerics.treepm_rcut_cells")),
      "numerics.treepm_rcut_cells");
  frozen.config.numerics.treepm_tree_opening_criterion = parseTreePmOpeningCriterion(requireString(
      entries,
      consumed,
      "numerics.treepm_tree_opening_criterion",
      defaultFor("numerics.treepm_tree_opening_criterion")));
  frozen.config.numerics.treepm_tree_opening_theta = parseFloating(
      requireString(
          entries,
          consumed,
          "numerics.treepm_tree_opening_theta",
          defaultFor("numerics.treepm_tree_opening_theta")),
      "numerics.treepm_tree_opening_theta");
  frozen.config.numerics.treepm_tree_relative_force_tolerance = parseFloating(
      requireString(
          entries,
          consumed,
          "numerics.treepm_tree_relative_force_tolerance",
          defaultFor("numerics.treepm_tree_relative_force_tolerance")),
      "numerics.treepm_tree_relative_force_tolerance");
  frozen.config.numerics.treepm_tree_relative_force_acceleration_floor = parseFloating(
      requireString(
          entries,
          consumed,
          "numerics.treepm_tree_relative_force_acceleration_floor",
          defaultFor("numerics.treepm_tree_relative_force_acceleration_floor")),
      "numerics.treepm_tree_relative_force_acceleration_floor");
  frozen.config.numerics.treepm_assignment_scheme = parseTreePmAssignmentScheme(requireString(
      entries,
      consumed,
      "numerics.treepm_assignment_scheme",
      defaultFor("numerics.treepm_assignment_scheme")));
  frozen.config.numerics.treepm_enable_window_deconvolution = parseBool(
      requireString(
          entries,
          consumed,
          "numerics.treepm_enable_window_deconvolution",
          defaultFor("numerics.treepm_enable_window_deconvolution")),
      "numerics.treepm_enable_window_deconvolution");
  frozen.config.numerics.treepm_update_cadence_steps = static_cast<int>(parseNumber<long>(
      requireString(
          entries,
          consumed,
          "numerics.treepm_update_cadence_steps",
          defaultFor("numerics.treepm_update_cadence_steps")),
      "numerics.treepm_update_cadence_steps"));
  frozen.config.numerics.treepm_pm_decomposition_mode = parsePmDecompositionMode(requireString(
      entries,
      consumed,
      "numerics.treepm_pm_decomposition_mode",
      defaultFor("numerics.treepm_pm_decomposition_mode")));
  frozen.config.numerics.treepm_tree_exchange_batch_bytes = parseNumber<std::uint64_t>(
      requireString(
          entries,
          consumed,
          "numerics.treepm_tree_exchange_batch_bytes",
          defaultFor("numerics.treepm_tree_exchange_batch_bytes")),
      "numerics.treepm_tree_exchange_batch_bytes");

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
  frozen.config.output.snapshot_interval_time_code = parseFloating(
      requireString(entries, consumed, "output.snapshot_interval_time_code", "0.0"),
      "output.snapshot_interval_time_code");
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
  frozen.config.parallel.decomposition_particle_count_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_particle_count_weight",
                    defaultFor("parallel.decomposition_particle_count_weight")),
      "parallel.decomposition_particle_count_weight");
  frozen.config.parallel.decomposition_gas_cell_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_gas_cell_weight",
                    defaultFor("parallel.decomposition_gas_cell_weight")),
      "parallel.decomposition_gas_cell_weight");
  frozen.config.parallel.decomposition_tree_interaction_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_tree_interaction_weight",
                    defaultFor("parallel.decomposition_tree_interaction_weight")),
      "parallel.decomposition_tree_interaction_weight");
  frozen.config.parallel.decomposition_pm_mesh_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_pm_mesh_weight",
                    defaultFor("parallel.decomposition_pm_mesh_weight")),
      "parallel.decomposition_pm_mesh_weight");
  frozen.config.parallel.decomposition_amr_patch_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_amr_patch_weight",
                    defaultFor("parallel.decomposition_amr_patch_weight")),
      "parallel.decomposition_amr_patch_weight");
  frozen.config.parallel.decomposition_active_fraction_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_active_fraction_weight",
                    defaultFor("parallel.decomposition_active_fraction_weight")),
      "parallel.decomposition_active_fraction_weight");
  frozen.config.parallel.decomposition_memory_pressure_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_memory_pressure_weight",
                    defaultFor("parallel.decomposition_memory_pressure_weight")),
      "parallel.decomposition_memory_pressure_weight");
  frozen.config.parallel.decomposition_gpu_occupancy_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_gpu_occupancy_weight",
                    defaultFor("parallel.decomposition_gpu_occupancy_weight")),
      "parallel.decomposition_gpu_occupancy_weight");
  frozen.config.parallel.decomposition_generic_work_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_generic_work_weight",
                    defaultFor("parallel.decomposition_generic_work_weight")),
      "parallel.decomposition_generic_work_weight");
  frozen.config.parallel.decomposition_runtime_rebalance_enabled = parseBool(
      requireString(entries, consumed, "parallel.decomposition_runtime_rebalance_enabled",
                    defaultFor("parallel.decomposition_runtime_rebalance_enabled")),
      "parallel.decomposition_runtime_rebalance_enabled");
  frozen.config.parallel.decomposition_debug_exact_ownership_audit = parseBool(
      requireString(entries, consumed, "parallel.decomposition_debug_exact_ownership_audit",
                    defaultFor("parallel.decomposition_debug_exact_ownership_audit")),
      "parallel.decomposition_debug_exact_ownership_audit");
  frozen.config.parallel.decomposition_rebalance_imbalance_trigger = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_rebalance_imbalance_trigger",
                    defaultFor("parallel.decomposition_rebalance_imbalance_trigger")),
      "parallel.decomposition_rebalance_imbalance_trigger");
  frozen.config.parallel.decomposition_rebalance_memory_trigger = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_rebalance_memory_trigger",
                    defaultFor("parallel.decomposition_rebalance_memory_trigger")),
      "parallel.decomposition_rebalance_memory_trigger");
  frozen.config.parallel.decomposition_rebalance_max_migrated_load_fraction = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_rebalance_max_migrated_load_fraction",
                    defaultFor("parallel.decomposition_rebalance_max_migrated_load_fraction")),
      "parallel.decomposition_rebalance_max_migrated_load_fraction");
  frozen.config.parallel.decomposition_measured_tree_pair_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_measured_tree_pair_weight",
                    defaultFor("parallel.decomposition_measured_tree_pair_weight")),
      "parallel.decomposition_measured_tree_pair_weight");
  frozen.config.parallel.decomposition_measured_pm_cell_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_measured_pm_cell_weight",
                    defaultFor("parallel.decomposition_measured_pm_cell_weight")),
      "parallel.decomposition_measured_pm_cell_weight");
  frozen.config.parallel.decomposition_measured_amr_cell_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_measured_amr_cell_weight",
                    defaultFor("parallel.decomposition_measured_amr_cell_weight")),
      "parallel.decomposition_measured_amr_cell_weight");
  frozen.config.parallel.decomposition_measured_hydro_face_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_measured_hydro_face_weight",
                    defaultFor("parallel.decomposition_measured_hydro_face_weight")),
      "parallel.decomposition_measured_hydro_face_weight");
  frozen.config.parallel.decomposition_measured_wall_ms_weight = parseNumber<double>(
      requireString(entries, consumed, "parallel.decomposition_measured_wall_ms_weight",
                    defaultFor("parallel.decomposition_measured_wall_ms_weight")),
      "parallel.decomposition_measured_wall_ms_weight");
  frozen.config.parallel.isolated_pm_root_workspace_limit_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.isolated_pm_root_workspace_limit_bytes",
                    defaultFor("parallel.isolated_pm_root_workspace_limit_bytes")),
      "parallel.isolated_pm_root_workspace_limit_bytes");
  frozen.config.parallel.zoom_high_res_allgather_limit_bytes = parseNumber<std::uint64_t>(
      requireString(entries, consumed, "parallel.zoom_high_res_allgather_limit_bytes",
                    defaultFor("parallel.zoom_high_res_allgather_limit_bytes")),
      "parallel.zoom_high_res_allgather_limit_bytes");

  frozen.config.analysis.enable_diagnostics = parseBool(
      requireString(entries, consumed, "analysis.enable_diagnostics", "true"),
      "analysis.enable_diagnostics");
  frozen.config.analysis.enable_halo_workflow = parseBool(
      requireString(entries, consumed, "analysis.enable_halo_workflow", "false"),
      "analysis.enable_halo_workflow");
  frozen.config.analysis.halo_on_the_fly = parseBool(
      requireString(entries, consumed, "analysis.halo_on_the_fly", "false"),
      "analysis.halo_on_the_fly");
  frozen.config.analysis.diagnostics_execution_policy = parseDiagnosticsExecutionPolicy(requireString(
      entries,
      consumed,
      "analysis.diagnostics_execution_policy",
      diagnosticsExecutionPolicyToString(frozen.config.analysis.diagnostics_execution_policy)));
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



  for (const auto& [key, entry] : entries) {
    UserConfigEntry user_entry;
    user_entry.canonical_key = key;
    user_entry.value = trim(entry.value);
    user_entry.source_key = key;
    user_entry.source_line = entry.line_number;
    frozen.user_config.entries.push_back(std::move(user_entry));
  }

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

  frozen.normalized_text = buildNormalizedText(frozen);
  frozen.provenance.config_hash = stableConfigHash(frozen.normalized_text);
  frozen.provenance.config_hash_hex = stableConfigHashHex(frozen.normalized_text);

  return frozen;
}

}  // namespace

ConfigError::ConfigError(const std::string& message) : std::runtime_error(message) {}

SimulationConfig makeUnvalidatedSimulationConfigForTests() {
  return SimulationConfig(SimulationConfig::ConstructionToken{});
}

FrozenConfig::FrozenConfig()
    : config(SimulationConfig::ConstructionToken{}) {}

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
  return normalizeValidateFreeze(entries, config_text, source_name, options);
}

DerivedRuntimeConfig deriveRuntimeConfig(const FrozenConfig& frozen_config) {
  const SimulationConfig& config = frozen_config.config;

  if (!(config.numerics.t_code_end > config.numerics.t_code_begin)) {
    throw ConfigError(
        "derived runtime config requires numerics.t_code_end > numerics.t_code_begin");
  }
  if (config.cosmology.box_size_x_mpc_comoving <= 0.0 ||
      config.cosmology.box_size_y_mpc_comoving <= 0.0 ||
      config.cosmology.box_size_z_mpc_comoving <= 0.0) {
    throw ConfigError("derived runtime config requires positive comoving box-size axes");
  }
  if (config.numerics.treepm_pm_grid_nx <= 0 ||
      config.numerics.treepm_pm_grid_ny <= 0 ||
      config.numerics.treepm_pm_grid_nz <= 0) {
    throw ConfigError("derived runtime config requires positive TreePM PM grid shape");
  }
  if (!(config.numerics.gravity_softening_kpc_comoving > 0.0)) {
    throw ConfigError("derived runtime config requires positive numerics.gravity_softening");
  }

  const auto choose_softening = [fallback = config.numerics.gravity_softening_kpc_comoving](
                                    double species_value) {
    return species_value > 0.0 ? species_value : fallback;
  };

  DerivedRuntimeConfig derived;
  derived.t_code_begin = config.numerics.t_code_begin;
  derived.t_code_end = config.numerics.t_code_end;
  derived.box_size_mpc_comoving = {
      config.cosmology.box_size_x_mpc_comoving,
      config.cosmology.box_size_y_mpc_comoving,
      config.cosmology.box_size_z_mpc_comoving};
  derived.treepm_pm_grid_shape = {
      config.numerics.treepm_pm_grid_nx,
      config.numerics.treepm_pm_grid_ny,
      config.numerics.treepm_pm_grid_nz};
  derived.gravity_softening_kpc_comoving_by_species = {
      choose_softening(config.numerics.gravity_softening_dark_matter_kpc_comoving),
      choose_softening(config.numerics.gravity_softening_gas_kpc_comoving),
      choose_softening(config.numerics.gravity_softening_star_kpc_comoving),
      choose_softening(config.numerics.gravity_softening_black_hole_kpc_comoving),
      choose_softening(config.numerics.gravity_softening_tracer_kpc_comoving)};
  derived.normalized_config_hash = frozen_config.provenance.config_hash;
  derived.normalized_config_hash_hex = frozen_config.provenance.config_hash_hex;
  return derived;
}

std::string serializeDerivedRuntimeConfig(const DerivedRuntimeConfig& derived_config) {
  std::ostringstream stream;
  stream << std::setprecision(std::numeric_limits<double>::max_digits10);
  stream << "t_code_begin=" << derived_config.t_code_begin << '\n';
  stream << "t_code_end=" << derived_config.t_code_end << '\n';
  stream << "box_size_x_mpc_comoving=" << derived_config.box_size_mpc_comoving[0] << '\n';
  stream << "box_size_y_mpc_comoving=" << derived_config.box_size_mpc_comoving[1] << '\n';
  stream << "box_size_z_mpc_comoving=" << derived_config.box_size_mpc_comoving[2] << '\n';
  stream << "treepm_pm_grid_nx=" << derived_config.treepm_pm_grid_shape[0] << '\n';
  stream << "treepm_pm_grid_ny=" << derived_config.treepm_pm_grid_shape[1] << '\n';
  stream << "treepm_pm_grid_nz=" << derived_config.treepm_pm_grid_shape[2] << '\n';
  stream << "gravity_softening_dark_matter_kpc_comoving="
         << derived_config.gravity_softening_kpc_comoving_by_species[0] << '\n';
  stream << "gravity_softening_gas_kpc_comoving="
         << derived_config.gravity_softening_kpc_comoving_by_species[1] << '\n';
  stream << "gravity_softening_star_kpc_comoving="
         << derived_config.gravity_softening_kpc_comoving_by_species[2] << '\n';
  stream << "gravity_softening_black_hole_kpc_comoving="
         << derived_config.gravity_softening_kpc_comoving_by_species[3] << '\n';
  stream << "gravity_softening_tracer_kpc_comoving="
         << derived_config.gravity_softening_kpc_comoving_by_species[4] << '\n';
  stream << "normalized_config_hash_hex=" << derived_config.normalized_config_hash_hex << '\n';
  return stream.str();
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

std::string initialConditionConventionToString(
    InitialConditionConvention convention) {
  switch (convention) {
    case InitialConditionConvention::kGenerated:
      return "generated";
    case InitialConditionConvention::kChuiCanonicalV1:
      return "chui_canonical_v1";
    case InitialConditionConvention::kGadgetArepoBridgeV1:
      return "gadget_arepo_bridge_v1";
    case InitialConditionConvention::kManifestV1:
      return "manifest_v1";
  }
  throw ConfigError(
      "unhandled InitialConditionConvention enum value during serialization");
}

std::string initialConditionSpeciesPolicyToString(
    InitialConditionSpeciesPolicy policy) {
  switch (policy) {
    case InitialConditionSpeciesPolicy::kReject:
      return "reject";
    case InitialConditionSpeciesPolicy::kDarkMatter:
      return "dark_matter";
    case InitialConditionSpeciesPolicy::kStar:
      return "star";
    case InitialConditionSpeciesPolicy::kBlackHole:
      return "black_hole";
    case InitialConditionSpeciesPolicy::kTracer:
      return "tracer";
  }
  throw ConfigError(
      "unhandled InitialConditionSpeciesPolicy enum value during serialization");
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

std::string integratorTimeVariableToString(IntegratorTimeVariable variable) {
  switch (variable) {
    case IntegratorTimeVariable::kScaleFactor:
      return "scale_factor";
    case IntegratorTimeVariable::kLogScaleFactor:
      return "ln_a";
    case IntegratorTimeVariable::kCodeTime:
      return "t_code";
    case IntegratorTimeVariable::kPhysicalTime:
      return "t_phys";
  }
  throw ConfigError("unhandled IntegratorTimeVariable enum value during serialization");
}

}  // namespace cosmosim::core
