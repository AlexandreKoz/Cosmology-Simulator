#include "cosmosim/io/ic_reader.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace cosmosim::io {
namespace {

[[nodiscard]] std::string escapeJson(std::string_view value) {
  std::string result;
  result.reserve(value.size());
  for (const char ch : value) {
    if (ch == '\\' || ch == '"') {
      result.push_back('\\');
    }
    result.push_back(ch);
  }
  return result;
}

[[nodiscard]] std::string_view dialectName(IcDialect dialect) {
  switch (dialect) {
    case IcDialect::kChuiCanonicalV1:
      return "chui_canonical_v1";
    case IcDialect::kGadgetArepoBridgeV1:
      return "gadget_arepo_bridge_v1";
  }
  throw std::invalid_argument("unknown IC dialect");
}

[[nodiscard]] std::string_view frameName(IcCoordinateFrame frame) {
  return frame == IcCoordinateFrame::kComoving ? "comoving" : "physical";
}

[[nodiscard]] std::string_view velocityConventionName(
    IcVelocityConvention convention) {
  switch (convention) {
    case IcVelocityConvention::kNotVelocity:
      return "not_velocity";
    case IcVelocityConvention::kPhysicalPeculiar:
      return "physical_peculiar";
    case IcVelocityConvention::kSqrtAScaledPeculiar:
      return "sqrt_a_scaled_peculiar";
    case IcVelocityConvention::kComovingCoordinateRate:
      return "comoving_coordinate_rate";
  }
  return "unknown";
}

[[nodiscard]] std::string_view semanticsName(IcFieldSemantics semantics) {
  switch (semantics) {
    case IcFieldSemantics::kCoordinate: return "coordinate";
    case IcFieldSemantics::kVelocity: return "velocity";
    case IcFieldSemantics::kExtensive: return "extensive";
    case IcFieldSemantics::kIntensive: return "intensive";
    case IcFieldSemantics::kSpecific: return "specific";
    case IcFieldSemantics::kIdentifier: return "identifier";
  }
  return "unknown";
}

[[nodiscard]] std::string_view dispositionName(
    IcFieldDisposition disposition) {
  switch (disposition) {
    case IcFieldDisposition::kConverted: return "converted";
    case IcFieldDisposition::kPreserved: return "preserved";
    case IcFieldDisposition::kDefaulted: return "defaulted";
    case IcFieldDisposition::kDropped: return "dropped";
    case IcFieldDisposition::kRejected: return "rejected";
  }
  return "unknown";
}

[[nodiscard]] std::string_view speciesPolicyName(IcSpeciesPolicy policy) {
  switch (policy) {
    case IcSpeciesPolicy::kReject: return "reject";
    case IcSpeciesPolicy::kGas: return "gas";
    case IcSpeciesPolicy::kDarkMatter: return "dark_matter";
    case IcSpeciesPolicy::kCollisionlessFamily2AsDarkMatter:
      return "family2_as_dark_matter";
    case IcSpeciesPolicy::kCollisionlessFamily3AsDarkMatter:
      return "family3_as_dark_matter";
    case IcSpeciesPolicy::kStar: return "star";
    case IcSpeciesPolicy::kBlackHole: return "black_hole";
    case IcSpeciesPolicy::kTracer: return "tracer";
  }
  return "unknown";
}

}  // namespace

void validateIcManifest(const IcManifest& manifest) {
  if (manifest.schema_version != 1U) {
    throw std::invalid_argument("IcManifest schema_version must be 1");
  }
  if (manifest.dialect_version.empty()) {
    throw std::invalid_argument("IcManifest dialect_version is required");
  }
  if (manifest.source_files.empty()) {
    throw std::invalid_argument("IcManifest requires at least one source file");
  }
  if (manifest.num_files_per_snapshot != manifest.source_files.size() ||
      manifest.num_part_this_file.size() != manifest.source_files.size()) {
    throw std::invalid_argument(
        "IcManifest file/count arrays must match NumFilesPerSnapshot");
  }
  if (!manifest.source_provenance_ids.empty() &&
      manifest.source_provenance_ids.size() != manifest.source_files.size()) {
    throw std::invalid_argument(
        "IcManifest provenance IDs must be empty or one-per-file");
  }
  std::unordered_set<std::string> source_paths;
  std::string previous_path;
  for (const auto& path : manifest.source_files) {
    const std::string normalized = path.lexically_normal().generic_string();
    if (normalized.empty() || !source_paths.insert(normalized).second) {
      throw std::invalid_argument(
          "IcManifest source files must be non-empty and unique");
    }
    if (!previous_path.empty() && normalized < previous_path) {
      throw std::invalid_argument(
          "IcManifest source files must use deterministic lexical ordering");
    }
    previous_path = normalized;
  }
  if (!(manifest.scale_factor > 0.0) ||
      !std::isfinite(manifest.scale_factor) ||
      !std::isfinite(manifest.redshift) ||
      !(manifest.hubble_param > 0.0) ||
      !std::isfinite(manifest.hubble_param) ||
      !std::isfinite(manifest.box_size) || manifest.box_size <= 0.0 ||
      !std::isfinite(manifest.omega_matter) ||
      !std::isfinite(manifest.omega_lambda)) {
    throw std::invalid_argument(
        "IcManifest cosmology, BoxSize, h, and start epoch must be finite and physical");
  }
  const double expected_redshift = 1.0 / manifest.scale_factor - 1.0;
  if (std::abs(expected_redshift - manifest.redshift) >
      1.0e-10 * std::max(1.0, std::abs(expected_redshift))) {
    throw std::invalid_argument(
        "IcManifest Time/Redshift values are inconsistent");
  }
  for (std::size_t type = 0; type < manifest.num_part_total.size(); ++type) {
    if ((manifest.num_part_total[type] >> 32U) !=
        manifest.num_part_total_high_word[type]) {
      throw std::invalid_argument(
          "IcManifest NumPart_Total_HighWord disagrees with canonical 64-bit total");
    }
    std::uint64_t summed_count = 0U;
    for (const auto& per_file : manifest.num_part_this_file) {
      if (summed_count >
          std::numeric_limits<std::uint64_t>::max() - per_file[type]) {
        throw std::overflow_error("IcManifest particle count overflow");
      }
      summed_count += per_file[type];
    }
    if (summed_count != manifest.num_part_total[type]) {
      throw std::invalid_argument(
          "IcManifest per-file counts do not sum to NumPart_Total");
    }
    if (manifest.num_part_total[type] > 0U &&
        manifest.species_policy[type] == IcSpeciesPolicy::kReject) {
      throw std::invalid_argument(
          "IcManifest contains a populated particle type whose species policy is reject");
    }
    if (!std::isfinite(manifest.mass_table[type]) ||
        manifest.mass_table[type] < 0.0) {
      throw std::invalid_argument(
          "IcManifest MassTable values must be finite and non-negative");
    }
  }
  std::unordered_set<std::string> field_paths;
  for (const IcFieldManifest& field : manifest.fields) {
    if (field.dataset_path.empty() || field.scalar_type.empty() ||
        !field_paths.insert(field.dataset_path).second) {
      throw std::invalid_argument(
          "IcManifest field paths/types must be non-empty and unique");
    }
    if (field.rank != field.dimensions.size() || field.rank == 0U ||
        field.dimensions.empty() || field.dimensions.front() != field.record_count) {
      throw std::invalid_argument(
          "IcManifest field rank/dimensions/count are inconsistent");
    }
    if (!(field.base_unit_to_si > 0.0) ||
        !std::isfinite(field.base_unit_to_si) ||
        !std::isfinite(field.hubble_exponent) ||
        !std::isfinite(field.scale_factor_exponent)) {
      throw std::invalid_argument(
          "IcManifest field conversion factors must be finite and explicit");
    }
    if (field.semantics == IcFieldSemantics::kVelocity &&
        field.velocity_convention == IcVelocityConvention::kNotVelocity) {
      throw std::invalid_argument(
          "IcManifest velocity field requires a velocity convention");
    }
  }
}

double icStoredToSiMultiplier(
    const IcFieldManifest& field,
    double hubble_param,
    double scale_factor) {
  if (!(field.base_unit_to_si > 0.0) || !(hubble_param > 0.0) ||
      !(scale_factor > 0.0)) {
    throw std::invalid_argument(
        "IC unit conversion requires positive base unit, h, and scale factor");
  }
  return field.base_unit_to_si *
      std::pow(hubble_param, field.hubble_exponent) *
      std::pow(scale_factor, field.scale_factor_exponent);
}

double icVelocityConventionMultiplier(
    IcVelocityConvention convention,
    double scale_factor) {
  if (!(scale_factor > 0.0)) {
    throw std::invalid_argument(
        "IC velocity conversion requires scale_factor > 0");
  }
  switch (convention) {
    case IcVelocityConvention::kNotVelocity:
    case IcVelocityConvention::kPhysicalPeculiar:
      return 1.0;
    case IcVelocityConvention::kSqrtAScaledPeculiar:
      return 1.0 / std::sqrt(scale_factor);
    case IcVelocityConvention::kComovingCoordinateRate:
      return scale_factor;
  }
  throw std::invalid_argument("unknown IC velocity convention");
}

std::string serializeIcManifestJson(const IcManifest& manifest) {
  validateIcManifest(manifest);
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  out << "{\n  \"schema_version\": " << manifest.schema_version
      << ",\n  \"dialect\": \"" << dialectName(manifest.dialect)
      << "\",\n  \"dialect_version\": \""
      << escapeJson(manifest.dialect_version) << "\",\n";
  out << "  \"num_files_per_snapshot\": "
      << manifest.num_files_per_snapshot << ",\n  \"source_files\": [";
  for (std::size_t i = 0; i < manifest.source_files.size(); ++i) {
    out << (i == 0U ? "" : ", ") << '"'
        << escapeJson(manifest.source_files[i].generic_string()) << '"';
  }
  out << "],\n  \"source_provenance_ids\": [";
  for (std::size_t i = 0; i < manifest.source_provenance_ids.size(); ++i) {
    out << (i == 0U ? "" : ", ") << '"'
        << escapeJson(manifest.source_provenance_ids[i]) << '"';
  }
  const auto write_u64_array = [&](std::string_view name, const auto& values) {
    out << "],\n  \"" << name << "\": [";
    for (std::size_t i = 0; i < values.size(); ++i) {
      out << (i == 0U ? "" : ", ") << values[i];
    }
  };
  write_u64_array("num_part_total", manifest.num_part_total);
  write_u64_array("num_part_total_high_word", manifest.num_part_total_high_word);
  write_u64_array("mass_table", manifest.mass_table);
  out << "],\n  \"num_part_this_file\": [";
  for (std::size_t file = 0; file < manifest.num_part_this_file.size(); ++file) {
    out << (file == 0U ? "" : ", ") << '[';
    for (std::size_t type = 0; type < 6U; ++type) {
      out << (type == 0U ? "" : ", ")
          << manifest.num_part_this_file[file][type];
    }
    out << ']';
  }
  out << "],\n  \"box_size\": " << manifest.box_size
      << ",\n  \"time\": " << manifest.scale_factor
      << ",\n  \"redshift\": " << manifest.redshift
      << ",\n  \"omega_matter\": " << manifest.omega_matter
      << ",\n  \"omega_lambda\": " << manifest.omega_lambda
      << ",\n  \"hubble_param\": " << manifest.hubble_param
      << ",\n  \"species_policy\": [";
  for (std::size_t type = 0; type < manifest.species_policy.size(); ++type) {
    out << (type == 0U ? "" : ", ") << '"'
        << speciesPolicyName(manifest.species_policy[type]) << '"';
  }
  out << "],\n  \"fields\": [";
  for (std::size_t i = 0; i < manifest.fields.size(); ++i) {
    const IcFieldManifest& field = manifest.fields[i];
    out << (i == 0U ? "\n" : ",\n")
        << "    {\"dataset_path\": \"" << escapeJson(field.dataset_path)
        << "\", \"scalar_type\": \"" << escapeJson(field.scalar_type)
        << "\", \"rank\": " << static_cast<unsigned int>(field.rank)
        << ", \"dimensions\": [";
    for (std::size_t d = 0; d < field.dimensions.size(); ++d) {
      out << (d == 0U ? "" : ", ") << field.dimensions[d];
    }
    out << "], \"record_count\": " << field.record_count
        << ", \"base_unit_to_si\": " << field.base_unit_to_si
        << ", \"hubble_exponent\": " << field.hubble_exponent
        << ", \"scale_factor_exponent\": " << field.scale_factor_exponent
        << ", \"coordinate_frame\": \"" << frameName(field.coordinate_frame)
        << "\", \"velocity_convention\": \""
        << velocityConventionName(field.velocity_convention)
        << "\", \"semantics\": \"" << semanticsName(field.semantics)
        << "\", \"disposition\": \""
        << dispositionName(field.disposition) << "\"}";
  }
  out << (manifest.fields.empty() ? "" : "\n") << "  ]\n}\n";
  return out.str();
}

void writeIcManifestJson(
    const IcManifest& manifest,
    const std::filesystem::path& output_path) {
  std::ofstream output(output_path, std::ios::binary | std::ios::trunc);
  if (!output) {
    throw std::runtime_error(
        "failed to open IC manifest output: " + output_path.string());
  }
  output << serializeIcManifestJson(manifest);
  if (!output) {
    throw std::runtime_error(
        "failed to write IC manifest output: " + output_path.string());
  }
}

}  // namespace cosmosim::io
