#include "cosmosim/core/simulation_state.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>

#include "core/internal/text_codec.hpp"

namespace cosmosim::core {
namespace {
[[nodiscard]] std::string trim(std::string value) {
  const auto first = value.find_first_not_of(" \t\n\r");
  if (first == std::string::npos) {
    return {};
  }
  const auto last = value.find_last_not_of(" \t\n\r");
  return value.substr(first, last - first + 1);
}

[[nodiscard]] std::uint64_t parseUint64(std::string_view value, const std::string& key) {
  std::uint64_t parsed = 0;
  const auto* begin = value.data();
  const auto* end = value.data() + value.size();
  const auto [ptr, ec] = std::from_chars(begin, end, parsed);
  if (ec != std::errc{} || ptr != end) {
    throw std::invalid_argument("StateMetadata.deserialize: invalid integer for key '" + key + "'");
  }
  return parsed;
}

[[nodiscard]] std::uint32_t parseUint32(std::string_view value, const std::string& key) {
  const auto parsed = parseUint64(value, key);
  if (parsed > std::numeric_limits<std::uint32_t>::max()) {
    throw std::invalid_argument("StateMetadata.deserialize: integer overflow for key '" + key + "'");
  }
  return static_cast<std::uint32_t>(parsed);
}

[[nodiscard]] double parseDouble(std::string_view value, const std::string& key) {
  std::string temp(value);
  char* end = nullptr;
  const double parsed = std::strtod(temp.c_str(), &end);
  if (end != temp.c_str() + static_cast<std::ptrdiff_t>(temp.size()) || !std::isfinite(parsed)) {
    throw std::invalid_argument(
        "StateMetadata.deserialize: invalid floating-point value for key '" + key + "'");
  }
  return parsed;
}
}  // namespace

std::string StateMetadata::serialize() const {
  if (schema_version != 2U && schema_version != 3U) {
    throw std::invalid_argument("StateMetadata.serialize: unsupported schema_version");
  }
  const auto encode_string = [this](std::string_view value) {
    if (schema_version >= 3U) {
      return internal::escapeTextLine(value);
    }
    if (value.find_first_of("\n\r") != std::string_view::npos) {
      throw std::invalid_argument(
          "StateMetadata.serialize: schema v2 cannot represent embedded line breaks; use schema v3");
    }
    return std::string(value);
  };

  std::ostringstream out;
  out << "schema_version=" << schema_version << '\n';
  out << "run_name=" << encode_string(run_name) << '\n';
  out << "normalized_config_hash=" << normalized_config_hash << '\n';
  out << "normalized_config_hash_hex=" << encode_string(normalized_config_hash_hex) << '\n';
  out << "step_index=" << step_index << '\n';
  out << "scale_factor=" << scale_factor << '\n';
  out << "snapshot_stem=" << encode_string(snapshot_stem) << '\n';
  out << "restart_stem=" << encode_string(restart_stem) << '\n';
  return out.str();
}

StateMetadata StateMetadata::deserialize(std::string_view text) {
  std::unordered_set<std::string> seen_keys;
  std::unordered_map<std::string, std::string> raw_values;
  std::istringstream in{std::string(text)};
  std::string line;

  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    const auto equal_pos = line.find('=');
    if (equal_pos == std::string::npos || equal_pos == 0U) {
      throw std::invalid_argument("StateMetadata.deserialize: malformed line '" + line + "'");
    }
    const std::string key = trim(line.substr(0, equal_pos));
    if (key.empty()) {
      throw std::invalid_argument("StateMetadata.deserialize: malformed line '" + line + "'");
    }
    if (!seen_keys.insert(key).second) {
      throw std::invalid_argument("StateMetadata.deserialize: duplicate key '" + key + "'");
    }
    raw_values.emplace(key, line.substr(equal_pos + 1U));
  }

  static const std::array<std::string_view, 8> k_required_keys = {
      "schema_version", "run_name", "normalized_config_hash", "normalized_config_hash_hex",
      "step_index", "scale_factor", "snapshot_stem", "restart_stem"};
  for (const std::string_view key : k_required_keys) {
    if (!seen_keys.contains(std::string(key))) {
      throw std::invalid_argument("StateMetadata.deserialize: missing required key '" + std::string(key) + "'");
    }
  }
  if (seen_keys.size() != k_required_keys.size()) {
    for (const auto& [key, value] : raw_values) {
      (void)value;
      if (std::find(k_required_keys.begin(), k_required_keys.end(), key) == k_required_keys.end()) {
        throw std::invalid_argument("StateMetadata.deserialize: unknown key '" + key + "'");
      }
    }
  }

  StateMetadata metadata;
  metadata.schema_version = parseUint32(trim(raw_values.at("schema_version")), "schema_version");
  if (metadata.schema_version != 2U && metadata.schema_version != 3U) {
    throw std::invalid_argument("StateMetadata.deserialize: unsupported schema_version");
  }
  const auto decode_string = [&metadata](const std::string& raw, const char* key) {
    if (metadata.schema_version >= 3U) {
      return internal::unescapeTextLineStrict(raw, std::string("StateMetadata.deserialize key '") + key + "'");
    }
    return trim(raw);
  };

  metadata.run_name = decode_string(raw_values.at("run_name"), "run_name");
  metadata.normalized_config_hash = parseUint64(
      trim(raw_values.at("normalized_config_hash")), "normalized_config_hash");
  metadata.normalized_config_hash_hex = decode_string(
      raw_values.at("normalized_config_hash_hex"), "normalized_config_hash_hex");
  metadata.step_index = parseUint64(trim(raw_values.at("step_index")), "step_index");
  metadata.scale_factor = parseDouble(trim(raw_values.at("scale_factor")), "scale_factor");
  metadata.snapshot_stem = decode_string(raw_values.at("snapshot_stem"), "snapshot_stem");
  metadata.restart_stem = decode_string(raw_values.at("restart_stem"), "restart_stem");

  if (!(metadata.scale_factor > 0.0)) {
    throw std::invalid_argument("StateMetadata.deserialize: scale_factor must be finite and positive");
  }
  return metadata;
}

bool ModuleSidecarBlock::isParticleIndexed() const noexcept {
  return particle_indexed || row_stride_bytes != 0U || !particle_id_by_row.empty();
}

bool ModuleSidecarBlock::requiresSpecies(ParticleSpecies species) const noexcept {
  const auto bit = static_cast<std::uint32_t>(species);
  if (bit >= 32U) {
    return false;
  }
  const std::uint32_t mask = required_species_mask | requirement.species_mask;
  return (mask & (1U << bit)) != 0U;
}

std::size_t ModuleSidecarBlock::rowCount() const noexcept {
  return particle_id_by_row.size();
}

std::span<const std::byte> ModuleSidecarBlock::rowPayload(std::size_t row) const {
  if (!isParticleIndexed()) {
    throw std::invalid_argument("ModuleSidecarBlock.rowPayload: block is not particle-indexed");
  }
  if (row >= particle_id_by_row.size()) {
    throw std::out_of_range("ModuleSidecarBlock.rowPayload: row out of range");
  }
  const std::size_t stride = static_cast<std::size_t>(row_stride_bytes);
  if (stride == 0U || row > std::numeric_limits<std::size_t>::max() / stride) {
    throw std::invalid_argument("ModuleSidecarBlock.rowPayload: particle-indexed row offset overflow");
  }
  const std::size_t begin = row * stride;
  if (begin > payload.size() || stride > payload.size() - begin) {
    throw std::invalid_argument("ModuleSidecarBlock.rowPayload: invalid particle-indexed row layout");
  }
  return std::span<const std::byte>(payload.data() + begin, stride);
}

namespace {

[[nodiscard]] std::uint32_t validModuleRequirementSpeciesMask() noexcept {
  return (1U << static_cast<std::uint32_t>(ParticleSpecies::kDarkMatter)) |
      (1U << static_cast<std::uint32_t>(ParticleSpecies::kGas)) |
      (1U << static_cast<std::uint32_t>(ParticleSpecies::kStar)) |
      (1U << static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) |
      (1U << static_cast<std::uint32_t>(ParticleSpecies::kTracer));
}

void validateModuleSidecarRequirement(const ModuleSidecarRequirement& requirement) {
  const std::uint32_t valid_species_mask = validModuleRequirementSpeciesMask();
  if ((requirement.species_mask & ~valid_species_mask) != 0U) {
    throw std::invalid_argument("ModuleSidecarRegistry.upsert: module sidecar requirement has invalid species bits");
  }
  switch (requirement.kind) {
    case ModuleSidecarRequirementKind::kSparse:
      if (requirement.species_mask != 0U || requirement.particle_flags_mask != 0U || requirement.threshold_code != 0.0) {
        throw std::invalid_argument("ModuleSidecarRegistry.upsert: sparse sidecar requirement must not carry predicate parameters");
      }
      break;
    case ModuleSidecarRequirementKind::kSpeciesMask:
      if (requirement.species_mask == 0U) {
        throw std::invalid_argument("ModuleSidecarRegistry.upsert: species-mask sidecar requirement needs a nonzero species mask");
      }
      if (requirement.particle_flags_mask != 0U || requirement.threshold_code != 0.0) {
        throw std::invalid_argument("ModuleSidecarRegistry.upsert: species-mask sidecar requirement cannot carry flag or threshold predicates");
      }
      break;
    case ModuleSidecarRequirementKind::kGasDensityAtLeast:
      if (requirement.species_mask != 0U || requirement.particle_flags_mask != 0U ||
          !std::isfinite(requirement.threshold_code) || requirement.threshold_code < 0.0) {
        throw std::invalid_argument("ModuleSidecarRegistry.upsert: gas-density sidecar requirement has invalid predicate parameters");
      }
      break;
    case ModuleSidecarRequirementKind::kBlackHoleAccretionAtLeast:
      if (requirement.species_mask != 0U || requirement.particle_flags_mask != 0U ||
          !std::isfinite(requirement.threshold_code) || requirement.threshold_code < 0.0) {
        throw std::invalid_argument("ModuleSidecarRegistry.upsert: black-hole-accretion sidecar requirement has invalid predicate parameters");
      }
      break;
    case ModuleSidecarRequirementKind::kParticleFlagMask:
      if (requirement.particle_flags_mask == 0U || requirement.species_mask != 0U || requirement.threshold_code != 0.0) {
        throw std::invalid_argument("ModuleSidecarRegistry.upsert: particle-flag sidecar requirement needs a nonzero flag mask only");
      }
      break;
    default:
      throw std::invalid_argument("ModuleSidecarRegistry.upsert: unknown module sidecar requirement kind");
  }
}

void validateModuleSidecarBlock(const ModuleSidecarBlock& block) {
  if (block.module_name.empty()) {
    throw std::invalid_argument("ModuleSidecarRegistry.upsert: module_name cannot be empty");
  }
  if (!block.isParticleIndexed()) {
    if (block.row_stride_bytes != 0U || !block.particle_id_by_row.empty()) {
      throw std::invalid_argument("ModuleSidecarRegistry.upsert: non-indexed block has particle row metadata");
    }
    if (block.required_species_mask != 0U || block.requirement.kind != ModuleSidecarRequirementKind::kSparse) {
      throw std::invalid_argument("ModuleSidecarRegistry.upsert: module sidecar coverage requirements are valid only for particle-indexed sidecars");
    }
    return;
  }

  if (!block.particle_indexed) {
    throw std::invalid_argument("ModuleSidecarRegistry.upsert: particle-indexed sidecar must set particle_indexed=true");
  }
  if (block.row_stride_bytes == 0U) {
    throw std::invalid_argument("ModuleSidecarRegistry.upsert: particle-indexed sidecar row_stride_bytes must be positive");
  }
  const std::uint32_t valid_species_mask = validModuleRequirementSpeciesMask();
  if ((block.required_species_mask & ~valid_species_mask) != 0U) {
    throw std::invalid_argument("ModuleSidecarRegistry.upsert: required_species_mask contains invalid species bits");
  }
  validateModuleSidecarRequirement(block.requirement);
  const std::size_t stride = static_cast<std::size_t>(block.row_stride_bytes);
  if (block.particle_id_by_row.size() > std::numeric_limits<std::size_t>::max() / stride ||
      block.payload.size() != stride * block.particle_id_by_row.size()) {
    throw std::invalid_argument("ModuleSidecarRegistry.upsert: particle-indexed payload size must equal row_stride_bytes * row_count");
  }
  std::unordered_set<std::uint64_t> seen;
  seen.reserve(block.particle_id_by_row.size());
  for (const std::uint64_t particle_id : block.particle_id_by_row) {
    if (particle_id == 0U) {
      throw std::invalid_argument("ModuleSidecarRegistry.upsert: zero particle_id is invalid in particle-indexed sidecar");
    }
    if (!seen.insert(particle_id).second) {
      throw std::invalid_argument("ModuleSidecarRegistry.upsert: duplicate particle_id in particle-indexed sidecar");
    }
  }
}

}  // namespace

void ModuleSidecarRegistry::upsert(ModuleSidecarBlock block) {
  validateModuleSidecarBlock(block);
  m_sidecars[block.module_name] = std::move(block);
}

void ModuleSidecarRegistry::clear() noexcept {
  m_sidecars.clear();
}

const ModuleSidecarBlock* ModuleSidecarRegistry::find(std::string_view module_name) const {
  const auto it = m_sidecars.find(std::string(module_name));
  if (it == m_sidecars.end()) {
    return nullptr;
  }
  return &it->second;
}

std::size_t ModuleSidecarRegistry::size() const noexcept { return m_sidecars.size(); }

std::vector<const ModuleSidecarBlock*> ModuleSidecarRegistry::blocksSortedByName() const {
  std::vector<const ModuleSidecarBlock*> ordered;
  ordered.reserve(m_sidecars.size());
  for (const auto& [name, block] : m_sidecars) {
    (void)name;
    ordered.push_back(&block);
  }
  std::sort(
      ordered.begin(),
      ordered.end(),
      [](const ModuleSidecarBlock* lhs, const ModuleSidecarBlock* rhs) {
        return lhs->module_name < rhs->module_name;
      });
  return ordered;
}

}  // namespace cosmosim::core
