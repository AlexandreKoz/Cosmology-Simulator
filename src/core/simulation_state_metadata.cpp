#include "cosmosim/core/simulation_state.hpp"

#include <charconv>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <unordered_set>
#include <stdexcept>

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
  if (end != temp.c_str() + static_cast<std::ptrdiff_t>(temp.size())) {
    throw std::invalid_argument(
        "StateMetadata.deserialize: invalid floating-point value for key '" + key + "'");
  }
  return parsed;
}
}  // namespace

std::string StateMetadata::serialize() const {
  std::ostringstream out;
  out << "schema_version=" << schema_version << '\n';
  out << "run_name=" << run_name << '\n';
  out << "normalized_config_hash=" << normalized_config_hash << '\n';
  out << "normalized_config_hash_hex=" << normalized_config_hash_hex << '\n';
  out << "step_index=" << step_index << '\n';
  out << "scale_factor=" << scale_factor << '\n';
  out << "snapshot_stem=" << snapshot_stem << '\n';
  out << "restart_stem=" << restart_stem << '\n';
  return out.str();
}

StateMetadata StateMetadata::deserialize(std::string_view text) {
  StateMetadata metadata;
  std::istringstream in{std::string(text)};
  std::string line;

  while (std::getline(in, line)) {
    const std::string trimmed = trim(line);
    if (trimmed.empty()) {
      continue;
    }

    const auto equal_pos = trimmed.find('=');
    if (equal_pos == std::string::npos || equal_pos == 0) {
      throw std::invalid_argument("StateMetadata.deserialize: malformed line '" + line + "'");
    }

    const std::string key = trim(trimmed.substr(0, equal_pos));
    const std::string value = trim(trimmed.substr(equal_pos + 1));
    if (key.empty()) {
      throw std::invalid_argument("StateMetadata.deserialize: malformed line '" + line + "'");
    }

    if (key == "schema_version") {
      metadata.schema_version = parseUint32(value, key);
    } else if (key == "run_name") {
      metadata.run_name = value;
    } else if (key == "normalized_config_hash") {
      metadata.normalized_config_hash = parseUint64(value, key);
    } else if (key == "normalized_config_hash_hex") {
      metadata.normalized_config_hash_hex = value;
    } else if (key == "step_index") {
      metadata.step_index = parseUint64(value, key);
    } else if (key == "scale_factor") {
      metadata.scale_factor = parseDouble(value, key);
    } else if (key == "snapshot_stem") {
      metadata.snapshot_stem = value;
    } else if (key == "restart_stem") {
      metadata.restart_stem = value;
    }
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
  const std::size_t begin = row * stride;
  if (stride == 0U || begin + stride > payload.size()) {
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
      if (requirement.species_mask != 0U || requirement.particle_flags_mask != 0U || requirement.threshold_code < 0.0) {
        throw std::invalid_argument("ModuleSidecarRegistry.upsert: gas-density sidecar requirement has invalid predicate parameters");
      }
      break;
    case ModuleSidecarRequirementKind::kBlackHoleAccretionAtLeast:
      if (requirement.species_mask != 0U || requirement.particle_flags_mask != 0U || requirement.threshold_code < 0.0) {
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
  if (block.payload.size() != stride * block.particle_id_by_row.size()) {
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
