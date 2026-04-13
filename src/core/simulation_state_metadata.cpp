#include "cosmosim/core/simulation_state.hpp"

#include <charconv>
#include <cstdlib>
#include <limits>
#include <sstream>
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
    if (equal_pos == std::string::npos || equal_pos == 0 || equal_pos == trimmed.size() - 1) {
      throw std::invalid_argument("StateMetadata.deserialize: malformed line '" + line + "'");
    }

    const std::string key = trim(trimmed.substr(0, equal_pos));
    const std::string value = trim(trimmed.substr(equal_pos + 1));

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

void ModuleSidecarRegistry::upsert(ModuleSidecarBlock block) {
  if (block.module_name.empty()) {
    throw std::invalid_argument("ModuleSidecarRegistry.upsert: module_name cannot be empty");
  }
  m_sidecars[block.module_name] = std::move(block);
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
