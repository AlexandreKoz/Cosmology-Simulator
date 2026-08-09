#include "cosmosim/core/provenance.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string_view>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>

#include "cosmosim/core/build_config.hpp"
#include "core/internal/sha256.hpp"
#include "core/internal/text_codec.hpp"

#include <cstdlib>
#if !defined(_WIN32)
#include <unistd.h>
#endif

namespace cosmosim::core {
namespace {

[[nodiscard]] std::string trim(const std::string& input) {
  auto begin = input.begin();
  while (begin != input.end() && std::isspace(static_cast<unsigned char>(*begin)) != 0) {
    ++begin;
  }

  auto end = input.end();
  while (end != begin && std::isspace(static_cast<unsigned char>(*(end - 1))) != 0) {
    --end;
  }

  return std::string(begin, end);
}

[[nodiscard]] std::uint64_t fnv1a64(const std::string& text) {
  constexpr std::uint64_t k_offset_basis = 14695981039346656037ull;
  constexpr std::uint64_t k_prime = 1099511628211ull;
  std::uint64_t hash = k_offset_basis;
  for (const char raw_c : text) {
    const auto c = static_cast<unsigned char>(raw_c);
    hash ^= static_cast<std::uint64_t>(c);
    hash *= k_prime;
  }
  return hash;
}

[[nodiscard]] std::string toHex(std::uint64_t value) {
  std::ostringstream stream;
  stream << std::hex << std::setfill('0') << std::setw(16) << value;
  return stream.str();
}

[[nodiscard]] std::string escapeMultilineV6(std::string_view text) {
  std::string out;
  out.reserve(text.size());
  for (const char c : text) {
    if (c == '\\') {
      out += "\\\\";
    } else if (c == '\n') {
      out += "\\n";
    } else {
      out.push_back(c);
    }
  }
  return out;
}

[[nodiscard]] std::string unescapeMultilineV6Strict(
    std::string_view text,
    std::string_view context) {
  std::string out;
  out.reserve(text.size());
  for (std::size_t i = 0; i < text.size(); ++i) {
    if (text[i] != '\\') {
      out.push_back(text[i]);
      continue;
    }
    if (i + 1U >= text.size()) {
      throw std::invalid_argument(std::string(context) + ": dangling escape delimiter");
    }
    const char next = text[++i];
    if (next == 'n') {
      out.push_back('\n');
    } else if (next == '\\') {
      out.push_back('\\');
    } else {
      throw std::invalid_argument(
          std::string(context) + ": unknown provenance_v6 escape sequence");
    }
  }
  return out;
}

[[nodiscard]] bool isSha256Hex(std::string_view value) {
  return value.size() == 64U && std::all_of(value.begin(), value.end(), [](char c) {
    const unsigned char uc = static_cast<unsigned char>(c);
    return (uc >= static_cast<unsigned char>('0') && uc <= static_cast<unsigned char>('9')) ||
        (uc >= static_cast<unsigned char>('a') && uc <= static_cast<unsigned char>('f'));
  });
}

[[nodiscard]] std::string collectCpuModel() {
#if defined(__linux__)
  std::ifstream input("/proc/cpuinfo");
  std::string line;
  while (std::getline(input, line)) {
    const auto pos = line.find(':');
    if (pos == std::string::npos) {
      continue;
    }
    const std::string key = trim(line.substr(0, pos));
    if (key == "model name" || key == "Hardware") {
      const std::string value = trim(line.substr(pos + 1U));
      if (!value.empty()) {
        return value;
      }
    }
  }
#elif defined(_WIN32)
  if (const char* value = std::getenv("PROCESSOR_IDENTIFIER"); value != nullptr && *value != '\0') {
    return value;
  }
#endif
  return "unavailable";
}


[[nodiscard]] std::uint32_t collectPhysicalCoreCount() {
#if defined(__linux__)
  std::ifstream input("/proc/cpuinfo");
  std::set<std::pair<std::string, std::string>> cores;
  std::string physical_id;
  std::string core_id;
  std::string line;
  auto commit = [&]() {
    if (!physical_id.empty() && !core_id.empty()) {
      cores.emplace(physical_id, core_id);
    }
    physical_id.clear();
    core_id.clear();
  };
  while (std::getline(input, line)) {
    if (line.empty()) {
      commit();
      continue;
    }
    const auto pos = line.find(':');
    if (pos == std::string::npos) {
      continue;
    }
    const std::string key = trim(line.substr(0, pos));
    const std::string value = trim(line.substr(pos + 1U));
    if (key == "physical id") {
      physical_id = value;
    } else if (key == "core id") {
      core_id = value;
    }
  }
  commit();
  if (!cores.empty() && cores.size() <= std::numeric_limits<std::uint32_t>::max()) {
    return static_cast<std::uint32_t>(cores.size());
  }
#endif
  return 0U;
}

[[nodiscard]] std::uint64_t collectSystemRamBytes() noexcept {
#if defined(__linux__) || defined(__APPLE__)
  const long pages = sysconf(_SC_PHYS_PAGES);
  const long page_size = sysconf(_SC_PAGESIZE);
  if (pages > 0 && page_size > 0) {
    const auto upages = static_cast<std::uint64_t>(pages);
    const auto usize = static_cast<std::uint64_t>(page_size);
    if (upages <= std::numeric_limits<std::uint64_t>::max() / usize) {
      return upages * usize;
    }
  }
#endif
  return 0U;
}

[[nodiscard]] std::string collectHostName() {
#if !defined(_WIN32)
  char name[256]{};
  if (gethostname(name, sizeof(name)) == 0) {
    name[sizeof(name) - 1U] = '\0';
    if (name[0] != '\0') {
      return name;
    }
  }
#else
  if (const char* value = std::getenv("COMPUTERNAME"); value != nullptr && *value != '\0') {
    return value;
  }
#endif
  return "unavailable";
}

[[nodiscard]] int parseIntStrict(std::string_view value, std::string_view key) {
  int parsed = 0;
  const auto [ptr, ec] = std::from_chars(value.data(), value.data() + value.size(), parsed);
  if (ec != std::errc{} || ptr != value.data() + value.size()) {
    throw std::invalid_argument("invalid integer provenance value for key '" + std::string(key) + "'");
  }
  return parsed;
}

[[nodiscard]] std::uint64_t parseUint64Strict(std::string_view value, std::string_view key) {
  std::uint64_t parsed = 0;
  const auto [ptr, ec] = std::from_chars(value.data(), value.data() + value.size(), parsed);
  if (ec != std::errc{} || ptr != value.data() + value.size()) {
    throw std::invalid_argument("invalid uint64 provenance value for key '" + std::string(key) + "'");
  }
  return parsed;
}

[[nodiscard]] std::uint32_t parseUint32Strict(std::string_view value, std::string_view key) {
  const std::uint64_t parsed = parseUint64Strict(value, key);
  if (parsed > std::numeric_limits<std::uint32_t>::max()) {
    throw std::invalid_argument("uint32 provenance value out of range for key '" + std::string(key) + "'");
  }
  return static_cast<std::uint32_t>(parsed);
}

[[nodiscard]] double parseDoubleStrict(std::string_view value, std::string_view key) {
  std::string text(value);
  char* end = nullptr;
  const double parsed = std::strtod(text.c_str(), &end);
  if (end != text.c_str() + static_cast<std::ptrdiff_t>(text.size()) || !std::isfinite(parsed)) {
    throw std::invalid_argument("invalid finite floating provenance value for key '" + std::string(key) + "'");
  }
  return parsed;
}

[[nodiscard]] bool parseBoolStrict(std::string_view value, std::string_view key) {
  if (value == "true") {
    return true;
  }
  if (value == "false") {
    return false;
  }
  throw std::invalid_argument("invalid boolean provenance value for key '" + std::string(key) + "'");
}

}  // namespace

std::string collectCompilerId() {
#if defined(__clang__)
  return "clang";
#elif defined(__GNUC__)
  return "gcc";
#elif defined(_MSC_VER)
  return "msvc";
#else
  return "unknown";
#endif
}

std::string collectCompilerVersion() {
#if defined(__clang__)
  return std::to_string(__clang_major__) + "." + std::to_string(__clang_minor__) + "." +
         std::to_string(__clang_patchlevel__);
#elif defined(__GNUC__)
  return std::to_string(__GNUC__) + "." + std::to_string(__GNUC_MINOR__) + "." +
         std::to_string(__GNUC_PATCHLEVEL__);
#elif defined(_MSC_VER)
  return std::to_string(_MSC_VER);
#else
  return "unknown";
#endif
}

std::string collectHardwareSummary() {
  const unsigned int threads = std::thread::hardware_concurrency();
  std::ostringstream stream;
  stream << "logical_threads=" << threads;
  return stream.str();
}

std::string utcTimestampNowIso8601() {
  const auto now = std::chrono::system_clock::now();
  const std::time_t now_time = std::chrono::system_clock::to_time_t(now);
  std::tm now_tm{};
#if defined(_WIN32)
  if (gmtime_s(&now_tm, &now_time) != 0) {
    throw std::runtime_error("failed to convert provenance timestamp to UTC");
  }
#else
  if (gmtime_r(&now_time, &now_tm) == nullptr) {
    throw std::runtime_error("failed to convert provenance timestamp to UTC");
  }
#endif

  std::ostringstream stream;
  stream << std::put_time(&now_tm, "%Y-%m-%dT%H:%M:%SZ");
  return stream.str();
}

std::uint64_t stableConfigHash(const std::string& normalized_config_text) {
  return fnv1a64(normalized_config_text);
}

std::string stableConfigHashHex(const std::string& normalized_config_text) {
  return toHex(stableConfigHash(normalized_config_text));
}

std::string strongConfigHashSha256Hex(std::string_view normalized_config_text) {
  return internal::sha256Hex(normalized_config_text);
}

ProvenanceRecord makeProvenanceRecord(
    const std::string& config_hash_hex,
    const std::string& git_sha,
    int rank,
    std::string_view normalized_config_text) {
  ProvenanceRecord record;
  record.git_sha = git_sha;
  record.compiler_id = collectCompilerId();
  record.compiler_version = collectCompilerVersion();
  record.build_preset = COSMOSIM_BUILD_PRESET;
  record.enabled_features = std::string("mpi=") + std::to_string(COSMOSIM_ENABLE_MPI) +
                            ",hdf5=" + std::to_string(COSMOSIM_ENABLE_HDF5) +
                            ",fftw=" + std::to_string(COSMOSIM_ENABLE_FFTW) +
                            ",cuda=" + std::to_string(COSMOSIM_ENABLE_CUDA) +
                            ",python=" + std::to_string(COSMOSIM_ENABLE_PYTHON) +
                            ",openmp=" + std::to_string(COSMOSIM_HAVE_OPENMP);
  record.config_hash_hex = config_hash_hex;
  record.normalized_config_hash_hex = config_hash_hex;
  record.integrity_digest_algorithm = "sha256";
  if (!normalized_config_text.empty()) {
    record.normalized_config_sha256_hex = strongConfigHashSha256Hex(normalized_config_text);
    record.normalized_config = std::string(normalized_config_text);
  }
  record.timestamp_utc = utcTimestampNowIso8601();
  record.hardware_summary = collectHardwareSummary();
  record.compiler_flags = COSMOSIM_CXX_FLAGS[0] == '\0' ? "unavailable" : COSMOSIM_CXX_FLAGS;
  record.cpu_model = collectCpuModel();
  record.logical_thread_count = std::thread::hardware_concurrency();
  record.physical_core_count = collectPhysicalCoreCount();
  record.system_ram_bytes = collectSystemRamBytes();
  record.host_name = collectHostName();
  record.mpi_summary = COSMOSIM_ENABLE_MPI ? "compiled" : "disabled";
  record.gpu_summary = COSMOSIM_ENABLE_CUDA ? "compiled_runtime_not_queried" : "disabled";
  record.author_rank = rank;
  return record;
}

std::string serializeProvenanceRecord(const ProvenanceRecord& record) {
  const bool is_v6 = record.schema_version == "provenance_v6";
  const bool is_v7 = record.schema_version == "provenance_v7";
  if (!is_v6 && !is_v7) {
    throw std::invalid_argument("serializeProvenanceRecord: unsupported provenance schema '" + record.schema_version + "'");
  }
  if (is_v7) {
    if (record.integrity_digest_algorithm != "sha256") {
      throw std::invalid_argument("serializeProvenanceRecord: provenance_v7 requires integrity_digest_algorithm=sha256");
    }
    if (record.normalized_config_sha256_hex != "unavailable" &&
        !isSha256Hex(record.normalized_config_sha256_hex)) {
      throw std::invalid_argument(
          "serializeProvenanceRecord: provenance_v7 normalized_config_sha256_hex must be 64 lowercase hex digits or unavailable");
    }
    if (!record.normalized_config.empty()) {
      if (record.normalized_config_sha256_hex == "unavailable" ||
          strongConfigHashSha256Hex(record.normalized_config) != record.normalized_config_sha256_hex) {
        throw std::invalid_argument(
            "serializeProvenanceRecord: provenance_v7 normalized_config SHA-256 mismatch");
      }
    }
    if (record.mpi_world_size <= 0 || record.mpi_node_local_rank < 0) {
      throw std::invalid_argument(
          "serializeProvenanceRecord: provenance_v7 MPI topology fields must be non-negative and world size positive");
    }
  }

  std::ostringstream stream;
  stream << std::setprecision(std::numeric_limits<double>::max_digits10);
  const auto write_string = [&](std::string_view key, std::string_view value) {
    stream << key << '=';
    if (is_v7) {
      stream << internal::escapeTextLine(value);
    } else {
      stream << value;
    }
    stream << '\n';
  };
  const auto write_multiline = [&](std::string_view key, std::string_view value) {
    stream << key << '=';
    if (is_v7) {
      stream << internal::escapeTextLine(value);
    } else {
      stream << escapeMultilineV6(value);
    }
    stream << '\n';
  };

  stream << "schema_version=" << record.schema_version << '\n';
  write_string("config_schema_name", record.config_schema_name);
  write_string("config_schema_version", record.config_schema_version);
  write_string("git_sha", record.git_sha);
  write_string("compiler_id", record.compiler_id);
  write_string("compiler_version", record.compiler_version);
  write_string("build_preset", record.build_preset);
  write_string("enabled_features", record.enabled_features);
  write_string("config_hash_hex", record.config_hash_hex);
  write_string("normalized_config_hash_hex", record.normalized_config_hash_hex);
  if (is_v7) {
    write_string("integrity_digest_algorithm", record.integrity_digest_algorithm);
    write_string("normalized_config_sha256_hex", record.normalized_config_sha256_hex);
  }
  write_multiline("raw_input_config", record.raw_input_config);
  write_multiline("normalized_config", record.normalized_config);
  write_multiline("derived_runtime_state", record.derived_runtime_state);
  write_string("timestamp_utc", record.timestamp_utc);
  write_string("hardware_summary", record.hardware_summary);
  if (is_v7) {
    write_string("compiler_flags", record.compiler_flags);
    write_string("cpu_model", record.cpu_model);
    stream << "logical_thread_count=" << record.logical_thread_count << '\n';
    stream << "physical_core_count=" << record.physical_core_count << '\n';
    stream << "system_ram_bytes=" << record.system_ram_bytes << '\n';
    write_string("host_name", record.host_name);
    write_string("gpu_summary", record.gpu_summary);
    write_string("cuda_runtime_version", record.cuda_runtime_version);
    write_string("cuda_driver_version", record.cuda_driver_version);
    write_string("mpi_summary", record.mpi_summary);
    stream << "mpi_world_size=" << record.mpi_world_size << '\n';
    stream << "mpi_node_local_rank=" << record.mpi_node_local_rank << '\n';
    write_string("deterministic_mode", record.deterministic_mode);
  }
  stream << "author_rank=" << record.author_rank << '\n';
  stream << "gravity_treepm_pm_grid=" << record.gravity_treepm_pm_grid << '\n';
  stream << "gravity_treepm_pm_grid_nx=" << record.gravity_treepm_pm_grid_nx << '\n';
  stream << "gravity_treepm_pm_grid_ny=" << record.gravity_treepm_pm_grid_ny << '\n';
  stream << "gravity_treepm_pm_grid_nz=" << record.gravity_treepm_pm_grid_nz << '\n';
  write_string("gravity_treepm_assignment_scheme", record.gravity_treepm_assignment_scheme);
  stream << "gravity_treepm_window_deconvolution="
         << (record.gravity_treepm_window_deconvolution ? "true" : "false") << '\n';
  stream << "gravity_treepm_asmth_cells=" << record.gravity_treepm_asmth_cells << '\n';
  stream << "gravity_treepm_rcut_cells=" << record.gravity_treepm_rcut_cells << '\n';
  stream << "gravity_treepm_mesh_spacing_mpc_comoving=" << record.gravity_treepm_mesh_spacing_mpc_comoving << '\n';
  stream << "gravity_treepm_mesh_spacing_x_mpc_comoving=" << record.gravity_treepm_mesh_spacing_x_mpc_comoving << '\n';
  stream << "gravity_treepm_mesh_spacing_y_mpc_comoving=" << record.gravity_treepm_mesh_spacing_y_mpc_comoving << '\n';
  stream << "gravity_treepm_mesh_spacing_z_mpc_comoving=" << record.gravity_treepm_mesh_spacing_z_mpc_comoving << '\n';
  stream << "gravity_treepm_split_scale_mpc_comoving=" << record.gravity_treepm_split_scale_mpc_comoving << '\n';
  stream << "gravity_treepm_cutoff_radius_mpc_comoving=" << record.gravity_treepm_cutoff_radius_mpc_comoving << '\n';
  stream << "gravity_treepm_update_cadence_steps=" << record.gravity_treepm_update_cadence_steps << '\n';
  write_string("gravity_treepm_tree_opening_criterion", record.gravity_treepm_tree_opening_criterion);
  stream << "gravity_treepm_tree_opening_theta=" << record.gravity_treepm_tree_opening_theta << '\n';
  stream << "gravity_treepm_tree_relative_force_tolerance=" << record.gravity_treepm_tree_relative_force_tolerance << '\n';
  stream << "gravity_treepm_tree_relative_force_acceleration_floor=" << record.gravity_treepm_tree_relative_force_acceleration_floor << '\n';
  write_string("gravity_treepm_pm_decomposition_mode", record.gravity_treepm_pm_decomposition_mode);
  stream << "gravity_treepm_tree_exchange_batch_bytes=" << record.gravity_treepm_tree_exchange_batch_bytes << '\n';
  write_string("gravity_softening_policy", record.gravity_softening_policy);
  write_string("gravity_softening_kernel", record.gravity_softening_kernel);
  stream << "gravity_softening_epsilon_kpc_comoving=" << record.gravity_softening_epsilon_kpc_comoving << '\n';
  write_string("gravity_pm_fft_backend", record.gravity_pm_fft_backend);
  stream << "gravity_treepm_decomposition_epoch=" << record.gravity_treepm_decomposition_epoch << '\n';
  stream << "gravity_treepm_restart_world_size=" << record.gravity_treepm_restart_world_size << '\n';
  write_string("gravity_treepm_restart_pm_grid", record.gravity_treepm_restart_pm_grid);
  write_string("gravity_treepm_restart_slab_signature", record.gravity_treepm_restart_slab_signature);
  stream << "gravity_treepm_restart_kick_opportunity=" << record.gravity_treepm_restart_kick_opportunity << '\n';
  stream << "gravity_treepm_restart_field_version=" << record.gravity_treepm_restart_field_version << '\n';
  write_string("gravity_treepm_long_range_restart_policy", record.gravity_treepm_long_range_restart_policy);
  write_string("zoom_long_range_strategy", record.zoom_long_range_strategy);
  stream << "zoom_region_center_x_mpc_comoving=" << record.zoom_region_center_x_mpc_comoving << '\n';
  stream << "zoom_region_center_y_mpc_comoving=" << record.zoom_region_center_y_mpc_comoving << '\n';
  stream << "zoom_region_center_z_mpc_comoving=" << record.zoom_region_center_z_mpc_comoving << '\n';
  stream << "zoom_region_radius_mpc_comoving=" << record.zoom_region_radius_mpc_comoving << '\n';
  write_string("zoom_focused_pm_grid", record.zoom_focused_pm_grid);
  stream << "zoom_contamination_radius_mpc_comoving=" << record.zoom_contamination_radius_mpc_comoving << '\n';
  return stream.str();
}

ProvenanceRecord deserializeProvenanceRecord(std::string_view text) {
  std::unordered_map<std::string, std::string> raw_values;
  std::unordered_set<std::string> seen_keys;
  std::istringstream input{std::string(text)};
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }
    const std::size_t pos = line.find('=');
    if (pos == std::string::npos || pos == 0U) {
      throw std::invalid_argument("malformed provenance line: '" + line + "'");
    }
    const std::string key = trim(line.substr(0, pos));
    if (key.empty() || !seen_keys.insert(key).second) {
      throw std::invalid_argument("duplicate or empty provenance key: '" + key + "'");
    }
    raw_values.emplace(key, line.substr(pos + 1U));
  }
  if (!raw_values.contains("schema_version")) {
    throw std::invalid_argument("missing required provenance key: 'schema_version'");
  }
  const std::string schema = trim(raw_values.at("schema_version"));
  const bool is_v6 = schema == "provenance_v6";
  const bool is_v7 = schema == "provenance_v7";
  if (!is_v6 && !is_v7) {
    throw std::invalid_argument("unsupported provenance schema: '" + schema + "'");
  }

  const auto string_value = [&](const std::string& raw, std::string_view key) {
    if (is_v7) {
      return internal::unescapeTextLineStrict(raw, std::string("provenance key '") + std::string(key) + "'");
    }
    return trim(raw);
  };
  const auto multiline_value = [&](const std::string& raw, std::string_view key) {
    if (is_v7) {
      return internal::unescapeTextLineStrict(raw, std::string("provenance key '") + std::string(key) + "'");
    }
    return unescapeMultilineV6Strict(trim(raw), std::string("provenance key '") + std::string(key) + "'");
  };

  ProvenanceRecord record;
  record.schema_version = schema;
  for (const auto& [key, raw] : raw_values) {
    if (key == "schema_version") {
      continue;
    }
    const std::string value = trim(raw);
    if (key == "config_schema_name") {
      record.config_schema_name = string_value(raw, key);
    } else if (key == "config_schema_version") {
      record.config_schema_version = string_value(raw, key);
    } else if (key == "git_sha") {
      record.git_sha = string_value(raw, key);
    } else if (key == "compiler_id") {
      record.compiler_id = string_value(raw, key);
    } else if (key == "compiler_version") {
      record.compiler_version = string_value(raw, key);
    } else if (key == "build_preset") {
      record.build_preset = string_value(raw, key);
    } else if (key == "enabled_features") {
      record.enabled_features = string_value(raw, key);
    } else if (key == "config_hash_hex") {
      record.config_hash_hex = string_value(raw, key);
    } else if (key == "normalized_config_hash_hex") {
      record.normalized_config_hash_hex = string_value(raw, key);
    } else if (key == "integrity_digest_algorithm" && is_v7) {
      record.integrity_digest_algorithm = string_value(raw, key);
    } else if (key == "normalized_config_sha256_hex" && is_v7) {
      record.normalized_config_sha256_hex = string_value(raw, key);
    } else if (key == "raw_input_config") {
      record.raw_input_config = multiline_value(raw, key);
    } else if (key == "normalized_config") {
      record.normalized_config = multiline_value(raw, key);
    } else if (key == "derived_runtime_state") {
      record.derived_runtime_state = multiline_value(raw, key);
    } else if (key == "timestamp_utc") {
      record.timestamp_utc = string_value(raw, key);
    } else if (key == "hardware_summary") {
      record.hardware_summary = string_value(raw, key);
    } else if (key == "compiler_flags" && is_v7) {
      record.compiler_flags = string_value(raw, key);
    } else if (key == "cpu_model" && is_v7) {
      record.cpu_model = string_value(raw, key);
    } else if (key == "logical_thread_count" && is_v7) {
      record.logical_thread_count = parseUint32Strict(value, key);
    } else if (key == "physical_core_count" && is_v7) {
      record.physical_core_count = parseUint32Strict(value, key);
    } else if (key == "system_ram_bytes" && is_v7) {
      record.system_ram_bytes = parseUint64Strict(value, key);
    } else if (key == "host_name" && is_v7) {
      record.host_name = string_value(raw, key);
    } else if (key == "gpu_summary" && is_v7) {
      record.gpu_summary = string_value(raw, key);
    } else if (key == "cuda_runtime_version" && is_v7) {
      record.cuda_runtime_version = string_value(raw, key);
    } else if (key == "cuda_driver_version" && is_v7) {
      record.cuda_driver_version = string_value(raw, key);
    } else if (key == "mpi_summary" && is_v7) {
      record.mpi_summary = string_value(raw, key);
    } else if (key == "mpi_world_size" && is_v7) {
      record.mpi_world_size = parseIntStrict(value, key);
    } else if (key == "mpi_node_local_rank" && is_v7) {
      record.mpi_node_local_rank = parseIntStrict(value, key);
    } else if (key == "deterministic_mode" && is_v7) {
      record.deterministic_mode = string_value(raw, key);
    } else if (key == "author_rank") {
      record.author_rank = parseIntStrict(value, key);
    } else if (key == "gravity_treepm_pm_grid") {
      record.gravity_treepm_pm_grid = parseIntStrict(value, key);
    } else if (key == "gravity_treepm_pm_grid_nx") {
      record.gravity_treepm_pm_grid_nx = parseIntStrict(value, key);
    } else if (key == "gravity_treepm_pm_grid_ny") {
      record.gravity_treepm_pm_grid_ny = parseIntStrict(value, key);
    } else if (key == "gravity_treepm_pm_grid_nz") {
      record.gravity_treepm_pm_grid_nz = parseIntStrict(value, key);
    } else if (key == "gravity_treepm_assignment_scheme") {
      record.gravity_treepm_assignment_scheme = string_value(raw, key);
    } else if (key == "gravity_treepm_window_deconvolution") {
      record.gravity_treepm_window_deconvolution = parseBoolStrict(value, key);
    } else if (key == "gravity_treepm_asmth_cells") {
      record.gravity_treepm_asmth_cells = parseDoubleStrict(value, key);
    } else if (key == "gravity_treepm_rcut_cells") {
      record.gravity_treepm_rcut_cells = parseDoubleStrict(value, key);
    } else if (key == "gravity_treepm_mesh_spacing_mpc_comoving") {
      record.gravity_treepm_mesh_spacing_mpc_comoving = parseDoubleStrict(value, key);
    } else if (key == "gravity_treepm_mesh_spacing_x_mpc_comoving") {
      record.gravity_treepm_mesh_spacing_x_mpc_comoving = parseDoubleStrict(value, key);
    } else if (key == "gravity_treepm_mesh_spacing_y_mpc_comoving") {
      record.gravity_treepm_mesh_spacing_y_mpc_comoving = parseDoubleStrict(value, key);
    } else if (key == "gravity_treepm_mesh_spacing_z_mpc_comoving") {
      record.gravity_treepm_mesh_spacing_z_mpc_comoving = parseDoubleStrict(value, key);
    } else if (key == "gravity_treepm_split_scale_mpc_comoving") {
      record.gravity_treepm_split_scale_mpc_comoving = parseDoubleStrict(value, key);
    } else if (key == "gravity_treepm_cutoff_radius_mpc_comoving") {
      record.gravity_treepm_cutoff_radius_mpc_comoving = parseDoubleStrict(value, key);
    } else if (key == "gravity_treepm_update_cadence_steps") {
      record.gravity_treepm_update_cadence_steps = parseIntStrict(value, key);
    } else if (key == "gravity_treepm_tree_opening_criterion") {
      record.gravity_treepm_tree_opening_criterion = string_value(raw, key);
    } else if (key == "gravity_treepm_tree_opening_theta") {
      record.gravity_treepm_tree_opening_theta = parseDoubleStrict(value, key);
    } else if (key == "gravity_treepm_tree_relative_force_tolerance") {
      record.gravity_treepm_tree_relative_force_tolerance = parseDoubleStrict(value, key);
    } else if (key == "gravity_treepm_tree_relative_force_acceleration_floor") {
      record.gravity_treepm_tree_relative_force_acceleration_floor = parseDoubleStrict(value, key);
    } else if (key == "gravity_treepm_pm_decomposition_mode") {
      record.gravity_treepm_pm_decomposition_mode = string_value(raw, key);
    } else if (key == "gravity_treepm_tree_exchange_batch_bytes") {
      record.gravity_treepm_tree_exchange_batch_bytes = parseUint64Strict(value, key);
    } else if (key == "gravity_softening_policy") {
      record.gravity_softening_policy = string_value(raw, key);
    } else if (key == "gravity_softening_kernel") {
      record.gravity_softening_kernel = string_value(raw, key);
    } else if (key == "gravity_softening_epsilon_kpc_comoving") {
      record.gravity_softening_epsilon_kpc_comoving = parseDoubleStrict(value, key);
    } else if (key == "gravity_pm_fft_backend") {
      record.gravity_pm_fft_backend = string_value(raw, key);
    } else if (key == "gravity_treepm_decomposition_epoch") {
      record.gravity_treepm_decomposition_epoch = parseUint64Strict(value, key);
    } else if (key == "gravity_treepm_restart_world_size") {
      record.gravity_treepm_restart_world_size = parseIntStrict(value, key);
    } else if (key == "gravity_treepm_restart_pm_grid") {
      record.gravity_treepm_restart_pm_grid = string_value(raw, key);
    } else if (key == "gravity_treepm_restart_slab_signature") {
      record.gravity_treepm_restart_slab_signature = string_value(raw, key);
    } else if (key == "gravity_treepm_restart_kick_opportunity") {
      record.gravity_treepm_restart_kick_opportunity = parseUint64Strict(value, key);
    } else if (key == "gravity_treepm_restart_field_version") {
      record.gravity_treepm_restart_field_version = parseUint64Strict(value, key);
    } else if (key == "gravity_treepm_long_range_restart_policy") {
      record.gravity_treepm_long_range_restart_policy = string_value(raw, key);
    } else if (key == "zoom_long_range_strategy") {
      record.zoom_long_range_strategy = string_value(raw, key);
    } else if (key == "zoom_region_center_x_mpc_comoving") {
      record.zoom_region_center_x_mpc_comoving = parseDoubleStrict(value, key);
    } else if (key == "zoom_region_center_y_mpc_comoving") {
      record.zoom_region_center_y_mpc_comoving = parseDoubleStrict(value, key);
    } else if (key == "zoom_region_center_z_mpc_comoving") {
      record.zoom_region_center_z_mpc_comoving = parseDoubleStrict(value, key);
    } else if (key == "zoom_region_radius_mpc_comoving") {
      record.zoom_region_radius_mpc_comoving = parseDoubleStrict(value, key);
    } else if (key == "zoom_focused_pm_grid") {
      record.zoom_focused_pm_grid = string_value(raw, key);
    } else if (key == "zoom_contamination_radius_mpc_comoving") {
      record.zoom_contamination_radius_mpc_comoving = parseDoubleStrict(value, key);
    } else {
      throw std::invalid_argument("unknown provenance key: '" + key + "'");
    }
  }

  const auto require_key = [&](std::string_view key) {
    if (!seen_keys.contains(std::string(key))) {
      throw std::invalid_argument("missing required provenance key: '" + std::string(key) + "'");
    }
  };
  for (const std::string_view required : {
           "schema_version", "config_schema_name", "config_schema_version", "config_hash_hex"}) {
    require_key(required);
  }
  if (is_v7) {
    for (const std::string_view required : {
             "integrity_digest_algorithm", "normalized_config_sha256_hex", "compiler_flags",
             "cpu_model", "logical_thread_count", "physical_core_count", "system_ram_bytes",
             "host_name", "gpu_summary", "cuda_runtime_version", "cuda_driver_version",
             "mpi_summary", "mpi_world_size", "mpi_node_local_rank", "deterministic_mode"}) {
      require_key(required);
    }
    if (record.integrity_digest_algorithm != "sha256") {
      throw std::invalid_argument("provenance_v7 requires integrity_digest_algorithm=sha256");
    }
    if (record.normalized_config_sha256_hex != "unavailable" &&
        !isSha256Hex(record.normalized_config_sha256_hex)) {
      throw std::invalid_argument("provenance_v7 normalized_config_sha256_hex must be 64 lowercase hex digits or unavailable");
    }
    if (!record.normalized_config.empty()) {
      if (record.normalized_config_sha256_hex == "unavailable") {
        throw std::invalid_argument("provenance_v7 with normalized_config requires a SHA-256 digest");
      }
      const std::string observed = strongConfigHashSha256Hex(record.normalized_config);
      if (observed != record.normalized_config_sha256_hex) {
        throw std::invalid_argument("provenance_v7 normalized_config SHA-256 mismatch");
      }
    }
    if (record.mpi_world_size <= 0 || record.mpi_node_local_rank < 0) {
      throw std::invalid_argument("provenance_v7 MPI topology fields must be non-negative and world size positive");
    }
  }

  if (record.gravity_treepm_pm_grid_nx == 0 && record.gravity_treepm_pm_grid > 0) {
    record.gravity_treepm_pm_grid_nx = record.gravity_treepm_pm_grid;
    record.gravity_treepm_pm_grid_ny = record.gravity_treepm_pm_grid;
    record.gravity_treepm_pm_grid_nz = record.gravity_treepm_pm_grid;
  }
  if (record.gravity_treepm_mesh_spacing_x_mpc_comoving == 0.0 &&
      record.gravity_treepm_mesh_spacing_mpc_comoving > 0.0) {
    record.gravity_treepm_mesh_spacing_x_mpc_comoving = record.gravity_treepm_mesh_spacing_mpc_comoving;
    record.gravity_treepm_mesh_spacing_y_mpc_comoving = record.gravity_treepm_mesh_spacing_mpc_comoving;
    record.gravity_treepm_mesh_spacing_z_mpc_comoving = record.gravity_treepm_mesh_spacing_mpc_comoving;
  }
  if (record.normalized_config_hash_hex.empty()) {
    record.normalized_config_hash_hex = record.config_hash_hex;
  }
  return record;
}

void writeProvenanceRecord(
    const ProvenanceRecord& record,
    const std::filesystem::path& run_directory,
    const std::string& file_name) {
  if (record.author_rank != 0) {
    return;
  }

  std::filesystem::create_directories(run_directory);
  const auto path = run_directory / file_name;
  auto part_path = path;
  part_path += ".part";
  std::ofstream output(part_path, std::ios::out | std::ios::trunc);
  if (!output) {
    throw std::runtime_error("failed to write provenance temporary record: " + part_path.string());
  }
  output << serializeProvenanceRecord(record);
  output.flush();
  if (!output) {
    throw std::runtime_error("failed while writing provenance record: " + part_path.string());
  }
  output.close();
  if (!output) {
    throw std::runtime_error("failed while closing provenance record: " + part_path.string());
  }
  std::error_code ec;
  std::filesystem::rename(part_path, path, ec);
  if (ec) {
    // Do not remove an already-published provenance record to force replacement.
    // On platforms without atomic replacement semantics, preserve the complete
    // prior record and fail closed instead of creating a publication gap.
    std::error_code remove_ec;
    std::filesystem::remove(part_path, remove_ec);
    throw std::runtime_error("failed to atomically finalize provenance record: " + ec.message());
  }
}

ProvenanceRecord readProvenanceRecord(
    const std::filesystem::path& run_directory,
    const std::string& file_name) {
  const auto path = run_directory / file_name;
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("failed to read provenance record: " + path.string());
  }

  std::ostringstream stream;
  stream << input.rdbuf();
  return deserializeProvenanceRecord(stream.str());
}

}  // namespace cosmosim::core
