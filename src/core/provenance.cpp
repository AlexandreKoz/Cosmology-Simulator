#include "cosmosim/core/provenance.hpp"

#include <algorithm>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string_view>
#include <stdexcept>
#include <thread>

#include "cosmosim/core/build_config.hpp"

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
  for (unsigned char c : text) {
    hash ^= c;
    hash *= k_prime;
  }
  return hash;
}

[[nodiscard]] std::string toHex(std::uint64_t value) {
  std::ostringstream stream;
  stream << std::hex << std::setfill('0') << std::setw(16) << value;
  return stream.str();
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
  gmtime_s(&now_tm, &now_time);
#else
  gmtime_r(&now_time, &now_tm);
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

ProvenanceRecord makeProvenanceRecord(
    const std::string& config_hash_hex,
    const std::string& git_sha,
    int rank) {
  ProvenanceRecord record;
  record.git_sha = git_sha;
  record.compiler_id = collectCompilerId();
  record.compiler_version = collectCompilerVersion();
  record.build_preset = COSMOSIM_BUILD_PRESET;
  record.enabled_features = std::string("mpi=") + std::to_string(COSMOSIM_ENABLE_MPI) +
                            ",hdf5=" + std::to_string(COSMOSIM_ENABLE_HDF5) +
                            ",fftw=" + std::to_string(COSMOSIM_ENABLE_FFTW) +
                            ",cuda=" + std::to_string(COSMOSIM_ENABLE_CUDA) +
                            ",python=" + std::to_string(COSMOSIM_ENABLE_PYTHON);
  record.config_hash_hex = config_hash_hex;
  record.timestamp_utc = utcTimestampNowIso8601();
  record.hardware_summary = collectHardwareSummary();
  record.author_rank = rank;
  return record;
}

std::string serializeProvenanceRecord(const ProvenanceRecord& record) {
  std::ostringstream stream;
  stream << "schema_version=" << record.schema_version << '\n';
  stream << "git_sha=" << record.git_sha << '\n';
  stream << "compiler_id=" << record.compiler_id << '\n';
  stream << "compiler_version=" << record.compiler_version << '\n';
  stream << "build_preset=" << record.build_preset << '\n';
  stream << "enabled_features=" << record.enabled_features << '\n';
  stream << "config_hash_hex=" << record.config_hash_hex << '\n';
  stream << "timestamp_utc=" << record.timestamp_utc << '\n';
  stream << "hardware_summary=" << record.hardware_summary << '\n';
  stream << "author_rank=" << record.author_rank << '\n';
  stream << "gravity_treepm_pm_grid=" << record.gravity_treepm_pm_grid << '\n';
  stream << "gravity_treepm_assignment_scheme=" << record.gravity_treepm_assignment_scheme << '\n';
  stream << "gravity_treepm_window_deconvolution="
         << (record.gravity_treepm_window_deconvolution ? "true" : "false") << '\n';
  stream << "gravity_treepm_asmth_cells=" << record.gravity_treepm_asmth_cells << '\n';
  stream << "gravity_treepm_rcut_cells=" << record.gravity_treepm_rcut_cells << '\n';
  stream << "gravity_treepm_mesh_spacing_mpc_comoving=" << record.gravity_treepm_mesh_spacing_mpc_comoving << '\n';
  stream << "gravity_treepm_split_scale_mpc_comoving=" << record.gravity_treepm_split_scale_mpc_comoving << '\n';
  stream << "gravity_treepm_cutoff_radius_mpc_comoving=" << record.gravity_treepm_cutoff_radius_mpc_comoving << '\n';
  stream << "gravity_treepm_update_cadence_steps=" << record.gravity_treepm_update_cadence_steps << '\n';
  stream << "gravity_treepm_pm_decomposition_mode=" << record.gravity_treepm_pm_decomposition_mode << '\n';
  stream << "gravity_treepm_tree_exchange_batch_bytes="
         << record.gravity_treepm_tree_exchange_batch_bytes << '\n';
  stream << "gravity_softening_policy=" << record.gravity_softening_policy << '\n';
  stream << "gravity_softening_kernel=" << record.gravity_softening_kernel << '\n';
  stream << "gravity_softening_epsilon_kpc_comoving=" << record.gravity_softening_epsilon_kpc_comoving << '\n';
  stream << "gravity_pm_fft_backend=" << record.gravity_pm_fft_backend << '\n';
  return stream.str();
}


ProvenanceRecord deserializeProvenanceRecord(std::string_view text) {
  ProvenanceRecord record;
  std::istringstream input{std::string(text)};
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }
    const std::size_t pos = line.find('=');
    if (pos == std::string::npos) {
      continue;
    }
    const std::string key = trim(line.substr(0, pos));
    const std::string value = trim(line.substr(pos + 1));
    if (key == "schema_version") {
      record.schema_version = value;
    } else if (key == "git_sha") {
      record.git_sha = value;
    } else if (key == "compiler_id") {
      record.compiler_id = value;
    } else if (key == "compiler_version") {
      record.compiler_version = value;
    } else if (key == "build_preset") {
      record.build_preset = value;
    } else if (key == "enabled_features") {
      record.enabled_features = value;
    } else if (key == "config_hash_hex") {
      record.config_hash_hex = value;
    } else if (key == "timestamp_utc") {
      record.timestamp_utc = value;
    } else if (key == "hardware_summary") {
      record.hardware_summary = value;
    } else if (key == "author_rank") {
      record.author_rank = std::stoi(value);
    } else if (key == "gravity_treepm_pm_grid") {
      record.gravity_treepm_pm_grid = std::stoi(value);
    } else if (key == "gravity_treepm_assignment_scheme") {
      record.gravity_treepm_assignment_scheme = value;
    } else if (key == "gravity_treepm_window_deconvolution") {
      record.gravity_treepm_window_deconvolution = value == "true";
    } else if (key == "gravity_treepm_asmth_cells") {
      record.gravity_treepm_asmth_cells = std::stod(value);
    } else if (key == "gravity_treepm_rcut_cells") {
      record.gravity_treepm_rcut_cells = std::stod(value);
    } else if (key == "gravity_treepm_mesh_spacing_mpc_comoving") {
      record.gravity_treepm_mesh_spacing_mpc_comoving = std::stod(value);
    } else if (key == "gravity_treepm_split_scale_mpc_comoving") {
      record.gravity_treepm_split_scale_mpc_comoving = std::stod(value);
    } else if (key == "gravity_treepm_cutoff_radius_mpc_comoving") {
      record.gravity_treepm_cutoff_radius_mpc_comoving = std::stod(value);
    } else if (key == "gravity_treepm_update_cadence_steps") {
      record.gravity_treepm_update_cadence_steps = std::stoi(value);
    } else if (key == "gravity_treepm_pm_decomposition_mode") {
      record.gravity_treepm_pm_decomposition_mode = value;
    } else if (key == "gravity_treepm_tree_exchange_batch_bytes") {
      record.gravity_treepm_tree_exchange_batch_bytes = static_cast<std::uint64_t>(std::stoull(value));
    } else if (key == "gravity_softening_policy") {
      record.gravity_softening_policy = value;
    } else if (key == "gravity_softening_kernel") {
      record.gravity_softening_kernel = value;
    } else if (key == "gravity_softening_epsilon_kpc_comoving") {
      record.gravity_softening_epsilon_kpc_comoving = std::stod(value);
    } else if (key == "gravity_pm_fft_backend") {
      record.gravity_pm_fft_backend = value;
    }
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
  std::ofstream output(path);
  if (!output) {
    throw std::runtime_error("failed to write provenance record: " + path.string());
  }
  output << serializeProvenanceRecord(record);
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
