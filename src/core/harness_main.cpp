#include <charconv>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

#include "cosmosim/cosmosim.hpp"
#include "cosmosim/core/build_config.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
void requireMpiSuccess(int result, const char* operation) {
  if (result != MPI_SUCCESS) {
    throw std::runtime_error(std::string("MPI failure in ") + operation);
  }
}
#endif

void printUsage(std::ostream& out, const char* argv0) {
  out << "Usage: " << argv0 << " <config.param.txt> [options]\n"
      << "\n"
      << "Native CHUI console options:\n"
      << "  --quiet                 Suppress routine [CHUI] informational records.\n"
      << "  --status-every N        Emit a productive-step heartbeat every N steps\n"
      << "                          (0 disables the step-count trigger; default 10).\n"
      << "  --status-seconds SEC    Emit a heartbeat after SEC seconds without console\n"
      << "                          output at a completed step boundary (0 disables; default 30).\n"
      << "  -h, --help              Show this help.\n"
      << "\n"
      << "The simulation configuration remains authoritative; these flags affect only\n"
      << "process-local console presentation.\n";
}

[[nodiscard]] std::uint64_t parseUnsigned(
    std::string_view text,
    std::string_view option_name) {
  std::uint64_t value = 0U;
  const char* begin = text.data();
  const char* end = text.data() + text.size();
  const auto result = std::from_chars(begin, end, value, 10);
  if (text.empty() || result.ec != std::errc{} || result.ptr != end) {
    throw std::invalid_argument(
        std::string(option_name) + " requires a non-negative integer, got '" +
        std::string(text) + "'");
  }
  return value;
}

[[nodiscard]] double parseNonNegativeDouble(
    std::string_view text,
    std::string_view option_name) {
  if (text.empty()) {
    throw std::invalid_argument(std::string(option_name) + " requires a value");
  }
  std::size_t parsed = 0U;
  double value = 0.0;
  try {
    value = std::stod(std::string(text), &parsed);
  } catch (const std::exception&) {
    throw std::invalid_argument(
        std::string(option_name) + " requires a non-negative finite number, got '" +
        std::string(text) + "'");
  }
  if (parsed != text.size() || !std::isfinite(value) || value < 0.0) {
    throw std::invalid_argument(
        std::string(option_name) + " requires a non-negative finite number, got '" +
        std::string(text) + "'");
  }
  return value;
}

struct HarnessCliOptions {
  std::filesystem::path config_path;
  cosmosim::workflows::RuntimeConsoleOptions console;
  bool help = false;
};

[[nodiscard]] HarnessCliOptions parseHarnessCli(int argc, char** argv) {
  HarnessCliOptions parsed;
  parsed.console.enabled = true;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if (arg == "-h" || arg == "--help") {
      parsed.help = true;
      continue;
    }
    if (arg == "--quiet") {
      parsed.console.quiet = true;
      continue;
    }
    if (arg == "--status-every") {
      if (i + 1 >= argc) {
        throw std::invalid_argument("--status-every requires N");
      }
      parsed.console.status_every_steps = parseUnsigned(argv[++i], "--status-every");
      continue;
    }
    if (arg.starts_with("--status-every=")) {
      parsed.console.status_every_steps = parseUnsigned(
          arg.substr(std::string_view("--status-every=").size()),
          "--status-every");
      continue;
    }
    if (arg == "--status-seconds") {
      if (i + 1 >= argc) {
        throw std::invalid_argument("--status-seconds requires SEC");
      }
      parsed.console.status_seconds = parseNonNegativeDouble(argv[++i], "--status-seconds");
      continue;
    }
    if (arg.starts_with("--status-seconds=")) {
      parsed.console.status_seconds = parseNonNegativeDouble(
          arg.substr(std::string_view("--status-seconds=").size()),
          "--status-seconds");
      continue;
    }
    if (!arg.empty() && arg.front() == '-') {
      throw std::invalid_argument("unknown cosmosim_harness option: " + std::string(arg));
    }
    if (!parsed.config_path.empty()) {
      throw std::invalid_argument(
          "cosmosim_harness accepts exactly one config path; unexpected positional argument: " +
          std::string(arg));
    }
    parsed.config_path = std::filesystem::path(arg);
  }

  if (!parsed.help && parsed.config_path.empty()) {
    throw std::invalid_argument("missing required <config.param.txt> path");
  }
  if (!parsed.config_path.empty()) {
    parsed.console.config_path = parsed.config_path.string();
  }
  return parsed;
}

class ExecutableMpiSession {
 public:
  ExecutableMpiSession(int* argc, char*** argv) {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    int finalized = 0;
    requireMpiSuccess(MPI_Finalized(&finalized), "MPI_Finalized");
    if (finalized != 0) {
      throw std::runtime_error("cosmosim_harness cannot start after MPI_Finalize has already completed");
    }

    int initialized = 0;
    requireMpiSuccess(MPI_Initialized(&initialized), "MPI_Initialized");
    if (initialized == 0) {
      int provided = MPI_THREAD_SINGLE;
      const int init_result = MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
      if (init_result != MPI_SUCCESS) {
        throw std::runtime_error("cosmosim_harness MPI_Init_thread failed while requesting MPI_THREAD_FUNNELED");
      }
      m_owns_finalize = true;
      m_thread_level = provided;
    } else {
      requireMpiSuccess(MPI_Query_thread(&m_thread_level), "MPI_Query_thread");
    }

    requireMpiSuccess(MPI_Comm_size(MPI_COMM_WORLD, &m_world_size), "MPI_Comm_size");
    requireMpiSuccess(MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank), "MPI_Comm_rank");
    if (m_thread_level < MPI_THREAD_FUNNELED) {
      std::ostringstream msg;
      msg << "cosmosim_harness MPI thread support is insufficient: expected>=MPI_THREAD_FUNNELED("
          << MPI_THREAD_FUNNELED << "), provided=" << m_thread_level << ", rank=" << m_world_rank << '/'
          << m_world_size;
      if (m_owns_finalize) {
        int finalized_after_init = 0;
        requireMpiSuccess(MPI_Finalized(&finalized_after_init), "MPI_Finalized");
        if (finalized_after_init == 0) {
          requireMpiSuccess(MPI_Finalize(), "MPI_Finalize");
        }
      }
      throw std::runtime_error(msg.str());
    }
#else
    (void)argc;
    (void)argv;
#endif
  }

  ExecutableMpiSession(const ExecutableMpiSession&) = delete;
  ExecutableMpiSession& operator=(const ExecutableMpiSession&) = delete;

  ~ExecutableMpiSession() noexcept {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    if (!m_owns_finalize) {
      return;
    }
    int finalized = 0;
    if (MPI_Finalized(&finalized) == MPI_SUCCESS && finalized == 0) {
      (void)MPI_Finalize();
    }
#endif
  }

  [[nodiscard]] int worldSize() const noexcept { return m_world_size; }
  [[nodiscard]] int worldRank() const noexcept { return m_world_rank; }
  [[nodiscard]] bool isRoot() const noexcept { return m_world_rank == 0; }

  [[nodiscard]] std::string rankPrefix() const {
    std::ostringstream out;
    out << "rank=" << m_world_rank << '/' << m_world_size;
    return out.str();
  }

  void abortDistributed(int error_code) const noexcept {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    if (m_world_size > 1) {
      (void)MPI_Abort(MPI_COMM_WORLD, error_code);
    }
#else
    (void)error_code;
#endif
  }

 private:
  bool m_owns_finalize = false;
  int m_world_size = 1;
  int m_world_rank = 0;
  int m_thread_level = 0;
};

[[nodiscard]] std::string currentRankPrefix() {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized != 0 && finalized == 0) {
    int world_size = 1;
    int world_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    std::ostringstream out;
    out << "rank=" << world_rank << '/' << world_size;
    return out.str();
  }
#endif
  return "rank=0/1";
}

void abortCurrentMpiWorldIfDistributed(int error_code) noexcept {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  int initialized = 0;
  int finalized = 0;
  if (MPI_Initialized(&initialized) == MPI_SUCCESS && MPI_Finalized(&finalized) == MPI_SUCCESS &&
      initialized != 0 && finalized == 0) {
    int world_size = 1;
    if (MPI_Comm_size(MPI_COMM_WORLD, &world_size) == MPI_SUCCESS && world_size > 1) {
      (void)MPI_Abort(MPI_COMM_WORLD, error_code);
    }
  }
#else
  (void)error_code;
#endif
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 3 && std::string_view(argv[1]) == "--print-mpi-ranks-expected") {
      const cosmosim::core::FrozenConfig frozen =
          cosmosim::core::loadFrozenConfigFromFile(std::filesystem::path(argv[2]), {});
      std::cout << frozen.config.parallel.mpi_ranks_expected << '\n';
      return 0;
    }
    ExecutableMpiSession mpi_session(&argc, &argv);

    HarnessCliOptions cli;
    try {
      cli = parseHarnessCli(argc, argv);
    } catch (const std::exception& ex) {
      if (mpi_session.isRoot()) {
        std::cerr << "[CHUI][FATAL] " << ex.what() << "\n\n";
        printUsage(std::cerr, argc > 0 ? argv[0] : "cosmosim_harness");
      }
      mpi_session.abortDistributed(2);
      return 2;
    }

    if (cli.help) {
      if (mpi_session.isRoot()) {
        printUsage(std::cout, argc > 0 ? argv[0] : "cosmosim_harness");
      }
      return 0;
    }

    try {
      const cosmosim::core::FrozenConfig frozen =
          cosmosim::core::loadFrozenConfigFromFile(cli.config_path, {});
      cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
      cosmosim::workflows::ReferenceWorkflowOptions workflow_options;
      workflow_options.console = cli.console;
      static_cast<void>(runner.run(workflow_options));
      return 0;
    } catch (const std::exception& ex) {
      std::cerr << "[CHUI][FATAL] " << mpi_session.rankPrefix()
                << " message=" << ex.what() << '\n';
      mpi_session.abortDistributed(1);
      return 1;
    }
  } catch (const std::exception& ex) {
    std::cerr << "[CHUI][FATAL] " << currentRankPrefix()
              << " message=" << ex.what() << '\n';
    abortCurrentMpiWorldIfDistributed(1);
    return 1;
  }
}
