#include <exception>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>

#include "cosmosim/cosmosim.hpp"
#include "cosmosim/core/build_config.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

int printUsage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " <config.param.txt>\n";
  return 2;
}

class ExecutableMpiSession {
 public:
  ExecutableMpiSession(int* argc, char*** argv) {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized != 0) {
      throw std::runtime_error("cosmosim_harness cannot start after MPI_Finalize has already completed");
    }

    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized == 0) {
      int provided = MPI_THREAD_SINGLE;
      const int init_result = MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
      if (init_result != MPI_SUCCESS) {
        throw std::runtime_error("cosmosim_harness MPI_Init_thread failed while requesting MPI_THREAD_FUNNELED");
      }
      m_owns_finalize = true;
      m_thread_level = provided;
    } else {
      MPI_Query_thread(&m_thread_level);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &m_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
    if (m_thread_level < MPI_THREAD_FUNNELED) {
      std::ostringstream msg;
      msg << "cosmosim_harness MPI thread support is insufficient: expected>=MPI_THREAD_FUNNELED("
          << MPI_THREAD_FUNNELED << "), provided=" << m_thread_level << ", rank=" << m_world_rank << '/'
          << m_world_size;
      if (m_owns_finalize) {
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (finalized == 0) {
          MPI_Finalize();
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

  ~ExecutableMpiSession() {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    if (!m_owns_finalize) {
      return;
    }
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized == 0) {
      MPI_Finalize();
    }
#endif
  }

  [[nodiscard]] int worldSize() const noexcept { return m_world_size; }
  [[nodiscard]] int worldRank() const noexcept { return m_world_rank; }

  [[nodiscard]] std::string rankPrefix() const {
    std::ostringstream out;
    out << "rank=" << m_world_rank << '/' << m_world_size;
    return out.str();
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

}  // namespace

int main(int argc, char** argv) {
  try {
    ExecutableMpiSession mpi_session(&argc, &argv);

    if (argc != 2) {
      return printUsage(argc > 0 ? argv[0] : "cosmosim_harness");
    }

    try {
      const std::filesystem::path config_path = argv[1];
      const cosmosim::core::FrozenConfig frozen = cosmosim::core::loadFrozenConfigFromFile(config_path, {});
      cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
      const cosmosim::workflows::ReferenceWorkflowReport report = runner.run();

      std::cout << cosmosim::core::projectName() << ' ' << cosmosim::core::versionString() << '\n';
      std::cout << "config=" << config_path.string() << '\n';
      std::cout << "run_directory=" << report.run_directory.string() << '\n';
      std::cout << "completed_steps=" << report.completed_steps << '\n';
      std::cout << "normalized_config=" << report.normalized_config_snapshot_path.string() << '\n';
      std::cout << "operational_report=" << report.operational_report_json_path.string() << '\n';
      if (!report.snapshot_path.empty()) {
        std::cout << "last_snapshot=" << report.snapshot_path.string() << '\n';
      }
      if (!report.restart_path.empty()) {
        std::cout << "last_restart=" << report.restart_path.string() << '\n';
      }
      return 0;
    } catch (const std::exception& ex) {
      std::cerr << "cosmosim runtime failed [" << mpi_session.rankPrefix() << "]: " << ex.what() << '\n';
      return 1;
    }
  } catch (const std::exception& ex) {
    std::cerr << "cosmosim runtime failed [" << currentRankPrefix() << "]: " << ex.what() << '\n';
    return 1;
  }
}
