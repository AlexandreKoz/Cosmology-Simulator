#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <sstream>
#include <string>

#include "cosmosim/cosmosim.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {
std::string buildConfigText(int mpi_ranks_expected, std::uint64_t isolated_limit_bytes) {
  std::stringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[mode]\n";
  stream << "mode = isolated_galaxy\n";
  stream << "ic_file = generated\n\n";
  stream << "[numerics]\n";
  stream << "time_begin_code = 0.0\n";
  stream << "time_end_code = 0.02\n";
  stream << "max_global_steps = 2\n";
  stream << "hierarchical_max_rung = 0\n";
  stream << "treepm_pm_grid = 12\n";
  stream << "treepm_asmth_cells = 1.25\n";
  stream << "treepm_rcut_cells = 4.5\n";
  stream << "treepm_enable_window_deconvolution = false\n";
  stream << "treepm_update_cadence_steps = 1\n\n";
  stream << "[output]\n";
  stream << "run_name = reference_workflow_distributed_isolated_pm_mpi\n";
  stream << "output_directory = integration_outputs\n";
  stream << "output_stem = snapshot\n";
  stream << "restart_stem = restart\n\n";
  stream << "[parallel]\n";
  stream << "mpi_ranks_expected = " << mpi_ranks_expected << "\n";
  stream << "omp_threads = 1\n";
  stream << "gpu_devices = 0\n";
  stream << "isolated_pm_root_workspace_limit_bytes = " << isolated_limit_bytes << "\n";
  return stream.str();
}
}

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
  int world_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size != 2) {
    MPI_Finalize();
    return 0;
  }
  bool guard_threw = false;
  std::string guard_diagnostic = "no exception";
  try {
    const auto guarded = cosmosim::core::loadFrozenConfigFromString(
        buildConfigText(world_size, 1024ULL),
        "test_reference_workflow_distributed_isolated_pm_mpi_guard");
    assert(guarded.config.parallel.isolated_pm_root_workspace_limit_bytes == 1024ULL);
    cosmosim::workflows::ReferenceWorkflowRunner guarded_runner(guarded);
    (void)guarded_runner.run(cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = false});
  } catch (const std::exception& err) {
    guard_diagnostic = err.what();
    guard_threw = guard_diagnostic.find("isolated/open PM root-gather guard exceeded") !=
        std::string::npos;
  }
  if (!guard_threw) {
    throw std::runtime_error(
        "distributed isolated PM workspace guard did not produce its contract diagnostic: " +
        guard_diagnostic);
  }

  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      buildConfigText(world_size, 64ULL * 1024ULL * 1024ULL),
      "test_reference_workflow_distributed_isolated_pm_mpi");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  const auto report = runner.run(cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = false});
  assert(report.completed_steps == 2);
  assert(report.world_size == world_size);
  assert(report.treepm_pm_grid_nx == 12);
  assert(report.treepm_long_range_refresh_count >= 2);
  MPI_Finalize();
#endif
  return 0;
}
