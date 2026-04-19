#include <cassert>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include "cosmosim/cosmosim.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

std::string buildConfigText(int cadence_steps, int mpi_ranks_expected) {
  std::stringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[mode]\n";
  stream << "mode = zoom_in\n";
  stream << "ic_file = generated\n";
  stream << "zoom_high_res_region = false\n\n";
  stream << "[numerics]\n";
  stream << "time_begin_code = 0.01\n";
  stream << "time_end_code = 0.0102\n";
  stream << "max_global_steps = 2\n";
  stream << "hierarchical_max_rung = 2\n";
  stream << "treepm_pm_grid = 16\n";
  stream << "treepm_asmth_cells = 1.5\n";
  stream << "treepm_rcut_cells = 4.5\n";
  stream << "treepm_update_cadence_steps = " << cadence_steps << "\n";
  stream << "treepm_tree_exchange_batch_bytes = 256\n\n";
  stream << "[output]\n";
  stream << "run_name = reference_workflow_distributed_treepm_mpi\n";
  stream << "output_directory = integration_outputs\n";
  stream << "output_stem = snapshot\n";
  stream << "restart_stem = restart\n\n";
  stream << "[parallel]\n";
  stream << "mpi_ranks_expected = " << mpi_ranks_expected << "\n";
  stream << "omp_threads = 1\n";
  stream << "gpu_devices = 0\n";
  return stream.str();
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    MPI_Finalize();
    return 0;
  }

  const std::string config_text = buildConfigText(/*cadence_steps=*/2, world_size);
  const cosmosim::core::FrozenConfig frozen =
      cosmosim::core::loadFrozenConfigFromString(config_text, "test_reference_workflow_distributed_treepm_mpi");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  const cosmosim::workflows::ReferenceWorkflowReport report =
      runner.run(cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = false});

  assert(report.completed_steps == 2);
  assert(report.treepm_update_cadence_steps == 2);
  assert(report.treepm_long_range_refresh_count == 2);
  assert(report.treepm_long_range_reuse_count == 2);
  assert(report.treepm_cadence_records.size() == 4);

  for (std::size_t i = 0; i < report.treepm_cadence_records.size(); ++i) {
    const auto& record = report.treepm_cadence_records[i];
    assert(record.gravity_kick_opportunity == i + 1U);
  }

  std::vector<std::uint64_t> gathered_opportunity(report.treepm_cadence_records.size(), 0ULL);
  std::vector<std::uint64_t> gathered_field_version(report.treepm_cadence_records.size(), 0ULL);
  std::vector<std::uint64_t> gathered_refresh_flag(report.treepm_cadence_records.size(), 0ULL);
  for (std::size_t i = 0; i < report.treepm_cadence_records.size(); ++i) {
    gathered_opportunity[i] = report.treepm_cadence_records[i].gravity_kick_opportunity;
    gathered_field_version[i] = report.treepm_cadence_records[i].field_version;
    gathered_refresh_flag[i] = report.treepm_cadence_records[i].refreshed_long_range_field ? 1ULL : 0ULL;
  }

  std::vector<std::uint64_t> reduced_opportunity(gathered_opportunity.size(), 0ULL);
  std::vector<std::uint64_t> reduced_field_version(gathered_field_version.size(), 0ULL);
  std::vector<std::uint64_t> reduced_refresh_flag(gathered_refresh_flag.size(), 0ULL);
  MPI_Allreduce(
      gathered_opportunity.data(),
      reduced_opportunity.data(),
      static_cast<int>(gathered_opportunity.size()),
      MPI_UINT64_T,
      MPI_SUM,
      MPI_COMM_WORLD);
  MPI_Allreduce(
      gathered_field_version.data(),
      reduced_field_version.data(),
      static_cast<int>(gathered_field_version.size()),
      MPI_UINT64_T,
      MPI_SUM,
      MPI_COMM_WORLD);
  MPI_Allreduce(
      gathered_refresh_flag.data(),
      reduced_refresh_flag.data(),
      static_cast<int>(gathered_refresh_flag.size()),
      MPI_UINT64_T,
      MPI_SUM,
      MPI_COMM_WORLD);

  for (std::size_t i = 0; i < report.treepm_cadence_records.size(); ++i) {
    assert(reduced_opportunity[i] == gathered_opportunity[i] * static_cast<std::uint64_t>(world_size));
    assert(reduced_field_version[i] == gathered_field_version[i] * static_cast<std::uint64_t>(world_size));
    assert(reduced_refresh_flag[i] == gathered_refresh_flag[i] * static_cast<std::uint64_t>(world_size));
  }

  const std::uint64_t local_digest = report.final_state_digest;
  std::uint64_t digest_xor = 0ULL;
  MPI_Allreduce(&local_digest, &digest_xor, 1, MPI_UINT64_T, MPI_BXOR, MPI_COMM_WORLD);
  assert(digest_xor == 0ULL);

  MPI_Finalize();
#endif
  return 0;
}
