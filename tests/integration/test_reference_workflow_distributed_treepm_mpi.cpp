#include <cassert>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include "cosmosim/cosmosim.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

[[nodiscard]] std::uint64_t xorRangeOneToN(std::uint64_t n) {
  switch (n & 3ULL) {
    case 0ULL: return n;
    case 1ULL: return 1ULL;
    case 2ULL: return n + 1ULL;
    default: return 0ULL;
  }
}

std::vector<cosmosim::core::TimeBinSchedulerIdentityRecord> schedulerSeedForRank(int world_rank) {
  std::vector<cosmosim::core::TimeBinSchedulerIdentityRecord> records;
  records.reserve(42U);
  for (std::uint64_t particle_id = 1ULL; particle_id <= 42ULL; ++particle_id) {
    records.push_back(cosmosim::core::TimeBinSchedulerIdentityRecord{
        .element_id = particle_id,
        .bin_index = static_cast<std::uint8_t>(world_rank == 0 ? 0U : 1U),
        .next_activation_tick = world_rank == 0 ? 0ULL : 2ULL,
        .pending_bin_index = cosmosim::core::HierarchicalTimeBinScheduler::k_unset_pending_bin,
    });
  }
  return records;
}

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
  stream << "decomposition_runtime_rebalance_enabled = true\n";
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
  cosmosim::workflows::ReferenceWorkflowOptions options;
  options.write_outputs = false;
  options.initial_particle_scheduler_identity_records = schedulerSeedForRank(world_rank);
  const cosmosim::workflows::ReferenceWorkflowReport report = runner.run(options);

  assert(report.completed_steps == 2);
  assert(report.world_size == world_size);
  assert(report.world_rank == world_rank);
  assert(report.global_particle_count == 42ULL);
  assert(report.global_cell_count == 6ULL);
  assert(report.global_particle_id_sum == (42ULL * 43ULL) / 2ULL);
  assert(report.global_particle_id_xor == xorRangeOneToN(42ULL));
  assert(report.local_particle_ids_unique);
  assert(report.global_particle_partition_identity_match);
  assert(report.treepm_update_cadence_steps == 2);
  assert(report.treepm_long_range_refresh_count == 2);
  assert(report.treepm_long_range_reuse_count == 1);
  assert(report.treepm_cadence_records.size() == 3U);
  const std::vector<std::uint64_t> expected_field_versions{1ULL, 1ULL, 2ULL};
  const std::vector<std::uint64_t> expected_field_built_steps{0ULL, 0ULL, 1ULL};
  const std::vector<bool> expected_refresh_flags{true, false, true};
  const std::vector<std::string> expected_stage_names{
      "gravity_kick_pre", "force_refresh", "force_refresh"};
  for (std::size_t i = 0; i < report.treepm_cadence_records.size(); ++i) {
    const auto& record = report.treepm_cadence_records[i];
    assert(record.gravity_kick_opportunity == i + 1U);
    assert(record.stage_name == expected_stage_names[i]);
    assert(record.field_version == expected_field_versions[i]);
    assert(record.field_built_step_index == expected_field_built_steps[i]);
    assert(record.refreshed_long_range_field == expected_refresh_flags[i]);
    assert(record.active_particles_kicked + record.inactive_particles_skipped == report.local_particle_count);
  }
  assert(frozen.config.parallel.decomposition_runtime_rebalance_enabled);


  std::vector<std::uint64_t> gathered_opportunity(report.treepm_cadence_records.size(), 0ULL);
  std::vector<std::uint64_t> gathered_field_version(report.treepm_cadence_records.size(), 0ULL);
  std::vector<std::uint64_t> gathered_last_refresh_opportunity(report.treepm_cadence_records.size(), 0ULL);
  std::vector<std::uint64_t> gathered_refresh_flag(report.treepm_cadence_records.size(), 0ULL);
  for (std::size_t i = 0; i < report.treepm_cadence_records.size(); ++i) {
    gathered_opportunity[i] = report.treepm_cadence_records[i].gravity_kick_opportunity;
    gathered_field_version[i] = report.treepm_cadence_records[i].field_version;
    gathered_last_refresh_opportunity[i] = report.treepm_cadence_records[i].last_refresh_opportunity;
    gathered_refresh_flag[i] = report.treepm_cadence_records[i].refreshed_long_range_field ? 1ULL : 0ULL;
  }

  std::vector<std::uint64_t> reduced_opportunity(gathered_opportunity.size(), 0ULL);
  std::vector<std::uint64_t> reduced_field_version(gathered_field_version.size(), 0ULL);
  std::vector<std::uint64_t> reduced_last_refresh_opportunity(gathered_last_refresh_opportunity.size(), 0ULL);
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
      gathered_last_refresh_opportunity.data(),
      reduced_last_refresh_opportunity.data(),
      static_cast<int>(gathered_last_refresh_opportunity.size()),
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
    assert(
        reduced_last_refresh_opportunity[i] ==
        gathered_last_refresh_opportunity[i] * static_cast<std::uint64_t>(world_size));
    assert(reduced_refresh_flag[i] == gathered_refresh_flag[i] * static_cast<std::uint64_t>(world_size));
  }

  std::vector<std::uint64_t> local_active_by_boundary(report.treepm_cadence_records.size(), 0ULL);
  for (std::size_t i = 0; i < report.treepm_cadence_records.size(); ++i) {
    local_active_by_boundary[i] = report.treepm_cadence_records[i].active_particles_kicked;
  }
  std::vector<std::uint64_t> active_by_rank_and_boundary(
      static_cast<std::size_t>(world_size) * local_active_by_boundary.size(), 0ULL);
  MPI_Allgather(
      local_active_by_boundary.data(),
      static_cast<int>(local_active_by_boundary.size()),
      MPI_UINT64_T,
      active_by_rank_and_boundary.data(),
      static_cast<int>(local_active_by_boundary.size()),
      MPI_UINT64_T,
      MPI_COMM_WORLD);
  bool saw_zero_local_work_with_remote_activity = false;
  for (std::size_t boundary = 0; boundary < local_active_by_boundary.size(); ++boundary) {
    const std::uint64_t rank_zero_active = active_by_rank_and_boundary[boundary];
    const std::uint64_t rank_one_active =
        active_by_rank_and_boundary[local_active_by_boundary.size() + boundary];
    if (rank_zero_active > 0ULL && rank_one_active == 0ULL) {
      saw_zero_local_work_with_remote_activity = true;
    }
  }
  assert(saw_zero_local_work_with_remote_activity);

  cosmosim::parallel::DistributedRestartState invalid_restart_state;
  invalid_restart_state.world_size = world_size;
  invalid_restart_state.pm_grid_nx = 16;
  invalid_restart_state.pm_grid_ny = 16;
  invalid_restart_state.pm_grid_nz = 16;
  invalid_restart_state.pm_decomposition_mode = "slab";
  invalid_restart_state.pm_update_cadence_steps = 0;
  invalid_restart_state.gravity_kick_opportunity = 1;
  invalid_restart_state.last_long_range_refresh_opportunity = 2;
  invalid_restart_state.long_range_field_version = 0;
  invalid_restart_state.pm_slab_begin_x_by_rank.resize(static_cast<std::size_t>(world_size), 0);
  invalid_restart_state.pm_slab_end_x_by_rank.resize(static_cast<std::size_t>(world_size), 0);
  for (int rank = 0; rank < world_size; ++rank) {
    const auto slab = cosmosim::parallel::pmOwnedXRangeForRank(16, world_size, rank);
    invalid_restart_state.pm_slab_begin_x_by_rank[static_cast<std::size_t>(rank)] = slab.begin_x;
    invalid_restart_state.pm_slab_end_x_by_rank[static_cast<std::size_t>(rank)] = slab.end_x;
  }
  const cosmosim::parallel::MpiContext runtime_context(
      /*is_enabled=*/world_size > 1,
      world_size,
      world_rank);
  const auto compatibility = cosmosim::parallel::evaluateDistributedRestartCompatibility(
      invalid_restart_state,
      cosmosim::parallel::buildDistributedExecutionTopology(
          /*global_nx=*/16,
          /*global_ny=*/16,
          /*global_nz=*/16,
          runtime_context,
          /*mpi_ranks_expected=*/world_size,
          /*configured_gpu_devices=*/0,
          /*cuda_runtime_available=*/false,
          /*visible_device_count=*/0,
          /*pm_decomposition_mode=*/"slab"));
  assert(!compatibility.compatible());
  assert(!compatibility.pm_cadence_steps_match);
  assert(!compatibility.gravity_kick_state_match);
  assert(!compatibility.long_range_field_state_match);

  const std::uint64_t local_digest = report.final_state_digest;
  std::vector<std::uint64_t> gathered_digest(static_cast<std::size_t>(world_size), 0ULL);
  MPI_Allgather(&local_digest, 1, MPI_UINT64_T, gathered_digest.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
  for (const std::uint64_t digest : gathered_digest) {
    assert(digest != 0ULL);
  }

  const std::uint64_t local_particle_count = report.local_particle_count;
  std::vector<std::uint64_t> gathered_particle_count(static_cast<std::size_t>(world_size), 0ULL);
  MPI_Allgather(&local_particle_count, 1, MPI_UINT64_T, gathered_particle_count.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
  std::uint64_t reduced_particle_count = 0ULL;
  MPI_Allreduce(&local_particle_count, &reduced_particle_count, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  assert(reduced_particle_count == report.global_particle_count);
  bool all_ranks_hold_full_state = true;
  for (const std::uint64_t count : gathered_particle_count) {
    all_ranks_hold_full_state = all_ranks_hold_full_state && (count == report.global_particle_count);
  }
  assert(!all_ranks_hold_full_state);

  const std::uint64_t local_cell_count = report.local_cell_count;
  std::uint64_t reduced_cell_count = 0ULL;
  MPI_Allreduce(&local_cell_count, &reduced_cell_count, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  assert(reduced_cell_count == report.global_cell_count);

  const std::uint64_t local_id_sum = report.local_particle_id_sum;
  std::uint64_t reduced_id_sum = 0ULL;
  MPI_Allreduce(&local_id_sum, &reduced_id_sum, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  assert(reduced_id_sum == report.global_particle_id_sum);

  const std::uint64_t local_id_xor = report.local_particle_id_xor;
  std::uint64_t reduced_id_xor = 0ULL;
  MPI_Allreduce(&local_id_xor, &reduced_id_xor, 1, MPI_UINT64_T, MPI_BXOR, MPI_COMM_WORLD);
  assert(reduced_id_xor == report.global_particle_id_xor);

  const cosmosim::parallel::LocalOwnershipIdentitySummary reduced_identity{
      .local_owned_count = reduced_particle_count,
      .local_particle_id_sum = reduced_id_sum,
      .local_particle_id_xor = reduced_id_xor,
      .local_particle_ids_unique = report.local_particle_ids_unique,
  };
  assert(cosmosim::parallel::partitionIdentityMatchesGeneratedSet(
      reduced_identity,
      /*expected_global_count=*/42ULL,
      /*expected_particle_id_sum=*/(42ULL * 43ULL) / 2ULL,
      /*expected_particle_id_xor=*/xorRangeOneToN(42ULL)));

  MPI_Finalize();
#endif
  return 0;
}
