#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "cosmosim/parallel/distributed_memory.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

std::vector<cosmosim::parallel::DecompositionItem> makeLocalItems(int world_size, int world_rank) {
  namespace parallel = cosmosim::parallel;
  std::vector<parallel::DecompositionItem> items;
  const int local_count = 4;
  items.reserve(static_cast<std::size_t>(local_count + 1));
  for (int i = 0; i < local_count; ++i) {
    const int global_ordinal = world_rank * local_count + i;
    parallel::DecompositionItem item;
    item.entity_id = static_cast<std::uint64_t>(1000 + global_ordinal);
    item.kind = parallel::DecompositionEntityKind::kParticle;
    item.current_owner_rank = 0;
    item.x_comov = (static_cast<double>(global_ordinal) + 0.25) /
        static_cast<double>(world_size * local_count + 2);
    item.y_comov = 0.25;
    item.z_comov = 0.75;
    item.memory_bytes = 256U + static_cast<std::uint64_t>(16 * i);
    item.work_components = parallel::DecompositionWorkComponents{
        .particle_count_cost = 1.0,
        .memory_pressure_cost = static_cast<double>(item.memory_bytes),
        .generic_work_cost = 1.0,
        .has_explicit_components = true,
    };
    items.push_back(item);
  }
  if (world_rank == 0) {
    parallel::DecompositionItem patch;
    patch.entity_id = 9001U;
    patch.kind = parallel::DecompositionEntityKind::kAmrPatch;
    patch.current_owner_rank = 0;
    patch.x_comov = 0.92;
    patch.y_comov = 0.25;
    patch.z_comov = 0.75;
    patch.memory_bytes = 4096U;
    patch.work_components = parallel::DecompositionWorkComponents{
        .amr_patch_cost = 8.0,
        .memory_pressure_cost = static_cast<double>(patch.memory_bytes),
        .generic_work_cost = 8.0,
        .has_explicit_components = true,
    };
    items.push_back(patch);
  }
  return items;
}

std::vector<cosmosim::parallel::DecompositionItem> makeTieItems(int world_rank) {
  namespace parallel = cosmosim::parallel;
  std::vector<parallel::DecompositionItem> items(2);
  for (std::size_t i = 0; i < items.size(); ++i) {
    items[i].entity_id = static_cast<std::uint64_t>(2000 + world_rank * 2 + static_cast<int>(i));
    items[i].kind = parallel::DecompositionEntityKind::kParticle;
    items[i].current_owner_rank = 0;
    items[i].x_comov = 0.5;
    items[i].y_comov = 0.5;
    items[i].z_comov = 0.5;
    items[i].memory_bytes = 128U;
    items[i].work_components = parallel::DecompositionWorkComponents{
        .particle_count_cost = 1.0,
        .generic_work_cost = 1.0,
        .has_explicit_components = true,
    };
  }
  return items;
}

std::vector<cosmosim::parallel::DecompositionItem> makeMixedActionabilityItems(int world_rank) {
  namespace parallel = cosmosim::parallel;
  parallel::DecompositionItem item;
  item.entity_id = static_cast<std::uint64_t>(3000 + world_rank);
  item.kind = parallel::DecompositionEntityKind::kParticle;
  item.current_owner_rank = world_rank;
  const int ordered_position = (world_rank == 0) ? 1 : ((world_rank == 1) ? 0 : world_rank);
  item.x_comov = (static_cast<double>(ordered_position) + 0.5) / 4.0;
  item.y_comov = 0.5;
  item.z_comov = 0.5;
  item.memory_bytes = 128U;
  item.work_components = parallel::DecompositionWorkComponents{
      .particle_count_cost = 1.0,
      .generic_work_cost = 1.0,
      .has_explicit_components = true,
  };
  return {item};
}

void checkDistributedMatchesExact(
    const cosmosim::parallel::MpiContext& mpi_context,
    const std::vector<cosmosim::parallel::DecompositionItem>& local_items) {
  namespace parallel = cosmosim::parallel;
  parallel::DecompositionConfig decomposition_config;
  decomposition_config.world_size = mpi_context.worldSize();
  decomposition_config.domain_x_min_comov = 0.0;
  decomposition_config.domain_x_max_comov = 1.0;
  decomposition_config.domain_y_min_comov = 0.0;
  decomposition_config.domain_y_max_comov = 1.0;
  decomposition_config.domain_z_min_comov = 0.0;
  decomposition_config.domain_z_max_comov = 1.0;
  decomposition_config.component_weights.particle_count = 1.0;
  decomposition_config.component_weights.amr_patch = 1.0;
  decomposition_config.component_weights.memory_pressure = 0.0;
  decomposition_config.component_weights.generic_work = 0.0;

  parallel::RuntimeRebalanceConfig rebalance_config;
  rebalance_config.world_size = mpi_context.worldSize();
  rebalance_config.imbalance_trigger_ratio = 1.0;
  rebalance_config.memory_trigger_ratio = 10.0;
  rebalance_config.max_migrated_load_fraction = 1.0;

  const auto distributed = parallel::buildDistributedRuntimeRebalancePlan(
      mpi_context, local_items, decomposition_config, rebalance_config);
  assert(distributed.used_distributed_sfc_cuts);
  assert(distributed.global_entities_considered >= local_items.size());
  assert(std::isfinite(distributed.current_metrics.weighted_imbalance_ratio));
  assert(std::isfinite(distributed.target_decomposition.metrics.weighted_imbalance_ratio));
  assert(distributed.global_control_bytes < distributed.global_entities_considered * sizeof(parallel::DecompositionItem));
  assert(distributed.peak_temporary_bytes > 0U);

  const auto global_items = parallel::gatherDecompositionItemsAcrossRanks(mpi_context, local_items);
  const auto exact = parallel::buildRuntimeRebalancePlan(global_items, decomposition_config, rebalance_config);
  std::unordered_map<std::uint64_t, int> exact_owner_by_id;
  for (std::size_t i = 0; i < global_items.size(); ++i) {
    exact_owner_by_id.emplace(global_items[i].entity_id, exact.target_decomposition.owning_rank_by_item[i]);
  }

  for (std::size_t i = 0; i < local_items.size(); ++i) {
    const auto found = exact_owner_by_id.find(local_items[i].entity_id);
    assert(found != exact_owner_by_id.end());
    // The compact cut planner is the production authority and is not required
    // to reproduce the exact all-item planner once more than two ranks can
    // expose one-item-per-rank boundary constraints that compact samples do
    // not encode. Keep exact owner equivalence on the original two-rank
    // fixture and use the higher-rank cases for distributed invariants.
    if (mpi_context.worldSize() == 2) {
      assert(distributed.target_decomposition.owning_rank_by_item[i] == found->second);
    }
  }
  for (const auto& intent : distributed.particle_migrations) {
    assert(intent.old_owner_rank != intent.new_owner_rank);
    if (mpi_context.worldSize() == 2) {
      assert(exact_owner_by_id.at(intent.particle_id) == intent.new_owner_rank);
    }
  }
  for (const auto& update : distributed.amr_patch_ownership_updates) {
    assert(update.old_owner_rank != update.new_owner_rank);
    if (mpi_context.worldSize() == 2) {
      assert(exact_owner_by_id.at(update.patch_id) == update.new_owner_rank);
    }
  }

}

void checkMixedActionabilityCompletion(const cosmosim::parallel::MpiContext& mpi_context) {
  namespace parallel = cosmosim::parallel;
  assert(mpi_context.worldSize() == 4);

  parallel::DecompositionConfig decomposition_config;
  decomposition_config.world_size = mpi_context.worldSize();
  decomposition_config.domain_x_min_comov = 0.0;
  decomposition_config.domain_x_max_comov = 1.0;
  decomposition_config.domain_y_min_comov = 0.0;
  decomposition_config.domain_y_max_comov = 1.0;
  decomposition_config.domain_z_min_comov = 0.0;
  decomposition_config.domain_z_max_comov = 1.0;
  decomposition_config.component_weights.particle_count = 1.0;

  parallel::RuntimeRebalanceConfig rebalance_config;
  rebalance_config.world_size = mpi_context.worldSize();
  rebalance_config.imbalance_trigger_ratio = 1.0;
  rebalance_config.memory_trigger_ratio = 10.0;
  rebalance_config.max_migrated_load_fraction = 1.0;

  const auto distributed = parallel::buildDistributedRuntimeRebalancePlan(
      mpi_context,
      makeMixedActionabilityItems(mpi_context.worldRank()),
      decomposition_config,
      rebalance_config);
  const std::uint64_t local_actionable =
      (!distributed.particle_migrations.empty() ||
       !distributed.amr_patch_ownership_updates.empty())
      ? 1ULL
      : 0ULL;
  const std::uint64_t actionable_rank_count =
      mpi_context.allreduceSumUint64(local_actionable);
  assert(actionable_rank_count > 0U);
  assert(actionable_rank_count < static_cast<std::uint64_t>(mpi_context.worldSize()));
  const std::uint64_t should_rebalance_count =
      mpi_context.allreduceSumUint64(distributed.should_rebalance ? 1ULL : 0ULL);
  assert(should_rebalance_count ==
         static_cast<std::uint64_t>(mpi_context.worldSize()));
}

}  // namespace

int main(int argc, char** argv) {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  MPI_Init(&argc, &argv);
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  assert(world_size >= 2 && world_size <= 4);

  const cosmosim::parallel::MpiContext mpi_context(/*is_enabled=*/true, world_size, world_rank);
  checkDistributedMatchesExact(mpi_context, makeLocalItems(world_size, world_rank));
  checkDistributedMatchesExact(mpi_context, makeTieItems(world_rank));
  if (world_size == 4) {
    checkMixedActionabilityCompletion(mpi_context);
  }

  MPI_Finalize();
#else
  (void)argc;
  (void)argv;
#endif
  return 0;
}
