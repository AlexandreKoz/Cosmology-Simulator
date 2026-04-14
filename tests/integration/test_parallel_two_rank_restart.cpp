#include <cassert>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

#include "cosmosim/core/profiling.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace {

struct MiniEntity {
  std::uint64_t id = 0;
  double x_comov = 0.0;
  double y_comov = 0.0;
  double z_comov = 0.0;
  double mass_code = 0.0;
};

[[nodiscard]] std::vector<cosmosim::parallel::DecompositionItem> makeItems(std::span<const MiniEntity> entities) {
  std::vector<cosmosim::parallel::DecompositionItem> items;
  items.reserve(entities.size());
  for (const auto& entity : entities) {
    items.push_back(cosmosim::parallel::DecompositionItem{
        .entity_id = entity.id,
        .kind = cosmosim::parallel::DecompositionEntityKind::kParticle,
        .x_comov = entity.x_comov,
        .y_comov = entity.y_comov,
        .z_comov = entity.z_comov,
        .work_units = 1.0,
        .memory_bytes = 64,
    });
  }
  return items;
}

[[nodiscard]] double rankLocalMassSum(
    std::span<const MiniEntity> entities,
    std::span<const int> owning_rank_by_item,
    int world_rank) {
  double local_sum = 0.0;
  for (std::size_t i = 0; i < entities.size(); ++i) {
    if (owning_rank_by_item[i] == world_rank) {
      local_sum += entities[i].mass_code;
    }
  }
  return local_sum;
}

void runOneRankVsTwoRankMassConservation() {
  std::vector<MiniEntity> entities(64);
  for (std::size_t i = 0; i < entities.size(); ++i) {
    entities[i].id = static_cast<std::uint64_t>(i);
    entities[i].x_comov = static_cast<double>(i % 8) / 8.0;
    entities[i].y_comov = static_cast<double>((i / 8) % 4) / 4.0;
    entities[i].z_comov = static_cast<double>(i % 5) / 5.0;
    entities[i].mass_code = 1.0 + static_cast<double>(i % 3) * 0.25;
  }

  const auto items = makeItems(entities);

  cosmosim::parallel::DecompositionConfig one_rank_cfg;
  one_rank_cfg.world_size = 1;
  const auto one_rank_plan = cosmosim::parallel::buildMortonSfcDecomposition(items, one_rank_cfg);

  cosmosim::parallel::DecompositionConfig two_rank_cfg;
  two_rank_cfg.world_size = 2;
  two_rank_cfg.work_weight = 1.0;
  two_rank_cfg.memory_weight = 1.0 / 256.0;
  const auto two_rank_plan = cosmosim::parallel::buildMortonSfcDecomposition(items, two_rank_cfg);

  const double one_rank_sum = rankLocalMassSum(entities, one_rank_plan.owning_rank_by_item, 0);
  const double two_rank_sum = rankLocalMassSum(entities, two_rank_plan.owning_rank_by_item, 0) +
                              rankLocalMassSum(entities, two_rank_plan.owning_rank_by_item, 1);

  assert(std::abs(one_rank_sum - two_rank_sum) < 1.0e-12);
  assert(two_rank_plan.metrics.weighted_imbalance_ratio < 1.35);

  // Ghost exchange scaffolding in the two-rank split.
  std::vector<int> ghost_owners(entities.size(), 0);
  for (std::size_t i = 0; i < ghost_owners.size(); ++i) {
    ghost_owners[i] = (i % 2 == 0) ? 0 : 1;
  }
  const auto plan_rank0 = cosmosim::parallel::buildGhostExchangePlan(0, ghost_owners, sizeof(double) * 3 + sizeof(std::uint64_t));
  std::vector<cosmosim::parallel::LocalGhostDescriptor> explicit_rank0_descriptors;
  explicit_rank0_descriptors.reserve(ghost_owners.size());
  for (const int owner_rank : ghost_owners) {
    explicit_rank0_descriptors.push_back(cosmosim::parallel::LocalGhostDescriptor{
        .residency = (owner_rank == 0) ? cosmosim::parallel::LocalIndexResidency::kOwned
                                       : cosmosim::parallel::LocalIndexResidency::kGhost,
        .owning_rank = owner_rank,
    });
  }
  const auto typed_plan_rank0 = cosmosim::parallel::buildGhostExchangePlan(
      0,
      explicit_rank0_descriptors,
      sizeof(double) * 3 + sizeof(std::uint64_t));
  assert(typed_plan_rank0.neighbor_ranks == plan_rank0.neighbor_ranks);
  assert(typed_plan_rank0.recv_local_indices_by_neighbor == plan_rank0.recv_local_indices_by_neighbor);
  assert(typed_plan_rank0.send_bytes == plan_rank0.send_bytes);
  assert(typed_plan_rank0.recv_bytes == plan_rank0.recv_bytes);

  cosmosim::core::ProfilerSession profiler(true);
  cosmosim::parallel::recordDistributedProfiling(
      &profiler,
      two_rank_plan.metrics,
      plan_rank0.send_bytes,
      plan_rank0.recv_bytes);
  assert(profiler.counters().count("parallel.ghost_exchange_send_bytes") == plan_rank0.send_bytes);

  const std::array<double, 2> rank_local_masses = {
      rankLocalMassSum(entities, two_rank_plan.owning_rank_by_item, 0),
      rankLocalMassSum(entities, two_rank_plan.owning_rank_by_item, 1),
  };
  const double deterministic_mass_sum =
      cosmosim::parallel::deterministicRankOrderedSum(rank_local_masses);
  const auto reduction_agreement = cosmosim::parallel::compareReductionAgreement(
      rank_local_masses,
      two_rank_sum);
  assert(std::abs(deterministic_mass_sum - two_rank_sum) < 1.0e-12);
  assert(reduction_agreement.absolute_error < 1.0e-12);

  const std::vector<cosmosim::parallel::RankConfigDigest> rank_digests = {
      {.world_rank = 0, .normalized_config_hash = 0x52abU, .mpi_ranks_expected = 2, .deterministic_reduction = true},
      {.world_rank = 1, .normalized_config_hash = 0x52abU, .mpi_ranks_expected = 2, .deterministic_reduction = true},
  };
  const auto consensus = cosmosim::parallel::evaluateRankConfigConsensus(rank_digests);
  assert(consensus.allConsistent());
}

void runRestartRoundTripWithTwoRankPlan() {
  std::vector<MiniEntity> entities(24);
  for (std::size_t i = 0; i < entities.size(); ++i) {
    entities[i].id = static_cast<std::uint64_t>(10 + i);
    entities[i].x_comov = static_cast<double>(i) / static_cast<double>(entities.size());
    entities[i].mass_code = 1.0;
  }

  const auto items = makeItems(entities);

  cosmosim::parallel::DecompositionConfig cfg;
  cfg.world_size = 2;
  const auto plan = cosmosim::parallel::buildMortonSfcDecomposition(items, cfg);

  cosmosim::parallel::DistributedRestartState state;
  state.schema_version = 1;
  state.decomposition_epoch = 12;
  state.world_size = 2;
  state.owning_rank_by_item = plan.owning_rank_by_item;

  const std::string encoded = state.serialize();
  const auto restored = cosmosim::parallel::DistributedRestartState::deserialize(encoded);

  assert(restored.schema_version == state.schema_version);
  assert(restored.decomposition_epoch == state.decomposition_epoch);
  assert(restored.world_size == state.world_size);
  assert(restored.owning_rank_by_item == state.owning_rank_by_item);
}

}  // namespace

int main() {
  runOneRankVsTwoRankMassConservation();
  runRestartRoundTripWithTwoRankPlan();
  return 0;
}
