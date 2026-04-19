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
  assert(typed_plan_rank0.outbound_transfers.size() == typed_plan_rank0.neighbor_ranks.size());
  assert(typed_plan_rank0.inbound_transfers.size() == typed_plan_rank0.neighbor_ranks.size());
  for (std::size_t i = 0; i < typed_plan_rank0.neighbor_ranks.size(); ++i) {
    assert(typed_plan_rank0.outbound_transfers[i].role ==
           cosmosim::parallel::GhostTransferRole::kOutboundSend);
    assert(typed_plan_rank0.inbound_transfers[i].role ==
           cosmosim::parallel::GhostTransferRole::kInboundReceive);
    assert(typed_plan_rank0.outbound_transfers[i].intent ==
           cosmosim::parallel::GhostTransferIntent::kGhostRefreshRequest);
    assert(typed_plan_rank0.inbound_transfers[i].intent ==
           cosmosim::parallel::GhostTransferIntent::kGhostRefreshReceiveStaging);
    assert(typed_plan_rank0.outbound_transfers[i].peer_rank == typed_plan_rank0.neighbor_ranks[i]);
    assert(typed_plan_rank0.inbound_transfers[i].peer_rank == typed_plan_rank0.neighbor_ranks[i]);
    assert(typed_plan_rank0.outbound_transfers[i].neighbor_slot == i);
    assert(typed_plan_rank0.inbound_transfers[i].neighbor_slot == i);
  }
  cosmosim::parallel::validateGhostExchangePlan(typed_plan_rank0);

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
  assert(std::abs(reduction_agreement.deterministic_baseline_sum - deterministic_mass_sum) < 1.0e-12);
  assert(std::abs(reduction_agreement.measured_sum - two_rank_sum) < 1.0e-12);
  assert(reduction_agreement.absolute_error < 1.0e-12);
  assert(cosmosim::parallel::satisfiesReductionAgreement(
      reduction_agreement,
      {.mode = cosmosim::parallel::ReductionAgreementMode::kAbsoluteOrRelative,
       .absolute_tolerance = 1.0e-12,
       .relative_tolerance = 1.0e-12}));

  const std::vector<cosmosim::parallel::RankConfigDigest> rank_digests = {
      {.world_rank = 0, .normalized_config_hash = 0x52abU, .mpi_ranks_expected = 2, .deterministic_reduction = true},
      {.world_rank = 1, .normalized_config_hash = 0x52abU, .mpi_ranks_expected = 2, .deterministic_reduction = true},
  };
  const auto consensus = cosmosim::parallel::evaluateRankConfigConsensus(rank_digests);
  assert(consensus.allConsistent());
  assert(consensus.mismatches.empty());

  const std::vector<cosmosim::parallel::RankConfigDigest> rank_digests_mismatch = {
      {.world_rank = 0, .normalized_config_hash = 0x52abU, .mpi_ranks_expected = 2, .deterministic_reduction = true},
      {.world_rank = 1, .normalized_config_hash = 0x52abU, .mpi_ranks_expected = 2, .deterministic_reduction = false},
  };
  const auto bad_consensus = cosmosim::parallel::evaluateRankConfigConsensus(rank_digests_mismatch);
  assert(!bad_consensus.allConsistent());
  assert(bad_consensus.mismatches.size() == 1);
  assert(bad_consensus.mismatches[0].property ==
         cosmosim::parallel::RankConfigMismatchProperty::kDeterministicReduction);
  assert(bad_consensus.mismatches[0].baseline_rank == 0);
  assert(bad_consensus.mismatches[0].rank == 1);
  assert(bad_consensus.mismatches[0].baseline_value == "true");
  assert(bad_consensus.mismatches[0].rank_value == "false");
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
  state.schema_version = 2;
  state.decomposition_epoch = 12;
  state.world_size = 2;
  state.pm_grid_nx = 16;
  state.pm_grid_ny = 16;
  state.pm_grid_nz = 16;
  state.pm_decomposition_mode = "slab";
  state.gravity_kick_opportunity = 7;
  state.pm_update_cadence_steps = 2;
  state.long_range_field_version = 3;
  state.last_long_range_refresh_opportunity = 6;
  state.long_range_field_built_step_index = 24;
  state.long_range_field_built_scale_factor = 0.5;
  state.long_range_restart_policy = "deterministic_rebuild";
  state.owning_rank_by_item = plan.owning_rank_by_item;
  state.pm_slab_begin_x_by_rank = {0, 8};
  state.pm_slab_end_x_by_rank = {8, 16};

  const std::string encoded = state.serialize();
  const auto restored = cosmosim::parallel::DistributedRestartState::deserialize(encoded);

  assert(restored.schema_version == state.schema_version);
  assert(restored.decomposition_epoch == state.decomposition_epoch);
  assert(restored.world_size == state.world_size);
  assert(restored.pm_grid_nx == state.pm_grid_nx);
  assert(restored.gravity_kick_opportunity == state.gravity_kick_opportunity);
  assert(restored.owning_rank_by_item == state.owning_rank_by_item);
  assert(restored.pm_slab_begin_x_by_rank == state.pm_slab_begin_x_by_rank);
}

}  // namespace

int main() {
  runOneRankVsTwoRankMassConservation();
  runRestartRoundTripWithTwoRankPlan();
  return 0;
}
