#include <cassert>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "cosmosim/parallel/distributed_memory.hpp"

namespace {

void testGhostPackUnpackRoundTrip() {
  cosmosim::parallel::GhostExchangeBufferSoA source;
  source.entity_id = {101, 202, 303, 404};
  source.density_code = {1.0, 2.0, 3.0, 4.0};
  source.velocity_x_code = {10.0, 20.0, 30.0, 40.0};
  source.pressure_code = {100.0, 200.0, 300.0, 400.0};

  const std::vector<std::uint32_t> packed_indices = {1, 3};

  cosmosim::parallel::GhostExchangeBuffer buffer;
  buffer.packFrom(source, packed_indices);
  assert(buffer.byteSize() > 0);

  cosmosim::parallel::GhostExchangeBufferSoA destination;
  buffer.unpackAppendTo(destination);

  assert(destination.entity_id.size() == 2);
  assert(destination.entity_id[0] == 202);
  assert(destination.entity_id[1] == 404);
  assert(std::abs(destination.density_code[0] - 2.0) < 1.0e-12);
  assert(std::abs(destination.velocity_x_code[1] - 40.0) < 1.0e-12);
  assert(std::abs(destination.pressure_code[1] - 400.0) < 1.0e-12);
}

void testMortonDecompositionInvariants() {
  std::vector<cosmosim::parallel::DecompositionItem> items(16);
  for (std::size_t i = 0; i < items.size(); ++i) {
    items[i].entity_id = static_cast<std::uint64_t>(i);
    items[i].x_comov = static_cast<double>(i % 4) / 4.0;
    items[i].y_comov = static_cast<double>((i / 4) % 2) / 2.0;
    items[i].z_comov = static_cast<double>(i % 3) / 3.0;
    items[i].work_units = 1.0 + static_cast<double>(i % 5);
    items[i].memory_bytes = 128U + static_cast<std::uint64_t>(i * 8);
    items[i].kind = (i % 2 == 0) ? cosmosim::parallel::DecompositionEntityKind::kParticle
                                 : cosmosim::parallel::DecompositionEntityKind::kAmrPatch;
  }

  cosmosim::parallel::DecompositionConfig config;
  config.world_size = 3;
  config.work_weight = 1.0;
  config.memory_weight = 1.0 / 1024.0;

  const cosmosim::parallel::DecompositionPlan plan =
      cosmosim::parallel::buildMortonSfcDecomposition(items, config);

  assert(plan.owning_rank_by_item.size() == items.size());
  assert(plan.sorted_indices.size() == items.size());
  assert(plan.ranges_by_rank.size() == 3);
  for (const int rank : plan.owning_rank_by_item) {
    assert(rank >= 0 && rank < config.world_size);
  }

  for (std::size_t rank = 0; rank < plan.ranges_by_rank.size(); ++rank) {
    const auto range = plan.ranges_by_rank[rank];
    assert(range.begin_sorted <= range.end_sorted);
    assert(range.end_sorted <= items.size());
    for (std::size_t s = range.begin_sorted; s < range.end_sorted; ++s) {
      assert(plan.owning_rank_by_item[plan.sorted_indices[s]] == static_cast<int>(rank));
    }
  }

  assert(plan.metrics.weighted_imbalance_ratio >= 1.0);
  assert(plan.metrics.memory_imbalance_ratio >= 1.0);
  assert(plan.metrics.total_memory_bytes > 0);
}

void testRestartStateRoundTrip() {
  cosmosim::parallel::DistributedRestartState in;
  in.schema_version = 1;
  in.decomposition_epoch = 7;
  in.world_size = 2;
  in.owning_rank_by_item = {0, 1, 1, 0};

  const std::string encoded = in.serialize();
  const auto out = cosmosim::parallel::DistributedRestartState::deserialize(encoded);

  assert(out.schema_version == in.schema_version);
  assert(out.decomposition_epoch == in.decomposition_epoch);
  assert(out.world_size == in.world_size);
  assert(out.owning_rank_by_item == in.owning_rank_by_item);
}

void testExplicitOwnedVsGhostContracts() {
  const std::vector<cosmosim::parallel::LocalGhostDescriptor> descriptors = {
      {.residency = cosmosim::parallel::LocalIndexResidency::kOwned, .owning_rank = 0},
      {.residency = cosmosim::parallel::LocalIndexResidency::kGhost, .owning_rank = 1},
      {.residency = cosmosim::parallel::LocalIndexResidency::kOwned, .owning_rank = 0},
      {.residency = cosmosim::parallel::LocalIndexResidency::kGhost, .owning_rank = 2},
  };

  const auto plan = cosmosim::parallel::buildGhostExchangePlan(0, descriptors, sizeof(double) * 3 + sizeof(std::uint64_t));
  assert(plan.neighbor_ranks.size() == 2);
  assert(plan.neighbor_ranks[0] == 1);
  assert(plan.neighbor_ranks[1] == 2);
  assert(plan.recv_local_indices_by_neighbor[0].size() == 1);
  assert(plan.recv_local_indices_by_neighbor[0][0] == 1);
  assert(plan.recv_local_indices_by_neighbor[1].size() == 1);
  assert(plan.recv_local_indices_by_neighbor[1][0] == 3);
  assert(plan.send_local_indices_by_neighbor == plan.recv_local_indices_by_neighbor);
  assert(plan.outbound_transfers.size() == 2);
  assert(plan.inbound_transfers.size() == 2);
  assert(plan.outbound_transfers[0].role == cosmosim::parallel::GhostTransferRole::kOutboundSend);
  assert(plan.inbound_transfers[0].role == cosmosim::parallel::GhostTransferRole::kInboundReceive);
  assert(plan.outbound_transfers[0].peer_rank == 1);
  assert(plan.inbound_transfers[1].peer_rank == 2);
  assert(plan.outbound_transfers[1].local_indices == plan.send_local_indices_by_neighbor[1]);
  assert(plan.inbound_transfers[1].local_indices == plan.recv_local_indices_by_neighbor[1]);
}

void testGhostTransferInvariantFailures() {
  {
    bool threw = false;
    try {
      const std::vector<cosmosim::parallel::LocalGhostDescriptor> descriptors = {
          {.residency = cosmosim::parallel::LocalIndexResidency::kOwned, .owning_rank = 1},
      };
      static_cast<void>(cosmosim::parallel::buildGhostExchangePlan(0, descriptors, sizeof(double)));
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }

  {
    bool threw = false;
    try {
      const std::vector<cosmosim::parallel::LocalGhostDescriptor> descriptors = {
          {.residency = cosmosim::parallel::LocalIndexResidency::kGhost, .owning_rank = 0},
      };
      static_cast<void>(cosmosim::parallel::buildGhostExchangePlan(0, descriptors, sizeof(double)));
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }

  {
    bool threw = false;
    try {
      const std::vector<cosmosim::parallel::LocalGhostDescriptor> descriptors = {
          {.residency = cosmosim::parallel::LocalIndexResidency::kGhost, .owning_rank = 1},
      };
      static_cast<void>(cosmosim::parallel::buildGhostExchangePlan(0, descriptors, 0));
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }
}

void testDeterministicReductionAgreement() {
  const std::vector<double> rank_values = {1.5, -2.0, 3.25, 10.0};
  const double deterministic = cosmosim::parallel::deterministicRankOrderedSum(rank_values);
  assert(std::abs(deterministic - 12.75) < 1.0e-15);

  const auto exact = cosmosim::parallel::compareReductionAgreement(rank_values, deterministic);
  assert(std::abs(exact.deterministic_baseline_sum - deterministic) < 1.0e-15);
  assert(std::abs(exact.measured_sum - deterministic) < 1.0e-15);
  assert(exact.absolute_error < 1.0e-15);
  assert(exact.relative_error < 1.0e-15);

  const auto perturbed = cosmosim::parallel::compareReductionAgreement(rank_values, deterministic + 1.0e-6);
  assert(perturbed.absolute_error > 0.0);
  assert(perturbed.relative_error > 0.0);
  assert(!cosmosim::parallel::satisfiesReductionAgreement(
      perturbed,
      {.absolute_tolerance = 1.0e-9, .relative_tolerance = 1.0e-9}));
  assert(cosmosim::parallel::satisfiesReductionAgreement(
      perturbed,
      {.absolute_tolerance = 2.0e-6, .relative_tolerance = 0.0}));
}

void testRankConfigConsensus() {
  const std::vector<cosmosim::parallel::RankConfigDigest> matching = {
      {.world_rank = 0, .normalized_config_hash = 0xabcU, .mpi_ranks_expected = 2, .deterministic_reduction = true},
      {.world_rank = 1, .normalized_config_hash = 0xabcU, .mpi_ranks_expected = 2, .deterministic_reduction = true},
  };
  const auto ok = cosmosim::parallel::evaluateRankConfigConsensus(matching);
  assert(ok.allConsistent());
  assert(ok.mismatched_ranks.empty());
  assert(ok.mismatches.empty());

  const std::vector<cosmosim::parallel::RankConfigDigest> mismatch = {
      {.world_rank = 0, .normalized_config_hash = 0xabcU, .mpi_ranks_expected = 2, .deterministic_reduction = true},
      {.world_rank = 1, .normalized_config_hash = 0xdefU, .mpi_ranks_expected = 4, .deterministic_reduction = false},
  };
  const auto bad = cosmosim::parallel::evaluateRankConfigConsensus(mismatch);
  assert(!bad.normalized_config_hash_match);
  assert(!bad.mpi_ranks_expected_match);
  assert(!bad.deterministic_reduction_match);
  assert(!bad.allConsistent());
  assert(bad.mismatched_ranks.size() == 1);
  assert(bad.mismatched_ranks[0] == 1);
  assert(bad.mismatches.size() == 3);
  assert(bad.mismatches[0].property == cosmosim::parallel::RankConfigMismatchProperty::kNormalizedConfigHash);
  assert(bad.mismatches[0].baseline_rank == 0);
  assert(bad.mismatches[0].rank == 1);
  assert(bad.mismatches[0].baseline_value == "2748");
  assert(bad.mismatches[0].rank_value == "3567");
  assert(bad.mismatches[1].property == cosmosim::parallel::RankConfigMismatchProperty::kMpiRanksExpected);
  assert(bad.mismatches[2].property == cosmosim::parallel::RankConfigMismatchProperty::kDeterministicReduction);
}

void testGhostBufferPayloadShapeValidation() {
  cosmosim::parallel::GhostExchangeBufferSoA source;
  source.entity_id = {1};
  source.density_code = {2.0};
  source.velocity_x_code = {3.0};
  source.pressure_code = {4.0};
  cosmosim::parallel::GhostExchangeBuffer buffer;
  const std::vector<std::uint32_t> packed_indices = {0};
  buffer.packFrom(source, packed_indices);

  cosmosim::parallel::GhostExchangeBufferSoA destination;
  buffer.unpackAppendTo(destination);
  assert(destination.entity_id.size() == 1);
}

}  // namespace

int main() {
  testGhostPackUnpackRoundTrip();
  testMortonDecompositionInvariants();
  testRestartStateRoundTrip();
  testExplicitOwnedVsGhostContracts();
  testGhostTransferInvariantFailures();
  testDeterministicReductionAgreement();
  testRankConfigConsensus();
  testGhostBufferPayloadShapeValidation();
  return 0;
}
