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
  config.owned_particle_weight = 1.0;
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
  std::uint64_t counted_particles = 0;
  for (const auto value : plan.metrics.owned_particles_by_rank) {
    counted_particles += value;
  }
  assert(counted_particles == 8);
}

void testGravityAwareDecompositionTracksInteractionCost() {
  std::vector<cosmosim::parallel::DecompositionItem> items(8);
  for (std::size_t i = 0; i < items.size(); ++i) {
    items[i].entity_id = static_cast<std::uint64_t>(i + 1);
    items[i].kind = cosmosim::parallel::DecompositionEntityKind::kParticle;
    items[i].x_comov = static_cast<double>(i) / 8.0;
    items[i].y_comov = 0.0;
    items[i].z_comov = 0.0;
    items[i].memory_bytes = 256;
    items[i].active_target_count_recent = (i < 4) ? 1 : 40;
    items[i].remote_tree_interactions_recent = (i < 4) ? 1 : 50;
  }

  cosmosim::parallel::DecompositionConfig baseline_config;
  baseline_config.world_size = 2;
  baseline_config.owned_particle_weight = 1.0;
  const auto baseline = cosmosim::parallel::buildMortonSfcDecomposition(items, baseline_config);

  cosmosim::parallel::DecompositionConfig gravity_config;
  gravity_config.world_size = 2;
  gravity_config.owned_particle_weight = 1.0;
  gravity_config.active_target_weight = 2.0;
  gravity_config.remote_tree_interaction_weight = 3.0;
  gravity_config.memory_weight = 1.0 / 4096.0;

  const auto plan = cosmosim::parallel::buildMortonSfcDecomposition(items, gravity_config);
  assert(plan.metrics.active_targets_by_rank.size() == 2);
  assert(plan.metrics.remote_tree_interactions_by_rank.size() == 2);
  assert(plan.metrics.active_targets_by_rank[0] != plan.metrics.active_targets_by_rank[1]);
  assert(plan.metrics.remote_tree_interactions_by_rank[0] !=
         plan.metrics.remote_tree_interactions_by_rank[1]);
  assert(plan.owning_rank_by_item != baseline.owning_rank_by_item);
  std::uint64_t total_active_targets = 0;
  std::uint64_t total_remote_interactions = 0;
  for (const auto value : plan.metrics.active_targets_by_rank) {
    total_active_targets += value;
  }
  for (const auto value : plan.metrics.remote_tree_interactions_by_rank) {
    total_remote_interactions += value;
  }
  assert(total_active_targets == (4U * 1U + 4U * 40U));
  assert(total_remote_interactions == (4U * 1U + 4U * 50U));
}

void testRestartStateRoundTrip() {
  cosmosim::parallel::DistributedRestartState in;
  in.schema_version = 2;
  in.decomposition_epoch = 7;
  in.world_size = 2;
  in.pm_grid_nx = 16;
  in.pm_grid_ny = 16;
  in.pm_grid_nz = 16;
  in.pm_decomposition_mode = "slab";
  in.gravity_kick_opportunity = 5;
  in.pm_update_cadence_steps = 2;
  in.long_range_field_version = 3;
  in.last_long_range_refresh_opportunity = 4;
  in.long_range_field_built_step_index = 12;
  in.long_range_field_built_scale_factor = 0.4;
  in.long_range_restart_policy = "deterministic_rebuild";
  in.owning_rank_by_item = {0, 1, 1, 0};
  in.pm_slab_begin_x_by_rank = {0, 8};
  in.pm_slab_end_x_by_rank = {8, 16};

  const std::string encoded = in.serialize();
  const auto out = cosmosim::parallel::DistributedRestartState::deserialize(encoded);

  assert(out.schema_version == in.schema_version);
  assert(out.decomposition_epoch == in.decomposition_epoch);
  assert(out.world_size == in.world_size);
  assert(out.pm_grid_nx == in.pm_grid_nx);
  assert(out.pm_grid_ny == in.pm_grid_ny);
  assert(out.pm_grid_nz == in.pm_grid_nz);
  assert(out.gravity_kick_opportunity == in.gravity_kick_opportunity);
  assert(out.long_range_field_version == in.long_range_field_version);
  assert(out.owning_rank_by_item == in.owning_rank_by_item);
  assert(out.pm_slab_begin_x_by_rank == in.pm_slab_begin_x_by_rank);
  assert(out.pm_slab_end_x_by_rank == in.pm_slab_end_x_by_rank);
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
  assert(plan.send_local_indices_by_neighbor[0].empty());
  assert(plan.send_local_indices_by_neighbor[1].empty());
  assert(plan.outbound_transfers.size() == 2);
  assert(plan.inbound_transfers.size() == 2);
  assert(plan.outbound_transfers[0].role == cosmosim::parallel::GhostTransferRole::kOutboundSend);
  assert(plan.inbound_transfers[0].role == cosmosim::parallel::GhostTransferRole::kInboundReceive);
  assert(plan.outbound_transfers[0].intent ==
         cosmosim::parallel::GhostTransferIntent::kGhostRefreshRequest);
  assert(plan.inbound_transfers[0].intent ==
         cosmosim::parallel::GhostTransferIntent::kGhostRefreshReceiveStaging);
  assert(plan.outbound_transfers[0].expected_post_transfer_residency ==
         cosmosim::parallel::LocalIndexResidency::kGhost);
  assert(plan.inbound_transfers[0].expected_post_transfer_residency ==
         cosmosim::parallel::LocalIndexResidency::kGhost);
  assert(plan.outbound_transfers[0].neighbor_slot == 0);
  assert(plan.inbound_transfers[1].neighbor_slot == 1);
  assert(plan.outbound_transfers[0].peer_rank == 1);
  assert(plan.inbound_transfers[1].peer_rank == 2);
  assert(plan.outbound_transfers[1].local_indices == plan.send_local_indices_by_neighbor[1]);
  assert(plan.inbound_transfers[1].local_indices == plan.recv_local_indices_by_neighbor[1]);
  assert(plan.send_bytes == 0);
  assert(plan.recv_bytes == 2 * (sizeof(double) * 3 + sizeof(std::uint64_t)));
  cosmosim::parallel::validateGhostExchangePlan(plan);
}

void testRestartStateRejectsMissingOrDuplicateOwnershipEntries() {
  {
    bool threw = false;
    try {
      (void)cosmosim::parallel::DistributedRestartState::deserialize(
          "schema_version=1\ndecomposition_epoch=7\nworld_size=2\nitem_count=3\nrank[0]=1\n");
    } catch (const std::runtime_error&) {
      threw = true;
    }
    assert(threw);
  }

  {
    bool threw = false;
    try {
      (void)cosmosim::parallel::DistributedRestartState::deserialize(
          "schema_version=1\ndecomposition_epoch=7\nworld_size=2\nitem_count=2\nrank[0]=1\nrank[0]=0\nrank[1]=1\n");
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }

  {
    bool threw = false;
    try {
      (void)cosmosim::parallel::DistributedRestartState::deserialize(
          "schema_version=2\ndecomposition_epoch=7\nworld_size=2\npm_grid_nx=8\npm_grid_ny=8\npm_grid_nz=8\n"
          "pm_decomposition_mode=slab\ngravity_kick_opportunity=3\npm_update_cadence_steps=2\n"
          "long_range_field_version=0\nlast_long_range_refresh_opportunity=2\n"
          "long_range_field_built_step_index=0\nlong_range_field_built_scale_factor=1\n"
          "long_range_restart_policy=deterministic_rebuild\nitem_count=2\nrank[0]=0\nrank[1]=1\n"
          "pm_slab_rank_count=2\npm_slab_begin_x[0]=0\npm_slab_end_x[0]=4\npm_slab_begin_x[1]=4\npm_slab_end_x[1]=8\n");
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }
}

void testDistributedRestartCompatibilityReporting() {
  cosmosim::parallel::DistributedRestartState restart;
  restart.schema_version = 2;
  restart.world_size = 2;
  restart.pm_grid_nx = 8;
  restart.pm_grid_ny = 8;
  restart.pm_grid_nz = 8;
  restart.pm_decomposition_mode = "slab";
  restart.pm_slab_begin_x_by_rank = {0, 4};
  restart.pm_slab_end_x_by_rank = {4, 8};

  cosmosim::parallel::DistributedExecutionTopology runtime;
  runtime.world_size = 2;
  runtime.world_rank = 1;
  runtime.pm_slab = cosmosim::parallel::makePmSlabLayout(8, 8, 8, 2, 1);
  const auto ok = cosmosim::parallel::evaluateDistributedRestartCompatibility(restart, runtime);
  assert(ok.compatible());
  assert(ok.mismatch_messages.empty());

  runtime.pm_slab = cosmosim::parallel::makePmSlabLayout(10, 8, 8, 2, 1);
  const auto bad = cosmosim::parallel::evaluateDistributedRestartCompatibility(restart, runtime);
  assert(!bad.compatible());
  assert(!bad.pm_grid_shape_match);
  assert(!bad.mismatch_messages.empty());

  restart.pm_update_cadence_steps = 0;
  runtime.pm_slab = cosmosim::parallel::makePmSlabLayout(8, 8, 8, 2, 1);
  const auto cadence_bad = cosmosim::parallel::evaluateDistributedRestartCompatibility(restart, runtime);
  assert(!cadence_bad.compatible());
  assert(!cadence_bad.pm_cadence_steps_match);
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
      {.mode = cosmosim::parallel::ReductionAgreementMode::kAbsoluteOrRelative,
       .absolute_tolerance = 2.0e-6,
       .relative_tolerance = 0.0}));
}

void testGhostExchangePlanValidationDriftRejection() {
  const std::vector<cosmosim::parallel::LocalGhostDescriptor> descriptors = {
      {.residency = cosmosim::parallel::LocalIndexResidency::kOwned, .owning_rank = 0},
      {.residency = cosmosim::parallel::LocalIndexResidency::kGhost, .owning_rank = 2},
      {.residency = cosmosim::parallel::LocalIndexResidency::kGhost, .owning_rank = 1},
  };
  auto plan = cosmosim::parallel::buildGhostExchangePlan(0, descriptors, sizeof(double));

  {
    bool threw = false;
    auto bad = plan;
    bad.outbound_transfers[0].role = cosmosim::parallel::GhostTransferRole::kInboundReceive;
    try {
      cosmosim::parallel::validateGhostExchangePlan(bad);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }

  {
    bool threw = false;
    auto bad = plan;
    bad.inbound_transfers[0].peer_rank = 99;
    try {
      cosmosim::parallel::validateGhostExchangePlan(bad);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }

  {
    bool threw = false;
    auto bad = plan;
    bad.outbound_transfers[0].local_indices.push_back(42);
    try {
      cosmosim::parallel::validateGhostExchangePlan(bad);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }
}

void testReductionAgreementPolicyModes() {
  const cosmosim::parallel::ReductionAgreement agreement{
      .deterministic_baseline_sum = 10.0,
      .measured_sum = 10.2,
      .absolute_error = 0.2,
      .relative_error = 0.02,
  };

  assert(cosmosim::parallel::satisfiesReductionAgreement(
      agreement,
      {.mode = cosmosim::parallel::ReductionAgreementMode::kAbsoluteOnly,
       .absolute_tolerance = 0.25,
       .relative_tolerance = 0.0}));
  assert(!cosmosim::parallel::satisfiesReductionAgreement(
      agreement,
      {.mode = cosmosim::parallel::ReductionAgreementMode::kAbsoluteOnly,
       .absolute_tolerance = 0.1,
       .relative_tolerance = 1.0}));
  assert(cosmosim::parallel::satisfiesReductionAgreement(
      agreement,
      {.mode = cosmosim::parallel::ReductionAgreementMode::kRelativeOnly,
       .absolute_tolerance = 0.0,
       .relative_tolerance = 0.03}));
  assert(!cosmosim::parallel::satisfiesReductionAgreement(
      agreement,
      {.mode = cosmosim::parallel::ReductionAgreementMode::kRelativeOnly,
       .absolute_tolerance = 1.0,
       .relative_tolerance = 0.01}));
  assert(cosmosim::parallel::satisfiesReductionAgreement(
      agreement,
      {.mode = cosmosim::parallel::ReductionAgreementMode::kAbsoluteAndRelative,
       .absolute_tolerance = 0.3,
       .relative_tolerance = 0.03}));
  assert(!cosmosim::parallel::satisfiesReductionAgreement(
      agreement,
      {.mode = cosmosim::parallel::ReductionAgreementMode::kAbsoluteAndRelative,
       .absolute_tolerance = 0.3,
       .relative_tolerance = 0.01}));
  assert(cosmosim::parallel::satisfiesReductionAgreement(
      agreement,
      {.mode = cosmosim::parallel::ReductionAgreementMode::kAbsoluteOrRelative,
       .absolute_tolerance = 0.1,
       .relative_tolerance = 0.03}));
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


void testMpiContextContractValidation() {
  const cosmosim::parallel::MpiContext runtime(true, 4, 2);
  runtime.validateExpectedWorldSizeOrThrow(4);
  assert(runtime.isEnabled());
  assert(!runtime.isRoot());
  assert(runtime.worldSize() == 4);
  assert(runtime.worldRank() == 2);

  bool threw = false;
  try {
    runtime.validateExpectedWorldSizeOrThrow(3);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}

void testDistributedExecutionTopologyCpuOnly() {
  const cosmosim::parallel::MpiContext runtime(true, 3, 1);
  const auto topology = cosmosim::parallel::buildDistributedExecutionTopology(
      12, 8, 4, runtime, /*mpi_ranks_expected=*/3, /*configured_gpu_devices=*/0, /*cuda_runtime_available=*/false, /*visible_device_count=*/0);
  assert(topology.mpi_enabled);
  assert(topology.world_size == 3);
  assert(topology.world_rank == 1);
  assert(topology.isDistributed());
  assert(!topology.usesCuda());
  assert(topology.pm_slab.owned_x.begin_x == 4);
  assert(topology.pm_slab.owned_x.end_x == 8);
}

void testDistributedExecutionTopologyCudaAssignment() {
  const cosmosim::parallel::MpiContext runtime(true, 4, 3);
  const auto topology = cosmosim::parallel::buildDistributedExecutionTopology(
      16, 8, 8, runtime, /*mpi_ranks_expected=*/4, /*configured_gpu_devices=*/2, /*cuda_runtime_available=*/true, /*visible_device_count=*/4);
  assert(topology.usesCuda());
  assert(topology.device_assignment.requested_device_count == 2);
  assert(topology.device_assignment.active_device_count == 2);
  assert(topology.device_assignment.assigned_device_index == 1);
  assert(topology.device_assignment.isValid());
}

void testDistributedExecutionTopologyRejectsInvalidGpuRequest() {
  bool threw = false;
  try {
    const cosmosim::parallel::MpiContext runtime(true, 2, 0);
    (void)cosmosim::parallel::buildDistributedExecutionTopology(
        8, 8, 8, runtime, /*mpi_ranks_expected=*/2, /*configured_gpu_devices=*/1, /*cuda_runtime_available=*/false, /*visible_device_count=*/0);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}

void testPmSlabUnevenPartitionOwnership() {
  const std::size_t nx = 10;
  const int world_size = 3;
  const auto r0 = cosmosim::parallel::pmOwnedXRangeForRank(nx, world_size, 0);
  const auto r1 = cosmosim::parallel::pmOwnedXRangeForRank(nx, world_size, 1);
  const auto r2 = cosmosim::parallel::pmOwnedXRangeForRank(nx, world_size, 2);

  assert(r0.begin_x == 0 && r0.end_x == 4);
  assert(r1.begin_x == 4 && r1.end_x == 7);
  assert(r2.begin_x == 7 && r2.end_x == 10);

  for (std::size_t x = 0; x < nx; ++x) {
    const int owner = cosmosim::parallel::pmOwnerRankForGlobalX(nx, world_size, x);
    if (x < 4) {
      assert(owner == 0);
    } else if (x < 7) {
      assert(owner == 1);
    } else {
      assert(owner == 2);
    }
  }
}

void testPmSlabLayoutRoundTripAndCellOwnership() {
  const auto layout = cosmosim::parallel::makePmSlabLayout(
      /*global_nx=*/10, /*global_ny=*/6, /*global_nz=*/4, /*world_size=*/3, /*world_rank=*/1);
  assert(layout.isValid());
  assert(layout.local_nx() == 3);
  assert(layout.localCellCount() == 3 * 6 * 4);
  assert(!layout.ownsFullDomain());
  assert(layout.owned_x.begin_x == 4);
  assert(layout.owned_x.end_x == 7);

  for (std::size_t global_x = layout.owned_x.begin_x; global_x < layout.owned_x.end_x; ++global_x) {
    const std::size_t local_x = layout.localXFromGlobal(global_x);
    assert(layout.globalXFromLocal(local_x) == global_x);
  }

  const std::size_t idx = layout.localLinearIndex(/*global_x=*/5, /*global_y=*/2, /*global_z=*/3);
  assert(idx == ((1 * 6 + 2) * 4 + 3));
  assert(layout.ownsGlobalCell(/*global_x=*/6, /*global_y=*/5, /*global_z=*/3));
  assert(!layout.ownsGlobalCell(/*global_x=*/7, /*global_y=*/1, /*global_z=*/1));

  bool threw = false;
  try {
    (void)layout.localXFromGlobal(8);
  } catch (const std::out_of_range&) {
    threw = true;
  }
  assert(threw);

  threw = false;
  try {
    (void)cosmosim::parallel::pmOwnerRankForGlobalCell(
        10, 6, 4, 3, /*global_x=*/9, /*global_y=*/5, /*global_z=*/3);
  } catch (...) {
    threw = true;
  }
  assert(!threw);
  assert(cosmosim::parallel::pmOwnerRankForGlobalCell(10, 6, 4, 3, 9, 5, 3) == 2);
}

}  // namespace

int main() {
  testGhostPackUnpackRoundTrip();
  testMortonDecompositionInvariants();
  testRestartStateRoundTrip();
  testGravityAwareDecompositionTracksInteractionCost();
  testExplicitOwnedVsGhostContracts();
  testRestartStateRejectsMissingOrDuplicateOwnershipEntries();
  testDistributedRestartCompatibilityReporting();
  testGhostTransferInvariantFailures();
  testDeterministicReductionAgreement();
  testGhostExchangePlanValidationDriftRejection();
  testReductionAgreementPolicyModes();
  testRankConfigConsensus();
  testGhostBufferPayloadShapeValidation();
  testMpiContextContractValidation();
  testDistributedExecutionTopologyCpuOnly();
  testDistributedExecutionTopologyCudaAssignment();
  testDistributedExecutionTopologyRejectsInvalidGpuRequest();
  testPmSlabUnevenPartitionOwnership();
  testPmSlabLayoutRoundTripAndCellOwnership();
  return 0;
}
