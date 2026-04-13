#include <cassert>
#include <cmath>
#include <cstdint>
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

}  // namespace

int main() {
  testGhostPackUnpackRoundTrip();
  testMortonDecompositionInvariants();
  testRestartStateRoundTrip();
  return 0;
}
