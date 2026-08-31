#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/parallel/distributed_memory.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
void setTestEnvironment(const char* name, const std::string& value) {
  const int status = ::setenv(name, value.c_str(), 1);
  assert(status == 0);
}

void clearTestEnvironment(const char* name) {
  const int status = ::unsetenv(name);
  assert(status == 0);
}

void requireAllRanksObservedFailure(
    const cosmosim::parallel::MpiContext& mpi_context,
    bool local_failed) {
  const std::uint64_t failed_count =
      mpi_context.allreduceSumUint64(local_failed ? 1U : 0U);
  assert(failed_count == static_cast<std::uint64_t>(mpi_context.worldSize()));
}

void checkForcedBoundedAllgather(
    const cosmosim::parallel::MpiContext& mpi_context) {
  // Eight bytes across 2-3 ranks produces a tiny per-rank slice and forces
  // several classic-MPI-safe rounds without large allocations.
  setTestEnvironment("COSMOSIM_MPI_TEST_TRANSPORT_LIMIT_BYTES", "8");
  std::vector<std::uint8_t> local(
      static_cast<std::size_t>(5 + mpi_context.worldRank()),
      static_cast<std::uint8_t>(20 + mpi_context.worldRank()));
  const auto gathered = mpi_context.allgatherBytesBounded(local);

  std::vector<std::uint8_t> expected;
  for (int rank = 0; rank < mpi_context.worldSize(); ++rank) {
    expected.insert(
        expected.end(),
        static_cast<std::size_t>(5 + rank),
        static_cast<std::uint8_t>(20 + rank));
  }
  assert(gathered == expected);
  clearTestEnvironment("COSMOSIM_MPI_TEST_TRANSPORT_LIMIT_BYTES");
}

void checkCollectivePreparationRejection(
    const cosmosim::parallel::MpiContext& mpi_context) {
  assert(mpi_context.worldSize() >= 2);
  setTestEnvironment("COSMOSIM_MPI_TEST_FAULT", "gather_post_count:1");
  const std::vector<std::uint8_t> local{
      static_cast<std::uint8_t>(mpi_context.worldRank() + 1)};
  bool failed = false;
  try {
    (void)mpi_context.gatherBytesToRoot(local, 0);
  } catch (const std::runtime_error& error) {
    const std::string message = error.what();
    failed = message.find("gatherBytesToRoot payload preparation") !=
             std::string::npos;
  }
  requireAllRanksObservedFailure(mpi_context, failed);
  clearTestEnvironment("COSMOSIM_MPI_TEST_FAULT");
  assert(MPI_Barrier(MPI_COMM_WORLD) == MPI_SUCCESS);
}

cosmosim::parallel::GhostExchangeBufferSoA makeAuthoritativeGhostState(
    int world_size,
    int world_rank,
    const cosmosim::parallel::GhostLayerEpoch& epoch) {
  cosmosim::parallel::GhostExchangeBufferSoA state;
  state.epoch = epoch;
  state.entity_id.push_back(static_cast<std::uint64_t>(1000 + world_rank));
  for (int peer = 0; peer < world_size; ++peer) {
    if (peer != world_rank) {
      state.entity_id.push_back(static_cast<std::uint64_t>(1000 + peer));
    }
  }
  const std::size_t count = state.entity_id.size();
  state.position_x_comoving.resize(count, static_cast<double>(world_rank));
  state.position_y_comoving.resize(count, 0.25);
  state.position_z_comoving.resize(count, 0.75);
  state.mass_code.resize(count, 1.0 + static_cast<double>(world_rank));
  state.density_code.resize(count, 2.0);
  state.velocity_x_code.resize(count, 0.0);
  state.velocity_y_code.resize(count, 0.0);
  state.velocity_z_code.resize(count, 0.0);
  state.pressure_code.resize(count, 3.0);
  state.internal_energy_code.resize(count, 4.0);
  return state;
}

std::vector<cosmosim::parallel::LocalGhostDescriptor> makeGhostDescriptors(
    int world_size,
    int world_rank,
    const cosmosim::parallel::GhostLayerEpoch& epoch) {
  using cosmosim::parallel::LocalGhostDescriptor;
  using cosmosim::parallel::LocalIndexResidency;
  std::vector<LocalGhostDescriptor> descriptors;
  descriptors.push_back(LocalGhostDescriptor{
      .residency = LocalIndexResidency::kOwned,
      .owning_rank = world_rank,
      .particle_id = static_cast<std::uint64_t>(1000 + world_rank),
      .epoch = epoch,
  });
  for (int peer = 0; peer < world_size; ++peer) {
    if (peer == world_rank) {
      continue;
    }
    descriptors.push_back(LocalGhostDescriptor{
        .residency = LocalIndexResidency::kGhost,
        .owning_rank = peer,
        .particle_id = static_cast<std::uint64_t>(1000 + peer),
        .epoch = epoch,
    });
  }
  return descriptors;
}

std::vector<cosmosim::parallel::LocalGhostDescriptor> makeSparseGhostDescriptors(
    int world_rank,
    const cosmosim::parallel::GhostLayerEpoch& epoch) {
  using cosmosim::parallel::LocalGhostDescriptor;
  using cosmosim::parallel::LocalIndexResidency;
  std::vector<LocalGhostDescriptor> descriptors;
  descriptors.push_back(LocalGhostDescriptor{
      .residency = LocalIndexResidency::kOwned,
      .owning_rank = world_rank,
      .particle_id = static_cast<std::uint64_t>(1000 + world_rank),
      .epoch = epoch,
  });
  if (world_rank == 0 || world_rank == 1) {
    const int peer_rank = 1 - world_rank;
    descriptors.push_back(LocalGhostDescriptor{
        .residency = LocalIndexResidency::kGhost,
        .owning_rank = peer_rank,
        .particle_id = static_cast<std::uint64_t>(1000 + peer_rank),
        .epoch = epoch,
    });
  }
  return descriptors;
}

cosmosim::parallel::GhostExchangeBufferSoA makeGhostStateForDescriptors(
    int world_rank,
    const cosmosim::parallel::GhostLayerEpoch& epoch,
    const std::vector<cosmosim::parallel::LocalGhostDescriptor>& descriptors) {
  cosmosim::parallel::GhostExchangeBufferSoA state;
  state.epoch = epoch;
  state.entity_id.reserve(descriptors.size());
  for (const auto& descriptor : descriptors) {
    state.entity_id.push_back(descriptor.particle_id);
  }
  const std::size_t count = state.entity_id.size();
  state.position_x_comoving.resize(count, static_cast<double>(world_rank));
  state.position_y_comoving.resize(count, 0.25);
  state.position_z_comoving.resize(count, 0.75);
  state.mass_code.resize(count, 1.0 + static_cast<double>(world_rank));
  state.density_code.resize(count, 2.0);
  state.velocity_x_code.resize(count, 0.0);
  state.velocity_y_code.resize(count, 0.0);
  state.velocity_z_code.resize(count, 0.0);
  state.pressure_code.resize(count, 3.0);
  state.internal_energy_code.resize(count, 4.0);
  return state;
}

void checkSparseZeroNeighborSuccess(
    const cosmosim::parallel::MpiContext& mpi_context) {
  if (mpi_context.worldSize() != 3) {
    return;
  }
  const cosmosim::parallel::GhostLayerEpoch epoch{
      .decomposition_epoch = 42U,
      .ghost_sync_epoch = 18U,
      .particle_index_generation = 10U,
  };
  const auto descriptors = makeSparseGhostDescriptors(
      mpi_context.worldRank(), epoch);
  const auto state = makeGhostStateForDescriptors(
      mpi_context.worldRank(), epoch, descriptors);

  const auto exchange =
      cosmosim::parallel::executeBlockingGhostRefreshExchangeFromDescriptors(
          mpi_context, descriptors, state, epoch);
  if (mpi_context.worldRank() == 2) {
    assert(exchange.plan.neighbor_ranks.empty());
    assert(exchange.result.received_ghosts.size() == 0U);
    assert(exchange.result.sent_bytes == 0U);
    assert(exchange.result.received_bytes == 0U);
  } else {
    const int peer_rank = 1 - mpi_context.worldRank();
    assert(exchange.plan.neighbor_ranks.size() == 1U);
    assert(exchange.plan.neighbor_ranks.front() == peer_rank);
    assert(exchange.result.received_ghosts.size() == 1U);
    assert(exchange.result.received_ghosts.entity_id.front() ==
           static_cast<std::uint64_t>(1000 + peer_rank));
    assert(exchange.result.sent_bytes > 0U);
    assert(exchange.result.received_bytes > 0U);
  }
  assert(MPI_Barrier(MPI_COMM_WORLD) == MPI_SUCCESS);
}

void checkSparseZeroNeighborPreparationRejection(
    const cosmosim::parallel::MpiContext& mpi_context) {
  if (mpi_context.worldSize() != 3) {
    return;
  }
  const cosmosim::parallel::GhostLayerEpoch epoch{
      .decomposition_epoch = 43U,
      .ghost_sync_epoch = 19U,
      .particle_index_generation = 11U,
  };
  const auto descriptors = makeSparseGhostDescriptors(
      mpi_context.worldRank(), epoch);
  const auto state = makeGhostStateForDescriptors(
      mpi_context.worldRank(), epoch, descriptors);

  // Rank 2 has no payload peer, but it must still remain in the MPI-world
  // protocol and observe rank 0's coordinated payload-preparation failure.
  setTestEnvironment(
      "COSMOSIM_MPI_TEST_FAULT", "ghost_payload_post_metadata:0");
  bool failed = false;
  try {
    (void)cosmosim::parallel::executeBlockingGhostRefreshExchangeFromDescriptors(
        mpi_context, descriptors, state, epoch);
  } catch (const std::runtime_error& error) {
    const std::string message = error.what();
    failed = message.find("ghost refresh sequential peer payload exchange") !=
             std::string::npos;
  }
  requireAllRanksObservedFailure(mpi_context, failed);
  clearTestEnvironment("COSMOSIM_MPI_TEST_FAULT");
  assert(MPI_Barrier(MPI_COMM_WORLD) == MPI_SUCCESS);
}

void checkSequentialPeerPreparationRejection(
    const cosmosim::parallel::MpiContext& mpi_context) {
  if (mpi_context.worldSize() < 3) {
    return;
  }
  const cosmosim::parallel::GhostLayerEpoch epoch{
      .decomposition_epoch = 41U,
      .ghost_sync_epoch = 17U,
      .particle_index_generation = 9U,
  };
  const auto descriptors = makeGhostDescriptors(
      mpi_context.worldSize(), mpi_context.worldRank(), epoch);
  const auto state = makeAuthoritativeGhostState(
      mpi_context.worldSize(), mpi_context.worldRank(), epoch);

  // Rank 0 rejects the first peer's payload receive allocation after metadata.
  // The implementation must still service its later peer(s) through readiness
  // handshakes, then fail coherently only after the sequential schedule is
  // drained. A post-failure Barrier proves no later rank was stranded.
  setTestEnvironment(
      "COSMOSIM_MPI_TEST_FAULT", "ghost_payload_post_metadata:0");
  bool failed = false;
  try {
    (void)cosmosim::parallel::executeBlockingGhostRefreshExchangeFromDescriptors(
        mpi_context, descriptors, state, epoch);
  } catch (const std::runtime_error& error) {
    const std::string message = error.what();
    failed = message.find("ghost refresh sequential peer payload exchange") !=
             std::string::npos;
  }
  requireAllRanksObservedFailure(mpi_context, failed);
  clearTestEnvironment("COSMOSIM_MPI_TEST_FAULT");
  assert(MPI_Barrier(MPI_COMM_WORLD) == MPI_SUCCESS);
}
#endif

}  // namespace

int main(int argc, char** argv) {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  assert(MPI_Init(&argc, &argv) == MPI_SUCCESS);
  int world_size = 1;
  int world_rank = 0;
  assert(MPI_Comm_size(MPI_COMM_WORLD, &world_size) == MPI_SUCCESS);
  assert(MPI_Comm_rank(MPI_COMM_WORLD, &world_rank) == MPI_SUCCESS);
  assert(world_size >= 2 && world_size <= 3);

  const cosmosim::parallel::MpiContext mpi_context(
      /*is_enabled=*/true, world_size, world_rank);
  checkForcedBoundedAllgather(mpi_context);
  checkCollectivePreparationRejection(mpi_context);
  checkSparseZeroNeighborSuccess(mpi_context);
  checkSparseZeroNeighborPreparationRejection(mpi_context);
  checkSequentialPeerPreparationRejection(mpi_context);

  assert(MPI_Finalize() == MPI_SUCCESS);
#else
  (void)argc;
  (void)argv;
#endif
  return 0;
}
