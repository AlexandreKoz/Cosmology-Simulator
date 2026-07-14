#include <cassert>
#include <cstdint>
#include <exception>
#include <span>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

#if COSMOSIM_ENABLE_MPI
[[nodiscard]] cosmosim::parallel::TreePseudoParticlePacket makeHierarchyRootPacket(
    int world_rank,
    std::uint64_t decomposition_epoch,
    std::uint64_t exchange_sequence) {
  const double center = 0.2 + 0.1 * static_cast<double>(world_rank);
  return cosmosim::parallel::TreePseudoParticlePacket{
      .descriptor = {
          .wire_version = 1U,
          .pseudo_particle_id = static_cast<std::uint64_t>(world_rank) + 1U,
          .source_rank = world_rank,
          .decomposition_epoch = decomposition_epoch,
          .force_epoch = 17U,
          .exchange_sequence = exchange_sequence,
          .derived_not_authoritative = true,
      },
      .mass_code = 1.0 + static_cast<double>(world_rank),
      .center_x_comoving = center,
      .center_y_comoving = center,
      .center_z_comoving = center,
      .min_x_comoving = center - 0.01,
      .max_x_comoving = center + 0.01,
      .min_y_comoving = center - 0.01,
      .max_y_comoving = center + 0.01,
      .min_z_comoving = center - 0.01,
      .max_z_comoving = center + 0.01,
      .source_count = 1U,
      .hierarchy_level = 0U,
      .local_node_index = 0U,
      .child_count = 0U,
      .is_leaf = 1U,
      .geometry_frame = 0U,
  };
}

void requireAllRanksRejected(bool local_rejected, int world_size) {
  const std::uint64_t local_vote = local_rejected ? 1U : 0U;
  std::uint64_t global_votes = 0U;
  MPI_Allreduce(&local_vote, &global_votes, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  assert(global_votes == static_cast<std::uint64_t>(world_size));
}

void testHierarchyExchangeFailureCoordination(int world_size, int world_rank) {
  constexpr std::uint64_t k_valid_exchange_sequence = 9001U;
  const cosmosim::parallel::MpiContext mpi_context(true, world_size, world_rank);
  auto local_packet = makeHierarchyRootPacket(
      world_rank, /*decomposition_epoch=*/41U, k_valid_exchange_sequence);
  const auto gathered = cosmosim::parallel::executeBlockingTreePseudoParticleHierarchyExchange(
      mpi_context,
      std::span<const cosmosim::parallel::TreePseudoParticlePacket>(&local_packet, 1U),
      k_valid_exchange_sequence);
  assert(gathered.size() == static_cast<std::size_t>(world_size));
  for (int rank = 0; rank < world_size; ++rank) {
    assert(gathered[static_cast<std::size_t>(rank)].descriptor.source_rank == rank);
    assert(gathered[static_cast<std::size_t>(rank)].descriptor.decomposition_epoch == 41U);
  }

  // Only rank zero advertises stale communicator metadata. The real communicator
  // must vote before any rank advances into count/payload collectives.
  bool communicator_mismatch_rejected = false;
  try {
    const cosmosim::parallel::MpiContext mismatched_context(
        true,
        world_size + (world_rank == 0 ? 1 : 0),
        world_rank);
    static_cast<void>(cosmosim::parallel::executeBlockingTreePseudoParticleHierarchyExchange(
        mismatched_context,
        std::span<const cosmosim::parallel::TreePseudoParticlePacket>(&local_packet, 1U),
        k_valid_exchange_sequence));
  } catch (const std::exception&) {
    communicator_mismatch_rejected = true;
  }
  requireAllRanksRejected(communicator_mismatch_rejected, world_size);

  // Every local packet is valid, but the gathered hierarchy mixes decomposition
  // epochs. Decode/validation must coordinate the failure before any rank returns.
  constexpr std::uint64_t k_mixed_epoch_exchange_sequence = 9002U;
  local_packet = makeHierarchyRootPacket(
      world_rank,
      /*decomposition_epoch=*/51U + static_cast<std::uint64_t>(world_rank),
      k_mixed_epoch_exchange_sequence);
  bool mixed_epoch_rejected = false;
  try {
    static_cast<void>(cosmosim::parallel::executeBlockingTreePseudoParticleHierarchyExchange(
        mpi_context,
        std::span<const cosmosim::parallel::TreePseudoParticlePacket>(&local_packet, 1U),
        k_mixed_epoch_exchange_sequence));
  } catch (const std::exception&) {
    mixed_epoch_rejected = true;
  }
  requireAllRanksRejected(mixed_epoch_rejected, world_size);
}
#endif

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size > 1) {
    testHierarchyExchangeFailureCoordination(world_size, world_rank);
  }
  MPI_Finalize();
#endif
  return 0;
}
