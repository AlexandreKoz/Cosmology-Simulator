#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <vector>

#include "cosmosim/parallel/distributed_memory.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

[[nodiscard]] double encodedPlaneValue(
    int owner_rank,
    std::size_t global_x,
    std::size_t y,
    std::size_t z) {
  return static_cast<double>(owner_rank * 1'000'000) + static_cast<double>(global_x * 10'000) +
      static_cast<double>(y * 100) + static_cast<double>(z);
}

void runPeriodicCase(
    std::size_t global_nx,
    int world_size,
    int world_rank,
    std::uint64_t exchange_sequence) {
  constexpr std::size_t k_global_ny = 3;
  constexpr std::size_t k_global_nz = 2;
  constexpr std::size_t k_halo_depth = 1;
  const cosmosim::parallel::PmSlabLayout layout = cosmosim::parallel::makePmSlabLayout(
      global_nx, k_global_ny, k_global_nz, world_size, world_rank);

  std::vector<double> local_field(layout.localCellCount(), 0.0);
  for (std::size_t local_x = 0; local_x < layout.local_nx(); ++local_x) {
    const std::size_t global_x = layout.globalXFromLocal(local_x);
    for (std::size_t y = 0; y < k_global_ny; ++y) {
      for (std::size_t z = 0; z < k_global_nz; ++z) {
        local_field[layout.localLinearIndex(global_x, y, z)] =
            encodedPlaneValue(world_rank, global_x, y, z);
      }
    }
  }

  const cosmosim::parallel::MpiContext mpi_context(/*is_enabled=*/true, world_size, world_rank);
  const cosmosim::parallel::PmSlabHaloExchangeResult exchange =
      cosmosim::parallel::executeBlockingPmSlabHaloExchange(
          mpi_context,
          layout,
          local_field,
          k_halo_depth,
          /*periodic_x=*/true,
          exchange_sequence);

  if (layout.local_nx() == 0) {
    assert(exchange.halo_depth_x == 0U);
    assert(exchange.left_peer_rank == -1);
    assert(exchange.right_peer_rank == -1);
    assert(exchange.left_halo.empty());
    assert(exchange.right_halo.empty());
    assert(exchange.sent_bytes == 0U);
    assert(exchange.received_bytes == 0U);
    return;
  }

  const std::size_t left_global_x =
      (layout.owned_x.begin_x + global_nx - 1U) % global_nx;
  const std::size_t right_global_x = layout.owned_x.end_x % global_nx;
  const int left_owner =
      cosmosim::parallel::pmOwnerRankForGlobalX(global_nx, world_size, left_global_x);
  const int right_owner =
      cosmosim::parallel::pmOwnerRankForGlobalX(global_nx, world_size, right_global_x);
  const std::size_t plane_size = k_global_ny * k_global_nz;

  assert(exchange.left_peer_rank == left_owner);
  assert(exchange.right_peer_rank == right_owner);
  assert(exchange.halo_depth_x == k_halo_depth);
  assert(exchange.left_halo.size() == plane_size);
  assert(exchange.right_halo.size() == plane_size);
  const std::uint64_t remote_side_count =
      static_cast<std::uint64_t>(left_owner != world_rank) +
      static_cast<std::uint64_t>(right_owner != world_rank);
  const std::uint64_t expected_bytes =
      remote_side_count * static_cast<std::uint64_t>(plane_size * sizeof(double));
  assert(exchange.sent_bytes == expected_bytes);
  assert(exchange.received_bytes == expected_bytes);

  for (std::size_t y = 0; y < k_global_ny; ++y) {
    for (std::size_t z = 0; z < k_global_nz; ++z) {
      const std::size_t plane_index = y * k_global_nz + z;
      assert(exchange.left_halo[plane_index] == encodedPlaneValue(left_owner, left_global_x, y, z));
      assert(exchange.right_halo[plane_index] == encodedPlaneValue(right_owner, right_global_x, y, z));
    }
  }
}

void requireCoordinatedRejection(bool local_rejected, int world_size) {
  int local_vote = local_rejected ? 1 : 0;
  int rejected_rank_count = 0;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  MPI_Allreduce(
      &local_vote,
      &rejected_rank_count,
      1,
      MPI_INT,
      MPI_SUM,
      MPI_COMM_WORLD);
#endif
  assert(rejected_rank_count == world_size);
}

void runDivergentControlRejectionCase(int world_size, int world_rank) {
  constexpr std::size_t k_global_ny = 2;
  constexpr std::size_t k_global_nz = 2;
  const cosmosim::parallel::PmSlabLayout layout =
      cosmosim::parallel::makePmSlabLayout(
          static_cast<std::size_t>(world_size) * 2U,
          k_global_ny,
          k_global_nz,
          world_size,
          world_rank);
  const std::vector<double> local_field(layout.localCellCount(), 0.0);
  const cosmosim::parallel::MpiContext mpi_context(
      /*is_enabled=*/true, world_size, world_rank);

  bool rejected = false;
  try {
    (void)cosmosim::parallel::executeBlockingPmSlabHaloExchange(
        mpi_context,
        layout,
        local_field,
        /*halo_depth_x=*/world_rank == 0 ? 0U : 1U,
        /*periodic_x=*/true,
        /*exchange_sequence=*/41U);
  } catch (const std::exception&) {
    rejected = true;
  }
  requireCoordinatedRejection(rejected, world_size);
}

void runRankLocalPreparationFailureCase(int world_size, int world_rank) {
  constexpr std::size_t k_global_ny = 2;
  constexpr std::size_t k_global_nz = 2;
  const cosmosim::parallel::PmSlabLayout layout =
      cosmosim::parallel::makePmSlabLayout(
          static_cast<std::size_t>(world_size) * 2U,
          k_global_ny,
          k_global_nz,
          world_size,
          world_rank);
  std::vector<double> local_field(layout.localCellCount(), 0.0);
  if (world_rank == 0) {
    local_field.pop_back();
  }
  const cosmosim::parallel::MpiContext mpi_context(
      /*is_enabled=*/true, world_size, world_rank);

  bool rejected = false;
  try {
    (void)cosmosim::parallel::executeBlockingPmSlabHaloExchange(
        mpi_context,
        layout,
        local_field,
        /*halo_depth_x=*/1U,
        /*periodic_x=*/true,
        /*exchange_sequence=*/43U);
  } catch (const std::exception&) {
    rejected = true;
  }
  requireCoordinatedRejection(rejected, world_size);
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

  // Every rank owns at least two planes. This exercises the ordinary periodic
  // ring, including the same-peer left/right topology at two ranks.
  runPeriodicCase(
      static_cast<std::size_t>(world_size) * 2U,
      world_size,
      world_rank,
      /*exchange_sequence=*/17U);
  MPI_Barrier(MPI_COMM_WORLD);

  // nx < ranks leaves one or more zero-width slabs. Neighbor lookup must skip
  // those ranks and preserve the periodic ring among actual plane owners.
  runPeriodicCase(
      static_cast<std::size_t>(world_size - 1),
      world_size,
      world_rank,
      /*exchange_sequence=*/23U);
  MPI_Barrier(MPI_COMM_WORLD);

  // With one global plane, rank zero is both periodic owners and all remaining
  // ranks are empty. No MPI payload bytes should be reported for the local copy.
  runPeriodicCase(
      /*global_nx=*/1U,
      world_size,
      world_rank,
      /*exchange_sequence=*/31U);
  MPI_Barrier(MPI_COMM_WORLD);

  // Divergent controls and rank-local preparation errors must fail on every
  // rank before any point-to-point traffic is posted.
  runDivergentControlRejectionCase(world_size, world_rank);
  MPI_Barrier(MPI_COMM_WORLD);
  runRankLocalPreparationFailureCase(world_size, world_rank);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
#else
  (void)argc;
  (void)argv;
#endif
  return 0;
}
