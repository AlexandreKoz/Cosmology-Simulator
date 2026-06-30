#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "cosmosim/parallel/distributed_memory.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

[[nodiscard]] double encodedPlaneValue(int rank, std::size_t local_x, std::size_t y, std::size_t z) {
  return static_cast<double>(rank * 1000) + static_cast<double>(local_x * 100) +
      static_cast<double>(y * 10) + static_cast<double>(z);
}

}  // namespace

int main(int argc, char** argv) {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  MPI_Init(&argc, &argv);

  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  assert(world_size == 2);

  constexpr std::size_t k_global_nx = 4;
  constexpr std::size_t k_global_ny = 3;
  constexpr std::size_t k_global_nz = 2;
  constexpr std::size_t k_halo_depth = 1;
  const cosmosim::parallel::PmSlabLayout layout = cosmosim::parallel::makePmSlabLayout(
      k_global_nx, k_global_ny, k_global_nz, world_size, world_rank);
  assert(layout.local_nx() == 2U);

  std::vector<double> local_field(layout.localCellCount(), 0.0);
  for (std::size_t local_x = 0; local_x < layout.local_nx(); ++local_x) {
    const std::size_t global_x = layout.globalXFromLocal(local_x);
    for (std::size_t y = 0; y < k_global_ny; ++y) {
      for (std::size_t z = 0; z < k_global_nz; ++z) {
        local_field[layout.localLinearIndex(global_x, y, z)] = encodedPlaneValue(world_rank, local_x, y, z);
      }
    }
  }

  const cosmosim::parallel::MpiContext mpi_context(/*is_enabled=*/true, world_size, world_rank);
  const cosmosim::parallel::PmSlabHaloExchangeResult exchange =
      cosmosim::parallel::executeBlockingPmSlabHaloExchange(
          mpi_context, layout, local_field, k_halo_depth, /*periodic_x=*/true, /*exchange_sequence=*/17);

  const int peer_rank = 1 - world_rank;
  const std::size_t plane_size = k_global_ny * k_global_nz;
  assert(exchange.left_peer_rank == peer_rank);
  assert(exchange.right_peer_rank == peer_rank);
  assert(exchange.left_peer_rank == exchange.right_peer_rank);
  assert(exchange.halo_depth_x == k_halo_depth);
  assert(exchange.left_halo.size() == plane_size);
  assert(exchange.right_halo.size() == plane_size);
  const std::uint64_t expected_bytes =
      static_cast<std::uint64_t>(2U * plane_size * sizeof(double));
  assert(exchange.sent_bytes == expected_bytes);
  assert(exchange.received_bytes == expected_bytes);

  for (std::size_t y = 0; y < k_global_ny; ++y) {
    for (std::size_t z = 0; z < k_global_nz; ++z) {
      const std::size_t plane_index = y * k_global_nz + z;
      assert(exchange.left_halo[plane_index] == encodedPlaneValue(peer_rank, /*local_x=*/1U, y, z));
      assert(exchange.right_halo[plane_index] == encodedPlaneValue(peer_rank, /*local_x=*/0U, y, z));
    }
  }

  MPI_Finalize();
#else
  (void)argc;
  (void)argv;
#endif
  return 0;
}
