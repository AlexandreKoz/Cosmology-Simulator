#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

cosmosim::parallel::AmrPatchPayloadRecord makePatch(int rank) {
  cosmosim::parallel::AmrPatchPayloadRecord record;
  record.patch_id = rank == 0 ? 9001U : 9002U;
  record.parent_patch_id = 0U;
  record.morton_key = static_cast<std::uint64_t>(rank + 1);
  record.owner_rank = rank;
  record.level = 0U;
  record.first_cell = 0U;
  record.cell_count = 2U;
  record.origin_x_comoving = rank == 0 ? 0.0 : 0.5;
  record.origin_y_comoving = 0.0;
  record.origin_z_comoving = 0.0;
  record.extent_x_comoving = 0.5;
  record.extent_y_comoving = 1.0;
  record.extent_z_comoving = 1.0;
  record.cell_dim_x = 2U;
  record.cell_dim_y = 1U;
  record.cell_dim_z = 1U;
  record.decomposition_epoch = 42U;
  record.cell_mass_sum_code = 2.0 + 0.25 * static_cast<double>(rank);
  record.gas_internal_energy_sum_code = 3.0 + 0.5 * static_cast<double>(rank);
  cosmosim::parallel::validateAmrPatchPayloadRecord(record);
  return record;
}

std::vector<cosmosim::parallel::AmrPatchCellPayloadRecord> makeCells(int rank) {
  std::vector<cosmosim::parallel::AmrPatchCellPayloadRecord> records;
  records.reserve(2U);
  for (std::uint32_t i = 0; i < 2U; ++i) {
    cosmosim::parallel::AmrPatchCellPayloadRecord record;
    record.patch_id = rank == 0 ? 9001U : 9002U;
    record.owner_rank = rank;
    record.local_cell_offset = i;
    record.patch_index = 0U;
    record.center_x_comoving = (rank == 0 ? 0.125 : 0.625) + 0.25 * static_cast<double>(i);
    record.center_y_comoving = 0.5;
    record.center_z_comoving = 0.5;
    record.mass_code = 1.0 + 0.1 * static_cast<double>(rank) + 0.01 * static_cast<double>(i);
    record.time_bin = i;
    record.gas_cell_id = static_cast<std::uint64_t>(8000 + rank * 10 + static_cast<int>(i));
    record.parent_particle_id = static_cast<std::uint64_t>(5000 + rank * 10 + static_cast<int>(i));
    record.velocity_x_peculiar = rank == 0 ? 0.03 : -0.02;
    record.velocity_y_peculiar = 0.01 * static_cast<double>(i + 1U);
    record.velocity_z_peculiar = -0.015 * static_cast<double>(rank + 1);
    record.density_code = 1.0 + 0.1 * static_cast<double>(rank) + 0.05 * static_cast<double>(i);
    record.pressure_code = 1.2 + 0.05 * static_cast<double>(rank) + 0.02 * static_cast<double>(i);
    record.internal_energy_code = 2.5;
    record.temperature_code = 1.0;
    record.sound_speed_code = 1.0;
    cosmosim::parallel::validateAmrPatchCellPayloadRecord(record);
    records.push_back(record);
  }
  return records;
}

cosmosim::parallel::AmrFluxRegisterPayloadRecord makeFlux(int rank, int peer) {
  cosmosim::parallel::AmrFluxRegisterPayloadRecord record;
  record.register_key = static_cast<std::uint64_t>(99000 + rank);
  record.coarse_patch_id = peer == 0 ? 9001U : 9002U;
  record.coarse_gas_cell_id = static_cast<std::uint64_t>(8000 + peer * 10);
  record.coarse_cell_index = 0U;
  record.level = 1U;
  record.axis = 0U;
  record.orientation = rank == 0 ? 1U : 0U;
  record.source_rank = rank;
  record.owner_rank = peer;
  record.gas_cell_identity_generation = 42U;
  record.patch_geometry_generation = 7U;
  record.coarse_mass_flux_code = 0.1;
  record.coarse_momentum_x_flux_code = 0.2;
  record.coarse_momentum_y_flux_code = 0.03;
  record.coarse_momentum_z_flux_code = -0.04;
  record.coarse_total_energy_flux_code = 0.5;
  record.fine_mass_flux_code = 0.09;
  record.fine_momentum_x_flux_code = 0.18;
  record.fine_momentum_y_flux_code = 0.025;
  record.fine_momentum_z_flux_code = -0.035;
  record.fine_total_energy_flux_code = 0.45;
  record.face_area_comov = 0.25;
  record.coarse_area_comov = 0.25;
  record.fine_area_comov = 0.25;
  record.dt_code = 0.001;
  record.coarse_face_count = 1U;
  record.fine_face_count = 1U;
  cosmosim::parallel::validateAmrFluxRegisterPayloadRecord(record);
  return record;
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI && COSMOSIM_ENABLE_HDF5
  MPI_Init(nullptr, nullptr);
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    MPI_Finalize();
    return 0;
  }

  const int peer_rank = 1 - world_rank;
  const cosmosim::parallel::MpiContext mpi_context;
  const std::vector<cosmosim::parallel::AmrPatchPayloadRecord> local_patches{makePatch(world_rank)};
  const std::vector<cosmosim::parallel::AmrPatchCellPayloadRecord> local_cells = makeCells(world_rank);
  const auto exchange = cosmosim::parallel::executeBlockingDirectedAmrPatchPayloadExchange(
      mpi_context,
      local_patches,
      local_cells,
      123U);

  if (exchange.patch_payloads_received.size() != 1U || exchange.patch_cell_payloads_received.size() != 2U) {
    std::cerr << "rank=" << world_rank
              << " remote_patches=" << exchange.patch_payloads_received.size()
              << " remote_cells=" << exchange.patch_cell_payloads_received.size()
              << " candidates=" << exchange.diagnostics.candidate_peer_count
              << " neighbors=" << exchange.diagnostics.neighbor_peer_count << '\n';
  }
  assert(exchange.diagnostics.candidate_peer_count == 1U);
  assert(exchange.diagnostics.neighbor_peer_count == 1U);
  assert(exchange.diagnostics.directed_patch_descriptor_records_sent == 1U);
  assert(exchange.diagnostics.directed_patch_descriptor_records_received == 1U);
  assert(exchange.diagnostics.directed_patch_cell_records_sent == 2U);
  assert(exchange.diagnostics.directed_patch_cell_records_received == 2U);
  assert(exchange.diagnostics.control_plane_bytes > 0U);
  assert(exchange.diagnostics.patch_descriptor_bytes > 0U);
  assert(exchange.diagnostics.patch_cell_payload_bytes > 0U);
  assert(exchange.diagnostics.remote_patch_ghost_count == 1U);
  assert(exchange.diagnostics.remote_interface_count == 1U);
  assert(exchange.patch_payloads_received.front().owner_rank == peer_rank);
  for (const auto& cell : exchange.patch_cell_payloads_received) {
    assert(cell.owner_rank == peer_rank);
    assert(cell.gas_cell_id != 0U);
  }

  const std::vector<cosmosim::parallel::AmrFluxRegisterPayloadRecord> outbound_flux{makeFlux(world_rank, peer_rank)};
  const auto inbound_flux = cosmosim::parallel::executeBlockingAmrFluxRegisterPayloadExchange(
      mpi_context,
      outbound_flux,
      124U);
  assert(inbound_flux.size() == 1U);
  assert(inbound_flux.front().owner_rank == world_rank);
  assert(inbound_flux.front().source_rank == peer_rank);

  MPI_Finalize();
#endif
  return 0;
}
