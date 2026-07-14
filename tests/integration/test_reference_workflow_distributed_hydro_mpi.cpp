#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>

#include "cosmosim/cosmosim.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

std::string configText(int world_size, std::string_view run_name) {
  std::ostringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[mode]\n";
  stream << "mode = cosmo_cube\n";
  stream << "ic_file = generated\n";
  stream << "hydro_boundary = periodic\n\n";
  stream << "[cosmology]\n";
  stream << "box_size_x = 2.0\n";
  stream << "box_size_y = 0.5\n";
  stream << "box_size_z = 0.5\n\n";
  stream << "[numerics]\n";
  stream << "time_begin_code = 0.01\n";
  stream << "time_end_code = 0.0102\n";
  stream << "max_global_steps = 2\n";
  stream << "hierarchical_max_rung = 0\n";
  // Match mesh spacing across the 4:1 x/y,z aspect ratio. The resulting
  // r_cut remains below half the shortest periodic axis.
  stream << "treepm_pm_grid_nx = 32\n";
  stream << "treepm_pm_grid_ny = 8\n";
  stream << "treepm_pm_grid_nz = 8\n";
  stream << "treepm_asmth_cells = 1.25\n";
  stream << "treepm_rcut_cells = 3.9\n";
  stream << "treepm_update_cadence_steps = 1\n\n";
  stream << "[output]\n";
  stream << "run_name = " << run_name << '\n';
  stream << "output_directory = integration_outputs\n";
  stream << "output_stem = snapshot\n";
  stream << "restart_stem = restart\n";
  stream << "snapshot_interval_steps = 1\n";
  stream << "write_restarts = true\n\n";
  stream << "[parallel]\n";
  stream << "mpi_ranks_expected = " << world_size << "\n";
  stream << "decomposition_runtime_rebalance_enabled = false\n";
  return stream.str();
}

cosmosim::core::SimulationState makeRankState(int world_rank) {
  (void)world_rank;
  namespace core = cosmosim::core;
  core::SimulationState state;
  state.resizeParticles(4);
  state.resizeCells(4);
  state.resizePatches(2);
  state.metadata.run_name = "distributed_hydro_mpi";
  state.metadata.snapshot_stem = "snapshot";
  state.metadata.restart_stem = "restart";

  state.species.count_by_species.fill(0U);
  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)] = 4U;
  for (std::uint32_t row = 0; row < 4U; ++row) {
    const std::uint32_t owner_rank = row < 2U ? 0U : 1U;
    const std::uint64_t particle_id = 1U + row;
    state.particle_sidecar.particle_id[row] = particle_id;
    state.particle_sidecar.sfc_key[row] = row;
    state.particle_sidecar.species_tag[row] = static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    state.particle_sidecar.owning_rank[row] = owner_rank;
    state.particles.position_x_comoving[row] = 0.25 + 0.5 * static_cast<double>(row);
    state.particles.position_y_comoving[row] = 0.25;
    state.particles.position_z_comoving[row] = 0.25;
    state.particles.mass_code[row] = 0.5;
    state.particles.time_bin[row] = 0U;
  }
  state.rebuildSpeciesIndex();

  for (std::uint32_t patch = 0; patch < 2U; ++patch) {
    state.patches.patch_id[patch] = 9001U + patch;
    state.patches.level[patch] = 0;
    state.patches.first_cell[patch] = patch * 2U;
    state.patches.cell_count[patch] = 2U;
    state.patches.owning_rank[patch] = patch;
    state.patches.origin_x_comoving[patch] = 0.0;
    state.patches.origin_y_comoving[patch] = 0.0;
    state.patches.origin_z_comoving[patch] = 0.0;
    state.patches.extent_x_comoving[patch] = 0.0;
    state.patches.extent_y_comoving[patch] = 0.0;
    state.patches.extent_z_comoving[patch] = 0.0;
    state.patches.cell_dim_x[patch] = 0U;
    state.patches.cell_dim_y[patch] = 0U;
    state.patches.cell_dim_z[patch] = 0U;
  }

  std::vector<core::GasCellIdentityRecord> records;
  records.reserve(4U);
  for (std::uint32_t row = 0; row < 4U; ++row) {
    const std::uint32_t owner_rank = row < 2U ? 0U : 1U;
    const double rho = 1.0 + 0.1 * static_cast<double>(owner_rank) + 0.03 * static_cast<double>(row % 2U);
    const double pressure = 1.0 + 0.2 * static_cast<double>(owner_rank) - 0.05 * static_cast<double>(row % 2U);
    const double vx = owner_rank == 0 ? 0.04 : -0.03;
    state.cells.center_x_comoving[row] = 0.25 + 0.5 * static_cast<double>(row);
    state.cells.center_y_comoving[row] = 0.25;
    state.cells.center_z_comoving[row] = 0.25;
    state.cells.mass_code[row] = rho * 0.125;
    state.cells.patch_index[row] = owner_rank;
    state.cells.time_bin[row] = 0U;
    state.gas_cells.density_code[row] = rho;
    state.gas_cells.pressure_code[row] = pressure;
    state.gas_cells.internal_energy_code[row] = pressure / ((5.0 / 3.0 - 1.0) * rho);
    state.gas_cells.temperature_code[row] = pressure / rho;
    state.gas_cells.sound_speed_code[row] = std::sqrt((5.0 / 3.0) * pressure / rho);
    state.gas_cells.velocity_x_peculiar[row] = vx;
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = 70001U + static_cast<std::uint64_t>(world_rank) * 10U + row,
        .parent_particle_id = state.particle_sidecar.particle_id[row],
        .owning_patch_id = state.patches.patch_id[owner_rank],
        .local_cell_row = row,
    });
  }
  state.replaceGasCellIdentityRecords(std::move(records));
  assert(state.validateOwnershipInvariants());
  return state;
}

double localMass(const cosmosim::core::SimulationState& state) {
  double mass = 0.0;
  for (double value : state.cells.mass_code) {
    mass += value;
  }
  return mass;
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

  const cosmosim::core::SimulationState initial_state = makeRankState(world_rank);
  const std::filesystem::path root =
      std::filesystem::temp_directory_path() / "cosmosim_reference_workflow_distributed_hydro_mpi";

  const auto direct_frozen = cosmosim::core::loadFrozenConfigFromString(
      configText(world_size, "distributed_hydro_mpi_direct"),
      "test_reference_workflow_distributed_hydro_mpi_direct");
  cosmosim::workflows::ReferenceWorkflowRunner direct_runner(direct_frozen);
  const auto direct_report = direct_runner.run(
      root / ("rank_" + std::to_string(world_rank) + "_direct"),
      cosmosim::workflows::ReferenceWorkflowOptions{
          .write_outputs = true,
          .initial_state_override = &initial_state,
          .max_steps_override = 2});
  assert(direct_report.restart_roundtrip_ok);
  if (direct_report.final_hydro_imported_mpi_ghosts != 1U) {
    std::cerr << "rank=" << world_rank
              << " local_cells=" << direct_report.local_cell_count
              << " global_cells=" << direct_report.global_cell_count
              << " imported=" << direct_report.final_hydro_imported_mpi_ghosts
              << " faces=" << direct_report.final_hydro_remote_interface_faces
              << " stale=" << direct_report.final_hydro_remote_stale_payloads
              << '\n';
  }
  assert(direct_report.final_hydro_imported_mpi_ghosts == 1U);
  assert(direct_report.final_hydro_remote_interface_faces == 1U);
  assert(direct_report.final_hydro_remote_stale_payloads == 0U);
  const auto direct_restart = cosmosim::io::readRestartCheckpointHdf5(direct_report.restart_path);
  const double direct_global_mass = cosmosim::parallel::MpiContext{}.allreduceSumDouble(localMass(direct_restart.state));

  const auto restart_frozen = cosmosim::core::loadFrozenConfigFromString(
      configText(world_size, "distributed_hydro_mpi_restart"),
      "test_reference_workflow_distributed_hydro_mpi_restart");
  cosmosim::workflows::ReferenceWorkflowRunner restart_runner(restart_frozen);
  const auto first_report = restart_runner.run(
      root / ("rank_" + std::to_string(world_rank) + "_first"),
      cosmosim::workflows::ReferenceWorkflowOptions{
          .write_outputs = true,
          .initial_state_override = &initial_state,
          .max_steps_override = 1});
  assert(first_report.restart_roundtrip_ok);
  const auto first_restart = cosmosim::io::readRestartCheckpointHdf5(first_report.restart_path);
  const auto resumed_report = restart_runner.run(
      root / ("rank_" + std::to_string(world_rank) + "_resumed"),
      cosmosim::workflows::ReferenceWorkflowOptions{
          .write_outputs = true,
          .restart_state_override = &first_restart,
          .max_steps_override = 1});
  assert(resumed_report.restart_roundtrip_ok);
  const auto resumed_restart = cosmosim::io::readRestartCheckpointHdf5(resumed_report.restart_path);
  const double resumed_global_mass = cosmosim::parallel::MpiContext{}.allreduceSumDouble(localMass(resumed_restart.state));
  assert(std::abs(resumed_global_mass - direct_global_mass) < 1.0e-10);

  MPI_Finalize();
#endif
  return 0;
}
