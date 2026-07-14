#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "cosmosim/cosmosim.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

std::string configText(int world_size) {
  std::ostringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[mode]\n";
  stream << "mode = cosmo_cube\n";
  stream << "ic_file = generated\n\n";
  stream << "[numerics]\n";
  stream << "time_begin_code = 0.01\n";
  stream << "time_end_code = 0.0103\n";
  stream << "max_global_steps = 2\n";
  stream << "hierarchical_max_rung = 0\n";
  stream << "treepm_pm_grid = 8\n";
  stream << "treepm_asmth_cells = 1.25\n";
  // This ownership test uses a deliberately small mesh; keep r_cut < L/2.
  stream << "treepm_rcut_cells = 3.9\n";
  stream << "treepm_update_cadence_steps = 1\n\n";
  stream << "[output]\n";
  stream << "run_name = distributed_gas_cell_migration_mpi\n";
  stream << "output_directory = integration_outputs\n";
  stream << "output_stem = snapshot\n";
  stream << "restart_stem = restart\n";
  stream << "snapshot_interval_steps = 1\n";
  stream << "write_restarts = true\n\n";
  stream << "[parallel]\n";
  stream << "mpi_ranks_expected = " << world_size << "\n";
  stream << "decomposition_runtime_rebalance_enabled = true\n";
  stream << "decomposition_rebalance_imbalance_trigger = 1.0\n";
  stream << "decomposition_rebalance_memory_trigger = 1.0\n";
  stream << "decomposition_rebalance_max_migrated_load_fraction = 1.0\n";
  return stream.str();
}

cosmosim::core::SimulationState makeState() {
  namespace core = cosmosim::core;
  core::SimulationState state;
  state.resizeParticles(2);
  state.resizeCells(3);
  state.resizePatches(2);
  state.metadata.run_name = "distributed_gas_cell_migration_mpi";
  state.metadata.snapshot_stem = "snapshot";
  state.metadata.restart_stem = "restart";

  state.species.count_by_species.fill(0U);
  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)] = 1U;
  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)] = 1U;

  state.particle_sidecar.particle_id = {5001U, 6001U};
  state.particle_sidecar.sfc_key = {10U, 90U};
  state.particle_sidecar.species_tag = {
      static_cast<std::uint32_t>(core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter)};
  state.particle_sidecar.owning_rank = {0U, 0U};
  state.particles.position_x_comoving = {0.75, 0.25};
  state.particles.position_y_comoving = {0.50, 0.50};
  state.particles.position_z_comoving = {0.50, 0.50};
  state.particles.mass_code = {2.0, 4.0};
  state.rebuildSpeciesIndex();

  state.patches.patch_id = {9001U, 9002U};
  state.patches.level = {0, 0};
  state.patches.first_cell = {0U, 1U};
  state.patches.cell_count = {1U, 2U};
  state.patches.owning_rank = {0U, 0U};
  state.patches.origin_x_comoving = {0.0, 0.65};
  state.patches.origin_y_comoving = {0.45, 0.45};
  state.patches.origin_z_comoving = {0.45, 0.45};
  state.patches.extent_x_comoving = {0.1, 0.2};
  state.patches.extent_y_comoving = {0.1, 0.1};
  state.patches.extent_z_comoving = {0.1, 0.1};
  state.patches.cell_dim_x = {1U, 2U};
  state.patches.cell_dim_y = {1U, 1U};
  state.patches.cell_dim_z = {1U, 1U};

  struct CellSeed {
    std::uint64_t gas_cell_id;
    std::optional<std::uint64_t> parent_particle_id;
    std::uint64_t patch_id;
    std::uint32_t patch_index;
    double x;
    double density;
  };
  const std::array<CellSeed, 3> cells{{
      {8001U, std::nullopt, 9001U, 0U, 0.05, 1.0},
      {8002U, 5001U, 9002U, 1U, 0.70, 1.2},
      {8003U, 5001U, 9002U, 1U, 0.80, 0.9},
  }};
  std::vector<core::GasCellIdentityRecord> records;
  records.reserve(cells.size());
  for (std::uint32_t row = 0; row < cells.size(); ++row) {
    const CellSeed& seed = cells[row];
    state.cells.center_x_comoving[row] = seed.x;
    state.cells.center_y_comoving[row] = 0.5;
    state.cells.center_z_comoving[row] = 0.5;
    state.cells.mass_code[row] = seed.density * 0.1;
    state.cells.patch_index[row] = seed.patch_index;
    state.gas_cells.density_code[row] = seed.density;
    state.gas_cells.pressure_code[row] = 1.0 + 0.1 * static_cast<double>(row);
    state.gas_cells.internal_energy_code[row] = 2.5;
    state.gas_cells.temperature_code[row] = 1.0;
    state.gas_cells.sound_speed_code[row] = 1.0;
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = seed.gas_cell_id,
        .parent_particle_id = seed.parent_particle_id,
        .owning_patch_id = seed.patch_id,
        .local_cell_row = row,
    });
  }
  state.replaceGasCellIdentityRecords(std::move(records));
  assert(state.validateOwnershipInvariants());
  return state;
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

  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      configText(world_size),
      "test_distributed_gas_cell_migration_mpi");
  const auto initial_state = makeState();
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  const std::filesystem::path root =
      std::filesystem::temp_directory_path() / "cosmosim_distributed_gas_cell_migration_mpi";
  const auto first_report = runner.run(
      root / ("rank_" + std::to_string(world_rank) + "_first"),
      cosmosim::workflows::ReferenceWorkflowOptions{
          .write_outputs = true,
          .initial_state_override = &initial_state,
          .max_steps_override = 1});
  assert(first_report.restart_roundtrip_ok);
  const auto first_restart = cosmosim::io::readRestartCheckpointHdf5(first_report.restart_path);
  first_restart.state.requireGasCellIdentityMapCoversDenseRows("distributed gas-cell MPI first restart");

  const bool has_parentless = first_restart.state.gas_cell_identity.findByGasCellId(8001U) != nullptr;
  const bool has_shared_parent =
      first_restart.state.gas_cell_identity.rowsForParentParticleId(5001U).size() == 2U;
  const std::uint64_t global_parentless = cosmosim::parallel::MpiContext{}.allreduceSumUint64(has_parentless ? 1U : 0U);
  const std::uint64_t global_shared_parent =
      cosmosim::parallel::MpiContext{}.allreduceSumUint64(has_shared_parent ? 1U : 0U);
  assert(global_parentless == 1U);
  assert(global_shared_parent == 1U);

  const auto resumed_report = runner.run(
      root / ("rank_" + std::to_string(world_rank) + "_resumed"),
      cosmosim::workflows::ReferenceWorkflowOptions{
          .write_outputs = true,
          .restart_state_override = &first_restart,
          .max_steps_override = 1});
  assert(resumed_report.restart_roundtrip_ok);
  const auto resumed_restart = cosmosim::io::readRestartCheckpointHdf5(resumed_report.restart_path);
  resumed_restart.state.requireGasCellIdentityMapCoversDenseRows("distributed gas-cell MPI resumed restart");

  const std::uint64_t local_cell_count = resumed_restart.state.cells.size();
  const std::uint64_t global_cell_count = cosmosim::parallel::MpiContext{}.allreduceSumUint64(local_cell_count);
  assert(global_cell_count == 3U);
  MPI_Finalize();
#endif
  return 0;
}
