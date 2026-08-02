#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
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

std::string configText(int world_size, std::string_view run_name) {
  std::ostringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[units]\nlength_unit = kpc\nmass_unit = msun\nvelocity_unit = km_s\n"
            "coordinate_frame = physical\n\n";
  stream << "[mode]\nmode = cosmo_cube\nic_file = generated\nhydro_boundary = periodic\n"
            "gravity_boundary = periodic\n\n";
  stream << "[cosmology]\nbox_size_x = 1.0\nbox_size_y = 1.0\nbox_size_z = 1.0\n\n";
  stream << "[numerics]\ntime_begin_code = 0.01\ntime_end_code = 0.010000001\n"
            "max_global_steps = 1\nhierarchical_max_rung = 0\ntreepm_pm_grid = 8\n"
            "treepm_asmth_cells = 1.0\ntreepm_rcut_cells = 3.0\n"
            "treepm_enable_window_deconvolution = false\ntreepm_update_cadence_steps = 1\n\n";
  stream << "[physics]\nenable_star_formation = true\n"
            "star_formation_model = effective_multiphase_tng_like\n"
            "sf_effective_t_star_at_threshold = 1.0e-9 code\n"
            "sf_effective_n_h_threshold = 0.13 cm^-3\n"
            "sf_star_particle_mass_policy = fixed\n"
            "sf_target_star_particle_mass_code = 100.0\n"
            "sf_min_star_particle_mass_code = 1.0\n"
            "sf_max_star_particle_mass_code = 1000000.0\n"
            "sf_max_spawn_particles_per_cell_step = 1\n"
            "sf_max_fractional_mass_conversion = 0.25\n"
            "sf_min_remaining_gas_fraction = 0.01\n"
            "sf_min_remaining_gas_mass_code = 1.0\n"
            "sf_stochastic_spawning = false\n"
            "sf_random_seed = 424242\n\n";
  stream << "[output]\nrun_name = " << run_name
         << "\noutput_directory = integration_outputs\noutput_stem = snapshot\n"
            "restart_stem = restart\nsnapshot_interval_steps = 1\nwrite_restarts = true\n\n";
  stream << "[parallel]\nmpi_ranks_expected = " << world_size
         << "\ndecomposition_runtime_rebalance_enabled = false\n"
            "deterministic_reduction = true\n";
  return stream.str();
}

cosmosim::core::SimulationState makeGlobalState() {
  namespace core = cosmosim::core;
  core::SimulationState state;
  state.resizeParticles(2);
  state.resizeCells(2);
  state.resizePatches(2);
  state.metadata.run_name = "star_formation_source_runtime_mpi";
  state.metadata.snapshot_stem = "snapshot";
  state.metadata.restart_stem = "restart";

  std::vector<core::GasCellIdentityRecord> records;
  for (std::uint32_t row = 0; row < 2U; ++row) {
    const std::uint32_t owner = row;
    const std::uint64_t particle_id = 1001U + row;
    const std::uint64_t gas_cell_id = 2001U + row;
    state.patches.patch_id[row] = 3001U + row;
    state.patches.level[row] = 0U;
    state.patches.first_cell[row] = row;
    state.patches.cell_count[row] = 1U;
    state.patches.owning_rank[row] = owner;
    state.patches.origin_x_comoving[row] = 0.5 * static_cast<double>(row);
    state.patches.origin_y_comoving[row] = 0.0;
    state.patches.origin_z_comoving[row] = 0.0;
    state.patches.extent_x_comoving[row] = 0.5;
    state.patches.extent_y_comoving[row] = 1.0;
    state.patches.extent_z_comoving[row] = 1.0;
    state.patches.cell_dim_x[row] = 1U;
    state.patches.cell_dim_y[row] = 1U;
    state.patches.cell_dim_z[row] = 1U;

    state.particles.position_x_comoving[row] = 0.25 + 0.5 * static_cast<double>(row);
    state.particles.position_y_comoving[row] = 0.5;
    state.particles.position_z_comoving[row] = 0.5;
    state.particles.mass_code[row] = 1.0e6;
    state.particle_sidecar.particle_id[row] = particle_id;
    state.particle_sidecar.sfc_key[row] = particle_id;
    state.particle_sidecar.species_tag[row] =
        static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    state.particle_sidecar.owning_rank[row] = owner;

    state.cells.center_x_comoving[row] = state.particles.position_x_comoving[row];
    state.cells.center_y_comoving[row] = 0.5;
    state.cells.center_z_comoving[row] = 0.5;
    state.cells.mass_code[row] = 1.0e6;
    state.cells.patch_index[row] = row;
    state.cells.time_bin[row] = 0U;
    state.gas_cells.gas_cell_id[row] = gas_cell_id;
    state.gas_cells.parent_particle_id[row] = particle_id;
    state.gas_cells.density_code[row] = row == 0U ? 1.0e15 : 1.0e-20;
    state.gas_cells.pressure_code[row] = 1.0e3;
    state.gas_cells.internal_energy_code[row] = 1.0e-3;
    state.gas_cells.temperature_code[row] = 100.0;
    state.gas_cells.sound_speed_code[row] = 1.0e-3;
    state.gas_cells.metal_mass_code[row] = 2.0e4;
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = gas_cell_id,
        .parent_particle_id = particle_id,
        .owning_patch_id = 3001U + row,
        .local_cell_row = row,
    });
  }
  state.species.count_by_species.fill(0U);
  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)] = 2U;
  state.rebuildSpeciesIndex();
  state.replaceGasCellIdentityRecords(std::move(records));
  assert(state.validateOwnershipInvariants());
  return state;
}

std::uint64_t localStarCount(const cosmosim::core::SimulationState& state) {
  return static_cast<std::uint64_t>(state.star_particles.size());
}

std::string sidecarText(const cosmosim::core::SimulationState& state, const std::string& name) {
  const auto* block = state.sidecars.find(name);
  if (block == nullptr) return {};
  return std::string(reinterpret_cast<const char*>(block->payload.data()), block->payload.size());
}

}  // namespace

int main(int argc, char** argv) {
#if COSMOSIM_ENABLE_MPI && COSMOSIM_ENABLE_HDF5
  MPI_Init(&argc, &argv);
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    MPI_Finalize();
    return 0;
  }
  const std::string mode = argc > 1 ? argv[1] : "source_runtime";
  const auto initial_state = makeGlobalState();
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      configText(world_size, "star_formation_source_runtime_mpi_" + mode),
      "test_star_formation_source_runtime_mpi");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  const auto root = std::filesystem::temp_directory_path() /
      ("cosmosim_star_formation_source_runtime_mpi_" + mode);
  const auto report = runner.run(
      root / ("rank_" + std::to_string(world_rank)),
      cosmosim::workflows::ReferenceWorkflowOptions{
          .dt_time_code = 1.0e-9,
          .write_outputs = true,
          .initial_state_override = &initial_state,
          .max_steps_override = 1});
  assert(report.restart_roundtrip_ok);
  const auto restart = cosmosim::io::readRestartCheckpointHdf5(report.restart_path);
  const std::uint64_t global_stars = cosmosim::parallel::MpiContext{}.allreduceSumUint64(
      localStarCount(restart.state));
  assert(global_stars == 1U);
  assert(restart.state.validateUniqueParticleIds());

  const std::string eos_sidecar = sidecarText(restart.state, "effective_multiphase_ism");
  assert(eos_sidecar.find("table_hash=") != std::string::npos);
  const std::uint64_t local_hash = cosmosim::core::stableConfigHash(eos_sidecar);
  const std::uint64_t hash_xor = cosmosim::parallel::MpiContext{}.allreduceXorUint64(local_hash);
  assert(hash_xor == 0U);

  if (mode == "restart") {
    const auto resumed = runner.run(
        root / ("rank_" + std::to_string(world_rank) + "_resumed"),
        cosmosim::workflows::ReferenceWorkflowOptions{
            .dt_time_code = 1.0e-9,
            .write_outputs = true,
            .restart_state_override = &restart,
            .max_steps_override = 1});
    assert(resumed.restart_roundtrip_ok);
  }
  MPI_Finalize();
#else
  (void)argc;
  (void)argv;
#endif
  return 0;
}
