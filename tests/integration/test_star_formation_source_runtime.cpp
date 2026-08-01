#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"

namespace {

[[nodiscard]] std::string buildConfig(std::string_view run_name) {
  std::ostringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[units]\n";
  stream << "length_unit = kpc\n";
  stream << "mass_unit = msun\n";
  stream << "velocity_unit = km_s\n";
  stream << "coordinate_frame = physical\n\n";
  stream << "[mode]\n";
  stream << "mode = isolated_galaxy\n";
  stream << "ic_file = generated\n";
  stream << "zoom_high_res_region = false\n";
  stream << "hydro_boundary = open\n";
  stream << "gravity_boundary = isolated_monopole\n\n";
  stream << "[numerics]\n";
  stream << "time_begin_code = 0.0\n";
  stream << "time_end_code = 0.000000001\n";
  stream << "max_global_steps = 1\n";
  stream << "hierarchical_max_rung = 0\n";
  stream << "treepm_pm_grid = 4\n";
  stream << "treepm_enable_window_deconvolution = false\n";
  stream << "treepm_update_cadence_steps = 1\n\n";
  stream << "[physics]\n";
  stream << "enable_star_formation = true\n";
  stream << "star_formation_model = adaptive_bound_jeans\n";
  stream << "sf_epsilon_ff = 1.0\n";
  stream << "sf_bound_alpha_vir_max = 1.0\n";
  stream << "sf_require_converging_flow = true\n";
  stream << "sf_collapse_timescale = minimum_free_fall_or_compression\n";
  stream << "sf_jeans_mass_floor_code = 100.0\n";
  stream << "sf_star_particle_mass_policy = fixed\n";
  stream << "sf_target_star_particle_mass_code = 100.0\n";
  stream << "sf_min_star_particle_mass_code = 1.0\n";
  stream << "sf_max_star_particle_mass_code = 1000000.0\n";
  stream << "sf_max_spawn_particles_per_cell_step = 8\n";
  stream << "sf_max_fractional_mass_conversion = 0.25\n";
  stream << "sf_min_remaining_gas_fraction = 0.01\n";
  stream << "sf_min_remaining_gas_mass_code = 1.0\n";
  stream << "sf_temperature_safety_ceiling_k = 0.0\n";
  stream << "sf_stochastic_spawning = false\n";
  stream << "sf_random_seed = 987654321\n\n";
  stream << "[output]\n";
  stream << "run_name = " << run_name << "\n";
  stream << "output_directory = integration_outputs\n";
  stream << "output_stem = snapshot\n";
  stream << "restart_stem = restart\n";
  stream << "snapshot_interval_steps = 1\n";
  stream << "write_restarts = true\n\n";
  stream << "[parallel]\n";
  stream << "mpi_ranks_expected = 1\n";
  stream << "omp_threads = 1\n";
  stream << "gpu_devices = 0\n";
  stream << "deterministic_reduction = true\n";
  return stream.str();
}

[[nodiscard]] cosmosim::core::SimulationState makeState(bool converging) {
  cosmosim::core::SimulationState state;
  state.resizeParticles(3);
  state.resizeCells(3);
  state.resizePatches(1);

  state.patches.patch_id[0] = 7001;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 3;
  state.patches.owning_rank[0] = 0;
  state.patches.origin_x_comoving[0] = 0.0;
  state.patches.origin_y_comoving[0] = 0.0;
  state.patches.origin_z_comoving[0] = 0.0;
  state.patches.extent_x_comoving[0] = 0.003;
  state.patches.extent_y_comoving[0] = 0.001;
  state.patches.extent_z_comoving[0] = 0.001;
  state.patches.cell_dim_x[0] = 3;
  state.patches.cell_dim_y[0] = 1;
  state.patches.cell_dim_z[0] = 1;

  const std::array<double, 3> velocity_x = converging
      ? std::array<double, 3>{1.0, 0.0, -1.0}
      : std::array<double, 3>{-1.0, 0.0, 1.0};
  std::vector<cosmosim::core::GasCellIdentityRecord> identity_records;
  identity_records.reserve(3);
  for (std::size_t row = 0; row < 3; ++row) {
    const std::uint64_t particle_id = 5001 + row;
    state.particles.position_x_comoving[row] = 0.0005 + 0.001 * row;
    state.particles.position_y_comoving[row] = 0.0005;
    state.particles.position_z_comoving[row] = 0.0005;
    state.particles.velocity_x_peculiar[row] = velocity_x[row];
    state.particles.mass_code[row] = 1.0e6;
    state.particle_sidecar.particle_id[row] = particle_id;
    state.particle_sidecar.sfc_key[row] = particle_id;
    state.particle_sidecar.species_tag[row] =
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
    state.particle_sidecar.owning_rank[row] = 0;

    state.cells.center_x_comoving[row] = state.particles.position_x_comoving[row];
    state.cells.center_y_comoving[row] = state.particles.position_y_comoving[row];
    state.cells.center_z_comoving[row] = state.particles.position_z_comoving[row];
    state.cells.mass_code[row] = 1.0e6;
    state.cells.patch_index[row] = 0;
    state.cells.time_bin[row] = 0;
    state.gas_cells.gas_cell_id[row] = 8001 + row;
    state.gas_cells.parent_particle_id[row] = particle_id;
    state.gas_cells.density_code[row] = 1.0e15;
    state.gas_cells.pressure_code[row] = 1.0e3;
    state.gas_cells.internal_energy_code[row] = 1.0e-3;
    state.gas_cells.temperature_code[row] = 100.0;
    state.gas_cells.sound_speed_code[row] = 1.0e-3;
    state.gas_cells.velocity_x_peculiar[row] = velocity_x[row];
    state.gas_cells.metal_mass_code[row] = 2.0e4;
    identity_records.push_back({
        .gas_cell_id = 8001 + row,
        .parent_particle_id = particle_id,
        .owning_patch_id = 7001,
        .local_cell_row = static_cast<std::uint32_t>(row),
    });
  }
  state.species.count_by_species.fill(0);
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)] = 3;
  state.rebuildSpeciesIndex();
  state.replaceGasCellIdentityRecords(std::move(identity_records));
  return state;
}

[[nodiscard]] cosmosim::workflows::ReferenceWorkflowReport runCase(
    const std::filesystem::path& output_root,
    std::string_view run_name,
    const cosmosim::core::SimulationState& state) {
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      buildConfig(run_name), "test_star_formation_source_runtime");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  return runner.run(
      output_root,
      cosmosim::workflows::ReferenceWorkflowOptions{
          .dt_time_code = 1.0e-9,
          .write_outputs = COSMOSIM_ENABLE_HDF5 != 0,
          .initial_state_override = &state,
          .max_steps_override = 1,
      });
}

}  // namespace

int main() {
  const auto output_root = std::filesystem::temp_directory_path() /
      "cosmosim_star_formation_source_runtime";
  std::filesystem::remove_all(output_root);

  const auto converging_state = makeState(true);
  const auto expanding_state = makeState(false);
  const auto converging_report = runCase(output_root, "sf_runtime_converging", converging_state);
  const auto expanding_report = runCase(output_root, "sf_runtime_expanding", expanding_state);

  assert(converging_report.completed_steps == 1);
  assert(expanding_report.completed_steps == 1);
  assert(converging_report.local_particle_count > 3);
  assert(expanding_report.local_particle_count == 3);

#if COSMOSIM_ENABLE_HDF5
  assert(converging_report.restart_roundtrip_ok);
  const auto restart = cosmosim::io::readRestartCheckpointHdf5(converging_report.restart_path);
  assert(restart.state.star_particles.size() == converging_report.local_particle_count - 3);
  double gas_particle_mass_code = 0.0;
  double stellar_particle_mass_code = 0.0;
  for (std::size_t particle_row = 0; particle_row < restart.state.particles.size(); ++particle_row) {
    const auto species = static_cast<cosmosim::core::ParticleSpecies>(
        restart.state.particle_sidecar.species_tag[particle_row]);
    if (species == cosmosim::core::ParticleSpecies::kGas) {
      gas_particle_mass_code += restart.state.particles.mass_code[particle_row];
    } else if (species == cosmosim::core::ParticleSpecies::kStar) {
      stellar_particle_mass_code += restart.state.particles.mass_code[particle_row];
    }
  }
  double gas_cell_mass_code = 0.0;
  for (const double mass_code : restart.state.cells.mass_code) {
    gas_cell_mass_code += mass_code;
  }
  assert(std::abs(gas_particle_mass_code - gas_cell_mass_code) < 1.0e-8);
  assert(stellar_particle_mass_code > 0.0);
  assert(std::abs(
      gas_particle_mass_code + stellar_particle_mass_code -
      (gas_cell_mass_code + stellar_particle_mass_code)) < 1.0e-8);
  for (std::size_t star_row = 0; star_row < restart.state.star_particles.size(); ++star_row) {
    assert(restart.state.star_particles.parent_gas_cell_id[star_row] >= 8001);
    assert(restart.state.star_particles.parent_gas_cell_id[star_row] <= 8003);
    assert(restart.state.star_particles.birth_key[star_row] != 0);
    assert(std::abs(restart.state.star_particles.metallicity_mass_fraction[star_row] - 0.02) < 1.0e-10);
  }
#endif

  std::filesystem::remove_all(output_root);
  return 0;
}
