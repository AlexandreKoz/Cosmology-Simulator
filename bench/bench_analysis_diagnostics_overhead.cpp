#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/core/time_integration.hpp"

int main() {
  constexpr std::size_t k_particle_count = 2000;
  constexpr std::size_t k_cell_count = 2000;
  constexpr int k_step_count = 16;

  cosmosim::core::SimulationConfig config;
  config.output.output_directory = "bench_outputs";
  config.output.run_name = "analysis_overhead";
  config.cosmology.box_size_mpc_comoving = 10.0;
  config.analysis.enable_diagnostics = true;
  config.analysis.run_health_interval_steps = 1;
  config.analysis.science_light_interval_steps = 1;
  config.analysis.science_heavy_interval_steps = 4;
  config.analysis.retention_bundle_count = 128;
  config.analysis.quicklook_grid_n = 16;
  config.analysis.power_spectrum_mesh_n = 8;
  config.analysis.power_spectrum_bin_count = 8;

  const std::filesystem::path output_root =
      std::filesystem::path(config.output.output_directory) / config.output.run_name;
  std::filesystem::remove_all(output_root);

  cosmosim::core::SimulationState state;
  state.resizeParticles(k_particle_count);
  state.resizeCells(k_cell_count);
  state.species.count_by_species = {1000, 600, 400, 0, 0};

  for (std::size_t i = 0; i < k_particle_count; ++i) {
    state.particle_sidecar.particle_id[i] = 100000 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(
        (i < 1000) ? cosmosim::core::ParticleSpecies::kDarkMatter
                   : (i < 1600) ? cosmosim::core::ParticleSpecies::kGas : cosmosim::core::ParticleSpecies::kStar);
    state.particles.mass_code[i] = 1.0;
    state.particles.position_x_comoving[i] = static_cast<double>(i % 97) / 97.0 * config.cosmology.box_size_mpc_comoving;
    state.particles.position_y_comoving[i] = static_cast<double>(i % 89) / 89.0 * config.cosmology.box_size_mpc_comoving;
    state.particles.position_z_comoving[i] = static_cast<double>(i % 83) / 83.0 * config.cosmology.box_size_mpc_comoving;
    state.particles.velocity_x_peculiar[i] = 0.01;
    state.particles.velocity_y_peculiar[i] = 0.02;
    state.particles.velocity_z_peculiar[i] = 0.03;
  }

  for (std::size_t i = 0; i < k_cell_count; ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i % 53) / 53.0 * config.cosmology.box_size_mpc_comoving;
    state.cells.center_y_comoving[i] = static_cast<double>(i % 47) / 47.0 * config.cosmology.box_size_mpc_comoving;
    state.cells.center_z_comoving[i] = static_cast<double>(i % 41) / 41.0 * config.cosmology.box_size_mpc_comoving;
    state.cells.mass_code[i] = 1.0;
    state.gas_cells.density_code[i] = 1.0 + static_cast<double>(i % 32) * 0.01;
  }

  state.star_particles.resize(400);
  for (std::size_t i = 0; i < state.star_particles.size(); ++i) {
    state.star_particles.formation_scale_factor[i] = 0.1 + 0.8 * static_cast<double>(i) / 400.0;
    state.star_particles.birth_mass_code[i] = 1.0;
  }

  state.rebuildSpeciesIndex();

  cosmosim::analysis::DiagnosticsCallback diagnostics(config);

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 1.0e-3;
  integrator_state.current_scale_factor = 0.5;

  std::vector<std::uint32_t> particle_indices(k_particle_count);
  std::vector<std::uint32_t> cell_indices(k_cell_count);
  for (std::uint32_t i = 0; i < k_particle_count; ++i) {
    particle_indices[i] = i;
  }
  for (std::uint32_t i = 0; i < k_cell_count; ++i) {
    cell_indices[i] = i;
  }

  const cosmosim::core::ActiveSetDescriptor active{
      .particle_indices = particle_indices,
      .cell_indices = cell_indices,
      .particles_are_subset = false,
      .cells_are_subset = false,
  };

  const auto t_start = std::chrono::steady_clock::now();
  for (int step = 0; step < k_step_count; ++step) {
    cosmosim::core::StepContext context{
        .state = state,
        .integrator_state = integrator_state,
        .active_set = active,
        .workspace = nullptr,
        .cosmology_background = nullptr,
        .mode_policy = nullptr,
        .stage = cosmosim::core::IntegrationStage::kAnalysisHooks,
    };
    diagnostics.onStage(context);
    ++integrator_state.step_index;
  }
  const auto t_end = std::chrono::steady_clock::now();

  const double total_ms =
      static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count()) / 1000.0;
  const auto timing = diagnostics.timing();

  std::cout << "build_type="
#ifdef NDEBUG
            << "Release";
#else
            << "Debug";
#endif
  std::cout << "\nthreads=1\nfeatures=cpu_serial\n";
  std::cout << "steps=" << k_step_count << "\n";
  std::cout << "particle_count=" << k_particle_count << "\n";
  std::cout << "cell_count=" << k_cell_count << "\n";
  std::cout << "total_wall_ms=" << total_ms << "\n";
  std::cout << "run_health_calls=" << timing.run_health_calls << "\n";
  std::cout << "light_calls=" << timing.light_calls << "\n";
  std::cout << "heavy_calls=" << timing.heavy_calls << "\n";
  std::cout << "run_health_ms=" << timing.cumulative_run_health_ms << "\n";
  std::cout << "light_ms=" << timing.cumulative_light_ms << "\n";
  std::cout << "heavy_ms=" << timing.cumulative_heavy_ms << "\n";
  std::cout << "ms_per_step=" << (total_ms / static_cast<double>(k_step_count)) << "\n";

  std::filesystem::remove_all(output_root);
  return 0;
}
