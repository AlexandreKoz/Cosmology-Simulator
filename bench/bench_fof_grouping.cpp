#include <chrono>
#include <cstdint>
#include <iostream>

#include "cosmosim/analysis/halo_workflow.hpp"

int main() {
  constexpr std::uint32_t k_particle_count = 4096;
  const auto setup_start = std::chrono::steady_clock::now();

  cosmosim::core::SimulationConfig config;
  config.output.run_name = "bench_fof";
  config.cosmology.box_size_mpc_comoving = 50.0;
  config.analysis.halo_fof_linking_length_factor = 0.2;
  config.analysis.halo_fof_min_group_size = 16;

  cosmosim::core::SimulationState state;
  state.resizeParticles(k_particle_count);
  state.species.count_by_species = {k_particle_count, 0, 0, 0, 0};
  state.metadata.normalized_config_hash = 42;

  for (std::uint32_t i = 0; i < k_particle_count; ++i) {
    state.particle_sidecar.particle_id[i] = 100000 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particles.mass_code[i] = 1.0;
    state.particles.position_x_comoving[i] = static_cast<double>(i % 64) / 64.0 * config.cosmology.box_size_mpc_comoving;
    state.particles.position_y_comoving[i] = static_cast<double>((i / 64) % 64) / 64.0 * config.cosmology.box_size_mpc_comoving;
    state.particles.position_z_comoving[i] = static_cast<double>(i % 17) / 17.0 * config.cosmology.box_size_mpc_comoving;
  }

  cosmosim::analysis::FofConfig fof_cfg;
  fof_cfg.linking_length_factor_mean_interparticle = config.analysis.halo_fof_linking_length_factor;
  fof_cfg.min_group_size = static_cast<std::uint64_t>(config.analysis.halo_fof_min_group_size);
  cosmosim::analysis::FofHaloFinder finder(fof_cfg);

  const auto setup_done = std::chrono::steady_clock::now();
  cosmosim::analysis::FofProfilingCounters profiling;
  const auto t0 = std::chrono::steady_clock::now();
  const auto catalog = finder.buildCatalog(state, config, 32, 0.5, &profiling);
  const auto t1 = std::chrono::steady_clock::now();

  const double setup_ms =
      static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(setup_done - setup_start).count()) /
      1000.0;
  const double steady_ms =
      static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) / 1000.0;

  const double particles_per_sec = (steady_ms > 0.0)
                                       ? static_cast<double>(profiling.candidate_particle_count) / (steady_ms * 1.0e-3)
                                       : 0.0;

  std::cout << "build_type="
#ifdef NDEBUG
            << "Release";
#else
            << "Debug";
#endif
  std::cout << "\nthreads=1\nfeatures=cpu_serial\n";
  std::cout << "hardware_hint=generic_ci_vm\n";
  std::cout << "particle_count=" << profiling.candidate_particle_count << "\n";
  std::cout << "pair_checks=" << profiling.pair_checks << "\n";
  std::cout << "pair_links=" << profiling.pair_links << "\n";
  std::cout << "halo_count=" << catalog.halos.size() << "\n";
  std::cout << "setup_cost_ms=" << setup_ms << "\n";
  std::cout << "steady_state_ms=" << steady_ms << "\n";
  std::cout << "throughput_particles_per_sec=" << particles_per_sec << "\n";
  return 0;
}
