#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"

int main() {
  constexpr std::size_t k_star_count = 1u << 20;
  constexpr std::size_t k_iterations = 6;

  cosmosim::core::SimulationState state;
  state.resizeParticles(k_star_count);
  state.star_particles.resize(k_star_count);
  for (std::size_t i = 0; i < k_star_count; ++i) {
    state.particles.mass_code[i] = 1.0;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar);
    state.star_particles.particle_index[i] = static_cast<std::uint32_t>(i);
    state.star_particles.formation_scale_factor[i] = 0.1 + 0.8 * static_cast<double>(i % 1024) / 1024.0;
    state.star_particles.birth_mass_code[i] = 1.0;
    state.star_particles.metallicity_mass_fraction[i] = 0.01;
  }

  std::vector<std::uint32_t> active(k_star_count);
  for (std::size_t i = 0; i < k_star_count; ++i) {
    active[i] = static_cast<std::uint32_t>(i);
  }

  cosmosim::physics::StellarEvolutionConfig config;
  cosmosim::physics::StellarEvolutionBookkeeper bookkeeper(
      config,
      cosmosim::physics::StellarEvolutionTable::makeBuiltinReference());

  const auto setup_start = std::chrono::steady_clock::now();
  auto warm = bookkeeper.apply(state, active, 0.95, 1.0e-4);
  const auto setup_end = std::chrono::steady_clock::now();

  double returned_mass = warm.counters.returned_mass_code;
  const auto steady_start = std::chrono::steady_clock::now();
  for (std::size_t iter = 0; iter < k_iterations; ++iter) {
    const auto report = bookkeeper.apply(state, active, 0.95, 1.0e-4);
    returned_mass += report.counters.returned_mass_code;
  }
  const auto steady_end = std::chrono::steady_clock::now();

  const double setup_ms = std::chrono::duration<double, std::milli>(setup_end - setup_start).count();
  const double steady_ms = std::chrono::duration<double, std::milli>(steady_end - steady_start).count();
  const double steady_s = std::max(steady_ms * 1.0e-3, 1.0e-12);
  const double star_updates = static_cast<double>(k_star_count * k_iterations);

  std::cout << "bench_stellar_evolution_bookkeeping"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " hardware=cpu"
            << " threads=1"
            << " features=stellar_evolution_bookkeeping"
            << " setup_ms=" << setup_ms
            << " steady_ms=" << steady_ms
            << " stars=" << k_star_count
            << " iterations=" << k_iterations
            << " star_updates_per_s=" << (star_updates / steady_s)
            << " effective_read_write_bandwidth_gb_s="
            << ((star_updates * 8.0 * sizeof(double)) / steady_s * 1.0e-9)
            << " returned_mass_code=" << returned_mass << '\n';

  return 0;
}
