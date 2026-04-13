#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>

#include "cosmosim/core/simulation_state.hpp"

int main() {
  cosmosim::core::SimulationState state;
  constexpr std::size_t k_particle_count = 200000;
  state.resizeParticles(k_particle_count);

  for (std::size_t i = 0; i < k_particle_count; ++i) {
    state.particle_sidecar.particle_id[i] = 1000000 + i;
    state.particle_sidecar.species_tag[i] =
        static_cast<std::uint32_t>((i % 10 == 0) ? cosmosim::core::ParticleSpecies::kStar
                                                 : cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particles.position_x_comoving[i] = static_cast<double>(i) * 0.001;
    state.particles.velocity_x_peculiar[i] = 0.01;
    state.particles.mass_code[i] = 1.0;
    state.particles.time_bin[i] = 1;
  }
  state.species.count_by_species = {180000, 0, 20000, 0, 0};

  state.star_particles.resize(20000);
  std::size_t star_local = 0;
  for (std::size_t i = 0; i < k_particle_count; ++i) {
    if (state.particle_sidecar.species_tag[i] ==
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar)) {
      state.star_particles.particle_index[star_local] = static_cast<std::uint32_t>(i);
      state.star_particles.formation_scale_factor[star_local] = 0.6;
      state.star_particles.birth_mass_code[star_local] = 1.2;
      state.star_particles.metallicity_mass_fraction[star_local] = 0.02;
      ++star_local;
    }
  }

  state.rebuildSpeciesIndex();

  volatile double hot_checksum = 0.0;
  const auto hot_begin = std::chrono::steady_clock::now();
  for (int iter = 0; iter < 20; ++iter) {
    for (std::size_t i = 0; i < k_particle_count; ++i) {
      hot_checksum += state.particles.position_x_comoving[i] + state.particles.mass_code[i];
    }
  }
  const auto hot_end = std::chrono::steady_clock::now();

  volatile double cold_checksum = 0.0;
  const auto cold_begin = std::chrono::steady_clock::now();
  for (int iter = 0; iter < 20; ++iter) {
    for (std::size_t i = 0; i < state.star_particles.size(); ++i) {
      cold_checksum += state.star_particles.formation_scale_factor[i] +
                       state.star_particles.metallicity_mass_fraction[i];
    }
  }
  const auto cold_end = std::chrono::steady_clock::now();

  const auto hot_us =
      std::chrono::duration_cast<std::chrono::microseconds>(hot_end - hot_begin).count();
  const auto cold_us =
      std::chrono::duration_cast<std::chrono::microseconds>(cold_end - cold_begin).count();

  std::cout << "gravity_hot_sweep_us=" << hot_us << '\n';
  std::cout << "species_cold_metadata_us=" << cold_us << '\n';
  std::cout << "checksums=" << hot_checksum << "," << cold_checksum << '\n';

  return 0;
}
