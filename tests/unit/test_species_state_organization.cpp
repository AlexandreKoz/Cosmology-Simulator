#include <cassert>
#include <cstddef>
#include <cstdint>

#include "cosmosim/core/simulation_state.hpp"

int main() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(6);

  const std::array<cosmosim::core::ParticleSpecies, 6> species = {
      cosmosim::core::ParticleSpecies::kDarkMatter,
      cosmosim::core::ParticleSpecies::kGas,
      cosmosim::core::ParticleSpecies::kStar,
      cosmosim::core::ParticleSpecies::kStar,
      cosmosim::core::ParticleSpecies::kBlackHole,
      cosmosim::core::ParticleSpecies::kTracer,
  };

  for (std::size_t i = 0; i < species.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 2000 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(species[i]);
    state.particle_sidecar.owning_rank[i] = 0;
    state.particles.position_x_comoving[i] = static_cast<double>(i);
    state.particles.position_y_comoving[i] = 0.0;
    state.particles.position_z_comoving[i] = 0.0;
    state.particles.velocity_x_peculiar[i] = 0.0;
    state.particles.velocity_y_peculiar[i] = 0.0;
    state.particles.velocity_z_peculiar[i] = 0.0;
    state.particles.mass_code[i] = 1.0;
    state.particles.time_bin[i] = 0;
  }

  state.species.count_by_species = {1, 1, 2, 1, 1};
  state.star_particles.resize(2);
  state.star_particles.particle_index = {2, 3};
  state.star_particles.formation_scale_factor = {0.5, 0.6};
  state.star_particles.birth_mass_code = {1.0, 1.0};
  state.star_particles.metallicity_mass_fraction = {0.01, 0.02};

  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 4;
  state.black_holes.subgrid_mass_code[0] = 5.0;
  state.black_holes.accretion_rate_code[0] = 0.2;
  state.black_holes.feedback_energy_code[0] = 1.0;

  state.tracers.resize(1);
  state.tracers.particle_index[0] = 5;
  state.tracers.parent_particle_id[0] = 2004;
  state.tracers.injection_step[0] = 11;

  state.rebuildSpeciesIndex();
  assert(state.particle_species_index.count(cosmosim::core::ParticleSpecies::kStar) == 2);
  assert(state.particle_species_index.globalIndex(cosmosim::core::ParticleSpecies::kStar, 1) == 3);
  assert(state.particle_species_index.localIndex(3) == 1);
  assert(state.validateOwnershipInvariants());
  assert(state.validateUniqueParticleIds());

  state.particle_sidecar.particle_id[5] = 2004;
  assert(!state.validateUniqueParticleIds());

  return 0;
}
