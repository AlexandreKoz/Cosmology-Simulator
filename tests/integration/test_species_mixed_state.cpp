#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>

#include "cosmosim/core/simulation_state.hpp"

int main() {

  cosmosim::core::SimulationState state;
  state.resizeParticles(7);
  state.resizeCells(2);
  state.resizePatches(1);

  state.species.count_by_species = {2, 1, 2, 1, 1};

  const std::array<std::uint32_t, 7> species_tags{0, 0, 1, 2, 2, 3, 4};
  for (std::size_t i = 0; i < species_tags.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 3000 + i;
    state.particle_sidecar.species_tag[i] = species_tags[i];
    state.particle_sidecar.owning_rank[i] = static_cast<std::uint32_t>(i % 2);
    state.particles.position_x_comoving[i] = static_cast<double>(i);
    state.particles.position_y_comoving[i] = static_cast<double>(i) + 1.0;
    state.particles.position_z_comoving[i] = static_cast<double>(i) + 2.0;
    state.particles.velocity_x_peculiar[i] = 0.1;
    state.particles.velocity_y_peculiar[i] = 0.2;
    state.particles.velocity_z_peculiar[i] = 0.3;
    state.particles.mass_code[i] = 1.0;
    state.particles.time_bin[i] = 2;
  }

  state.star_particles.resize(2);
  state.star_particles.particle_index = {3, 4};
  state.star_particles.formation_scale_factor = {0.7, 0.71};
  state.star_particles.birth_mass_code = {0.8, 0.9};
  state.star_particles.metallicity_mass_fraction = {0.015, 0.018};

  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 5;
  state.black_holes.subgrid_mass_code[0] = 7.0;
  state.black_holes.accretion_rate_code[0] = 0.05;
  state.black_holes.feedback_energy_code[0] = 0.2;

  state.tracers.resize(1);
  state.tracers.particle_index[0] = 6;
  state.tracers.parent_particle_id[0] = 3005;
  state.tracers.injection_step[0] = 15;

  state.patches.patch_id[0] = 77;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 2;

  for (std::size_t i = 0; i < 2; ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.center_y_comoving[i] = static_cast<double>(i);
    state.cells.center_z_comoving[i] = static_cast<double>(i);
    state.cells.mass_code[i] = 3.0;
    state.cells.time_bin[i] = 0;
    state.cells.patch_index[i] = 0;
    state.gas_cells.density_code[i] = 1.0 + static_cast<double>(i);
    state.gas_cells.pressure_code[i] = 2.0 + static_cast<double>(i);
    state.gas_cells.internal_energy_code[i] = 3.0 + static_cast<double>(i);
    state.gas_cells.temperature_code[i] = 100.0;
    state.gas_cells.sound_speed_code[i] = 1.0;
    state.gas_cells.recon_gradient_x[i] = 0.0;
    state.gas_cells.recon_gradient_y[i] = 0.0;
    state.gas_cells.recon_gradient_z[i] = 0.0;
  }

  state.rebuildSpeciesIndex();
  assert(state.validateOwnershipInvariants());

  const auto star_packet = state.packSpeciesTransferPacket(cosmosim::core::ParticleSpecies::kStar);
  assert(star_packet.particle_id.size() == 2);
  assert(star_packet.particle_id[0] == 3003);

  const auto tracer_packet = state.packSpeciesTransferPacket(cosmosim::core::ParticleSpecies::kTracer);
  assert(tracer_packet.particle_id.size() == 1);
  assert(tracer_packet.particle_id[0] == 3006);

  return 0;
}
