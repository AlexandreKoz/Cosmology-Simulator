#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

#include "cosmosim/core/simulation_state.hpp"

int main() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(5);
  state.resizeCells(3);
  state.resizePatches(1);

  for (std::size_t i = 0; i < 5; ++i) {
    state.particle_sidecar.particle_id[i] = 1000 + i;
    state.particle_sidecar.species_tag[i] =
        static_cast<std::uint32_t>((i < 2) ? cosmosim::core::ParticleSpecies::kDarkMatter
                                           : cosmosim::core::ParticleSpecies::kStar);
    state.particles.position_x_comoving[i] = static_cast<double>(i);
    state.particles.position_y_comoving[i] = static_cast<double>(i + 1);
    state.particles.position_z_comoving[i] = static_cast<double>(i + 2);
    state.particles.velocity_x_peculiar[i] = 0.1;
    state.particles.velocity_y_peculiar[i] = 0.2;
    state.particles.velocity_z_peculiar[i] = 0.3;
    state.particles.mass_code[i] = 1.0;
    state.particles.time_bin[i] = static_cast<std::uint8_t>(i % 3);
  }
  state.species.count_by_species = {2, 0, 3, 0, 0};

  state.star_particles.resize(3);
  for (std::size_t i = 0; i < 3; ++i) {
    state.star_particles.particle_index[i] = static_cast<std::uint32_t>(i + 2);
    state.star_particles.formation_scale_factor[i] = 0.5;
    state.star_particles.birth_mass_code[i] = 1.1;
    state.star_particles.metallicity_mass_fraction[i] = 0.02;
  }

  state.patches.patch_id[0] = 42;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 3;

  for (std::size_t i = 0; i < 3; ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.center_y_comoving[i] = static_cast<double>(i) + 0.5;
    state.cells.center_z_comoving[i] = static_cast<double>(i) + 1.0;
    state.cells.mass_code[i] = 3.0;
    state.cells.time_bin[i] = 1;
    state.cells.patch_index[i] = 0;

    state.gas_cells.density_code[i] = 10.0 + static_cast<double>(i);
    state.gas_cells.pressure_code[i] = 1.0 + static_cast<double>(i);
    state.gas_cells.internal_energy_code[i] = 0.1 + static_cast<double>(i);
    state.gas_cells.temperature_code[i] = 100.0 + static_cast<double>(i);
    state.gas_cells.sound_speed_code[i] = 2.0;
    state.gas_cells.recon_gradient_x[i] = 0.01;
    state.gas_cells.recon_gradient_y[i] = 0.02;
    state.gas_cells.recon_gradient_z[i] = 0.03;
  }

  state.rebuildSpeciesIndex();
  assert(state.validateOwnershipInvariants());
  assert(state.validateUniqueParticleIds());

  state.metadata.run_name = "unit_state";
  state.metadata.normalized_config_hash = 1234;
  state.metadata.normalized_config_hash_hex = "0x4d2";
  state.metadata.step_index = 17;
  state.metadata.scale_factor = 0.5;

  const std::string serialized = state.metadata.serialize();
  const auto parsed = cosmosim::core::StateMetadata::deserialize(serialized);
  assert(parsed.run_name == state.metadata.run_name);
  assert(parsed.normalized_config_hash == state.metadata.normalized_config_hash);
  assert(parsed.step_index == state.metadata.step_index);
  assert(parsed.scale_factor == state.metadata.scale_factor);

  const auto packet = state.packSpeciesTransferPacket(cosmosim::core::ParticleSpecies::kStar);
  assert(packet.particle_id.size() == 3);
  assert(packet.particle_id[0] == 1002);

  cosmosim::core::TransientStepWorkspace workspace;
  const std::array<std::uint32_t, 2> particle_indices{1, 3};
  auto particle_view = cosmosim::core::buildParticleActiveView(state, particle_indices, workspace);
  assert(particle_view.size() == 2);
  assert(particle_view.particle_id[0] == 1001);
  assert(particle_view.position_x_comoving[1] == 3.0);

  const std::array<std::uint32_t, 1> cell_indices{2};
  auto cell_view = cosmosim::core::buildCellActiveView(state, cell_indices, workspace);
  assert(cell_view.size() == 1);
  assert(cell_view.density_code[0] == 12.0);

  bool threw = false;
  try {
    const std::array<std::uint32_t, 1> bad_indices{8};
    (void)cosmosim::core::buildParticleActiveView(state, bad_indices, workspace);
  } catch (const std::out_of_range&) {
    threw = true;
  }
  assert(threw);

  return 0;
}
