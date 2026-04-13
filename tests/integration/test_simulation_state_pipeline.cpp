#include <cassert>
#include <cstddef>
#include <cstdint>

#include "cosmosim/cosmosim.hpp"

int main() {
  cosmosim::core::SimulationState state;
  state.metadata.run_name = "integration_pipeline";
  state.resizeParticles(8);
  state.resizeCells(4);
  state.resizePatches(1);

  state.species.count_by_species = {2, 2, 2, 1, 1};

  for (std::size_t i = 0; i < 8; ++i) {
    state.particle_sidecar.particle_id[i] = 100000 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(
        (i < 2) ? cosmosim::core::ParticleSpecies::kDarkMatter
                : (i < 4) ? cosmosim::core::ParticleSpecies::kGas
                          : (i < 6) ? cosmosim::core::ParticleSpecies::kStar
                                    : (i == 6) ? cosmosim::core::ParticleSpecies::kBlackHole
                                               : cosmosim::core::ParticleSpecies::kTracer);
    state.particles.position_x_comoving[i] = static_cast<double>(i);
    state.particles.position_y_comoving[i] = static_cast<double>(i) * 2.0;
    state.particles.position_z_comoving[i] = static_cast<double>(i) * 3.0;
    state.particles.velocity_x_peculiar[i] = 0.0;
    state.particles.velocity_y_peculiar[i] = 0.0;
    state.particles.velocity_z_peculiar[i] = 0.0;
    state.particles.mass_code[i] = 1.0;
    state.particles.time_bin[i] = 1;
  }

  state.star_particles.resize(2);
  state.star_particles.particle_index = {4, 5};
  state.star_particles.formation_scale_factor = {0.4, 0.45};
  state.star_particles.birth_mass_code = {1.2, 1.3};
  state.star_particles.metallicity_mass_fraction = {0.01, 0.02};

  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 6;
  state.black_holes.subgrid_mass_code[0] = 10.0;
  state.black_holes.accretion_rate_code[0] = 0.1;
  state.black_holes.feedback_energy_code[0] = 5.0;

  state.tracers.resize(1);
  state.tracers.particle_index[0] = 7;
  state.tracers.parent_particle_id[0] = 100006;
  state.tracers.injection_step[0] = 20;

  state.patches.patch_id[0] = 9000;
  state.patches.level[0] = 1;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 4;

  for (std::size_t i = 0; i < 4; ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i) + 0.25;
    state.cells.center_y_comoving[i] = static_cast<double>(i) + 0.5;
    state.cells.center_z_comoving[i] = static_cast<double>(i) + 0.75;
    state.cells.mass_code[i] = 2.0 + static_cast<double>(i);
    state.cells.time_bin[i] = 1;
    state.cells.patch_index[i] = 0;

    state.gas_cells.density_code[i] = 1.0 + static_cast<double>(i);
    state.gas_cells.pressure_code[i] = 2.0 + static_cast<double>(i);
    state.gas_cells.internal_energy_code[i] = 3.0 + static_cast<double>(i);
    state.gas_cells.temperature_code[i] = 100.0 + static_cast<double>(i);
    state.gas_cells.sound_speed_code[i] = 1.2;
    state.gas_cells.recon_gradient_x[i] = 0.01;
    state.gas_cells.recon_gradient_y[i] = 0.01;
    state.gas_cells.recon_gradient_z[i] = 0.01;
  }

  state.rebuildSpeciesIndex();
  assert(state.validateOwnershipInvariants());

  cosmosim::core::ActiveIndexSet active_set;
  active_set.particle_indices = {0, 2, 4, 6};
  active_set.cell_indices = {1, 3};

  cosmosim::core::TransientStepWorkspace workspace;
  const auto particle_view =
      cosmosim::core::buildParticleActiveView(state, active_set.particle_indices, workspace);
  const auto cell_view = cosmosim::core::buildCellActiveView(state, active_set.cell_indices, workspace);

  double checksum = 0.0;
  for (std::size_t i = 0; i < particle_view.size(); ++i) {
    checksum += particle_view.position_x_comoving[i] + particle_view.position_y_comoving[i] +
                particle_view.position_z_comoving[i];
  }
  for (std::size_t i = 0; i < cell_view.size(); ++i) {
    checksum += cell_view.mass_code[i] + cell_view.density_code[i] + cell_view.pressure_code[i];
  }

  const auto bh_packet = state.packSpeciesTransferPacket(cosmosim::core::ParticleSpecies::kBlackHole);
  checksum += static_cast<double>(bh_packet.particle_id[0]);

  assert(checksum > 100000.0);
  return 0;
}
