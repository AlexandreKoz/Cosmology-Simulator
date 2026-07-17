#include "workflows/internal/output_verification.hpp"

#include <cstddef>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::workflows::internal {
namespace {

[[nodiscard]] bool moduleSidecarsEqualForRestart(
    const core::ModuleSidecarRegistry& lhs,
    const core::ModuleSidecarRegistry& rhs) {
  const auto lhs_blocks = lhs.blocksSortedByName();
  const auto rhs_blocks = rhs.blocksSortedByName();
  if (lhs_blocks.size() != rhs_blocks.size()) {
    return false;
  }
  for (std::size_t i = 0; i < lhs_blocks.size(); ++i) {
    if (lhs_blocks[i]->module_name != rhs_blocks[i]->module_name ||
        lhs_blocks[i]->schema_version != rhs_blocks[i]->schema_version ||
        lhs_blocks[i]->payload != rhs_blocks[i]->payload) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] bool restartRuntimeStateExactlyEquivalentImpl(
    const core::SimulationState& restored,
    const core::SimulationState& reference) {
  return restored.metadata.serialize() == reference.metadata.serialize() &&
      restored.particles.position_x_comoving == reference.particles.position_x_comoving &&
      restored.particles.position_y_comoving == reference.particles.position_y_comoving &&
      restored.particles.position_z_comoving == reference.particles.position_z_comoving &&
      restored.particles.velocity_x_peculiar == reference.particles.velocity_x_peculiar &&
      restored.particles.velocity_y_peculiar == reference.particles.velocity_y_peculiar &&
      restored.particles.velocity_z_peculiar == reference.particles.velocity_z_peculiar &&
      restored.particles.mass_code == reference.particles.mass_code &&
      restored.particles.time_bin == reference.particles.time_bin &&
      restored.particle_sidecar.particle_id == reference.particle_sidecar.particle_id &&
      restored.particle_sidecar.sfc_key == reference.particle_sidecar.sfc_key &&
      restored.particle_sidecar.species_tag == reference.particle_sidecar.species_tag &&
      restored.particle_sidecar.particle_flags == reference.particle_sidecar.particle_flags &&
      restored.particle_sidecar.owning_rank == reference.particle_sidecar.owning_rank &&
      restored.particle_sidecar.gravity_softening_comoving == reference.particle_sidecar.gravity_softening_comoving &&
      restored.particle_sidecar.has_gravity_softening_override == reference.particle_sidecar.has_gravity_softening_override &&
      restored.cells.center_x_comoving == reference.cells.center_x_comoving &&
      restored.cells.center_y_comoving == reference.cells.center_y_comoving &&
      restored.cells.center_z_comoving == reference.cells.center_z_comoving &&
      restored.cells.mass_code == reference.cells.mass_code &&
      restored.cells.time_bin == reference.cells.time_bin &&
      restored.cells.patch_index == reference.cells.patch_index &&
      restored.gas_cells.gas_cell_id == reference.gas_cells.gas_cell_id &&
      restored.gas_cells.parent_particle_id == reference.gas_cells.parent_particle_id &&
      restored.gas_cells.density_code == reference.gas_cells.density_code &&
      restored.gas_cells.pressure_code == reference.gas_cells.pressure_code &&
      restored.gas_cells.internal_energy_code == reference.gas_cells.internal_energy_code &&
      restored.gas_cells.temperature_code == reference.gas_cells.temperature_code &&
      restored.gas_cells.sound_speed_code == reference.gas_cells.sound_speed_code &&
      restored.patches.patch_id == reference.patches.patch_id &&
      restored.patches.level == reference.patches.level &&
      restored.patches.first_cell == reference.patches.first_cell &&
      restored.patches.cell_count == reference.patches.cell_count &&
      restored.patches.owning_rank == reference.patches.owning_rank &&
      restored.star_particles.particle_index == reference.star_particles.particle_index &&
      restored.star_particles.formation_scale_factor == reference.star_particles.formation_scale_factor &&
      restored.star_particles.birth_mass_code == reference.star_particles.birth_mass_code &&
      restored.star_particles.metallicity_mass_fraction == reference.star_particles.metallicity_mass_fraction &&
      restored.star_particles.stellar_age_years_last == reference.star_particles.stellar_age_years_last &&
      restored.star_particles.stellar_returned_mass_cumulative_code == reference.star_particles.stellar_returned_mass_cumulative_code &&
      restored.star_particles.stellar_returned_metals_cumulative_code == reference.star_particles.stellar_returned_metals_cumulative_code &&
      restored.star_particles.stellar_feedback_energy_cumulative_erg == reference.star_particles.stellar_feedback_energy_cumulative_erg &&
      restored.star_particles.stellar_returned_mass_channel_cumulative_code == reference.star_particles.stellar_returned_mass_channel_cumulative_code &&
      restored.star_particles.stellar_returned_metals_channel_cumulative_code == reference.star_particles.stellar_returned_metals_channel_cumulative_code &&
      restored.star_particles.stellar_feedback_energy_channel_cumulative_erg == reference.star_particles.stellar_feedback_energy_channel_cumulative_erg &&
      restored.black_holes.particle_index == reference.black_holes.particle_index &&
      restored.black_holes.host_cell_index == reference.black_holes.host_cell_index &&
      restored.black_holes.subgrid_mass_code == reference.black_holes.subgrid_mass_code &&
      restored.black_holes.accretion_rate_code == reference.black_holes.accretion_rate_code &&
      restored.black_holes.feedback_energy_code == reference.black_holes.feedback_energy_code &&
      restored.black_holes.eddington_ratio == reference.black_holes.eddington_ratio &&
      restored.black_holes.cumulative_accreted_mass_code == reference.black_holes.cumulative_accreted_mass_code &&
      restored.black_holes.cumulative_feedback_energy_code == reference.black_holes.cumulative_feedback_energy_code &&
      restored.black_holes.duty_cycle_active_time_code == reference.black_holes.duty_cycle_active_time_code &&
      restored.black_holes.duty_cycle_total_time_code == reference.black_holes.duty_cycle_total_time_code &&
      restored.tracers.particle_index == reference.tracers.particle_index &&
      restored.tracers.parent_particle_id == reference.tracers.parent_particle_id &&
      restored.tracers.injection_step == reference.tracers.injection_step &&
      restored.tracers.host_cell_index == reference.tracers.host_cell_index &&
      restored.tracers.mass_fraction_of_host == reference.tracers.mass_fraction_of_host &&
      restored.tracers.last_host_mass_code == reference.tracers.last_host_mass_code &&
      restored.tracers.cumulative_exchanged_mass_code == reference.tracers.cumulative_exchanged_mass_code &&
      restored.species.count_by_species == reference.species.count_by_species &&
      moduleSidecarsEqualForRestart(restored.sidecars, reference.sidecars);
}

}  // namespace

bool restartRuntimeStateExactlyEquivalent(
    const core::SimulationState& restored,
    const core::SimulationState& reference) {
  return restartRuntimeStateExactlyEquivalentImpl(restored, reference);
}

}  // namespace cosmosim::workflows::internal
