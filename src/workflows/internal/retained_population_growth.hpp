#pragma once
#include "cosmosim/core/retained_capacity_transaction.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::workflows::internal {
inline void planParticlePopulationGrowth(core::RetainedCapacityTransaction& plan,
                                          core::SimulationState& state,
                                          std::size_t particle_count,
                                          std::size_t species_count,
                                          core::ParticleSpecies species) {
  (void)core::checkedLocalCount(particle_count, core::kMaxLocalParticleCount,
      "particle", "source population growth");
  plan.add(state.particles.position_x_comoving, particle_count);
  plan.add(state.particles.position_y_comoving, particle_count);
  plan.add(state.particles.position_z_comoving, particle_count);
  plan.add(state.particles.velocity_x_peculiar, particle_count);
  plan.add(state.particles.velocity_y_peculiar, particle_count);
  plan.add(state.particles.velocity_z_peculiar, particle_count);
  plan.add(state.particles.mass_code, particle_count);
  plan.add(state.particles.time_bin, particle_count);
  plan.add(state.particle_sidecar.particle_id, particle_count);
  plan.add(state.particle_sidecar.sfc_key, particle_count);
  plan.add(state.particle_sidecar.species_tag, particle_count);
  plan.add(state.particle_sidecar.particle_flags, particle_count);
  plan.add(state.particle_sidecar.owning_rank, particle_count);
  plan.add(state.particle_sidecar.last_drift_time_code, particle_count);
  plan.add(state.particle_sidecar.last_drift_scale_factor, particle_count);
  if (!state.particle_sidecar.gravity_softening_comoving.empty() ||
      !state.particle_sidecar.has_gravity_softening_override.empty()) {
    plan.add(state.particle_sidecar.gravity_softening_comoving, particle_count);
  }
  if (!state.particle_sidecar.has_gravity_softening_override.empty()) {
    plan.add(state.particle_sidecar.has_gravity_softening_override, particle_count);
  }
  plan.add(state.particle_species_index.local_index_by_global, particle_count);
  plan.add(state.particle_species_index.global_index_by_species[
      core::particleSpeciesIndex(species)], species_count);
  if (species == core::ParticleSpecies::kStar) {
    plan.add(state.star_particles.particle_index, species_count);
    plan.add(state.star_particles.formation_scale_factor, species_count);
    plan.add(state.star_particles.birth_mass_code, species_count);
    plan.add(state.star_particles.metallicity_mass_fraction, species_count);
    plan.add(state.star_particles.birth_key, species_count);
    plan.add(state.star_particles.parent_gas_cell_id, species_count);
    plan.add(state.star_particles.birth_tick, species_count);
    plan.add(state.star_particles.birth_ordinal, species_count);
    plan.add(state.star_particles.stellar_age_years_last, species_count);
    plan.add(state.star_particles.stellar_returned_mass_cumulative_code, species_count);
    plan.add(state.star_particles.stellar_returned_metals_cumulative_code, species_count);
    plan.add(state.star_particles.stellar_newly_synthesized_metals_cumulative_code, species_count);
    plan.add(state.star_particles.stellar_feedback_energy_cumulative_erg, species_count);
    plan.add(state.star_particles.enrichment_carry_mass_code, species_count);
    plan.add(state.star_particles.enrichment_carry_metals_code, species_count);
    plan.add(state.star_particles.enrichment_carry_feedback_energy_erg, species_count);
    plan.add(state.star_particles.enrichment_carry_momentum_code, species_count);
    plan.add(state.star_particles.stellar_deposited_mass_cumulative_code, species_count);
    plan.add(state.star_particles.stellar_deposited_metals_cumulative_code, species_count);
    plan.add(state.star_particles.stellar_deposited_feedback_energy_cumulative_erg, species_count);
    for (auto& lane : state.star_particles.stellar_returned_mass_channel_cumulative_code) plan.add(lane, species_count);
    for (auto& lane : state.star_particles.stellar_returned_metals_channel_cumulative_code) plan.add(lane, species_count);
    for (auto& lane : state.star_particles.stellar_feedback_energy_channel_cumulative_erg) plan.add(lane, species_count);
  } else if (species == core::ParticleSpecies::kBlackHole) {
    plan.add(state.black_holes.particle_index, species_count);
    plan.add(state.black_holes.host_cell_index, species_count);
    plan.add(state.black_holes.subgrid_mass_code, species_count);
    plan.add(state.black_holes.accretion_rate_code, species_count);
    plan.add(state.black_holes.feedback_energy_code, species_count);
    plan.add(state.black_holes.eddington_ratio, species_count);
    plan.add(state.black_holes.cumulative_accreted_mass_code, species_count);
    plan.add(state.black_holes.cumulative_feedback_energy_code, species_count);
    plan.add(state.black_holes.duty_cycle_active_time_code, species_count);
    plan.add(state.black_holes.duty_cycle_total_time_code, species_count);
  } else {
    throw std::invalid_argument("source population growth requires star or black-hole species");
  }
}
}  // namespace cosmosim::workflows::internal
