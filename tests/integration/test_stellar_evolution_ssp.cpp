#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"

int main() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.particles.mass_code[0] = 1.0;
  state.particle_sidecar.species_tag[0] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar);
  state.resizeCells(0);

  state.star_particles.resize(1);
  state.star_particles.particle_index[0] = 0;
  state.star_particles.formation_scale_factor[0] = 0.2;
  state.star_particles.birth_mass_code[0] = 1.0;
  state.star_particles.metallicity_mass_fraction[0] = 0.01;

  cosmosim::physics::StellarEvolutionConfig config;
  cosmosim::physics::StellarEvolutionBookkeeper bookkeeper(
      config,
      cosmosim::physics::StellarEvolutionTable::makeBuiltinReference());

  std::vector<std::uint32_t> active{0};
  const auto report = bookkeeper.apply(state, active, 0.25, 0.01);

  assert(report.counters.evolved_stars == 1);
  assert(report.counters.returned_mass_code > 0.0);
  assert(report.counters.returned_metals_code > 0.0);

  const double mass_old = report.budgets[0].mass_old_code;
  const double mass_new = report.budgets[0].mass_new_code;
  const double returned = report.budgets[0].interval.returned_mass_code;
  const double remnant_change = report.budgets[0].interval.remnant_change_code;
  assert(std::abs(mass_old - (mass_new + returned + remnant_change)) < 1.0e-12);

  assert(state.star_particles.stellar_returned_mass_cumulative_code[0] == returned);
  assert(state.star_particles.stellar_age_years_last[0] > 0.0);

  return 0;
}
