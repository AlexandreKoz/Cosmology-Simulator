#include <cassert>
#include <cmath>
#include <string>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"

int main() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.particles.mass_code[0] = 1.0;
  state.particle_sidecar.species_tag[0] =
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar);
  state.resizeCells(0);

  state.star_particles.resize(1);
  state.star_particles.particle_index[0] = 0;
  state.star_particles.formation_scale_factor[0] = 0.2;
  state.star_particles.birth_mass_code[0] = 1.0;
  state.star_particles.metallicity_mass_fraction[0] = 0.01;

  cosmosim::physics::StellarEvolutionConfig config;
  const auto table = cosmosim::physics::StellarEvolutionTable::loadFromTextFile(
      std::string(COSMOSIM_SOURCE_DIR) +
      "/resources/stellar_evolution/test_synthetic_v2.txt");
  cosmosim::physics::StellarEvolutionBookkeeper bookkeeper(config, table);

  std::vector<std::uint32_t> active{0};
  const auto evaluated = bookkeeper.evaluateElapsedYears(state, active, 1.0e8);
  assert(evaluated.counters.evolved_stars == 1);
  assert(evaluated.counters.returned_mass_code > 0.0);
  assert(evaluated.counters.returned_metals_code > 0.0);
  assert(state.particles.mass_code[0] == 1.0);

  bookkeeper.commitBudgets(state, evaluated);
  const double mass_old = evaluated.budgets[0].mass_old_code;
  const double mass_new = evaluated.budgets[0].mass_new_code;
  const double returned = evaluated.budgets[0].interval.returned_mass_code;
  const double remnant_change = evaluated.budgets[0].interval.remnant_change_code;
  assert(std::abs(mass_old - (mass_new + returned + remnant_change)) < 1.0e-12);
  assert(state.star_particles.stellar_returned_mass_cumulative_code[0] == returned);
  assert(state.star_particles.stellar_newly_synthesized_metals_cumulative_code[0] > 0.0);
  assert(state.star_particles.stellar_age_years_last[0] == 1.0e8);
  return 0;
}
