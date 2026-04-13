#include <cassert>

#include "cosmosim/cosmosim.hpp"

int main() {
  cosmosim::core::SimulationConfig config;
  config.output.run_name = "integration_ic_generated";
  config.units.length_unit = "kpc";

  const cosmosim::io::IcReadResult result = cosmosim::io::convertGeneratedIsolatedIcToState(config, 6);
  const auto& state = result.state;

  assert(state.metadata.run_name == "integration_ic_generated");
  assert(state.particles.size() == 42);
  assert(state.validateOwnershipInvariants());
  assert(state.validateUniqueParticleIds());

  const auto dm_indices =
      state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kDarkMatter);
  const auto gas_indices = state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kGas);
  assert(dm_indices.size() == 36);
  assert(gas_indices.size() == 6);

  double checksum = 0.0;
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    checksum += state.particles.position_x_comoving[i];
    checksum += state.particles.mass_code[i] * 1.0e-20;
  }
  assert(checksum > 0.0);

  return 0;
}
