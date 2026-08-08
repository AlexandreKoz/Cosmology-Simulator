#include <cassert>

#include "cosmosim/analysis/halo_workflow.hpp"

namespace {

cosmosim::core::SimulationState makeTwoGroupState() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(6);
  state.species.count_by_species = {6, 0, 0, 0, 0};

  for (std::size_t i = 0; i < 6; ++i) {
    state.particle_sidecar.particle_id[i] = 100 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particles.mass_code[i] = 1.0;
    state.particles.velocity_x_peculiar[i] = 0.0;
    state.particles.velocity_y_peculiar[i] = 0.0;
    state.particles.velocity_z_peculiar[i] = 0.0;
  }

  // Group A
  state.particles.position_x_comoving[0] = 0.10;
  state.particles.position_y_comoving[0] = 0.10;
  state.particles.position_z_comoving[0] = 0.10;
  state.particles.position_x_comoving[1] = 0.11;
  state.particles.position_y_comoving[1] = 0.10;
  state.particles.position_z_comoving[1] = 0.10;
  state.particles.position_x_comoving[2] = 0.12;
  state.particles.position_y_comoving[2] = 0.10;
  state.particles.position_z_comoving[2] = 0.10;

  // Group B
  state.particles.position_x_comoving[3] = 0.80;
  state.particles.position_y_comoving[3] = 0.80;
  state.particles.position_z_comoving[3] = 0.80;
  state.particles.position_x_comoving[4] = 0.81;
  state.particles.position_y_comoving[4] = 0.80;
  state.particles.position_z_comoving[4] = 0.80;
  state.particles.position_x_comoving[5] = 0.82;
  state.particles.position_y_comoving[5] = 0.80;
  state.particles.position_z_comoving[5] = 0.80;

  return state;
}

void testFofFindsTwoGroups() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.output.run_name = "unit_halo";
  config.cosmology.box_size_mpc_comoving = 1.0;

  cosmosim::analysis::FofConfig fof;
  fof.linking_length_factor_mean_interparticle = 0.30;
  fof.min_group_size = 2;

  cosmosim::analysis::FofHaloFinder finder(fof);
  const auto state = makeTwoGroupState();

  cosmosim::analysis::FofProfilingCounters profiling;
  const auto catalog = finder.buildCatalog(state, config, 7, 0.5, &profiling);
  const auto reference = finder.buildCatalogReferenceAllPairsFromView(
      cosmosim::analysis::buildHaloParticleView(state), config, 7, 0.5);

  assert(catalog.halos.size() == 2);
  assert(catalog.halo_id_by_particle == reference.halo_id_by_particle);
  assert(catalog.halos.size() == reference.halos.size());
  assert(profiling.candidate_particle_count == 6);
  assert(profiling.pair_checks <= 15);
  assert(profiling.neighbor_search == cosmosim::analysis::FofNeighborSearch::kSpatialHash);

  std::uint64_t assigned_count = 0;
  for (const auto halo_id : catalog.halo_id_by_particle) {
    if (halo_id != cosmosim::analysis::k_unbound_halo_id) {
      ++assigned_count;
    }
  }
  assert(assigned_count == 6);
  assert(catalog.halos[0].particle_count == 3);
  assert(catalog.halos[1].particle_count == 3);
}

void testPeriodicBoundaryGroupMatchesReference() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.output.run_name = "unit_halo_periodic";
  config.cosmology.box_size_mpc_comoving = 1.0;
  cosmosim::core::SimulationState state;
  state.resizeParticles(3);
  for (std::size_t i=0;i<3;++i) {
    state.particle_sidecar.particle_id[i]=200+i;
    state.particle_sidecar.species_tag[i]=static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particles.mass_code[i]=1.0;
  }
  state.particles.position_x_comoving[0]=0.99;
  state.particles.position_x_comoving[1]=0.01;
  state.particles.position_x_comoving[2]=0.03;
  for (std::size_t i=0;i<3;++i) {
    state.particles.position_y_comoving[i]=0.5;
    state.particles.position_z_comoving[i]=0.5;
  }
  cosmosim::analysis::FofConfig fof;
  fof.linking_length_factor_mean_interparticle=0.35;
  fof.min_group_size=2;
  cosmosim::analysis::FofHaloFinder finder(fof);
  const auto view=cosmosim::analysis::buildHaloParticleView(state);
  const auto spatial=finder.buildCatalogFromView(view,config,3,1.0);
  const auto reference=finder.buildCatalogReferenceAllPairsFromView(view,config,3,1.0);
  assert(spatial.halo_id_by_particle==reference.halo_id_by_particle);
  assert(spatial.halos.size()==1);
  const double x=spatial.halos.front().center_of_mass_comov[0];
  assert(x < 0.1 || x > 0.9);
  assert(spatial.subhalo_candidates.empty());
}

}  // namespace

int main() {
  testFofFindsTwoGroups();
  testPeriodicBoundaryGroupMatchesReference();
  return 0;
}
