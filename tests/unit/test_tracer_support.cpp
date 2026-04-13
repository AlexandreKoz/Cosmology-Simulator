#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/physics/tracer_support.hpp"

namespace {

void testInjectionAndAssociation() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(1);
  state.particle_sidecar.species_tag[0] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kTracer);
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kTracer)] = 1;
  state.cells.mass_code[0] = 8.0;

  cosmosim::physics::TracerConfig config;
  config.enabled = true;
  config.track_mass = true;
  cosmosim::physics::TracerModel model(config);

  cosmosim::physics::TracerInjectionRequest request;
  request.tracer_particle_index = 0;
  request.host_cell_index = 0;
  request.parent_particle_id = 42;
  request.injection_step = 7;
  request.injected_mass_code = 2.0;

  model.inject(state, request);
  assert(state.tracers.size() == 1);
  assert(state.tracers.particle_index[0] == 0);
  assert(state.tracers.host_cell_index[0] == 0);
  assert(state.tracers.parent_particle_id[0] == 42);
  assert(state.tracers.injection_step[0] == 7);
  assert(std::abs(state.tracers.mass_fraction_of_host[0] - 0.25) < 1.0e-14);
  assert(std::abs(state.particles.mass_code[0] - 2.0) < 1.0e-14);
  assert(state.validateOwnershipInvariants());
}

void testConservativeHostMassScaling() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(1);
  state.particle_sidecar.species_tag[0] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kTracer);
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kTracer)] = 1;
  state.cells.mass_code[0] = 4.0;

  cosmosim::physics::TracerConfig config;
  config.enabled = true;
  config.track_mass = true;
  config.min_host_mass_code = 0.0;
  cosmosim::physics::TracerModel model(config);

  model.inject(state, cosmosim::physics::TracerInjectionRequest{
                          .tracer_particle_index = 0,
                          .host_cell_index = 0,
                          .parent_particle_id = 99,
                          .injection_step = 0,
                          .injected_mass_code = 1.0,
                      });

  state.cells.mass_code[0] = 10.0;
  const std::array<std::uint32_t, 1> active_cells{0};
  const auto counters = model.updateMassFromHostCells(state, active_cells);
  assert(counters.updated_tracers == 1);
  assert(std::abs(state.particles.mass_code[0] - 2.5) < 1.0e-12);
  assert(std::abs(state.tracers.cumulative_exchanged_mass_code[0] - 1.5) < 1.0e-12);
}

void testStageCallbackHook() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(1);
  state.particle_sidecar.species_tag[0] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kTracer);
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kTracer)] = 1;
  state.cells.mass_code[0] = 5.0;

  cosmosim::physics::TracerConfig config;
  config.enabled = true;
  config.track_mass = true;
  cosmosim::physics::TracerModel model(config);
  model.inject(state, cosmosim::physics::TracerInjectionRequest{
                          .tracer_particle_index = 0,
                          .host_cell_index = 0,
                          .parent_particle_id = 101,
                          .injection_step = 0,
                          .injected_mass_code = 1.0,
                      });

  cosmosim::physics::TracerCallback callback(std::move(model));
  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 1.0;
  const std::array<std::uint32_t, 1> active_cells{0};
  cosmosim::core::ActiveSetDescriptor active_set{
      .cell_indices = active_cells,
      .cells_are_subset = true,
  };
  state.cells.mass_code[0] = 15.0;

  cosmosim::core::StepContext context{
      .state = state,
      .integrator_state = integrator_state,
      .active_set = active_set,
      .stage = cosmosim::core::IntegrationStage::kSourceTerms,
  };

  callback.onStage(context);
#if COSMOSIM_ENABLE_TRACERS
  assert(callback.lastUpdateCounters().updated_tracers == 1);
  assert(std::abs(state.particles.mass_code[0] - 3.0) < 1.0e-12);
#else
  assert(callback.lastUpdateCounters().updated_tracers == 0);
#endif
}

}  // namespace

int main() {
  testInjectionAndAssociation();
  testConservativeHostMassScaling();
  testStageCallbackHook();
  return 0;
}
