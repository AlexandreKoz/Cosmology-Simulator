#include <array>
#include <cassert>
#include <cmath>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/physics/tracer_support.hpp"

namespace {

void testTracerMassTracksHostAcrossOrchestratorStep() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(2);
  state.particle_sidecar.species_tag[0] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kTracer);
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kTracer)] = 1;

  state.cells.mass_code[0] = 6.0;
  state.cells.mass_code[1] = 5.0;

  cosmosim::physics::TracerConfig config;
  config.enabled = true;
  config.track_mass = true;
  cosmosim::physics::TracerModel model(config);
  model.inject(state, cosmosim::physics::TracerInjectionRequest{
                          .tracer_particle_index = 0,
                          .host_cell_index = 0,
                          .parent_particle_id = 5001,
                          .injection_step = 0,
                          .injected_mass_code = 3.0,
                      });

  cosmosim::physics::TracerCallback callback(std::move(model));
  cosmosim::core::StepOrchestrator orchestrator;
  orchestrator.registerCallback(callback);

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 0.1;

  // Mimic conservative host-cell exchange from hydro/source terms before tracer update.
  state.cells.mass_code[0] = 4.0;
  state.cells.mass_code[1] = 7.0;

  const std::array<std::uint32_t, 2> active_cells{0, 1};
  cosmosim::core::ActiveSetDescriptor active_set{
      .cell_indices = active_cells,
      .cells_are_subset = true,
  };

  orchestrator.executeSingleStep(state, integrator_state, active_set, nullptr, nullptr, nullptr);

#if COSMOSIM_ENABLE_TRACERS
  assert(std::abs(state.particles.mass_code[0] - 2.0) < 1.0e-12);
  assert(callback.lastUpdateCounters().updated_tracers == 1);
#else
  assert(callback.lastUpdateCounters().updated_tracers == 0);
#endif
}

}  // namespace

int main() {
  testTracerMassTracksHostAcrossOrchestratorStep();
  return 0;
}
