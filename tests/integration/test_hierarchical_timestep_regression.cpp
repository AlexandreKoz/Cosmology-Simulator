#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace {

bool throwsWithContext(const std::function<void()>& action, const std::string& required_context) {
  try {
    action();
  } catch (const std::exception& ex) {
    return std::string(ex.what()).find(required_context) != std::string::npos;
  }
  return false;
}

void testMappingRegressionReference() {
  const cosmosim::core::TimeStepLimits limits{
      .min_dt_time_code = 0.0625,
      .max_dt_time_code = 1.0,
      .max_bin = 4,
  };

  const auto coarse = cosmosim::core::mapDtToTimeBin(0.999, limits);
  const auto medium = cosmosim::core::mapDtToTimeBin(0.2, limits);
  const auto fine = cosmosim::core::mapDtToTimeBin(0.0625, limits);

  assert(coarse.bin_index == 3);
  assert(medium.bin_index == 1);
  assert(fine.bin_index == 0);

  assert(cosmosim::core::binIndexToDt(coarse.bin_index, limits) == 0.5);
  assert(cosmosim::core::binIndexToDt(medium.bin_index, limits) == 0.125);
}

void testHydroCflRejectsTooLargeAcceptedDtBeforeMutation() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(1);
  state.particle_sidecar.particle_id[0] = 1001;
  state.particle_sidecar.species_tag[0] =
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
  state.species.count_by_species = {};
  state.species.count_by_species[
      static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)] = 1;
  state.rebuildSpeciesIndex();
  state.refreshGasCellIdentityFromParticleOrder();
  state.gas_cells.density_code[0] = 1.0;
  state.gas_cells.pressure_code[0] = 1.0;
  state.gas_cells.internal_energy_code[0] = 1.5;
  state.gas_cells.sound_speed_code[0] = 1.0;
  state.cells.mass_code[0] = 1.0;
  state.particles.mass_code[0] = 1.0;
  state.particles.velocity_x_peculiar[0] = 100.0;
  state.particles.velocity_y_peculiar[0] = 0.0;
  state.particles.velocity_z_peculiar[0] = 0.0;

  const double density_before = state.gas_cells.density_code[0];
  const double pressure_before = state.gas_cells.pressure_code[0];
  const double mass_before = state.cells.mass_code[0];
  const cosmosim::core::DirectionalCflTimeStepInput cfl_input{
      .cell_width_axis_code = {0.25, 0.25, 0.25},
      .velocity_axis_code = {
          state.particles.velocity_x_peculiar[0],
          state.particles.velocity_y_peculiar[0],
          state.particles.velocity_z_peculiar[0]},
      .sound_speed_code = state.gas_cells.sound_speed_code[0],
  };
  const double accepted_dt_time_code = 0.1;
  const auto diagnostics = cosmosim::core::makeHydroCflDiagnostics(
      0,
      cfl_input,
      0.4,
      accepted_dt_time_code,
      state.gas_cells.gas_cell_id[0],
      std::nullopt,
      std::nullopt);
  assert(diagnostics.proposed_dt_time_code < accepted_dt_time_code);
  assert(throwsWithContext(
      [&]() { cosmosim::core::assertHydroCflStable(diagnostics); },
      "hydro CFL violation"));

  assert(state.gas_cells.density_code[0] == density_before);
  assert(state.gas_cells.pressure_code[0] == pressure_before);
  assert(state.cells.mass_code[0] == mass_before);
}

void testCoarseDisplacementBound() {
  using namespace cosmosim::core;
  ComovingDisplacementTimeStepInput in{.mesh_or_split_length_comoving_code = 0.125,
      .velocity_magnitude_peculiar_code = 2.0, .scale_free_acceleration_magnitude_code = 0.0,
      .scale_factor = 0.25};
  assert(computeComovingDisplacementTimeStep(in) == 0.015625);
  in.velocity_magnitude_peculiar_code = 0.0;
  in.scale_free_acceleration_magnitude_code = 4.0;
  const double expected = std::sqrt(2.0*0.125*std::pow(0.25,3)/4.0);
  assert(std::abs(computeComovingDisplacementTimeStep(in)-expected) < 1e-15);
  for (const double a : {0.01,0.25,1.0}) {
    in.scale_factor = a;
    in.velocity_magnitude_peculiar_code = 3.0;
    const double dt = computeComovingDisplacementTimeStep(in);
    const double displacement = 3.0/a*dt + 0.5*4.0/std::pow(a,3)*dt*dt;
    assert(std::abs(displacement-in.mesh_or_split_length_comoving_code) < 1e-14);
  }
  in.velocity_magnitude_peculiar_code = 0.0;
  in.scale_free_acceleration_magnitude_code = 0.0;
  assert(std::isinf(computeComovingDisplacementTimeStep(in)));
  in.scale_factor = 0.0;
  assert(throwsWithContext([&] { (void)computeComovingDisplacementTimeStep(in); }, "displacement"));
}

void testProductionHierarchicalDispatcherFreeDrift() {
  using namespace cosmosim::core;
  SimulationState state;
  state.resizeParticles(3);
  for (std::size_t row=0; row<3; ++row) {
    state.particle_sidecar.particle_id[row] = row+1U;
    state.particle_sidecar.species_tag[row] = static_cast<std::uint32_t>(ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[row] = 0U;
    state.particles.mass_code[row] = 1.0;
    state.particles.velocity_x_peculiar[row] = static_cast<double>(row+1U);
  }
  state.rebuildSpeciesIndex();
  state.updateAllParticleDriftEpoch(0.0,0.25);
  assert(state.compactHomogeneousDmoMetadata(0U));
  CosmologyBackgroundConfig cfg;
  cfg.omega_matter=1.0;cfg.omega_lambda=0.0;cfg.omega_radiation=0.0;cfg.omega_curvature=0.0;
  LambdaCdmBackground background(cfg);
  IntegratorState integrator;
  integrator.current_time_code=0.0;integrator.current_scale_factor=0.25;
  integrator.time_si_per_code=1.0/background.hubble0Si();integrator.dt_time_code=0.01;
  HierarchicalTimeBinScheduler particles(2),cells(2);
  particles.reset(3,0,0);cells.reset(0,0,0);
  particles.setElementBin(1,1,0);particles.setElementBin(2,2,0);
  TransientStepWorkspace workspace;
  workspace.prepareGravityParticleIndexScratch(3);
  workspace.gravity_particle_index_scratch={0,1,2};
  StepOrchestrator orchestrator;
  const auto generation=state.gravitySourceGeneration();
  unsigned drifts=0,refreshes=0,endpoint_outputs=0;
  orchestrator.executeHierarchicalBlockWithDispatcher(state,integrator,particles,cells,
      [&](StepContext& context,bool safe_output) {
        assert(context.particle_scheduler==&particles);
        assert(context.gas_cell_scheduler==&cells);
        if (context.stage==IntegrationStage::kDrift) {
          assert(context.active_set.particle_indices.size()==3U);
          for (const auto row:context.active_set.particle_indices) {
            state.particles.position_x_comoving[row] +=
                context.timeline_step.drift_factor_code*state.particles.velocity_x_peculiar[row];
            state.particles.velocity_x_peculiar[row] *= context.timeline_step.hubble_drag_factor;
          }
          ++drifts;
        } else if (context.stage==IntegrationStage::kForceRefresh) {
          ++refreshes;
          const auto tick=particles.currentTick();
          for (const auto row:context.active_set.particle_indices) {
            assert(tick%particles.binPeriodTicks(particles.binIndex(row))==0U);
            assert(state.particleLastDriftTimeCode(row)==context.timeline_step.time_end_code);
            assert(state.particleLastDriftScaleFactor(row)==context.timeline_step.scale_factor_end);
          }
          assert(context.active_set.particle_indices.size()==(tick==4U?3U:tick==2U?2U:1U));
          assert(context.hierarchical_kdk.include_long_range_force==(tick==4U));
          assert(state.gravitySourceGeneration()==generation+tick);
        } else if (context.stage==IntegrationStage::kGravityKickPre ||
                   context.stage==IntegrationStage::kGravityKickPost) {
          for (const auto row:context.active_set.particle_indices) {
            assert(context.hierarchical_kdk.tree_kick_factor_code[particles.binIndex(row)]>0.0);
          }
          // Zero force: kick does not change the peculiar momentum.
        } else if (context.stage==IntegrationStage::kOutputCheck) {
          if (safe_output) {++endpoint_outputs;assert(context.boundary.restart_safe);}
          else assert(!context.boundary.output_safe);
        }
      },&background,workspace,nullptr);
  assert(drifts==4U && refreshes==4U && endpoint_outputs==1U);
  const double a_end=std::pow(std::pow(0.25,1.5)+1.5*0.04,2.0/3.0);
  assert(std::abs(integrator.current_scale_factor-a_end)<1e-7);
  for (std::size_t row=0;row<3;++row) {
    const double u0=static_cast<double>(row+1U);
    assert(std::abs(integrator.current_scale_factor*state.particles.velocity_x_peculiar[row]-0.25*u0)<1e-14);
    const double x=2.0*0.25*u0*(1.0/std::sqrt(0.25)-1.0/std::sqrt(a_end));
    assert(std::abs(state.particles.position_x_comoving[row]-x)<2e-6);
  }
  assert(particles.currentTick()==4U && cells.currentTick()==4U);
  assert(integrator.step_index==1U && !integrator.inside_kdk_step && integrator.last_completed_restart_safe);
  const auto saved=particles.exportPersistentState();
  for (const auto tick:saved.next_activation_tick) assert(tick==4U);
  for (const auto flag:saved.active_flag) assert(flag==0U);
  // A duplicated/missing source row must fail before any next block mutates.
  workspace.gravity_particle_index_scratch={0,0,2};
  assert(throwsWithContext([&] {
    orchestrator.executeHierarchicalBlockWithDispatcher(state,integrator,particles,cells,
        [](StepContext&,bool){},&background,workspace,nullptr);
  },"exact identity row set"));
  assert(integrator.step_index==1U && !integrator.inside_kdk_step);
}

}  // namespace

int main() {
  testCoarseDisplacementBound();
  testProductionHierarchicalDispatcherFreeDrift();
  testMappingRegressionReference();
  testHydroCflRejectsTooLargeAcceptedDtBeforeMutation();
  return 0;
}
