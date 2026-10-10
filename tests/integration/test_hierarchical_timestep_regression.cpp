#include "../support/test_require.hpp"
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

  CHUI_REQUIRE(coarse.bin_index == 3);
  CHUI_REQUIRE(medium.bin_index == 1);
  CHUI_REQUIRE(fine.bin_index == 0);

  CHUI_REQUIRE(cosmosim::core::binIndexToDt(coarse.bin_index, limits) == 0.5);
  CHUI_REQUIRE(cosmosim::core::binIndexToDt(medium.bin_index, limits) == 0.125);
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
  CHUI_REQUIRE(diagnostics.proposed_dt_time_code < accepted_dt_time_code);
  CHUI_REQUIRE(throwsWithContext(
      [&]() { cosmosim::core::assertHydroCflStable(diagnostics); },
      "hydro CFL violation"));

  CHUI_REQUIRE(state.gas_cells.density_code[0] == density_before);
  CHUI_REQUIRE(state.gas_cells.pressure_code[0] == pressure_before);
  CHUI_REQUIRE(state.cells.mass_code[0] == mass_before);
}

void testCoarseDisplacementBound() {
  using namespace cosmosim::core;
  ComovingDisplacementTimeStepInput in{.mesh_or_split_length_comoving_code = 0.125,
      .velocity_magnitude_peculiar_code = 2.0, .scale_free_acceleration_magnitude_code = 0.0,
      .scale_factor = 0.25};
  CHUI_REQUIRE(computeComovingDisplacementTimeStep(in) == 0.015625);
  in.velocity_magnitude_peculiar_code = 0.0;
  in.scale_free_acceleration_magnitude_code = 4.0;
  const double expected = std::sqrt(2.0*0.125*std::pow(0.25,3)/4.0);
  CHUI_REQUIRE(std::abs(computeComovingDisplacementTimeStep(in)-expected) < 1e-15);
  for (const double a : {0.01,0.25,1.0}) {
    in.scale_factor = a;
    in.velocity_magnitude_peculiar_code = 3.0;
    const double dt = computeComovingDisplacementTimeStep(in);
    const double displacement = 3.0/a*dt + 0.5*4.0/std::pow(a,3)*dt*dt;
    CHUI_REQUIRE(std::abs(displacement-in.mesh_or_split_length_comoving_code) < 1e-14);
  }
  in.velocity_magnitude_peculiar_code = 0.0;
  in.scale_free_acceleration_magnitude_code = 0.0;
  CHUI_REQUIRE(std::isinf(computeComovingDisplacementTimeStep(in)));
  in.scale_factor = 0.0;
  CHUI_REQUIRE(throwsWithContext([&] { (void)computeComovingDisplacementTimeStep(in); }, "displacement"));
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
  CHUI_REQUIRE(state.compactHomogeneousDmoMetadata(0U));
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
        CHUI_REQUIRE(context.particle_scheduler==&particles);
        CHUI_REQUIRE(context.gas_cell_scheduler==&cells);
        if (context.stage==IntegrationStage::kDrift) {
          CHUI_REQUIRE(context.active_set.particle_indices.size()==3U);
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
            CHUI_REQUIRE(tick%particles.binPeriodTicks(particles.binIndex(row))==0U);
            CHUI_REQUIRE(state.particleLastDriftTimeCode(row)==context.timeline_step.time_end_code);
            CHUI_REQUIRE(state.particleLastDriftScaleFactor(row)==context.timeline_step.scale_factor_end);
          }
          CHUI_REQUIRE(context.active_set.particle_indices.size()==(tick==4U?3U:tick==2U?2U:1U));
          CHUI_REQUIRE(context.hierarchical_kdk.include_long_range_force==(tick==4U));
          CHUI_REQUIRE(state.gravitySourceGeneration()==generation+tick);
        } else if (context.stage==IntegrationStage::kGravityKickPre ||
                   context.stage==IntegrationStage::kGravityKickPost) {
          for (const auto row:context.active_set.particle_indices) {
            CHUI_REQUIRE(context.hierarchical_kdk.tree_kick_factor_code[particles.binIndex(row)]>0.0);
          }
          // Zero force: kick does not change the peculiar momentum.
        } else if (context.stage==IntegrationStage::kOutputCheck) {
          if (safe_output) {++endpoint_outputs;CHUI_REQUIRE(context.boundary.restart_safe);}
          else CHUI_REQUIRE(!context.boundary.output_safe);
        }
      },&background,workspace,nullptr);
  CHUI_REQUIRE(drifts==4U && refreshes==4U && endpoint_outputs==1U);
  const double a_end=std::pow(std::pow(0.25,1.5)+1.5*0.04,2.0/3.0);
  CHUI_REQUIRE(std::abs(integrator.current_scale_factor-a_end)<1e-7);
  for (std::size_t row=0;row<3;++row) {
    const double u0=static_cast<double>(row+1U);
    CHUI_REQUIRE(std::abs(integrator.current_scale_factor*state.particles.velocity_x_peculiar[row]-0.25*u0)<1e-14);
    const double x=2.0*0.25*u0*(1.0/std::sqrt(0.25)-1.0/std::sqrt(a_end));
    CHUI_REQUIRE(std::abs(state.particles.position_x_comoving[row]-x)<2e-6);
  }
  CHUI_REQUIRE(particles.currentTick()==4U && cells.currentTick()==4U);
  CHUI_REQUIRE(integrator.step_index==1U && !integrator.inside_kdk_step && integrator.last_completed_restart_safe);
  const auto saved=particles.exportPersistentState();
  for (const auto tick:saved.next_activation_tick) CHUI_REQUIRE(tick==4U);
  for (const auto flag:saved.active_flag) CHUI_REQUIRE(flag==0U);
  // A duplicated/missing source row must fail before any next block mutates.
  workspace.gravity_particle_index_scratch={0,0,2};
  CHUI_REQUIRE(throwsWithContext([&] {
    orchestrator.executeHierarchicalBlockWithDispatcher(state,integrator,particles,cells,
        [](StepContext&,bool){},&background,workspace,nullptr);
  },"exact identity row set"));
  CHUI_REQUIRE(integrator.step_index==1U && !integrator.inside_kdk_step);
}

void testProductionHierarchicalDispatcherRejectsMissingAndStalePmEpoch() {
  using namespace cosmosim::core;
  for (const bool missing_event : {false, true}) {
    SimulationState state;
    state.resizeParticles(1);
    state.particle_sidecar.particle_id[0] = 77;
    state.particle_sidecar.species_tag[0] = static_cast<std::uint32_t>(ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[0] = 0;
    state.particles.mass_code[0] = 1.0;
    state.rebuildSpeciesIndex();
    state.updateAllParticleDriftEpoch(0.0, 0.5);
    CHUI_REQUIRE(state.compactHomogeneousDmoMetadata(0));
    IntegratorState time;
    time.current_scale_factor = 0.5;
    time.dt_time_code = 0.01;
    time.pm_refresh_enabled = true;
    time.pm_sync_state.reset(1);
    HierarchicalTimeBinScheduler particles(1), cells(1);
    particles.reset(1, 0, 0); cells.reset(0, 0, 0);
    TransientStepWorkspace workspace;
    workspace.prepareGravityParticleIndexScratch(1);
    workspace.gravity_particle_index_scratch = {0};
    StepOrchestrator orchestrator;
    bool injected = false;
    const auto dispatcher = [&](StepContext& context, bool) {
      if (context.stage == IntegrationStage::kForceRefresh && context.hierarchical_kdk.synchronization_end) {
        injected = true;
        auto& d = context.pm_refresh_directive;
        d.has_sync_event = !missing_event;
        d.solver_executed = true;
        d.refresh_long_range_field = true;
        d.gravity_kick_opportunity = 1;
        d.field_version = 1;
        d.last_refresh_opportunity = 1;
        d.field_built_step_index = time.step_index;
        d.field_built_scale_factor = 0.25; // Old field, evaluated at a=0.5.
      }
    };
    CHUI_REQUIRE(throwsWithContext([&] {
      orchestrator.executeHierarchicalBlockWithDispatcher(state, time, particles, cells,
          dispatcher, nullptr, workspace, nullptr);
    }, missing_event ? "omitted synchronized refresh event" : "stale field epoch"));
    CHUI_REQUIRE(injected);
    CHUI_REQUIRE(time.step_index == 0U && !time.last_completed_restart_safe);
    CHUI_REQUIRE(time.pm_sync_state.fieldVersion() == 0U);
  }
}

}  // namespace

int main() {
  testProductionHierarchicalDispatcherRejectsMissingAndStalePmEpoch();
  testCoarseDisplacementBound();
  testProductionHierarchicalDispatcherFreeDrift();
  testMappingRegressionReference();
  testHydroCflRejectsTooLargeAcceptedDtBeforeMutation();
  return 0;
}
