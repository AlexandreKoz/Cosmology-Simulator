#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "cosmosim/core/time_integration.hpp"

namespace {

constexpr double k_tolerance = 1.0e-12;

class StageRecorder final : public cosmosim::core::IntegrationCallback {
 public:
  std::string_view callbackName() const override { return "stage_recorder"; }

  void onStage(cosmosim::core::StepContext& context) override { observed_stages.push_back(context.stage); }

  std::vector<cosmosim::core::IntegrationStage> observed_stages;
};

void testKickDriftKickOrdering() {
  cosmosim::core::SimulationState state;
  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 1.0;

  StageRecorder recorder;
  cosmosim::core::StepOrchestrator orchestrator;
  orchestrator.registerCallback(recorder);

  const cosmosim::core::ActiveSetDescriptor active_set{};
  orchestrator.executeSingleStep(state, integrator_state, active_set, nullptr, nullptr);

  const auto expected = cosmosim::core::StageScheduler::kickDriftKickOrder();
  assert(cosmosim::core::isCanonicalIntegrationStageOrder(expected));
  assert(recorder.observed_stages.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    assert(recorder.observed_stages[i] == expected[i]);
  }
  assert(cosmosim::core::isCanonicalIntegrationStageOrder(recorder.observed_stages));
}

void testCosmologyHelpers() {
  cosmosim::core::CosmologyBackgroundConfig cfg;
  cfg.hubble_param = 0.7;
  cfg.omega_matter = 1.0;
  cfg.omega_lambda = 0.0;
  cfg.omega_radiation = 0.0;
  cfg.omega_curvature = 0.0;

  cosmosim::core::LambdaCdmBackground background(cfg);

  const double a0 = 0.5;
  const double a1 = 1.0;
  const double drift = cosmosim::core::computeComovingDriftFactor(background, a0, a1, 2000);

  // For Einstein-de Sitter: H(a)=H0*a^{-3/2} => integral da/(a^2 H(a)) = 2*(sqrt(a1)-sqrt(a0))/H0.
  const double h0 = background.hubble0Si();
  const double expected = 2.0 * (std::sqrt(a1) - std::sqrt(a0)) / h0;
  assert(std::abs(drift - expected) / expected < 5.0e-4);

  const double drag = cosmosim::core::computeHubbleDragFactor(a0, a1);
  assert(std::abs(drag - 0.5) < k_tolerance);

  const double delta_a = 1.0e-4;
  const double dt = cosmosim::core::estimateDeltaTimeFromScaleFactorStep(background, a0, delta_a);
  const double a_next = cosmosim::core::advanceScaleFactorEuler(background, a0, dt);
  assert(std::abs(a_next - (a0 + delta_a)) / delta_a < 1.0e-12);
}

void testActiveSubsetDetection() {
  std::vector<std::uint32_t> active_particles = {0, 2, 4};
  std::vector<std::uint32_t> active_cells = {1, 3};

  cosmosim::core::ActiveSetDescriptor active_set{
      .particle_indices = active_particles,
      .cell_indices = active_cells,
      .particles_are_subset = true,
      .cells_are_subset = true,
  };

  assert(active_set.hasParticleSubset(10));
  assert(active_set.hasCellSubset(8));
  assert(!active_set.hasParticleSubset(3));
  assert(!active_set.hasCellSubset(2));
}

void testTimeBinMappingAndCriteria() {
  const cosmosim::core::TimeStepLimits limits{
      .min_dt_time_code = 0.125,
      .max_dt_time_code = 1.0,
      .max_bin = 3,
  };

  const auto mapped = cosmosim::core::mapDtToTimeBin(0.5, limits);
  assert(mapped.bin_index == 2);
  assert(!mapped.clipped_to_min);
  assert(!mapped.clipped_to_max);

  const auto clipped_min = cosmosim::core::mapDtToTimeBin(0.001, limits);
  assert(clipped_min.bin_index == 0);
  assert(clipped_min.clipped_to_min);

  const auto clipped_max = cosmosim::core::mapDtToTimeBin(8.0, limits);
  assert(clipped_max.bin_index == 3);
  assert(clipped_max.clipped_to_max);

  const cosmosim::core::CflTimeStepInput cfl_input{
      .cell_width_code = 0.25,
      .flow_speed_code = 1.0,
      .sound_speed_code = 0.5,
  };
  const double dt_cfl = cosmosim::core::computeCflTimeStep(cfl_input, 0.8);
  assert(std::abs(dt_cfl - (0.8 * 0.25 / 1.5)) < k_tolerance);

  const cosmosim::core::GravityTimeStepInput grav_input{
      .softening_length_code = 0.01,
      .acceleration_magnitude_code = 4.0,
  };
  const double dt_grav = cosmosim::core::computeGravityTimeStep(grav_input, 0.4);
  assert(std::abs(dt_grav - 0.4 * std::sqrt(0.01 / 4.0)) < k_tolerance);

  cosmosim::core::TimeStepCriteriaRegistry registry;
  registry.registerCflHook([](std::uint32_t index) { return index == 0 ? 0.2 : 0.4; });
  registry.registerGravityHook([](std::uint32_t index) { return index == 0 ? 0.3 : 0.1; });
  registry.registerSourceHook([](std::uint32_t) { return std::numeric_limits<double>::infinity(); });
  registry.registerUserClampHook([](std::uint32_t) { return 0.15; });

  const double dt0 = cosmosim::core::combineTimeStepCriteria(0, registry.hooks(), 0.5);
  const double dt1 = cosmosim::core::combineTimeStepCriteria(1, registry.hooks(), 0.5);
  assert(std::abs(dt0 - 0.15) < k_tolerance);
  assert(std::abs(dt1 - 0.1) < k_tolerance);
}

void testHierarchicalSchedulerTransitions() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(4, 2, 0);

  auto active = scheduler.beginSubstep();
  assert(active.size() == 4);

  scheduler.requestBinTransition(0, 0);
  scheduler.requestBinTransition(1, 3);
  scheduler.endSubstep();

  const auto& hot = scheduler.hotMetadata();
  assert(hot.bin_index[0] == 0);
  assert(hot.bin_index[1] == 3);

  active = scheduler.beginSubstep();
  assert(active.size() == 1);
  assert(active[0] == 0);

  scheduler.requestBinTransition(0, 3);
  scheduler.endSubstep();
  const auto illegal_before = scheduler.diagnostics().illegal_transition_attempts;
  assert(illegal_before >= 1);

  while (scheduler.currentTick() < 8) {
    scheduler.beginSubstep();
    scheduler.endSubstep();
  }

  active = scheduler.beginSubstep();
  assert(!active.empty() && active[0] == 0);
  scheduler.requestBinTransition(0, 3);
  scheduler.endSubstep();

  assert(scheduler.hotMetadata().bin_index[0] == 3);
  assert(scheduler.diagnostics().demoted_elements >= 2);
}

}  // namespace

int main() {
  testKickDriftKickOrdering();
  testCosmologyHelpers();
  testActiveSubsetDetection();
  testTimeBinMappingAndCriteria();
  testHierarchicalSchedulerTransitions();
  return 0;
}
