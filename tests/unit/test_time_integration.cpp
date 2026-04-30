#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/simulation_state.hpp"

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

void syncStateTimeBinsFromScheduler(
    const cosmosim::core::HierarchicalTimeBinScheduler& scheduler,
    cosmosim::core::SimulationState& state) {
  cosmosim::core::syncTimeBinMirrorsFromScheduler(
      scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticlesAndCells);
}

void testTimestepBinAuthorityInvariant() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(4);
  state.resizeCells(2);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particles.time_bin[i] = 0;
  }
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.time_bin[i] = 0;
  }

  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(4, 0, 0);
  scheduler.setElementBin(1, 1, scheduler.currentTick());
  scheduler.setElementBin(2, 2, scheduler.currentTick());
  syncStateTimeBinsFromScheduler(scheduler, state);

  assert(cosmosim::core::timeBinMirrorsMatchScheduler(scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticlesAndCells));
  cosmosim::core::debugAssertTimeBinMirrorAuthorityInvariant(scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticlesAndCells);

  // Unauthorized non-owner mutation of mirror lanes is detectable.
  state.particles.time_bin[2] = 0;
  assert(!cosmosim::core::timeBinMirrorsMatchScheduler(scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticlesAndCells));
  bool threw = false;
  try {
    cosmosim::core::debugAssertTimeBinMirrorAuthorityInvariant(scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticlesAndCells);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);

  // Active-set authority remains scheduler-owned despite mirror tampering.
  const auto active = scheduler.beginSubstep();
  assert(active.size() == 4);
  assert(active[0] == 0 && active[1] == 1 && active[2] == 2 && active[3] == 3);
  scheduler.endSubstep();

  // Authorized path re-synchronizes mirrors.
  syncStateTimeBinsFromScheduler(scheduler, state);
  assert(cosmosim::core::timeBinMirrorsMatchScheduler(scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticlesAndCells));
}

void testTimestepBinReassignmentAndRestartRoundTrip() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(5, 1, 0);

  scheduler.beginSubstep();
  scheduler.requestBinTransition(0, 0);
  scheduler.requestBinTransition(3, 2);
  scheduler.endSubstep();

  while (scheduler.currentTick() < 4) {
    scheduler.beginSubstep();
    scheduler.endSubstep();
  }

  scheduler.beginSubstep();
  scheduler.requestBinTransition(0, 2);
  scheduler.endSubstep();

  const auto saved = scheduler.exportPersistentState();
  cosmosim::core::HierarchicalTimeBinScheduler restored(saved.max_bin);
  restored.importPersistentState(saved);
  const auto loaded = restored.exportPersistentState();
  assert(loaded.current_tick == saved.current_tick);
  assert(loaded.max_bin == saved.max_bin);
  assert(loaded.bin_index == saved.bin_index);
  assert(loaded.next_activation_tick == saved.next_activation_tick);
  assert(loaded.active_flag == saved.active_flag);
  assert(loaded.pending_bin_index == saved.pending_bin_index);
}

void testTimestepBinReorderIdentitySurvival() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(5);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 100 + i;
    state.particle_sidecar.sfc_key[i] = static_cast<std::uint64_t>(9 - i);
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.particle_flags[i] = 0;
    state.particle_sidecar.owning_rank[i] = 0;
  }
  state.species.count_by_species = {5, 0, 0, 0, 0};

  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(5, 0, 0);
  scheduler.setElementBin(0, 2, scheduler.currentTick());
  scheduler.setElementBin(1, 0, scheduler.currentTick());
  scheduler.setElementBin(2, 1, scheduler.currentTick());
  scheduler.setElementBin(3, 3, scheduler.currentTick());
  scheduler.setElementBin(4, 1, scheduler.currentTick());
  syncStateTimeBinsFromScheduler(scheduler, state);

  std::vector<std::uint8_t> original_bin_by_id(5, 0);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    original_bin_by_id[i] = state.particles.time_bin[i];
  }

  const auto reorder_map = cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kBySfcKey);
  cosmosim::core::reorderParticles(state, reorder_map);

  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const std::uint64_t id = state.particle_sidecar.particle_id[i];
    const std::size_t old_index = static_cast<std::size_t>(id - 100);
    assert(state.particles.time_bin[i] == original_bin_by_id[old_index]);
  }
}

std::vector<std::uint32_t> buildCompetingActiveFromStateMirror(
    const cosmosim::core::SimulationState& state,
    std::uint8_t active_bin) {
  std::vector<std::uint32_t> active;
  for (std::uint32_t i = 0; i < state.particles.size(); ++i) {
    if (state.particles.time_bin[i] == active_bin) {
      active.push_back(i);
    }
  }
  return active;
}

void testActiveSetAuthority() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(6, 2, 0);
  scheduler.setElementBin(0, 0, scheduler.currentTick());
  scheduler.setElementBin(1, 1, scheduler.currentTick());
  scheduler.setElementBin(2, 2, scheduler.currentTick());
  scheduler.setElementBin(3, 1, scheduler.currentTick());
  scheduler.setElementBin(4, 0, scheduler.currentTick());
  scheduler.setElementBin(5, 2, scheduler.currentTick());

  const auto active0 = scheduler.beginSubstep();
  const std::vector<std::uint32_t> expected0{0, 1, 2, 3, 4, 5};
  assert(std::vector<std::uint32_t>(active0.begin(), active0.end()) == expected0);

  cosmosim::core::SimulationState state;
  state.resizeParticles(6);
  const auto descriptor = cosmosim::core::makeSchedulerActiveSetDescriptor(scheduler, state, active0);
  assert(descriptor.particles_from_scheduler);
  assert(descriptor.has_generation_metadata);
  assert(descriptor.source_scheduler_tick == scheduler.currentTick());
  cosmosim::core::debugAssertActiveSetDescriptorFresh(descriptor, state, scheduler);

  bool stale_scheduler_threw = false;
  try {
    cosmosim::core::debugAssertActiveSetDescriptorFresh(descriptor, state, scheduler.currentTick() + 1U);
  } catch (const std::runtime_error&) {
    stale_scheduler_threw = true;
  }
  assert(stale_scheduler_threw);

  state.bumpParticleIndexGeneration();
  bool stale_descriptor_threw = false;
  try {
    cosmosim::core::debugAssertActiveSetDescriptorFresh(descriptor, state);
  } catch (const std::runtime_error&) {
    stale_descriptor_threw = true;
  }
  assert(stale_descriptor_threw);

  scheduler.requestBinTransition(2, 0);

  const auto active1 = scheduler.beginSubstep();
  const std::vector<std::uint32_t> expected1{0, 2, 4};
  assert(std::vector<std::uint32_t>(active1.begin(), active1.end()) == expected1);
  scheduler.endSubstep();
}

void testActiveSetNoCompetingBuilders() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(6);
  state.resizeCells(2);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particles.time_bin[i] = 0;
  }

  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(6, 0, 0);
  scheduler.setElementBin(0, 0, scheduler.currentTick());
  scheduler.setElementBin(1, 1, scheduler.currentTick());
  scheduler.setElementBin(2, 2, scheduler.currentTick());
  scheduler.setElementBin(3, 0, scheduler.currentTick());
  scheduler.setElementBin(4, 1, scheduler.currentTick());
  scheduler.setElementBin(5, 2, scheduler.currentTick());
  syncStateTimeBinsFromScheduler(scheduler, state);

  const auto authoritative = scheduler.beginSubstep();
  assert(!authoritative.empty());

  // Simulate stale mirror mutation; competing local builders should diverge and be detectable.
  state.particles.time_bin[5] = 0;
  const auto competing_after = buildCompetingActiveFromStateMirror(state, 0);
  assert(competing_after != std::vector<std::uint32_t>(authoritative.begin(), authoritative.end()));
  bool threw = false;
  try {
    cosmosim::core::debugAssertTimeBinMirrorAuthorityInvariant(scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticlesAndCells);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
  scheduler.endSubstep();
}

}  // namespace

int main() {
  testKickDriftKickOrdering();
  testCosmologyHelpers();
  testActiveSubsetDetection();
  testTimeBinMappingAndCriteria();
  testHierarchicalSchedulerTransitions();
  testTimestepBinAuthorityInvariant();
  testTimestepBinReassignmentAndRestartRoundTrip();
  testTimestepBinReorderIdentitySurvival();
  testActiveSetAuthority();
  testActiveSetNoCompetingBuilders();
  return 0;
}
