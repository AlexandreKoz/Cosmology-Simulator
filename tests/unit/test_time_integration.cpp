#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace {

constexpr double k_tolerance = 1.0e-12;

bool throwsWithContext(const std::function<void()>& action, const std::string& required_context) {
  try {
    action();
  } catch (const std::exception& ex) {
    return std::string(ex.what()).find(required_context) != std::string::npos;
  }
  return false;
}


void initializeParticleIdsAndGasCells(
    cosmosim::core::SimulationState& state,
    const std::vector<std::uint32_t>& gas_particle_rows) {
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 1000 + i;
    state.particle_sidecar.sfc_key[i] = i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.particle_flags[i] = 0;
    state.particle_sidecar.owning_rank[i] = 0;
  }
  state.species.count_by_species = {};
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] = state.particles.size();
  for (const auto row : gas_particle_rows) {
    state.particle_sidecar.species_tag[row] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
    --state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)];
    ++state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)];
  }
  state.rebuildSpeciesIndex();
  state.refreshGasCellIdentityFromParticleOrder();
}

std::vector<std::uint64_t> particleIdsForRows(
    const cosmosim::core::SimulationState& state,
    std::span<const std::uint32_t> rows) {
  std::vector<std::uint64_t> ids;
  ids.reserve(rows.size());
  for (const auto row : rows) {
    ids.push_back(state.particle_sidecar.particle_id[row]);
  }
  return ids;
}

std::vector<std::uint64_t> activeParticleIdsAtCurrentTick(
    const cosmosim::core::HierarchicalTimeBinScheduler& scheduler,
    const cosmosim::core::SimulationState& state) {
  const auto persistent = scheduler.exportPersistentState();
  cosmosim::core::HierarchicalTimeBinScheduler probe(persistent.max_bin);
  probe.importPersistentState(persistent);
  const auto active = probe.beginSubstep();
  auto ids = particleIdsForRows(state, active);
  std::sort(ids.begin(), ids.end());
  return ids;
}


class SingleStageRecorder final : public cosmosim::core::IntegrationCallback {
 public:
  SingleStageRecorder(std::string_view name, cosmosim::core::IntegrationStage stage, std::vector<std::string>* order)
      : name(name), stage(stage), order(order) {
    contract.stage = stage;
  }

  std::string_view callbackName() const override { return name; }
  std::span<const cosmosim::core::IntegrationStage> integrationStages() const override { return {&stage, 1}; }
  std::span<const cosmosim::core::StageContract> stageContracts() const override { return {&contract, 1}; }

  void onStage(cosmosim::core::StepContext& context) override {
    assert(context.stage == stage);
    ++invocations;
    observed.push_back(context.stage);
    if (order != nullptr) {
      order->push_back(std::string(name));
    }
  }

  std::string_view name;
  cosmosim::core::IntegrationStage stage;
  cosmosim::core::StageContract contract{
      .stage = cosmosim::core::IntegrationStage::kDrift,
      .required_inputs = cosmosim::core::StageDataDomain::kParticles,
      .mutated_state = cosmosim::core::StageDataDomain::kDiagnostics,
      .produced_outputs = cosmosim::core::StageDataDomain::kDiagnostics,
      .allowed_side_effects = cosmosim::core::StageDataDomain::kDiagnostics,
      .sync_requirements = cosmosim::core::StageSyncRequirement::kLocalOnly,
      .active_set_family = cosmosim::core::StageActiveSetFamily::kActiveParticles,
      .restart_safety = cosmosim::core::StageSafety::kSafe,
      .output_safety = cosmosim::core::StageSafety::kSafe,
      .owner = cosmosim::core::StageSubsystem::kCore,
  };
  std::vector<std::string>* order = nullptr;
  int invocations = 0;
  std::vector<cosmosim::core::IntegrationStage> observed;
};

class StageRecorder final : public cosmosim::core::IntegrationCallback {
 public:
  std::string_view callbackName() const override { return "stage_recorder"; }
  std::span<const cosmosim::core::IntegrationStage> integrationStages() const override {
    static constexpr std::array stages{
        cosmosim::core::IntegrationStage::kGravityKickPre,
        cosmosim::core::IntegrationStage::kDrift,
        cosmosim::core::IntegrationStage::kForceRefresh,
        cosmosim::core::IntegrationStage::kHydroUpdate,
        cosmosim::core::IntegrationStage::kSourceTerms,
        cosmosim::core::IntegrationStage::kGravityKickPost,
        cosmosim::core::IntegrationStage::kAnalysisHooks,
        cosmosim::core::IntegrationStage::kOutputCheck,
    };
    return stages;
  }
  std::span<const cosmosim::core::StageContract> stageContracts() const override { return contracts; }

  void onStage(cosmosim::core::StepContext& context) override { observed_stages.push_back(context.stage); }

  std::vector<cosmosim::core::IntegrationStage> observed_stages;
  std::array<cosmosim::core::StageContract, 8> contracts{{
      {.stage = cosmosim::core::IntegrationStage::kGravityKickPre, .required_inputs = cosmosim::core::StageDataDomain::kParticles, .mutated_state = cosmosim::core::StageDataDomain::kPmField, .produced_outputs = cosmosim::core::StageDataDomain::kPmField, .allowed_side_effects = cosmosim::core::StageDataDomain::kPmField, .sync_requirements = cosmosim::core::StageSyncRequirement::kGlobal, .active_set_family = cosmosim::core::StageActiveSetFamily::kAllParticles, .restart_safety = cosmosim::core::StageSafety::kSafe, .output_safety = cosmosim::core::StageSafety::kSafe, .owner = cosmosim::core::StageSubsystem::kGravity},
      {.stage = cosmosim::core::IntegrationStage::kDrift, .required_inputs = cosmosim::core::StageDataDomain::kParticles, .mutated_state = cosmosim::core::StageDataDomain::kParticles, .produced_outputs = cosmosim::core::StageDataDomain::kParticles, .allowed_side_effects = cosmosim::core::StageDataDomain::kNone, .sync_requirements = cosmosim::core::StageSyncRequirement::kLocalOnly, .active_set_family = cosmosim::core::StageActiveSetFamily::kActiveParticles, .restart_safety = cosmosim::core::StageSafety::kSafe, .output_safety = cosmosim::core::StageSafety::kSafe, .owner = cosmosim::core::StageSubsystem::kCore},
      {.stage = cosmosim::core::IntegrationStage::kForceRefresh, .required_inputs = cosmosim::core::StageDataDomain::kPmField, .mutated_state = cosmosim::core::StageDataDomain::kPmField, .produced_outputs = cosmosim::core::StageDataDomain::kPmField, .allowed_side_effects = cosmosim::core::StageDataDomain::kPmField, .sync_requirements = cosmosim::core::StageSyncRequirement::kForceEvaluation, .active_set_family = cosmosim::core::StageActiveSetFamily::kAllParticles, .restart_safety = cosmosim::core::StageSafety::kSafe, .output_safety = cosmosim::core::StageSafety::kSafe, .owner = cosmosim::core::StageSubsystem::kGravity},
      {.stage = cosmosim::core::IntegrationStage::kHydroUpdate, .required_inputs = cosmosim::core::StageDataDomain::kGasCells, .mutated_state = cosmosim::core::StageDataDomain::kGasCells, .produced_outputs = cosmosim::core::StageDataDomain::kGasCells, .allowed_side_effects = cosmosim::core::StageDataDomain::kNone, .sync_requirements = cosmosim::core::StageSyncRequirement::kLocalOnly, .active_set_family = cosmosim::core::StageActiveSetFamily::kGasCells, .restart_safety = cosmosim::core::StageSafety::kSafe, .output_safety = cosmosim::core::StageSafety::kSafe, .owner = cosmosim::core::StageSubsystem::kHydro},
      {.stage = cosmosim::core::IntegrationStage::kSourceTerms, .required_inputs = cosmosim::core::StageDataDomain::kParticles, .mutated_state = cosmosim::core::StageDataDomain::kParticles, .produced_outputs = cosmosim::core::StageDataDomain::kParticles, .allowed_side_effects = cosmosim::core::StageDataDomain::kDiagnostics, .sync_requirements = cosmosim::core::StageSyncRequirement::kLocalOnly, .active_set_family = cosmosim::core::StageActiveSetFamily::kActiveParticles, .restart_safety = cosmosim::core::StageSafety::kSafe, .output_safety = cosmosim::core::StageSafety::kSafe, .owner = cosmosim::core::StageSubsystem::kSources},
      {.stage = cosmosim::core::IntegrationStage::kGravityKickPost, .required_inputs = cosmosim::core::StageDataDomain::kParticles, .mutated_state = cosmosim::core::StageDataDomain::kParticles, .produced_outputs = cosmosim::core::StageDataDomain::kParticles, .allowed_side_effects = cosmosim::core::StageDataDomain::kNone, .sync_requirements = cosmosim::core::StageSyncRequirement::kGlobal, .active_set_family = cosmosim::core::StageActiveSetFamily::kAllParticles, .restart_safety = cosmosim::core::StageSafety::kSafe, .output_safety = cosmosim::core::StageSafety::kSafe, .owner = cosmosim::core::StageSubsystem::kGravity},
      {.stage = cosmosim::core::IntegrationStage::kAnalysisHooks, .required_inputs = cosmosim::core::StageDataDomain::kDiagnostics, .mutated_state = cosmosim::core::StageDataDomain::kDiagnostics, .produced_outputs = cosmosim::core::StageDataDomain::kDiagnostics, .allowed_side_effects = cosmosim::core::StageDataDomain::kDiagnostics, .sync_requirements = cosmosim::core::StageSyncRequirement::kLocalOnly, .active_set_family = cosmosim::core::StageActiveSetFamily::kNone, .restart_safety = cosmosim::core::StageSafety::kSafe, .output_safety = cosmosim::core::StageSafety::kSafe, .owner = cosmosim::core::StageSubsystem::kAnalysis},
      {.stage = cosmosim::core::IntegrationStage::kOutputCheck, .required_inputs = cosmosim::core::StageDataDomain::kOutputState, .mutated_state = cosmosim::core::StageDataDomain::kOutputState, .produced_outputs = cosmosim::core::StageDataDomain::kOutputState, .allowed_side_effects = cosmosim::core::StageDataDomain::kOutputState, .sync_requirements = cosmosim::core::StageSyncRequirement::kGlobal, .active_set_family = cosmosim::core::StageActiveSetFamily::kOutputState, .restart_safety = cosmosim::core::StageSafety::kSafe, .output_safety = cosmosim::core::StageSafety::kSafe, .owner = cosmosim::core::StageSubsystem::kOutput},
  }};
};

class PmDirectiveRecorder final : public cosmosim::core::IntegrationCallback {
 public:
  std::string_view callbackName() const override { return "pm_directive_recorder"; }
  std::span<const cosmosim::core::IntegrationStage> integrationStages() const override {
    static constexpr std::array stages{
        cosmosim::core::IntegrationStage::kGravityKickPre,
        cosmosim::core::IntegrationStage::kForceRefresh,
    };
    return stages;
  }
  std::span<const cosmosim::core::StageContract> stageContracts() const override { return contracts; }
  void onStage(cosmosim::core::StepContext& context) override {
    if (context.stage == cosmosim::core::IntegrationStage::kGravityKickPre) {
      kick_pre = context.pm_refresh_directive;
      if (context.pm_refresh_directive.has_sync_event) {
        context.pm_refresh_directive.solver_executed = true;
      }
    } else if (context.stage == cosmosim::core::IntegrationStage::kForceRefresh) {
      force_refresh = context.pm_refresh_directive;
      if (context.pm_refresh_directive.has_sync_event) {
        context.pm_refresh_directive.solver_executed = true;
      }
    }
  }

  std::array<cosmosim::core::StageContract, 2> contracts{{
      {.stage = cosmosim::core::IntegrationStage::kGravityKickPre,
       .required_inputs = cosmosim::core::StageDataDomain::kParticles,
       .mutated_state = cosmosim::core::StageDataDomain::kPmField | cosmosim::core::StageDataDomain::kDiagnostics,
       .produced_outputs = cosmosim::core::StageDataDomain::kPmField | cosmosim::core::StageDataDomain::kDiagnostics,
       .allowed_side_effects = cosmosim::core::StageDataDomain::kDiagnostics,
       .sync_requirements = cosmosim::core::StageSyncRequirement::kGlobal,
       .active_set_family = cosmosim::core::StageActiveSetFamily::kAllParticles,
       .restart_safety = cosmosim::core::StageSafety::kSafe,
       .output_safety = cosmosim::core::StageSafety::kSafe,
       .owner = cosmosim::core::StageSubsystem::kGravity},
      {.stage = cosmosim::core::IntegrationStage::kForceRefresh,
       .required_inputs = cosmosim::core::StageDataDomain::kParticles | cosmosim::core::StageDataDomain::kPmField,
       .mutated_state = cosmosim::core::StageDataDomain::kPmField | cosmosim::core::StageDataDomain::kTreeState | cosmosim::core::StageDataDomain::kDiagnostics,
       .produced_outputs = cosmosim::core::StageDataDomain::kPmField | cosmosim::core::StageDataDomain::kTreeState | cosmosim::core::StageDataDomain::kDiagnostics,
       .allowed_side_effects = cosmosim::core::StageDataDomain::kDiagnostics,
       .sync_requirements = cosmosim::core::StageSyncRequirement::kForceEvaluation,
       .active_set_family = cosmosim::core::StageActiveSetFamily::kActiveParticles,
       .restart_safety = cosmosim::core::StageSafety::kSafe,
       .output_safety = cosmosim::core::StageSafety::kSafe,
       .owner = cosmosim::core::StageSubsystem::kGravity},
  }};
  cosmosim::core::PmRefreshDirective kick_pre{};
  cosmosim::core::PmRefreshDirective force_refresh{};
};

void testKickDriftKickOrdering() {
  cosmosim::core::SimulationState state;
  cosmosim::core::IntegratorState integrator_state;
  cosmosim::core::TransientStepWorkspace workspace;
  integrator_state.dt_time_code = 1.0;

  StageRecorder recorder;
  cosmosim::core::StepOrchestrator orchestrator;
  orchestrator.registerCallback(recorder);

  const cosmosim::core::ActiveSetDescriptor active_set{};
  orchestrator.executeSingleStep(state, integrator_state, active_set, nullptr, &workspace);
  orchestrator.executeOutputBoundary(
      state,
      integrator_state,
      nullptr,
      cosmosim::core::StepBoundaryKind::kGlobalSynchronizationPoint);

  const auto expected = cosmosim::core::StageScheduler::kickDriftKickOrder();
  assert(cosmosim::core::isCanonicalIntegrationStageOrder(expected));
  assert(recorder.observed_stages.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    assert(recorder.observed_stages[i] == expected[i]);
  }
  assert(cosmosim::core::isCanonicalIntegrationStageOrder(recorder.observed_stages));
}

void testPmRefreshDirectiveCapturesReasonAndForceEvalTime() {
  cosmosim::core::SimulationState state;
  cosmosim::core::IntegratorState integrator_state;
  cosmosim::core::TransientStepWorkspace workspace;
  integrator_state.dt_time_code = 0.25;
  integrator_state.pm_refresh_enabled = true;
  integrator_state.pm_long_range_field_valid = false;
  integrator_state.current_scale_factor = 1.0;

  PmDirectiveRecorder recorder;
  cosmosim::core::StepOrchestrator orchestrator;
  orchestrator.registerCallback(recorder);
  orchestrator.executeSingleStep(state, integrator_state, cosmosim::core::ActiveSetDescriptor{}, nullptr, &workspace);

  assert(recorder.kick_pre.reason == cosmosim::core::PmRefreshDirective::Reason::kInitialForceBootstrap);
  assert(recorder.kick_pre.sync_stage == cosmosim::core::PmSyncStage::kInitialLongRangeBootstrap);
  assert(recorder.kick_pre.force_evaluation_scale_factor > 0.0);
  assert(recorder.force_refresh.reason == cosmosim::core::PmRefreshDirective::Reason::kScheduledForceRefreshStage);
  assert(recorder.force_refresh.sync_stage == cosmosim::core::PmSyncStage::kScheduledLongRangeRefresh);
  assert(recorder.force_refresh.force_refresh_surface);
  assert(recorder.force_refresh.force_evaluation_scale_factor >= recorder.kick_pre.force_evaluation_scale_factor);
}


void testLocalForceRefreshDoesNotIssuePmSyncEvent() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(4);

  std::vector<std::uint32_t> active_particles{0, 2};
  cosmosim::core::ActiveSetDescriptor local_active{
      .particle_indices = active_particles,
      .particles_are_subset = true,
      .particles_from_scheduler = true,
      .has_generation_metadata = true,
      .source_particle_index_generation = state.particleIndexGeneration(),
      .source_cell_index_generation = state.cellIndexGeneration(),
      .source_scheduler_tick = 7,
  };

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 0.125;
  integrator_state.pm_refresh_enabled = true;
  integrator_state.pm_long_range_field_valid = true;
  integrator_state.pm_sync_state.importPersistentState(cosmosim::core::PmSynchronizationPersistentState{
      .cadence_steps = 1,
      .gravity_kick_opportunity = 3,
      .last_refresh_opportunity = 3,
      .field_version = 2,
      .last_refresh_step_index = 6,
      .last_refresh_scale_factor = 1.0,
  });

  PmDirectiveRecorder recorder;
  cosmosim::core::TransientStepWorkspace workspace;
  cosmosim::core::StepOrchestrator orchestrator;
  orchestrator.registerCallback(recorder);
  orchestrator.executeSingleStep(
      state,
      integrator_state,
      local_active,
      nullptr,
      &workspace,
      nullptr,
      nullptr,
      std::uint64_t{7},
      cosmosim::core::StepBoundaryKind::kPmRefreshPoint);

  assert(recorder.force_refresh.force_refresh_surface);
  assert(recorder.force_refresh.requires_predicted_inactive_sources);
  assert(!recorder.force_refresh.cadence_opportunity_allowed);
  assert(!recorder.force_refresh.has_sync_event);
  assert(recorder.force_refresh.sync_stage == cosmosim::core::PmSyncStage::kNone);
  assert(!recorder.force_refresh.refresh_long_range_field);
}


void testStageBoundDispatch() {
  cosmosim::core::SimulationState state;
  cosmosim::core::IntegratorState integrator_state;
  cosmosim::core::TransientStepWorkspace workspace;
  integrator_state.dt_time_code = 1.0;

  std::vector<std::string> order;
  SingleStageRecorder drift_first("drift_first", cosmosim::core::IntegrationStage::kDrift, &order);
  SingleStageRecorder source("source", cosmosim::core::IntegrationStage::kSourceTerms, &order);
  SingleStageRecorder drift_second("drift_second", cosmosim::core::IntegrationStage::kDrift, &order);

  cosmosim::core::StepOrchestrator orchestrator;
  orchestrator.registerCallback(drift_first);
  orchestrator.registerCallback(source);
  orchestrator.registerCallback(drift_second);

  assert(orchestrator.callbackCount() == 3U);
  assert(orchestrator.handlersFor(cosmosim::core::IntegrationStage::kDrift).size() == 2U);
  assert(orchestrator.contractsFor(cosmosim::core::IntegrationStage::kDrift).size() == 2U);
  assert(orchestrator.handlersFor(cosmosim::core::IntegrationStage::kHydroUpdate).empty());

  orchestrator.executeSingleStep(state, integrator_state, {}, nullptr, &workspace);

  assert(drift_first.invocations == 1);
  assert(drift_second.invocations == 1);
  assert(source.invocations == 1);
  assert(drift_first.observed.front() == cosmosim::core::IntegrationStage::kDrift);
  assert(drift_second.observed.front() == cosmosim::core::IntegrationStage::kDrift);
  assert(source.observed.front() == cosmosim::core::IntegrationStage::kSourceTerms);
  assert((order == std::vector<std::string>{"drift_first", "drift_second", "source"}));
}


void testRegisterCallbackRejectsExtraUnregisteredStageContract() {
  class MismatchCallback final : public cosmosim::core::IntegrationCallback {
   public:
    std::string_view callbackName() const override { return "mismatch_callback"; }
    std::span<const cosmosim::core::IntegrationStage> integrationStages() const override {
      static constexpr std::array stages{cosmosim::core::IntegrationStage::kDrift};
      return stages;
    }
    std::span<const cosmosim::core::StageContract> stageContracts() const override { return contracts; }
    void onStage(cosmosim::core::StepContext&) override {}

    std::array<cosmosim::core::StageContract, 2> contracts{{
        {.stage = cosmosim::core::IntegrationStage::kDrift, .restart_safety = cosmosim::core::StageSafety::kSafe, .output_safety = cosmosim::core::StageSafety::kSafe},
        {.stage = cosmosim::core::IntegrationStage::kSourceTerms, .restart_safety = cosmosim::core::StageSafety::kSafe, .output_safety = cosmosim::core::StageSafety::kSafe},
    }};
  } callback;

  cosmosim::core::StepOrchestrator orchestrator;
  assert(throwsWithContext(
      [&]() { orchestrator.registerCallback(callback); },
      "must declare exactly one executable contract for each registered stage"));
}

void testRegisterCallbackRejectsDuplicateStageContracts() {
  class DuplicateCallback final : public cosmosim::core::IntegrationCallback {
   public:
    std::string_view callbackName() const override { return "duplicate_contract_callback"; }
    std::span<const cosmosim::core::IntegrationStage> integrationStages() const override {
      static constexpr std::array stages{
          cosmosim::core::IntegrationStage::kDrift,
          cosmosim::core::IntegrationStage::kSourceTerms,
      };
      return stages;
    }
    std::span<const cosmosim::core::StageContract> stageContracts() const override { return contracts; }
    void onStage(cosmosim::core::StepContext&) override {}

    std::array<cosmosim::core::StageContract, 2> contracts{{
        {.stage = cosmosim::core::IntegrationStage::kDrift,
         .required_inputs = cosmosim::core::StageDataDomain::kParticles,
         .mutated_state = cosmosim::core::StageDataDomain::kParticles,
         .produced_outputs = cosmosim::core::StageDataDomain::kParticles,
         .sync_requirements = cosmosim::core::StageSyncRequirement::kLocalOnly,
         .active_set_family = cosmosim::core::StageActiveSetFamily::kActiveParticles,
         .restart_safety = cosmosim::core::StageSafety::kSafe,
         .output_safety = cosmosim::core::StageSafety::kSafe,
         .owner = cosmosim::core::StageSubsystem::kCore},
        {.stage = cosmosim::core::IntegrationStage::kDrift,
         .required_inputs = cosmosim::core::StageDataDomain::kParticles,
         .mutated_state = cosmosim::core::StageDataDomain::kParticles,
         .produced_outputs = cosmosim::core::StageDataDomain::kParticles,
         .sync_requirements = cosmosim::core::StageSyncRequirement::kLocalOnly,
         .active_set_family = cosmosim::core::StageActiveSetFamily::kActiveParticles,
         .restart_safety = cosmosim::core::StageSafety::kSafe,
         .output_safety = cosmosim::core::StageSafety::kSafe,
         .owner = cosmosim::core::StageSubsystem::kCore},
    }};
  } callback;

  cosmosim::core::StepOrchestrator orchestrator;
  assert(throwsWithContext(
      [&]() { orchestrator.registerCallback(callback); },
      "duplicate executable contracts"));
}


void testRegisterCallbackRejectsIncompleteExecutableContract() {
  class IncompleteContractCallback final : public cosmosim::core::IntegrationCallback {
   public:
    std::string_view callbackName() const override { return "incomplete_contract_callback"; }
    std::span<const cosmosim::core::IntegrationStage> integrationStages() const override {
      static constexpr std::array stages{cosmosim::core::IntegrationStage::kDrift};
      return stages;
    }
    std::span<const cosmosim::core::StageContract> stageContracts() const override { return contracts; }
    void onStage(cosmosim::core::StepContext&) override {}

    std::array<cosmosim::core::StageContract, 1> contracts{{
        {.stage = cosmosim::core::IntegrationStage::kDrift,
         .restart_safety = cosmosim::core::StageSafety::kSafe,
         .output_safety = cosmosim::core::StageSafety::kSafe,
         .owner = cosmosim::core::StageSubsystem::kCore},
    }};
  } callback;

  cosmosim::core::StepOrchestrator orchestrator;
  assert(throwsWithContext(
      [&]() { orchestrator.registerCallback(callback); },
      "must declare required_inputs"));
}

void testLocalPmSyncEventIsHardRejectedAfterHandlerMutation() {
  class IllegalPmSyncCallback final : public cosmosim::core::IntegrationCallback {
   public:
    std::string_view callbackName() const override { return "illegal_pm_sync_callback"; }
    std::span<const cosmosim::core::IntegrationStage> integrationStages() const override {
      static constexpr std::array stages{cosmosim::core::IntegrationStage::kForceRefresh};
      return stages;
    }
    std::span<const cosmosim::core::StageContract> stageContracts() const override { return contracts; }
    void onStage(cosmosim::core::StepContext& context) override {
      context.pm_refresh_directive.has_sync_event = true;
      context.pm_refresh_directive.refresh_long_range_field = true;
      context.pm_refresh_directive.cadence_opportunity_allowed = true;
      context.pm_refresh_directive.sync_stage = cosmosim::core::PmSyncStage::kScheduledLongRangeRefresh;
      context.pm_refresh_directive.force_evaluation_scale_factor = 1.0;
      context.pm_refresh_directive.solver_executed = true;
    }

    std::array<cosmosim::core::StageContract, 1> contracts{{
        {.stage = cosmosim::core::IntegrationStage::kForceRefresh,
         .required_inputs = cosmosim::core::StageDataDomain::kParticles | cosmosim::core::StageDataDomain::kPmField,
         .mutated_state = cosmosim::core::StageDataDomain::kPmField | cosmosim::core::StageDataDomain::kTreeState | cosmosim::core::StageDataDomain::kDiagnostics,
         .produced_outputs = cosmosim::core::StageDataDomain::kPmField | cosmosim::core::StageDataDomain::kTreeState | cosmosim::core::StageDataDomain::kDiagnostics,
         .allowed_side_effects = cosmosim::core::StageDataDomain::kDiagnostics,
         .sync_requirements = cosmosim::core::StageSyncRequirement::kForceEvaluation,
         .active_set_family = cosmosim::core::StageActiveSetFamily::kActiveParticles,
         .restart_safety = cosmosim::core::StageSafety::kUnsafe,
         .output_safety = cosmosim::core::StageSafety::kUnsafe,
         .owner = cosmosim::core::StageSubsystem::kGravity},
    }};
  } callback;

  cosmosim::core::SimulationState state;
  state.resizeParticles(4);
  std::vector<std::uint32_t> active_particles{0, 2};
  cosmosim::core::ActiveSetDescriptor local_active{
      .particle_indices = active_particles,
      .particles_are_subset = true,
      .particles_from_scheduler = true,
      .has_generation_metadata = true,
      .source_particle_index_generation = state.particleIndexGeneration(),
      .source_cell_index_generation = state.cellIndexGeneration(),
      .source_scheduler_tick = 9,
  };

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 0.125;
  integrator_state.pm_refresh_enabled = true;
  integrator_state.pm_long_range_field_valid = true;
  integrator_state.pm_sync_state.importPersistentState(cosmosim::core::PmSynchronizationPersistentState{
      .cadence_steps = 1,
      .gravity_kick_opportunity = 4,
      .last_refresh_opportunity = 4,
      .field_version = 3,
      .last_refresh_step_index = 8,
      .last_refresh_scale_factor = 1.0,
  });

  cosmosim::core::StepOrchestrator orchestrator;
  cosmosim::core::TransientStepWorkspace workspace;
  orchestrator.registerCallback(callback);
  assert(throwsWithContext(
      [&]() {
        orchestrator.executeSingleStep(
            state,
            integrator_state,
            local_active,
            nullptr,
            &workspace,
            nullptr,
            nullptr,
            std::uint64_t{9},
            cosmosim::core::StepBoundaryKind::kPmRefreshPoint);
      },
      "PM sync event reached an illegal local integration boundary"));
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
  const double a_next = cosmosim::core::advanceScaleFactorByCosmicTime(background, a0, dt, 128);
  assert(std::abs(a_next - (a0 + delta_a)) / delta_a < 2.0e-8);

  const cosmosim::core::CosmologicalTimeline timeline(&background);
  const auto step = timeline.prepareStep(0.0, a0, dt);
  assert(step.cosmological);
  assert(step.scale_factor_end > a0);
  assert(step.redshift_end < step.redshift_begin);
  assert(step.drift_factor_code > 0.0);
  assert(step.first_kick_factor_code > 0.0);
  assert(step.second_kick_factor_code > 0.0);
  assert(step.first_hubble_drag_factor < 1.0);
  assert(step.second_hubble_drag_factor < 1.0);
  assert(std::abs(step.first_hubble_drag_factor * step.second_hubble_drag_factor - step.hubble_drag_factor) < 1.0e-12);
  assert(step.first_kick_factor_code < step.drift_factor_code);
  assert(step.second_kick_factor_code < step.drift_factor_code);
  const double combined_force_response =
      step.second_hubble_drag_factor * step.first_kick_factor_code + step.second_kick_factor_code;
  assert(std::abs(combined_force_response -
                  cosmosim::core::computeComovingKickFactor(background, a0, step.scale_factor_end, 128)) /
         combined_force_response < 1.0e-5);
}


void testCosmologicalTimelineConvertsSiIntegralsToCodeTime() {
  cosmosim::core::CosmologyBackgroundConfig cfg;
  cfg.hubble_param = 0.7;
  cfg.omega_matter = 1.0;
  cfg.omega_lambda = 0.0;
  cfg.omega_radiation = 0.0;
  cfg.omega_curvature = 0.0;
  const cosmosim::core::LambdaCdmBackground background(cfg);

  const double a0 = 0.5;
  const double dt_si = cosmosim::core::estimateDeltaTimeFromScaleFactorStep(background, a0, 1.0e-4);
  const double time_si_per_code = 10.0;
  const cosmosim::core::CosmologicalTimeline timeline(&background, time_si_per_code);
  const auto step = timeline.prepareStep(0.0, a0, dt_si / time_si_per_code);

  assert(std::abs(step.dt_time_si - dt_si) / dt_si < 1.0e-14);
  assert(std::abs(step.hubble_begin_code - background.hubbleSi(a0) * time_si_per_code) <
         1.0e-15 * background.hubbleSi(a0) * time_si_per_code);
  assert(std::abs(step.drift_factor_code * time_si_per_code -
                  cosmosim::core::computeComovingDriftFactor(background, a0, step.scale_factor_end, 128)) /
         (step.drift_factor_code * time_si_per_code) < 1.0e-5);
  assert(std::abs(step.first_hubble_drag_factor * step.second_hubble_drag_factor -
                  cosmosim::core::computeHubbleDragFactor(a0, step.scale_factor_end)) < 1.0e-12);
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

  const cosmosim::core::DirectionalCflTimeStepInput directional_cfl{
      .cell_width_axis_code = {0.5, 0.125, 0.25},
      .velocity_axis_code = {18.0, 1.0, 4.0},
      .sound_speed_code = 5.0,
  };
  const double dt_directional = cosmosim::core::computeDirectionalCflTimeStep(directional_cfl, 0.4);
  assert(std::abs(dt_directional - (0.4 * 0.125 / 6.0)) < k_tolerance);
  const auto cfl_diagnostics = cosmosim::core::makeHydroCflDiagnostics(
      7,
      directional_cfl,
      0.4,
      0.005,
      9007,
      std::uint64_t{42},
      std::uint32_t{3});
  assert(cfl_diagnostics.local_row == 7);
  assert(cfl_diagnostics.gas_cell_id == 9007);
  assert(cfl_diagnostics.patch_id == 42);
  assert(cfl_diagnostics.patch_row == 3);
  assert(cfl_diagnostics.limiting_axis == 1);
  cosmosim::core::assertHydroCflStable(cfl_diagnostics);
  assert(throwsWithContext(
      [&]() {
        const auto unstable = cosmosim::core::makeHydroCflDiagnostics(
            7,
            directional_cfl,
            0.4,
            0.02,
            9007,
            std::uint64_t{42},
            std::uint32_t{3});
        cosmosim::core::assertHydroCflStable(unstable);
      },
      "hydro CFL violation"));

  const cosmosim::core::GravityTimeStepInput grav_input{
      .softening_length_code = 0.01,
      .acceleration_magnitude_code = 4.0,
  };
  const double dt_grav = cosmosim::core::computeGravityTimeStep(grav_input, 0.4);
  assert(std::abs(dt_grav - 0.4 * std::sqrt(0.01 / 4.0)) < k_tolerance);
  const cosmosim::core::ComovingGravityTimeStepInput comoving_grav_input{
      .softening_length_comoving_code = 0.01,
      .scale_free_acceleration_magnitude_code = 4.0,
      .scale_factor = 0.5,
  };
  const double dt_comoving_grav =
      cosmosim::core::computeComovingGravityTimeStep(comoving_grav_input, 0.4);
  assert(std::abs(dt_comoving_grav - 0.4 * std::sqrt(0.125 * 0.01 / 4.0)) <
      k_tolerance);
  assert(throwsWithContext(
      [&]() {
        (void)cosmosim::core::computeComovingGravityTimeStep(
            {.softening_length_comoving_code = 0.01,
             .scale_free_acceleration_magnitude_code = 4.0,
             .scale_factor = 0.0},
            0.4);
      },
      "scale_factor"));

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

  scheduler.submitCandidateBin(0, 0, cosmosim::core::TimeStepCandidateSource::kUserClamp);
  scheduler.submitCandidateBin(1, 3, cosmosim::core::TimeStepCandidateSource::kUserClamp);
  scheduler.endSubstep();

  const auto& hot = scheduler.hotMetadata();
  assert(hot.bin_index[0] == 0);
  assert(hot.bin_index[1] == 3);

  active = scheduler.beginSubstep();
  assert(active.size() == 1);
  assert(active[0] == 0);

  scheduler.endSubstep();

  while (scheduler.currentTick() < 8) {
    scheduler.beginSubstep();
    scheduler.endSubstep();
  }

  active = scheduler.beginSubstep();
  assert(!active.empty() && active[0] == 0);
  scheduler.submitCandidateBin(0, 3, cosmosim::core::TimeStepCandidateSource::kUserClamp);
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
  initializeParticleIdsAndGasCells(state, {0, 1});
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
  scheduler.submitCandidateBin(0, 0, cosmosim::core::TimeStepCandidateSource::kUserClamp);
  scheduler.submitCandidateBin(3, 2, cosmosim::core::TimeStepCandidateSource::kUserClamp);
  scheduler.endSubstep();

  while (scheduler.currentTick() < 4) {
    scheduler.beginSubstep();
    scheduler.endSubstep();
  }

  scheduler.beginSubstep();
  scheduler.submitCandidateBin(0, 2, cosmosim::core::TimeStepCandidateSource::kUserClamp);
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

  scheduler.submitCandidateBin(2, 0, cosmosim::core::TimeStepCandidateSource::kUserClamp);
  scheduler.endSubstep();

  const auto active1 = scheduler.beginSubstep();
  const std::vector<std::uint32_t> expected1{0, 2, 4};
  assert(std::vector<std::uint32_t>(active1.begin(), active1.end()) == expected1);
  scheduler.endSubstep();
}

void testActiveSetNoCompetingBuilders() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(6);
  state.resizeCells(2);
  initializeParticleIdsAndGasCells(state, {0, 1});
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


void testHydroGravityCandidateReconciliation() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(3, 3, 0);
  const cosmosim::core::TimeStepLimits limits{
      .min_dt_time_code = 0.125,
      .max_dt_time_code = 1.0,
      .max_bin = 3,
  };

  scheduler.submitCandidateTimeStep(0, 1.0, limits, cosmosim::core::TimeStepCandidateSource::kHydroCfl, "hydro_coarse");
  scheduler.submitCandidateTimeStep(0, 0.125, limits, cosmosim::core::TimeStepCandidateSource::kGravityAcceleration, "gravity_fine");
  scheduler.submitCandidateTimeStep(1, 0.5, limits, cosmosim::core::TimeStepCandidateSource::kHydroCfl, "hydro_mid");

  const auto reconciliation = scheduler.reconcileCandidateTransitions();
  assert(reconciliation.submitted_candidates == 3);
  assert(reconciliation.elements_with_candidates == 2);
  assert(reconciliation.committed_transition_requests == 2);

  scheduler.beginSubstep();
  scheduler.endSubstep();
  const auto& hot = scheduler.hotMetadata();
  assert(hot.bin_index[0] == 0);
  assert(hot.bin_index[1] == 2);
  assert(hot.bin_index[2] == 3);
}

void testHighMachHydroCflCandidateLimitsBin() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(4);
  scheduler.reset(1, 4, 0);
  const cosmosim::core::TimeStepLimits limits{
      .min_dt_time_code = 0.01,
      .max_dt_time_code = 0.16,
      .max_bin = 4,
  };
  const cosmosim::core::DirectionalCflTimeStepInput high_mach_cell{
      .cell_width_axis_code = {1.0, 1.0, 0.5},
      .velocity_axis_code = {2.0, 3.0, 49.0},
      .sound_speed_code = 1.0,
  };
  const double hydro_dt = cosmosim::core::computeDirectionalCflTimeStep(high_mach_cell, 0.4);
  assert(std::abs(hydro_dt - 0.004) < k_tolerance);

  scheduler.submitCandidateTimeStep(
      0,
      hydro_dt,
      limits,
      cosmosim::core::TimeStepCandidateSource::kHydroCfl,
      "high_mach_cell_hydro_cfl");
  scheduler.submitCandidateTimeStep(
      0,
      0.08,
      limits,
      cosmosim::core::TimeStepCandidateSource::kGravityAcceleration,
      "gravity_looser");
  const auto reconciliation = scheduler.reconcileCandidateTransitions();
  assert(reconciliation.submitted_candidates == 2);
  assert(reconciliation.limiting_candidates_by_source[
             static_cast<std::size_t>(cosmosim::core::TimeStepCandidateSource::kHydroCfl)] == 1);
  assert(reconciliation.dominant_limiting_source == cosmosim::core::TimeStepCandidateSource::kHydroCfl);
  scheduler.beginSubstep();
  scheduler.endSubstep();
  assert(scheduler.hotMetadata().bin_index[0] == 0);
}

void testOrchestratorRequiresSchedulerProvenanceWhenTickIsExpected() {
  cosmosim::core::SimulationState state;
  cosmosim::core::TransientStepWorkspace workspace;
  state.resizeParticles(1);
  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 1.0;
  std::vector<std::uint32_t> active_particles{0};
  cosmosim::core::ActiveSetDescriptor untrusted{
      .particle_indices = active_particles,
      .particles_are_subset = true,
  };

  cosmosim::core::StepOrchestrator orchestrator;
  bool threw = false;
  try {
    orchestrator.executeSingleStep(state, integrator_state, untrusted, nullptr, &workspace, nullptr, nullptr, 0);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}

void testOrchestratorRejectsSchedulerActiveSetWithoutTick() {
  cosmosim::core::SimulationState state;
  cosmosim::core::TransientStepWorkspace workspace;
  state.resizeParticles(2);
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(1);
  scheduler.reset(2, 0, 0);
  const auto active = scheduler.beginSubstep();
  const auto active_set = cosmosim::core::makeSchedulerActiveSetDescriptor(scheduler, state, active);

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 1.0;
  cosmosim::core::StepOrchestrator orchestrator;
  assert(throwsWithContext(
      [&]() { orchestrator.executeSingleStep(state, integrator_state, active_set, nullptr, &workspace); },
      "scheduler-derived active sets require explicit scheduler tick"));

  orchestrator.executeSchedulerSubstep(state, integrator_state, scheduler, active, {}, nullptr, &workspace);
  scheduler.endSubstep();
}

void testSchedulerBackedTimeBinReorderIgnoresStaleMirrors() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(4);
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(4, 0, 0);
  scheduler.setElementBin(0, 3, 0);
  scheduler.setElementBin(1, 1, 0);
  scheduler.setElementBin(2, 2, 0);
  scheduler.setElementBin(3, 0, 0);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particles.time_bin[i] = 0;  // Deliberately stale mirror: production reorder must not consume it.
  }

  bool legacy_mirror_reorder_threw = false;
  try {
    (void)cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kByTimeBin);
  } catch (const std::invalid_argument&) {
    legacy_mirror_reorder_threw = true;
  }
  assert(legacy_mirror_reorder_threw);

  const auto reorder = cosmosim::core::buildParticleReorderMapByScheduler(state, scheduler);
  assert(reorder.isConsistent(state.particles.size()));
  assert(reorder.new_to_old_index[0] == 3);
  assert(reorder.new_to_old_index[1] == 1);
  assert(reorder.new_to_old_index[2] == 2);
  assert(reorder.new_to_old_index[3] == 0);
}

void testPmSynchronizationPersistentRoundTrip() {
  cosmosim::core::PmSynchronizationState pm_sync;
  pm_sync.reset(2);
  const auto first = pm_sync.registerKickOpportunity(7, 0.25, false);
  pm_sync.commitRefresh(first);
  (void)pm_sync.registerKickOpportunity(8, 0.5, true);
  const auto pending = pm_sync.registerKickOpportunity(9, 0.75, true);
  assert(pending.refresh_long_range_field);

  const auto saved = pm_sync.exportPersistentState();
  cosmosim::core::PmSynchronizationState restored;
  restored.importPersistentState(saved);
  assert(restored.cadenceSteps() == 2);
  assert(restored.gravityKickOpportunity() == pm_sync.gravityKickOpportunity());
  assert(restored.fieldVersion() == pm_sync.fieldVersion());
  assert(restored.refreshCommitPending());
  restored.commitRefresh(pending);
  assert(restored.fieldVersion() == 2);
  assert(restored.lastRefreshOpportunity() == pending.last_refresh_opportunity);
}

void testInvalidEarlyActivationTrap() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(2, 2, 0);
  auto state = scheduler.exportPersistentState();
  state.current_tick = 1;
  state.next_activation_tick[0] = 0;
  cosmosim::core::HierarchicalTimeBinScheduler restored(state.max_bin);
  assert(throwsWithContext([&]() { restored.importPersistentState(state); }, "element=0"));
}

void testInvalidBinJumpTrap() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(1, 0, 0);
  scheduler.beginSubstep();
  scheduler.endSubstep();
  scheduler.beginSubstep();
  scheduler.submitCandidateBin(0, 3, cosmosim::core::TimeStepCandidateSource::kUserClamp, "illegal_jump_test");
  assert(throwsWithContext([&]() { scheduler.endSubstep(); }, "illegal_jump_test"));
}

void testSkippedPmSyncTrap() {
  cosmosim::core::PmSynchronizationState pm_sync;
  pm_sync.reset(2);
  const auto event = pm_sync.registerKickOpportunity(4, 0.5, false);
  assert(event.refresh_long_range_field);
  assert(throwsWithContext([&]() { (void)pm_sync.registerKickOpportunity(5, 0.6, true); }, "previous refresh event"));
}

void testStaleMirrorUseTrap() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(2);
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(1);
  scheduler.reset(2, 0, 0);
  scheduler.setElementBin(1, 1, scheduler.currentTick());
  cosmosim::core::syncTimeBinMirrorsFromScheduler(scheduler, state);
  state.particles.time_bin[1] = 0;
  assert(throwsWithContext(
      [&]() {
        cosmosim::core::debugAssertTimeBinMirrorAuthorityInvariant(
            scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticles);
      },
      "time_bin mirror"));
}

void testInvalidRestartTimestepStateTrap() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(1, 1, 0);
  auto state = scheduler.exportPersistentState();
  state.active_flag[0] = 1;
  cosmosim::core::HierarchicalTimeBinScheduler restored(state.max_bin);
  assert(throwsWithContext([&]() { restored.importPersistentState(state); }, "active flag"));
}

void testActiveSetMismatchTrap() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(1);
  scheduler.reset(2, 0, 0);
  const auto active = scheduler.beginSubstep();
  (void)active;
  cosmosim::core::SimulationState state;
  state.resizeParticles(2);
  std::vector<std::uint32_t> mismatched{1, 0};
  cosmosim::core::ActiveSetDescriptor descriptor{
      .particle_indices = mismatched,
      .particles_are_subset = false,
      .particles_from_scheduler = true,
      .has_generation_metadata = true,
      .source_particle_index_generation = state.particleIndexGeneration(),
      .source_cell_index_generation = state.cellIndexGeneration(),
      .source_scheduler_tick = scheduler.currentTick(),
  };
  assert(throwsWithContext(
      [&]() { cosmosim::core::debugAssertActiveSetDescriptorFresh(descriptor, state, scheduler); },
      "do not match scheduler active set"));
  scheduler.endSubstep();
}

void testLocalGlobalSyncBoundaryViolationTrap() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(1, 0, 0);
  scheduler.beginSubstep();
  scheduler.endSubstep();
  scheduler.beginSubstep();
  scheduler.submitCandidateBin(0, 2, cosmosim::core::TimeStepCandidateSource::kUserClamp, "local_global_boundary_test");
  assert(throwsWithContext([&]() { scheduler.endSubstep(); }, "synchronization boundary"));
}


void testSchedulerReorderRemapPreservesActiveParticleIds() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(5);
  initializeParticleIdsAndGasCells(state, {});
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.sfc_key[i] = static_cast<std::uint64_t>(50 - i);
  }

  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(5, 1, 0);
  scheduler.setElementBin(1, 0, scheduler.currentTick());
  scheduler.setElementBin(3, 0, scheduler.currentTick());
  scheduler.setElementBin(4, 2, scheduler.currentTick());
  const std::vector<std::uint64_t> old_ids(state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end());
  const auto before_active_ids = activeParticleIdsAtCurrentTick(scheduler, state);

  const auto reorder = cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kBySfcKey);
  cosmosim::core::reorderParticles(state, reorder);
  const std::vector<std::uint64_t> new_ids(state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end());
  cosmosim::core::remapSchedulerByParticleId(scheduler, old_ids, new_ids);
  cosmosim::core::syncTimeBinMirrorsFromScheduler(scheduler, state);

  const auto after_active_ids = activeParticleIdsAtCurrentTick(scheduler, state);
  assert(after_active_ids == before_active_ids);
}

void testSchedulerCompactionRemapPreservesActiveParticleIds() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(6);
  initializeParticleIdsAndGasCells(state, {});
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(6, 2, 0);
  scheduler.setElementBin(0, 0, scheduler.currentTick());
  scheduler.setElementBin(2, 0, scheduler.currentTick());
  scheduler.setElementBin(5, 1, scheduler.currentTick());
  const std::vector<std::uint64_t> old_ids(state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end());
  const auto before_active_ids = activeParticleIdsAtCurrentTick(scheduler, state);

  cosmosim::core::ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.outbound_local_indices = {1, 4};
  state.commitParticleMigration(commit);
  const std::vector<std::uint64_t> compacted_ids(state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end());
  cosmosim::core::remapSchedulerByParticleId(scheduler, old_ids, compacted_ids);

  std::vector<std::uint64_t> expected;
  for (const auto id : before_active_ids) {
    if (id != 1001 && id != 1004) {
      expected.push_back(id);
    }
  }
  assert(activeParticleIdsAtCurrentTick(scheduler, state) == expected);
}

void testGasCellMirrorSyncUsesParentParticleIdentity() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(4);
  state.resizeCells(2);
  initializeParticleIdsAndGasCells(state, {1, 3});
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(4, 0, 0);
  scheduler.setElementBin(0, 3, scheduler.currentTick());
  scheduler.setElementBin(1, 2, scheduler.currentTick());
  scheduler.setElementBin(2, 1, scheduler.currentTick());
  scheduler.setElementBin(3, 1, scheduler.currentTick());

  cosmosim::core::syncTimeBinMirrorsFromScheduler(
      scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticlesAndCells);
  assert(state.cells.time_bin[0] == 2);
  assert(state.cells.time_bin[1] == 1);
  assert(cosmosim::core::timeBinMirrorsMatchScheduler(
      scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticlesAndCells));
}

void testMigrationRequiresSchedulerIdentityRecords() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(2, 0, 0);
  scheduler.setElementBin(0, 1, scheduler.currentTick());
  std::vector<cosmosim::core::TimeBinSchedulerIdentityRecord> only_kept{
      {.element_id = 1000, .bin_index = 1, .next_activation_tick = 0, .pending_bin_index = cosmosim::core::HierarchicalTimeBinScheduler::k_unset_pending_bin}};
  const std::vector<std::uint64_t> destination_ids{1000, 2000};
  assert(throwsWithContext(
      [&]() { cosmosim::core::rebuildSchedulerFromParticleIdentityRecords(scheduler, only_kept, destination_ids); },
      "time_bin mirror is not scheduler authority"));
}

void testRestartActiveIdEquivalenceWithPendingTransitions() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(4);
  initializeParticleIdsAndGasCells(state, {});
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(4, 1, 0);
  scheduler.setElementBin(0, 0, scheduler.currentTick());
  scheduler.setElementBin(1, 2, scheduler.currentTick());
  scheduler.setElementBin(2, 3, scheduler.currentTick());
  scheduler.setElementBin(3, 1, scheduler.currentTick());
  scheduler.submitCandidateBin(3, 0, cosmosim::core::TimeStepCandidateSource::kUserClamp, "restart_pending");
  (void)scheduler.reconcileCandidateTransitions();

  const auto saved = scheduler.exportPersistentState();
  cosmosim::core::HierarchicalTimeBinScheduler restored(saved.max_bin);
  restored.importPersistentState(saved);
  assert(activeParticleIdsAtCurrentTick(restored, state) == activeParticleIdsAtCurrentTick(scheduler, state));
}


void testOutputBoundaryRequiresSafeContracts() {
  cosmosim::core::SimulationState state;
  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 1.0;

  class UnsafeOutputCallback final : public cosmosim::core::IntegrationCallback {
   public:
    std::string_view callbackName() const override { return "unsafe_output_boundary"; }
    std::span<const cosmosim::core::IntegrationStage> integrationStages() const override {
      static constexpr std::array stages{cosmosim::core::IntegrationStage::kOutputCheck};
      return stages;
    }
    std::span<const cosmosim::core::StageContract> stageContracts() const override { return contracts; }
    void onStage(cosmosim::core::StepContext&) override {}

    std::array<cosmosim::core::StageContract, 1> contracts{{
        {.stage = cosmosim::core::IntegrationStage::kOutputCheck,
         .required_inputs = cosmosim::core::StageDataDomain::kOutputState,
         .mutated_state = cosmosim::core::StageDataDomain::kOutputState,
         .produced_outputs = cosmosim::core::StageDataDomain::kOutputState,
         .allowed_side_effects = cosmosim::core::StageDataDomain::kOutputState,
         .sync_requirements = cosmosim::core::StageSyncRequirement::kGlobal,
         .active_set_family = cosmosim::core::StageActiveSetFamily::kOutputState,
         .restart_safety = cosmosim::core::StageSafety::kUnsafe,
         .output_safety = cosmosim::core::StageSafety::kSafe,
         .owner = cosmosim::core::StageSubsystem::kOutput},
    }};
  } callback;

  cosmosim::core::StepOrchestrator orchestrator;
  assert(throwsWithContext(
      [&]() { orchestrator.registerCallback(callback); },
      "must be restart-safe and output-safe"));
}

void testBoundarySafetyClassification() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(4);
  std::vector<std::uint32_t> active = {0, 2};
  cosmosim::core::ActiveSetDescriptor local_active{
      .particle_indices = active,
      .particles_are_subset = true,
  };
  const auto local_boundary = cosmosim::core::classifyStepBoundary(state, local_active, true);
  assert(local_boundary.kind == cosmosim::core::StepBoundaryKind::kLocalActiveBinStep);
  assert(!local_boundary.restart_safe);
  assert(!local_boundary.output_safe);
  assert(!local_boundary.pm_refresh_allowed);

  local_active.has_global_synchronization_metadata = true;
  local_active.globally_complete_active_set = false;
  const auto coordinated_local_boundary =
      cosmosim::core::classifyStepBoundary(state, local_active, true);
  assert(coordinated_local_boundary.kind ==
         cosmosim::core::StepBoundaryKind::kLocalActiveBinStep);
  assert(!coordinated_local_boundary.restart_safe);
  assert(!coordinated_local_boundary.output_safe);
  assert(coordinated_local_boundary.pm_refresh_allowed);

  cosmosim::core::IntegratorState unsafe_state;
  unsafe_state.current_boundary_kind = local_boundary.kind;
  unsafe_state.last_completed_boundary_kind = local_boundary.kind;
  unsafe_state.last_completed_restart_safe = false;
  const auto local_restart_decision = cosmosim::core::evaluateRestartBoundary(unsafe_state, 7);
  assert(!local_restart_decision.restart_safe);
  assert(local_restart_decision.local_substep_active);
  assert(local_restart_decision.diagnostic.find("local_substep_active=true") != std::string::npos);
  assert(local_restart_decision.diagnostic.find("scheduler_tick=7") != std::string::npos);
  assert(throwsWithContext([&]() { cosmosim::core::assertCanWriteCheckpointAtBoundary(unsafe_state, 7); }, "local active-bin substep restart is not represented"));

  cosmosim::core::IntegratorState half_step_state;
  half_step_state.inside_kdk_step = true;
  half_step_state.current_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
  half_step_state.last_completed_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
  half_step_state.last_completed_restart_safe = true;
  const auto half_step_decision = cosmosim::core::evaluateRestartBoundary(half_step_state);
  assert(!half_step_decision.restart_safe);
  assert(half_step_decision.diagnostic.find("inside_kdk_step=true") != std::string::npos);
  assert(throwsWithContext([&]() { cosmosim::core::assertCanWriteCheckpointAtBoundary(half_step_state); }, "half-step restart is not represented"));

  cosmosim::core::ActiveSetDescriptor global_active{};
  const auto global_boundary = cosmosim::core::classifyStepBoundary(state, global_active, false);
  assert(global_boundary.kind == cosmosim::core::StepBoundaryKind::kGlobalSynchronizationPoint);
  assert(global_boundary.restart_safe);

  const auto checkpoint_boundary = cosmosim::core::classifyStepBoundary(
      state,
      global_active,
      false,
      cosmosim::core::StepBoundaryKind::kCheckpointPoint);
  assert(checkpoint_boundary.kind == cosmosim::core::StepBoundaryKind::kCheckpointPoint);
  assert(checkpoint_boundary.restart_safe);
  assert(checkpoint_boundary.output_safe);

  cosmosim::core::IntegratorState checkpoint_state;
  checkpoint_state.current_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
  checkpoint_state.last_completed_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
  checkpoint_state.last_completed_restart_safe = true;
  const auto checkpoint_decision = cosmosim::core::evaluateRestartBoundary(checkpoint_state, 8);
  assert(checkpoint_decision.restart_safe);
  assert(cosmosim::core::canWriteRestart(checkpoint_state, 8));
  assert(checkpoint_decision.diagnostic.find("restart boundary safe") != std::string::npos);

  const auto pm_boundary = cosmosim::core::classifyStepBoundary(
      state,
      global_active,
      false,
      cosmosim::core::StepBoundaryKind::kPmRefreshPoint);
  assert(pm_boundary.kind == cosmosim::core::StepBoundaryKind::kPmRefreshPoint);
  assert(!pm_boundary.restart_safe);
  assert(!pm_boundary.output_safe);
  assert(pm_boundary.pm_refresh_allowed);
}

void testPmSynchronizationCadencePreservesRefreshBoundaries() {
  cosmosim::core::PmSynchronizationState pm_sync;
  pm_sync.reset(3);

  const auto first = pm_sync.registerKickOpportunity(10, 0.5, false);
  assert(first.gravity_kick_opportunity == 1);
  assert(first.refresh_long_range_field);
  assert(first.field_version == 1);
  pm_sync.commitRefresh(first);
  assert(pm_sync.fieldVersion() == 1);
  assert(pm_sync.lastRefreshOpportunity() == 1);

  const auto reuse1 = pm_sync.registerKickOpportunity(11, 0.6, true);
  assert(reuse1.gravity_kick_opportunity == 2);
  assert(!reuse1.refresh_long_range_field);
  assert(reuse1.field_version == 1);
  const auto reuse2 = pm_sync.registerKickOpportunity(12, 0.7, true);
  assert(reuse2.gravity_kick_opportunity == 3);
  assert(!reuse2.refresh_long_range_field);
  assert(reuse2.last_refresh_opportunity == 1);

  const auto refresh = pm_sync.registerKickOpportunity(13, 0.8, true);
  assert(refresh.gravity_kick_opportunity == 4);
  assert(refresh.refresh_long_range_field);
  assert(refresh.field_version == 2);
  assert(refresh.last_refresh_opportunity == 4);
  assert(refresh.field_built_step_index == 13);
  pm_sync.commitRefresh(refresh);
  assert(pm_sync.fieldVersion() == 2);
  assert(pm_sync.lastRefreshOpportunity() == 4);
  assert(std::abs(pm_sync.lastRefreshScaleFactor() - 0.8) < k_tolerance);
}

}  // namespace

int main() {
  testKickDriftKickOrdering();
  testPmRefreshDirectiveCapturesReasonAndForceEvalTime();
  testLocalForceRefreshDoesNotIssuePmSyncEvent();
  testStageBoundDispatch();
  testRegisterCallbackRejectsExtraUnregisteredStageContract();
  testRegisterCallbackRejectsDuplicateStageContracts();
  testRegisterCallbackRejectsIncompleteExecutableContract();
  testLocalPmSyncEventIsHardRejectedAfterHandlerMutation();
  testCosmologyHelpers();
  testCosmologicalTimelineConvertsSiIntegralsToCodeTime();
  testActiveSubsetDetection();
  testTimeBinMappingAndCriteria();
  testHierarchicalSchedulerTransitions();
  testTimestepBinAuthorityInvariant();
  testTimestepBinReassignmentAndRestartRoundTrip();
  testTimestepBinReorderIdentitySurvival();
  testSchedulerReorderRemapPreservesActiveParticleIds();
  testSchedulerCompactionRemapPreservesActiveParticleIds();
  testGasCellMirrorSyncUsesParentParticleIdentity();
  testMigrationRequiresSchedulerIdentityRecords();
  testRestartActiveIdEquivalenceWithPendingTransitions();
  testActiveSetAuthority();
  testActiveSetNoCompetingBuilders();
  testHydroGravityCandidateReconciliation();
  testHighMachHydroCflCandidateLimitsBin();
  testOrchestratorRequiresSchedulerProvenanceWhenTickIsExpected();
  testOrchestratorRejectsSchedulerActiveSetWithoutTick();
  testSchedulerBackedTimeBinReorderIgnoresStaleMirrors();
  testPmSynchronizationPersistentRoundTrip();
  testInvalidEarlyActivationTrap();
  testInvalidBinJumpTrap();
  testSkippedPmSyncTrap();
  testStaleMirrorUseTrap();
  testInvalidRestartTimestepStateTrap();
  testActiveSetMismatchTrap();
  testLocalGlobalSyncBoundaryViolationTrap();
  testOutputBoundaryRequiresSafeContracts();
  testBoundarySafetyClassification();
  testPmSynchronizationCadencePreservesRefreshBoundaries();
  return 0;
}
