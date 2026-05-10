#include "cosmosim/core/time_integration.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace cosmosim::core {
namespace {

constexpr std::array<IntegrationStage, 7> k_kick_drift_kick_order = {
    IntegrationStage::kGravityKickPre,
    IntegrationStage::kDrift,
    IntegrationStage::kHydroUpdate,
    IntegrationStage::kSourceTerms,
    IntegrationStage::kGravityKickPost,
    IntegrationStage::kAnalysisHooks,
    IntegrationStage::kOutputCheck,
};

[[nodiscard]] double midpointIntegrateDriftLike(
    const LambdaCdmBackground& background,
    double scale_factor_begin,
    double scale_factor_end,
    std::uint32_t midpoint_samples) {
  if (scale_factor_begin <= 0.0 || scale_factor_end <= 0.0) {
    throw std::invalid_argument("scale factors must be positive for drift-like integrals");
  }
  if (scale_factor_end < scale_factor_begin) {
    throw std::invalid_argument("scale_factor_end must be >= scale_factor_begin");
  }

  const std::uint32_t samples = std::max<std::uint32_t>(midpoint_samples, 1U);
  const double delta_a = (scale_factor_end - scale_factor_begin) / static_cast<double>(samples);
  if (delta_a == 0.0) {
    return 0.0;
  }

  double accum = 0.0;
  for (std::uint32_t i = 0; i < samples; ++i) {
    const double a_mid = scale_factor_begin + (static_cast<double>(i) + 0.5) * delta_a;
    accum += 1.0 / (a_mid * a_mid * background.hubbleSi(a_mid));
  }
  return accum * delta_a;
}

[[nodiscard]] std::uint64_t powerOfTwo(std::uint8_t exponent) {
  if (exponent >= 63) {
    throw std::invalid_argument("bin exponent too large for 64-bit tick period");
  }
  return 1ULL << exponent;
}

[[nodiscard]] double finitePositiveOrInf(double value) {
  if (!std::isfinite(value) || value <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  return value;
}

[[nodiscard]] std::string timeBinContextMessage(
    std::string_view invariant,
    std::string_view source_label,
    std::uint64_t current_tick,
    std::uint32_t element_index,
    std::uint8_t bin,
    std::uint64_t next_activation_tick,
    std::uint64_t integrator_step = 0,
    double integrator_time = 0.0,
    std::string_view pm_sync_state = "n/a") {
  std::ostringstream out;
  out << "time-bin invariant failure: " << invariant
      << "; source=" << source_label
      << "; element=" << element_index
      << "; current_tick=" << current_tick
      << "; bin=" << static_cast<unsigned>(bin)
      << "; next_activation_tick=" << next_activation_tick
      << "; integrator_step=" << integrator_step
      << "; integrator_time=" << integrator_time
      << "; pm_sync=" << pm_sync_state;
  return out.str();
}

[[nodiscard]] std::string schedulerContextMessage(
    std::string_view invariant,
    std::string_view source_label,
    std::uint64_t current_tick,
    std::uint64_t integrator_step = 0,
    double integrator_time = 0.0,
    std::string_view pm_sync_state = "n/a") {
  std::ostringstream out;
  out << "scheduler invariant failure: " << invariant
      << "; source=" << source_label
      << "; current_tick=" << current_tick
      << "; integrator_step=" << integrator_step
      << "; integrator_time=" << integrator_time
      << "; pm_sync=" << pm_sync_state;
  return out.str();
}

[[nodiscard]] std::string pmContextMessage(
    std::string_view invariant,
    std::string_view source_label,
    const PmSyncEvent& event,
    std::uint64_t cadence_steps,
    std::uint64_t current_field_version) {
  std::ostringstream out;
  out << "PM synchronization invariant failure: " << invariant
      << "; source=" << source_label
      << "; opportunity=" << event.gravity_kick_opportunity
      << "; cadence_steps=" << cadence_steps
      << "; refresh=" << (event.refresh_long_range_field ? "true" : "false")
      << "; event_field_version=" << event.field_version
      << "; current_field_version=" << current_field_version
      << "; last_refresh_opportunity=" << event.last_refresh_opportunity
      << "; field_built_step=" << event.field_built_step_index
      << "; field_built_scale_factor=" << event.field_built_scale_factor;
  return out.str();
}

}  // namespace

std::string_view timeStepCandidateSourceName(TimeStepCandidateSource source) {
  switch (source) {
    case TimeStepCandidateSource::kHydroCfl:
      return "hydro_cfl";
    case TimeStepCandidateSource::kGravityAcceleration:
      return "gravity_acceleration";
    case TimeStepCandidateSource::kSourceTerm:
      return "source_term";
    case TimeStepCandidateSource::kUserClamp:
      return "user_clamp";
  }
  return "unknown";
}

void PmSynchronizationState::reset(std::uint64_t cadence_steps) {
  if (cadence_steps == 0) {
    throw std::invalid_argument("PM synchronization cadence must be >= 1");
  }
  m_cadence_steps = cadence_steps;
  m_gravity_kick_opportunity = 0;
  m_last_refresh_opportunity = 0;
  m_field_version = 0;
  m_last_refresh_step_index = 0;
  m_last_refresh_scale_factor = 1.0;
  m_refresh_commit_pending = false;
  m_pending_refresh_opportunity = 0;
  m_pending_refresh_field_version = 0;
}

PmSyncEvent PmSynchronizationState::registerKickOpportunity(
    std::uint64_t step_index,
    double scale_factor,
    bool has_long_range_field) {
  if (m_refresh_commit_pending) {
    PmSyncEvent pending{
        .gravity_kick_opportunity = m_pending_refresh_opportunity,
        .refresh_long_range_field = true,
        .field_version = m_pending_refresh_field_version,
        .last_refresh_opportunity = m_pending_refresh_opportunity,
        .field_built_step_index = m_last_refresh_step_index,
        .field_built_scale_factor = m_last_refresh_scale_factor,
    };
    throw std::runtime_error(pmContextMessage(
        "previous refresh event was not committed before next kick opportunity",
        "PmSynchronizationState::registerKickOpportunity",
        pending,
        m_cadence_steps,
        m_field_version));
  }
  ++m_gravity_kick_opportunity;
  const bool refresh = !has_long_range_field ||
      ((m_gravity_kick_opportunity - m_last_refresh_opportunity) >= m_cadence_steps);
  PmSyncEvent event{
      .gravity_kick_opportunity = m_gravity_kick_opportunity,
      .refresh_long_range_field = refresh,
      .field_version = refresh ? (m_field_version + 1U) : m_field_version,
      .last_refresh_opportunity = refresh ? m_gravity_kick_opportunity : m_last_refresh_opportunity,
      .field_built_step_index = refresh ? step_index : m_last_refresh_step_index,
      .field_built_scale_factor = refresh ? scale_factor : m_last_refresh_scale_factor,
  };
  if (event.refresh_long_range_field) {
    m_refresh_commit_pending = true;
    m_pending_refresh_opportunity = event.gravity_kick_opportunity;
    m_pending_refresh_field_version = event.field_version;
  }
  return event;
}

void PmSynchronizationState::commitRefresh(const PmSyncEvent& event) {
  if (!event.refresh_long_range_field) {
    return;
  }
  if (!m_refresh_commit_pending || event.gravity_kick_opportunity != m_pending_refresh_opportunity) {
    throw std::runtime_error(pmContextMessage(
        "refresh event does not match pending kick opportunity",
        "PmSynchronizationState::commitRefresh",
        event,
        m_cadence_steps,
        m_field_version));
  }
  if (event.gravity_kick_opportunity != m_gravity_kick_opportunity) {
    throw std::runtime_error(pmContextMessage(
        "event does not match current kick opportunity",
        "PmSynchronizationState::commitRefresh",
        event,
        m_cadence_steps,
        m_field_version));
  }
  if (event.field_version != m_field_version + 1U || event.field_version != m_pending_refresh_field_version) {
    throw std::runtime_error(pmContextMessage(
        "field version transition is illegal",
        "PmSynchronizationState::commitRefresh",
        event,
        m_cadence_steps,
        m_field_version));
  }
  m_field_version = event.field_version;
  m_last_refresh_opportunity = event.last_refresh_opportunity;
  m_last_refresh_step_index = event.field_built_step_index;
  m_last_refresh_scale_factor = event.field_built_scale_factor;
  m_refresh_commit_pending = false;
  m_pending_refresh_opportunity = 0;
  m_pending_refresh_field_version = 0;
}

std::string_view integrationStageName(IntegrationStage stage) {
  switch (stage) {
    case IntegrationStage::kGravityKickPre:
      return "gravity_kick_pre";
    case IntegrationStage::kDrift:
      return "drift";
    case IntegrationStage::kHydroUpdate:
      return "hydro_update";
    case IntegrationStage::kSourceTerms:
      return "source_terms";
    case IntegrationStage::kGravityKickPost:
      return "gravity_kick_post";
    case IntegrationStage::kAnalysisHooks:
      return "analysis_hooks";
    case IntegrationStage::kOutputCheck:
      return "output_check";
  }
  return "unknown";
}

bool ActiveSetDescriptor::hasParticleSubset(std::size_t total_particle_count) const noexcept {
  return particles_are_subset && particle_indices.size() < total_particle_count;
}

bool ActiveSetDescriptor::hasCellSubset(std::size_t total_cell_count) const noexcept {
  return cells_are_subset && cell_indices.size() < total_cell_count;
}

ActiveSetDescriptor makeSchedulerActiveSetDescriptor(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint32_t> active_cell_indices) {
  ActiveSetDescriptor descriptor{
      .particle_indices = active_particle_indices,
      .cell_indices = active_cell_indices,
      .particles_are_subset = active_particle_indices.size() != state.particles.size(),
      .cells_are_subset = active_cell_indices.size() != state.cells.size(),
      .particles_from_scheduler = true,
      .cells_from_scheduler = !active_cell_indices.empty(),
      .has_generation_metadata = true,
      .source_particle_index_generation = state.particleIndexGeneration(),
      .source_cell_index_generation = state.cellIndexGeneration(),
      .source_scheduler_tick = scheduler.currentTick(),
  };
  debugAssertActiveSetDescriptorFresh(descriptor, state, scheduler);
  return descriptor;
}

void debugAssertActiveSetDescriptorFresh(
    const ActiveSetDescriptor& active_set,
    const SimulationState& state) {
  if (active_set.has_generation_metadata &&
      active_set.source_particle_index_generation != state.particleIndexGeneration()) {
    throw std::runtime_error(schedulerContextMessage(
        "ActiveSetDescriptor particle generation is stale",
        "debugAssertActiveSetDescriptorFresh",
        active_set.source_scheduler_tick));
  }
  if (active_set.has_generation_metadata &&
      active_set.source_cell_index_generation != state.cellIndexGeneration()) {
    throw std::runtime_error(schedulerContextMessage(
        "ActiveSetDescriptor cell generation is stale",
        "debugAssertActiveSetDescriptorFresh",
        active_set.source_scheduler_tick));
  }
  for (const std::uint32_t particle_index : active_set.particle_indices) {
    if (particle_index >= state.particles.size()) {
      throw std::out_of_range(timeBinContextMessage(
          "ActiveSetDescriptor contains stale particle index",
          "debugAssertActiveSetDescriptorFresh",
          active_set.source_scheduler_tick,
          particle_index,
          0,
          active_set.source_scheduler_tick));
    }
  }
  for (const std::uint32_t cell_index : active_set.cell_indices) {
    if (cell_index >= state.cells.size()) {
      throw std::out_of_range(timeBinContextMessage(
          "ActiveSetDescriptor contains stale cell index",
          "debugAssertActiveSetDescriptorFresh",
          active_set.source_scheduler_tick,
          cell_index,
          0,
          active_set.source_scheduler_tick));
    }
  }
}

void debugAssertActiveSetDescriptorFresh(
    const ActiveSetDescriptor& active_set,
    const SimulationState& state,
    std::uint64_t expected_scheduler_tick) {
  debugAssertActiveSetDescriptorFresh(active_set, state);
  if (!active_set.has_generation_metadata || !active_set.particles_from_scheduler) {
    throw std::runtime_error(schedulerContextMessage(
        "ActiveSetDescriptor must carry scheduler provenance before solver callbacks",
        "debugAssertActiveSetDescriptorFresh",
        expected_scheduler_tick));
  }
  if ((active_set.particles_from_scheduler || active_set.cells_from_scheduler) &&
      active_set.source_scheduler_tick != expected_scheduler_tick) {
    throw std::runtime_error(schedulerContextMessage(
        "ActiveSetDescriptor scheduler tick is stale",
        "debugAssertActiveSetDescriptorFresh",
        expected_scheduler_tick));
  }
}

void debugAssertActiveSetDescriptorFresh(
    const ActiveSetDescriptor& active_set,
    const SimulationState& state,
    const HierarchicalTimeBinScheduler& scheduler) {
  debugAssertActiveSetDescriptorFresh(active_set, state, scheduler.currentTick());
  if (active_set.particles_from_scheduler) {
    const auto active = scheduler.activeElements();
    if (active.size() != active_set.particle_indices.size() ||
        !std::equal(active.begin(), active.end(), active_set.particle_indices.begin())) {
      throw std::runtime_error(schedulerContextMessage(
          "ActiveSetDescriptor particle indices do not match scheduler active set",
          "debugAssertActiveSetDescriptorFresh",
          scheduler.currentTick()));
    }
  }
}

std::vector<IntegrationStage> StageScheduler::schedule(
    const IntegratorState& integrator_state,
    const ActiveSetDescriptor& /*active_set*/) const {
  if (integrator_state.scheme != TimeStepScheme::kKickDriftKick) {
    throw std::invalid_argument("unsupported timestep scheme");
  }

  return std::vector<IntegrationStage>(k_kick_drift_kick_order.begin(), k_kick_drift_kick_order.end());
}

std::span<const IntegrationStage> StageScheduler::kickDriftKickOrder() { return k_kick_drift_kick_order; }

bool isCanonicalIntegrationStageOrder(std::span<const IntegrationStage> ordered_stages) {
  return ordered_stages.size() == k_kick_drift_kick_order.size() &&
         std::equal(
             ordered_stages.begin(),
             ordered_stages.end(),
             k_kick_drift_kick_order.begin(),
             k_kick_drift_kick_order.end());
}

StepOrchestrator::StepOrchestrator(StageScheduler scheduler) : m_scheduler(std::move(scheduler)) {}

void StepOrchestrator::registerCallback(IntegrationCallback& callback) { m_callbacks.push_back(&callback); }

std::size_t StepOrchestrator::callbackCount() const noexcept { return m_callbacks.size(); }

void StepOrchestrator::executeSingleStep(
    SimulationState& state,
    IntegratorState& integrator_state,
    ActiveSetDescriptor active_set,
    const LambdaCdmBackground* cosmology_background,
    TransientStepWorkspace* workspace,
    const ModePolicy* mode_policy,
    ProfilerSession* profiler_session,
    std::optional<std::uint64_t> expected_scheduler_tick) const {
  if (integrator_state.dt_time_code <= 0.0) {
    throw std::invalid_argument("dt_time_code must be positive");
  }
  if (expected_scheduler_tick.has_value()) {
    debugAssertActiveSetDescriptorFresh(active_set, state, *expected_scheduler_tick);
  } else {
    debugAssertActiveSetDescriptorFresh(active_set, state);
  }

  COSMOSIM_PROFILE_SCOPE(profiler_session, "step_orchestrator.execute_single_step");

  if (profiler_session != nullptr) {
    profiler_session->counters().setCount(
        "active_particles",
        static_cast<std::uint64_t>(active_set.particle_indices.size()));
    profiler_session->counters().setCount(
        "active_cells",
        static_cast<std::uint64_t>(active_set.cell_indices.size()));
    profiler_session->counters().addCount("step_invocations", 1);
  }

  StepContext context{
      .state = state,
      .integrator_state = integrator_state,
      .active_set = active_set,
      .workspace = workspace,
      .cosmology_background = cosmology_background,
      .mode_policy = mode_policy,
      .profiler_session = profiler_session,
      .stage = IntegrationStage::kGravityKickPre,
  };

  const auto ordered_stages = m_scheduler.schedule(integrator_state, active_set);
  if (!isCanonicalIntegrationStageOrder(ordered_stages)) {
    throw std::runtime_error("stage scheduler order deviates from canonical kick-drift-kick contract");
  }

  for (const auto stage : ordered_stages) {
    context.stage = stage;
    const std::string stage_name = "stage." + std::string(integrationStageName(stage));
    COSMOSIM_PROFILE_SCOPE(profiler_session, stage_name);
    if (profiler_session != nullptr) {
      profiler_session->counters().addCount(stage_name + ".invocations", 1);
    }

    for (auto* callback : m_callbacks) {
      const std::string callback_phase = "callback." + std::string(callback->callbackName());
      COSMOSIM_PROFILE_SCOPE(profiler_session, callback_phase);
      callback->onStage(context);
      if (profiler_session != nullptr) {
        profiler_session->counters().addCount(callback_phase + ".invocations", 1);
      }
    }
  }

  integrator_state.current_time_code += integrator_state.dt_time_code;
  if (cosmology_background != nullptr) {
    integrator_state.current_scale_factor = advanceScaleFactorEuler(
        *cosmology_background,
        integrator_state.current_scale_factor,
        integrator_state.dt_time_code);
  }
  ++integrator_state.step_index;
}

void TimeStepCriteriaRegistry::registerCflHook(CriteriaHook hook) { m_hooks.cfl_hook = std::move(hook); }

void TimeStepCriteriaRegistry::registerGravityHook(CriteriaHook hook) { m_hooks.gravity_hook = std::move(hook); }

void TimeStepCriteriaRegistry::registerSourceHook(CriteriaHook hook) { m_hooks.source_hook = std::move(hook); }

void TimeStepCriteriaRegistry::registerUserClampHook(CriteriaHook hook) {
  m_hooks.user_clamp_hook = std::move(hook);
}

const TimeStepCriteriaHooks& TimeStepCriteriaRegistry::hooks() const noexcept { return m_hooks; }

HierarchicalTimeBinScheduler::HierarchicalTimeBinScheduler(std::uint8_t max_bin) : m_max_bin(max_bin) {}

void HierarchicalTimeBinScheduler::reset(
    std::uint32_t element_count,
    std::uint8_t initial_bin,
    std::uint64_t start_tick) {
  m_current_tick = start_tick;

  const std::uint8_t clamped_bin = clampBin(initial_bin);
  m_hot.bin_index.assign(element_count, clamped_bin);
  m_hot.next_activation_tick.assign(element_count, m_current_tick);
  m_hot.active_flag.assign(element_count, 0);
  m_hot.pending_bin_index.assign(element_count, k_unset_pending_bin);
  m_candidate_bin_index.assign(element_count, k_unset_pending_bin);
  m_candidate_label.assign(element_count, {});
  m_last_reconciliation = {};
  m_substep_open = false;

  m_position_in_bin.resize(element_count, 0);
  m_elements_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, {});
  for (auto& bin_members : m_elements_by_bin) {
    bin_members.reserve(element_count / (m_max_bin + 1U) + 1U);
  }

  for (std::uint32_t element = 0; element < element_count; ++element) {
    m_position_in_bin[element] = m_elements_by_bin[clamped_bin].size();
    m_elements_by_bin[clamped_bin].push_back(element);
  }

  m_active_elements.clear();
  m_diagnostics = {};
  m_diagnostics.occupancy_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, 0U);
  m_diagnostics.active_count_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, 0U);
  m_diagnostics.occupancy_by_bin[clamped_bin] = element_count;
}

void HierarchicalTimeBinScheduler::setElementBin(
    std::uint32_t element_index,
    std::uint8_t bin_index,
    std::uint64_t current_tick) {
  if (element_index >= m_hot.size()) {
    throw std::out_of_range("element_index out of range");
  }

  const std::uint8_t clamped_bin = clampBin(bin_index);
  const std::uint8_t old_bin = m_hot.bin_index[element_index];
  if (old_bin == clamped_bin) {
    const std::uint64_t period_ticks = binPeriodTicks(clamped_bin);
    m_hot.next_activation_tick[element_index] =
        (current_tick % period_ticks == 0) ? current_tick : ((current_tick / period_ticks) * period_ticks + period_ticks);
    return;
  }

  eraseFromBin(element_index, old_bin);
  insertIntoBin(element_index, clamped_bin);

  m_hot.bin_index[element_index] = clamped_bin;
  const std::uint64_t period_ticks = binPeriodTicks(clamped_bin);
  m_hot.next_activation_tick[element_index] =
      (current_tick % period_ticks == 0) ? current_tick : ((current_tick / period_ticks) * period_ticks + period_ticks);
}

void HierarchicalTimeBinScheduler::submitCandidateTimeStep(
    std::uint32_t element_index,
    double dt_time_code,
    const TimeStepLimits& limits,
    TimeStepCandidateSource source,
    std::string_view label) {
  const auto mapped = mapDtToTimeBin(dt_time_code, limits);
  submitCandidateBin(element_index, mapped.bin_index, source, label);
  if (mapped.clipped_to_min) {
    ++m_last_reconciliation.clipped_to_min_dt;
  }
  if (mapped.clipped_to_max) {
    ++m_last_reconciliation.clipped_to_max_dt;
  }
}

void HierarchicalTimeBinScheduler::submitCandidateBin(
    std::uint32_t element_index,
    std::uint8_t target_bin,
    TimeStepCandidateSource source,
    std::string_view label) {
  if (element_index >= m_hot.size()) {
    throw std::out_of_range("element_index out of range");
  }
  if (m_candidate_bin_index.size() != m_hot.size()) {
    m_candidate_bin_index.assign(m_hot.size(), k_unset_pending_bin);
    m_candidate_label.assign(m_hot.size(), {});
  }
  const std::uint8_t clamped = clampBin(target_bin);
  ++m_last_reconciliation.submitted_candidates;
  if (m_candidate_bin_index[element_index] == k_unset_pending_bin || clamped < m_candidate_bin_index[element_index]) {
    m_candidate_bin_index[element_index] = clamped;
    m_candidate_label[element_index] = label.empty()
        ? std::string(timeStepCandidateSourceName(source))
        : std::string(label);
  }
}

TimeStepReconciliationResult HierarchicalTimeBinScheduler::reconcileCandidateTransitions() {
  if (m_candidate_bin_index.size() != m_hot.size()) {
    m_candidate_bin_index.assign(m_hot.size(), k_unset_pending_bin);
    m_candidate_label.assign(m_hot.size(), {});
  }
  TimeStepReconciliationResult result = m_last_reconciliation;
  for (std::uint32_t element = 0; element < m_candidate_bin_index.size(); ++element) {
    const std::uint8_t candidate = m_candidate_bin_index[element];
    if (candidate == k_unset_pending_bin) {
      continue;
    }
    ++result.elements_with_candidates;
    validateTransitionRequest(
        element,
        candidate,
        m_candidate_label[element].empty() ? "HierarchicalTimeBinScheduler::reconcileCandidateTransitions" : m_candidate_label[element]);
    requestBinTransition(element, candidate);
    ++result.committed_transition_requests;
    m_candidate_bin_index[element] = k_unset_pending_bin;
    m_candidate_label[element].clear();
  }
  m_last_reconciliation = {};
  return result;
}

void HierarchicalTimeBinScheduler::requestBinTransition(
    std::uint32_t element_index,
    std::uint8_t target_bin) {
  if (element_index >= m_hot.size()) {
    throw std::out_of_range("element_index out of range");
  }
  validateTransitionRequest(element_index, target_bin, "HierarchicalTimeBinScheduler::requestBinTransition");
  m_hot.pending_bin_index[element_index] = clampBin(target_bin);
}

std::span<const std::uint32_t> HierarchicalTimeBinScheduler::activeElements() const noexcept {
  return m_active_elements;
}

std::uint64_t HierarchicalTimeBinScheduler::currentTick() const noexcept { return m_current_tick; }

std::uint8_t HierarchicalTimeBinScheduler::maxBin() const noexcept { return m_max_bin; }

std::uint32_t HierarchicalTimeBinScheduler::elementCount() const noexcept {
  return static_cast<std::uint32_t>(m_hot.size());
}

bool HierarchicalTimeBinScheduler::isBinActiveAtTick(std::uint8_t bin_index, std::uint64_t tick) const {
  const std::uint8_t clamped_bin = clampBin(bin_index);
  return tick % binPeriodTicks(clamped_bin) == 0;
}

std::uint64_t HierarchicalTimeBinScheduler::binPeriodTicks(std::uint8_t bin_index) const {
  return powerOfTwo(clampBin(bin_index));
}

std::span<const std::uint32_t> HierarchicalTimeBinScheduler::beginSubstep() {
  if (m_substep_open) {
    throw std::runtime_error(schedulerContextMessage(
        "beginSubstep called while a substep is already open",
        "HierarchicalTimeBinScheduler::beginSubstep",
        m_current_tick));
  }
  validateInternalState("HierarchicalTimeBinScheduler::beginSubstep");
  rebuildActiveSet();
  m_substep_open = true;
  validateInternalState("HierarchicalTimeBinScheduler::beginSubstep.active_set_created");
  return m_active_elements;
}

void HierarchicalTimeBinScheduler::endSubstep() {
  if (!m_substep_open) {
    throw std::runtime_error(schedulerContextMessage(
        "endSubstep called without an open substep",
        "HierarchicalTimeBinScheduler::endSubstep",
        m_current_tick));
  }
  (void)reconcileCandidateTransitions();
  applyPendingTransitions();
  for (const std::uint32_t element : m_active_elements) {
    const std::uint8_t bin = m_hot.bin_index[element];
    const std::uint64_t period_ticks = binPeriodTicks(bin);
    m_hot.next_activation_tick[element] = m_current_tick + period_ticks;
    m_hot.active_flag[element] = 0;
  }
  if (m_current_tick == std::numeric_limits<std::uint64_t>::max()) {
    throw std::runtime_error(schedulerContextMessage(
        "global integer time cannot advance monotonically without overflowing",
        "HierarchicalTimeBinScheduler::endSubstep",
        m_current_tick));
  }
  ++m_current_tick;
  m_substep_open = false;
  validateInternalState("HierarchicalTimeBinScheduler::endSubstep.committed");
}

const TimeBinHotMetadata& HierarchicalTimeBinScheduler::hotMetadata() const noexcept { return m_hot; }

const TimeBinDiagnostics& HierarchicalTimeBinScheduler::diagnostics() const noexcept { return m_diagnostics; }

TimeBinPersistentState HierarchicalTimeBinScheduler::exportPersistentState() const {
  validateInternalState("HierarchicalTimeBinScheduler::exportPersistentState");
  TimeBinPersistentState persistent_state;
  persistent_state.current_tick = m_current_tick;
  persistent_state.max_bin = m_max_bin;
  persistent_state.bin_index = m_hot.bin_index;
  persistent_state.next_activation_tick = m_hot.next_activation_tick;
  persistent_state.active_flag = m_hot.active_flag;
  persistent_state.pending_bin_index = m_hot.pending_bin_index;
  return persistent_state;
}

void HierarchicalTimeBinScheduler::importPersistentState(const TimeBinPersistentState& persistent_state) {
  if (persistent_state.bin_index.size() != persistent_state.next_activation_tick.size() ||
      persistent_state.bin_index.size() != persistent_state.active_flag.size() ||
      persistent_state.bin_index.size() != persistent_state.pending_bin_index.size()) {
    throw std::invalid_argument("TimeBinPersistentState arrays must have matching sizes");
  }

  for (std::size_t element_index = 0; element_index < persistent_state.bin_index.size(); ++element_index) {
    if (persistent_state.bin_index[element_index] > persistent_state.max_bin) {
      throw std::invalid_argument("TimeBinPersistentState bin_index exceeds max_bin");
    }
    if (persistent_state.active_flag[element_index] > 1U) {
      throw std::invalid_argument("TimeBinPersistentState active_flag must be 0 or 1");
    }
    const bool pending_is_unset = persistent_state.pending_bin_index[element_index] == k_unset_pending_bin;
    if (!pending_is_unset && persistent_state.pending_bin_index[element_index] > persistent_state.max_bin) {
      throw std::invalid_argument("TimeBinPersistentState pending_bin_index exceeds max_bin");
    }
  }

  m_current_tick = persistent_state.current_tick;
  m_max_bin = persistent_state.max_bin;
  m_hot.bin_index = persistent_state.bin_index;
  m_hot.next_activation_tick = persistent_state.next_activation_tick;
  m_hot.active_flag = persistent_state.active_flag;
  m_hot.pending_bin_index = persistent_state.pending_bin_index;
  m_candidate_bin_index.assign(m_hot.bin_index.size(), k_unset_pending_bin);
  m_candidate_label.assign(m_hot.bin_index.size(), {});
  m_last_reconciliation = {};

  m_elements_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, {});
  m_position_in_bin.assign(m_hot.bin_index.size(), 0);

  for (std::size_t element_index = 0; element_index < m_hot.bin_index.size(); ++element_index) {
    const std::uint8_t bin = clampBin(m_hot.bin_index[element_index]);
    m_hot.bin_index[element_index] = bin;
    auto& members = m_elements_by_bin[bin];
    m_position_in_bin[element_index] = members.size();
    members.push_back(static_cast<std::uint32_t>(element_index));
  }

  m_active_elements.clear();
  m_substep_open = false;
  m_diagnostics = {};
  m_diagnostics.occupancy_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, 0U);
  m_diagnostics.active_count_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, 0U);
  for (std::size_t bin = 0; bin < m_elements_by_bin.size(); ++bin) {
    m_diagnostics.occupancy_by_bin[bin] = static_cast<std::uint32_t>(m_elements_by_bin[bin].size());
  }
  validateInternalState("HierarchicalTimeBinScheduler::importPersistentState");
}

void HierarchicalTimeBinScheduler::validateInternalState(std::string_view source_label) const {
  if (m_hot.bin_index.size() != m_hot.next_activation_tick.size() ||
      m_hot.bin_index.size() != m_hot.active_flag.size() ||
      m_hot.bin_index.size() != m_hot.pending_bin_index.size()) {
    throw std::runtime_error(schedulerContextMessage(
        "hot metadata arrays have mismatched sizes",
        source_label,
        m_current_tick));
  }
  if (m_position_in_bin.size() != m_hot.size()) {
    throw std::runtime_error(schedulerContextMessage(
        "bin-position mirror size does not match hot metadata",
        source_label,
        m_current_tick));
  }
  for (std::uint32_t element = 0; element < m_hot.size(); ++element) {
    const std::uint8_t bin = m_hot.bin_index[element];
    const std::uint64_t next_tick = m_hot.next_activation_tick[element];
    if (bin > m_max_bin) {
      throw std::runtime_error(timeBinContextMessage(
          "bin exceeds scheduler max_bin",
          source_label,
          m_current_tick,
          element,
          bin,
          next_tick));
    }
    if (m_hot.active_flag[element] > 1U) {
      throw std::runtime_error(timeBinContextMessage(
          "active flag is not boolean",
          source_label,
          m_current_tick,
          element,
          bin,
          next_tick));
    }
    if (!m_substep_open && m_hot.active_flag[element] != 0U) {
      throw std::runtime_error(timeBinContextMessage(
          "stale active flag outside an open substep",
          source_label,
          m_current_tick,
          element,
          bin,
          next_tick));
    }
    const std::uint64_t period = binPeriodTicks(bin);
    if (next_tick < m_current_tick) {
      throw std::runtime_error(timeBinContextMessage(
          "next activation tick is behind the global scheduler tick",
          source_label,
          m_current_tick,
          element,
          bin,
          next_tick));
    }
    if (next_tick > m_current_tick && next_tick % period != 0U) {
      throw std::runtime_error(timeBinContextMessage(
          "future next activation tick is not aligned to bin period",
          source_label,
          m_current_tick,
          element,
          bin,
          next_tick));
    }
    const std::uint8_t pending = m_hot.pending_bin_index[element];
    if (pending != k_unset_pending_bin && pending > m_max_bin) {
      throw std::runtime_error(timeBinContextMessage(
          "pending bin exceeds scheduler max_bin",
          source_label,
          m_current_tick,
          element,
          bin,
          next_tick));
    }
  }
}

void HierarchicalTimeBinScheduler::validateTransitionRequest(
    std::uint32_t element_index,
    std::uint8_t target_bin,
    std::string_view source_label) const {
  if (element_index >= m_hot.size()) {
    throw std::out_of_range("element_index out of range");
  }
  const std::uint8_t old_bin = m_hot.bin_index[element_index];
  const std::uint8_t new_bin = clampBin(target_bin);
  const std::uint64_t next_tick = m_hot.next_activation_tick[element_index];
  const bool element_is_currently_active =
      next_tick == m_current_tick && isBinActiveAtTick(old_bin, m_current_tick);
  const std::uint64_t new_period = binPeriodTicks(new_bin);
  if (element_is_currently_active && m_current_tick % new_period != 0U) {
    throw std::runtime_error(timeBinContextMessage(
        "bin transition violates current synchronization boundary",
        source_label,
        m_current_tick,
        element_index,
        old_bin,
        next_tick));
  }
}

std::uint8_t HierarchicalTimeBinScheduler::clampBin(std::uint8_t requested) const noexcept {
  return std::min(requested, m_max_bin);
}

void HierarchicalTimeBinScheduler::eraseFromBin(std::uint32_t element_index, std::uint8_t bin_index) {
  auto& members = m_elements_by_bin[bin_index];
  const std::size_t remove_pos = m_position_in_bin[element_index];

  if (remove_pos >= members.size()) {
    throw std::runtime_error("bin position metadata corrupted");
  }

  const std::uint32_t last_element = members.back();
  members[remove_pos] = last_element;
  m_position_in_bin[last_element] = remove_pos;
  members.pop_back();
}

void HierarchicalTimeBinScheduler::insertIntoBin(std::uint32_t element_index, std::uint8_t bin_index) {
  auto& members = m_elements_by_bin[bin_index];
  m_position_in_bin[element_index] = members.size();
  members.push_back(element_index);
}

void HierarchicalTimeBinScheduler::rebuildActiveSet() {
  m_active_elements.clear();

  if (m_diagnostics.occupancy_by_bin.size() != m_elements_by_bin.size()) {
    m_diagnostics.occupancy_by_bin.assign(m_elements_by_bin.size(), 0U);
  }
  if (m_diagnostics.active_count_by_bin.size() != m_elements_by_bin.size()) {
    m_diagnostics.active_count_by_bin.assign(m_elements_by_bin.size(), 0U);
  }

  std::fill(m_diagnostics.active_count_by_bin.begin(), m_diagnostics.active_count_by_bin.end(), 0U);
  std::fill(m_diagnostics.occupancy_by_bin.begin(), m_diagnostics.occupancy_by_bin.end(), 0U);

  for (std::size_t bin = 0; bin < m_elements_by_bin.size(); ++bin) {
    const auto& members = m_elements_by_bin[bin];
    m_diagnostics.occupancy_by_bin[bin] = static_cast<std::uint32_t>(members.size());

    if (!isBinActiveAtTick(static_cast<std::uint8_t>(bin), m_current_tick)) {
      continue;
    }

    for (const std::uint32_t element : members) {
      if (m_hot.next_activation_tick[element] != m_current_tick) {
        continue;
      }
      m_active_elements.push_back(element);
      m_hot.active_flag[element] = 1U;
      ++m_diagnostics.active_count_by_bin[bin];
    }
  }

  std::sort(m_active_elements.begin(), m_active_elements.end());

  m_diagnostics.active_elements = static_cast<std::uint32_t>(m_active_elements.size());
  const auto total_elements = static_cast<double>(m_hot.size());
  m_diagnostics.active_fraction = total_elements > 0.0
      ? static_cast<double>(m_diagnostics.active_elements) / total_elements
      : 0.0;

  m_diagnostics.most_active_bin = 0;
  std::uint32_t best = 0;
  for (std::uint8_t bin = 0; bin <= m_max_bin; ++bin) {
    if (m_diagnostics.active_count_by_bin[bin] > best) {
      best = m_diagnostics.active_count_by_bin[bin];
      m_diagnostics.most_active_bin = bin;
    }
  }

  m_diagnostics.collapse_candidates = 0;
  if (m_hot.size() > 0) {
    const auto finest_occupancy = m_diagnostics.occupancy_by_bin[0];
    if (finest_occupancy * 4U >= m_hot.size() * 3U) {
      m_diagnostics.collapse_candidates = finest_occupancy;
    }
  }
}

void HierarchicalTimeBinScheduler::applyPendingTransitions() {
  for (const std::uint32_t element : m_active_elements) {
    const std::uint8_t pending = m_hot.pending_bin_index[element];
    if (pending == k_unset_pending_bin) {
      continue;
    }

    const std::uint8_t old_bin = m_hot.bin_index[element];
    const std::uint8_t new_bin = clampBin(pending);
    if (new_bin == old_bin) {
      m_hot.pending_bin_index[element] = k_unset_pending_bin;
      continue;
    }

    const std::uint64_t new_period = binPeriodTicks(new_bin);
    if (m_current_tick % new_period != 0) {
      ++m_diagnostics.illegal_transition_attempts;
      throw std::runtime_error(timeBinContextMessage(
          "pending bin transition violates current synchronization boundary",
          "HierarchicalTimeBinScheduler::applyPendingTransitions",
          m_current_tick,
          element,
          old_bin,
          m_hot.next_activation_tick[element]));
    }

    eraseFromBin(element, old_bin);
    insertIntoBin(element, new_bin);
    m_hot.bin_index[element] = new_bin;

    if (new_bin < old_bin) {
      ++m_diagnostics.promoted_elements;
    } else {
      ++m_diagnostics.demoted_elements;
    }

    m_hot.pending_bin_index[element] = k_unset_pending_bin;
  }
}

TimeBinMappingResult mapDtToTimeBin(double dt_time_code, const TimeStepLimits& limits) {
  if (limits.min_dt_time_code <= 0.0 || limits.max_dt_time_code <= 0.0) {
    throw std::invalid_argument("timestep limits must be positive");
  }
  if (limits.max_dt_time_code < limits.min_dt_time_code) {
    throw std::invalid_argument("max_dt_time_code must be >= min_dt_time_code");
  }

  TimeBinMappingResult result{};
  const double clamped_dt = std::clamp(dt_time_code, limits.min_dt_time_code, limits.max_dt_time_code);
  result.clipped_to_min = clamped_dt == limits.min_dt_time_code && dt_time_code < limits.min_dt_time_code;
  result.clipped_to_max = clamped_dt == limits.max_dt_time_code && dt_time_code > limits.max_dt_time_code;

  std::uint8_t mapped_bin = 0;
  for (std::uint8_t bin = 0; bin <= limits.max_bin; ++bin) {
    const double dt_for_bin = binIndexToDt(bin, limits);
    if (dt_for_bin <= clamped_dt) {
      mapped_bin = bin;
    } else {
      break;
    }
  }
  result.bin_index = mapped_bin;
  return result;
}

double binIndexToDt(std::uint8_t bin_index, const TimeStepLimits& limits) {
  return limits.min_dt_time_code * static_cast<double>(powerOfTwo(std::min(bin_index, limits.max_bin)));
}

void syncTimeBinMirrorsFromScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    SimulationState& state,
    TimeBinMirrorDomain domain) {
  const auto persistent = scheduler.exportPersistentState();
  auto sync_particles = [&]() {
    if (persistent.bin_index.size() < state.particles.size()) {
      throw std::invalid_argument("syncTimeBinMirrorsFromScheduler: scheduler lacks particle bin entries");
    }
    for (std::size_t i = 0; i < state.particles.size(); ++i) {
      state.particles.time_bin[i] = persistent.bin_index[i];
    }
  };
  auto sync_cells = [&]() {
    if (persistent.bin_index.size() < state.cells.size()) {
      throw std::invalid_argument("syncTimeBinMirrorsFromScheduler: scheduler lacks cell bin entries");
    }
    for (std::size_t i = 0; i < state.cells.size(); ++i) {
      state.cells.time_bin[i] = persistent.bin_index[i];
    }
  };

  switch (domain) {
    case TimeBinMirrorDomain::kParticles:
      sync_particles();
      break;
    case TimeBinMirrorDomain::kCells:
      sync_cells();
      break;
    case TimeBinMirrorDomain::kParticlesAndCells:
      sync_particles();
      sync_cells();
      break;
  }
}

bool timeBinMirrorsMatchScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    TimeBinMirrorDomain domain) {
  const auto& hot = scheduler.hotMetadata();
  auto particles_match = [&]() {
    if (hot.bin_index.size() < state.particles.size()) {
      return false;
    }
    for (std::size_t i = 0; i < state.particles.size(); ++i) {
      if (state.particles.time_bin[i] != hot.bin_index[i]) {
        return false;
      }
    }
    return true;
  };
  auto cells_match = [&]() {
    if (hot.bin_index.size() < state.cells.size()) {
      return false;
    }
    for (std::size_t i = 0; i < state.cells.size(); ++i) {
      if (state.cells.time_bin[i] != hot.bin_index[i]) {
        return false;
      }
    }
    return true;
  };

  switch (domain) {
    case TimeBinMirrorDomain::kParticles:
      return particles_match();
    case TimeBinMirrorDomain::kCells:
      return cells_match();
    case TimeBinMirrorDomain::kParticlesAndCells:
      return particles_match() && cells_match();
  }
  return false;
}

void debugAssertTimeBinMirrorAuthorityInvariant(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    TimeBinMirrorDomain domain) {
  const auto& hot = scheduler.hotMetadata();
  const auto throw_particle = [&](std::uint32_t index) {
    throw std::runtime_error(timeBinContextMessage(
        "time_bin mirror authority invariant violated for particle mirror",
        "debugAssertTimeBinMirrorAuthorityInvariant",
        scheduler.currentTick(),
        index,
        index < hot.bin_index.size() ? hot.bin_index[index] : 0,
        index < hot.next_activation_tick.size() ? hot.next_activation_tick[index] : scheduler.currentTick()));
  };
  const auto throw_cell = [&](std::uint32_t index) {
    throw std::runtime_error(timeBinContextMessage(
        "time_bin mirror authority invariant violated for cell mirror",
        "debugAssertTimeBinMirrorAuthorityInvariant",
        scheduler.currentTick(),
        index,
        index < hot.bin_index.size() ? hot.bin_index[index] : 0,
        index < hot.next_activation_tick.size() ? hot.next_activation_tick[index] : scheduler.currentTick()));
  };
  if ((domain == TimeBinMirrorDomain::kParticles || domain == TimeBinMirrorDomain::kParticlesAndCells) &&
      hot.bin_index.size() < state.particles.size()) {
    throw_particle(static_cast<std::uint32_t>(hot.bin_index.size()));
  }
  if ((domain == TimeBinMirrorDomain::kCells || domain == TimeBinMirrorDomain::kParticlesAndCells) &&
      hot.bin_index.size() < state.cells.size()) {
    throw_cell(static_cast<std::uint32_t>(hot.bin_index.size()));
  }
  if (domain == TimeBinMirrorDomain::kParticles || domain == TimeBinMirrorDomain::kParticlesAndCells) {
    for (std::uint32_t i = 0; i < state.particles.size(); ++i) {
      if (state.particles.time_bin[i] != hot.bin_index[i]) {
        throw_particle(i);
      }
    }
  }
  if (domain == TimeBinMirrorDomain::kCells || domain == TimeBinMirrorDomain::kParticlesAndCells) {
    for (std::uint32_t i = 0; i < state.cells.size(); ++i) {
      if (state.cells.time_bin[i] != hot.bin_index[i]) {
        throw_cell(i);
      }
    }
  }
}

double computeCflTimeStep(const CflTimeStepInput& input, double c_cfl) {
  if (input.cell_width_code <= 0.0) {
    throw std::invalid_argument("cell_width_code must be positive");
  }
  if (c_cfl <= 0.0) {
    throw std::invalid_argument("c_cfl must be positive");
  }

  const double denom = std::abs(input.flow_speed_code) + std::max(input.sound_speed_code, 0.0);
  if (denom <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }

  return c_cfl * (input.cell_width_code / denom);
}

double computeGravityTimeStep(const GravityTimeStepInput& input, double eta) {
  if (input.softening_length_code <= 0.0) {
    throw std::invalid_argument("softening_length_code must be positive");
  }
  if (eta <= 0.0) {
    throw std::invalid_argument("eta must be positive");
  }

  const double accel = std::abs(input.acceleration_magnitude_code);
  if (accel == 0.0) {
    return std::numeric_limits<double>::infinity();
  }

  return eta * std::sqrt(input.softening_length_code / accel);
}

double combineTimeStepCriteria(
    std::uint32_t element_index,
    const TimeStepCriteriaHooks& hooks,
    double fallback_dt_time_code) {
  if (fallback_dt_time_code <= 0.0) {
    throw std::invalid_argument("fallback_dt_time_code must be positive");
  }

  double dt = fallback_dt_time_code;
  if (hooks.cfl_hook) {
    dt = std::min(dt, finitePositiveOrInf(hooks.cfl_hook(element_index)));
  }
  if (hooks.gravity_hook) {
    dt = std::min(dt, finitePositiveOrInf(hooks.gravity_hook(element_index)));
  }
  if (hooks.source_hook) {
    dt = std::min(dt, finitePositiveOrInf(hooks.source_hook(element_index)));
  }
  if (hooks.user_clamp_hook) {
    dt = std::min(dt, finitePositiveOrInf(hooks.user_clamp_hook(element_index)));
  }

  return dt;
}

double computeScaleFactorRate(const LambdaCdmBackground& background, double scale_factor) {
  if (scale_factor <= 0.0) {
    throw std::invalid_argument("scale_factor must be positive");
  }
  return scale_factor * background.hubbleSi(scale_factor);
}

double advanceScaleFactorEuler(
    const LambdaCdmBackground& background,
    double scale_factor,
    double dt_time_code) {
  if (dt_time_code < 0.0) {
    throw std::invalid_argument("dt_time_code must be non-negative");
  }
  return scale_factor + dt_time_code * computeScaleFactorRate(background, scale_factor);
}

double estimateDeltaTimeFromScaleFactorStep(
    const LambdaCdmBackground& background,
    double scale_factor,
    double delta_scale_factor) {
  if (delta_scale_factor < 0.0) {
    throw std::invalid_argument("delta_scale_factor must be non-negative");
  }

  const double rate = computeScaleFactorRate(background, scale_factor);
  if (rate <= 0.0) {
    throw std::invalid_argument("scale-factor rate must be positive");
  }
  return delta_scale_factor / rate;
}

double computeComovingDriftFactor(
    const LambdaCdmBackground& background,
    double scale_factor_begin,
    double scale_factor_end,
    std::uint32_t midpoint_samples) {
  return midpointIntegrateDriftLike(background, scale_factor_begin, scale_factor_end, midpoint_samples);
}

double computeComovingKickFactor(
    const LambdaCdmBackground& background,
    double scale_factor_begin,
    double scale_factor_end,
    std::uint32_t midpoint_samples) {
  return midpointIntegrateDriftLike(background, scale_factor_begin, scale_factor_end, midpoint_samples);
}

double computeHubbleDragFactor(double scale_factor_begin, double scale_factor_end) {
  if (scale_factor_begin <= 0.0 || scale_factor_end <= 0.0) {
    throw std::invalid_argument("scale factors must be positive");
  }
  if (scale_factor_end < scale_factor_begin) {
    throw std::invalid_argument("scale_factor_end must be >= scale_factor_begin");
  }

  return scale_factor_begin / scale_factor_end;
}

}  // namespace cosmosim::core
