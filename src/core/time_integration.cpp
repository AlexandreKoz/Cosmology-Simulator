#include "cosmosim/core/time_integration.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
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

}  // namespace

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
    ProfilerSession* profiler_session) const {
  if (integrator_state.dt_time_code <= 0.0) {
    throw std::invalid_argument("dt_time_code must be positive");
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

void HierarchicalTimeBinScheduler::requestBinTransition(
    std::uint32_t element_index,
    std::uint8_t target_bin) {
  if (element_index >= m_hot.size()) {
    throw std::out_of_range("element_index out of range");
  }
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
  rebuildActiveSet();
  return m_active_elements;
}

void HierarchicalTimeBinScheduler::endSubstep() {
  applyPendingTransitions();
  for (const std::uint32_t element : m_active_elements) {
    const std::uint8_t bin = m_hot.bin_index[element];
    const std::uint64_t period_ticks = binPeriodTicks(bin);
    m_hot.next_activation_tick[element] = m_current_tick + period_ticks;
    m_hot.active_flag[element] = 0;
  }
  ++m_current_tick;
}

const TimeBinHotMetadata& HierarchicalTimeBinScheduler::hotMetadata() const noexcept { return m_hot; }

const TimeBinDiagnostics& HierarchicalTimeBinScheduler::diagnostics() const noexcept { return m_diagnostics; }

TimeBinPersistentState HierarchicalTimeBinScheduler::exportPersistentState() const {
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

  m_current_tick = persistent_state.current_tick;
  m_max_bin = persistent_state.max_bin;
  m_hot.bin_index = persistent_state.bin_index;
  m_hot.next_activation_tick = persistent_state.next_activation_tick;
  m_hot.active_flag = persistent_state.active_flag;
  m_hot.pending_bin_index = persistent_state.pending_bin_index;

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
  m_diagnostics = {};
  m_diagnostics.occupancy_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, 0U);
  m_diagnostics.active_count_by_bin.assign(static_cast<std::size_t>(m_max_bin) + 1U, 0U);
  for (std::size_t bin = 0; bin < m_elements_by_bin.size(); ++bin) {
    m_diagnostics.occupancy_by_bin[bin] = static_cast<std::uint32_t>(m_elements_by_bin[bin].size());
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
      continue;
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
  const auto persistent = scheduler.exportPersistentState();
  auto particles_match = [&]() {
    if (persistent.bin_index.size() < state.particles.size()) {
      return false;
    }
    for (std::size_t i = 0; i < state.particles.size(); ++i) {
      if (state.particles.time_bin[i] != persistent.bin_index[i]) {
        return false;
      }
    }
    return true;
  };
  auto cells_match = [&]() {
    if (persistent.bin_index.size() < state.cells.size()) {
      return false;
    }
    for (std::size_t i = 0; i < state.cells.size(); ++i) {
      if (state.cells.time_bin[i] != persistent.bin_index[i]) {
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
  if (!timeBinMirrorsMatchScheduler(scheduler, state, domain)) {
    throw std::runtime_error(
        "time-bin mirror authority invariant violated: state mirrors diverged from scheduler authority");
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
