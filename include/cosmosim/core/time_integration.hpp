#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::core {

// Explicit step stage contract shared by gravity, hydro, physics, analysis, and I/O modules.
enum class IntegrationStage : std::uint8_t {
  kGravityKickPre = 0,
  kDrift = 1,
  kHydroUpdate = 2,
  kSourceTerms = 3,
  kGravityKickPost = 4,
  kAnalysisHooks = 5,
  kOutputCheck = 6,
};

[[nodiscard]] std::string_view integrationStageName(IntegrationStage stage);

// Baseline stepping family; hierarchical bins can reuse the same stage contract.
enum class TimeStepScheme : std::uint8_t {
  kKickDriftKick = 0,
};

// Hierarchical stepping metadata kept outside particle arrays for auditable ownership.
struct TimeBinContext {
  bool hierarchical_enabled = false;
  std::uint8_t active_bin = 0;
  std::uint8_t max_bin = 0;
};

// Persistent integrator state tracked by the orchestrator.
struct IntegratorState {
  double current_time_code = 0.0;
  double current_scale_factor = 1.0;
  double dt_time_code = 0.0;
  std::uint64_t step_index = 0;
  TimeStepScheme scheme = TimeStepScheme::kKickDriftKick;
  TimeBinContext time_bins;
};

// Explicit compact active-set descriptor with optional subset spans.
struct ActiveSetDescriptor {
  std::span<const std::uint32_t> particle_indices;
  std::span<const std::uint32_t> cell_indices;
  bool particles_are_subset = false;
  bool cells_are_subset = false;
  // Runtime-truth metadata: active descriptors are derived products of the
  // scheduler/state generation that created them, never independent active-set
  // authority. Stale descriptors must fail before solver kernels consume them.
  bool particles_from_scheduler = false;
  bool cells_from_scheduler = false;
  bool has_generation_metadata = false;
  std::uint64_t source_particle_index_generation = 0;
  std::uint64_t source_cell_index_generation = 0;
  std::uint64_t source_scheduler_tick = 0;

  [[nodiscard]] bool hasParticleSubset(std::size_t total_particle_count) const noexcept;
  [[nodiscard]] bool hasCellSubset(std::size_t total_cell_count) const noexcept;
};

// Per-stage execution context passed to callbacks without mutating interface contracts.
struct StepContext {
  SimulationState& state;
  IntegratorState& integrator_state;
  ActiveSetDescriptor active_set;
  TransientStepWorkspace* workspace = nullptr;
  const LambdaCdmBackground* cosmology_background = nullptr;
  const ModePolicy* mode_policy = nullptr;
  ProfilerSession* profiler_session = nullptr;
  IntegrationStage stage = IntegrationStage::kGravityKickPre;
};

// Callback interface implemented by gravity, hydro, source, analysis, and output modules.
class IntegrationCallback {
 public:
  virtual ~IntegrationCallback() = default;

  [[nodiscard]] virtual std::string_view callbackName() const = 0;
  virtual void onStage(StepContext& context) = 0;
};

// Stage scheduler isolates ordering from solver implementation details.
class StageScheduler {
 public:
  [[nodiscard]] std::vector<IntegrationStage> schedule(
      const IntegratorState& integrator_state,
      const ActiveSetDescriptor& active_set) const;

  [[nodiscard]] static std::span<const IntegrationStage> kickDriftKickOrder();
};

[[nodiscard]] bool isCanonicalIntegrationStageOrder(std::span<const IntegrationStage> ordered_stages);

// Single authoritative step orchestrator for current baseline stepping.
class StepOrchestrator {
 public:
  explicit StepOrchestrator(StageScheduler scheduler = {});

  void registerCallback(IntegrationCallback& callback);
  [[nodiscard]] std::size_t callbackCount() const noexcept;

  void executeSingleStep(
      SimulationState& state,
      IntegratorState& integrator_state,
      ActiveSetDescriptor active_set,
      const LambdaCdmBackground* cosmology_background,
      TransientStepWorkspace* workspace = nullptr,
      const ModePolicy* mode_policy = nullptr,
      ProfilerSession* profiler_session = nullptr,
      std::optional<std::uint64_t> expected_scheduler_tick = std::nullopt) const;

 private:
  StageScheduler m_scheduler;
  std::vector<IntegrationCallback*> m_callbacks;
};

// Hot metadata sidecar for element-local time-bin ownership and scheduling state.
struct TimeBinHotMetadata {
  std::vector<std::uint8_t> bin_index;
  std::vector<std::uint64_t> next_activation_tick;
  std::vector<std::uint8_t> active_flag;
  std::vector<std::uint8_t> pending_bin_index;

  [[nodiscard]] std::size_t size() const noexcept { return bin_index.size(); }
};

// Diagnostics emitted by hierarchical scheduler for auditing and pathological collapse detection.
struct TimeBinDiagnostics {
  std::vector<std::uint32_t> occupancy_by_bin;
  std::vector<std::uint32_t> active_count_by_bin;
  std::uint32_t active_elements = 0;
  std::uint32_t promoted_elements = 0;
  std::uint32_t demoted_elements = 0;
  std::uint32_t clipped_to_min_dt = 0;
  std::uint32_t clipped_to_max_dt = 0;
  std::uint32_t illegal_transition_attempts = 0;
  std::uint32_t collapse_candidates = 0;
  double active_fraction = 0.0;
  std::uint8_t most_active_bin = 0;
};

struct TimeBinMappingResult {
  std::uint8_t bin_index = 0;
  bool clipped_to_min = false;
  bool clipped_to_max = false;
};

// Auditable labels attached to solver timestep criteria before scheduler reconciliation.
enum class TimeStepCandidateSource : std::uint8_t {
  kHydroCfl = 0,
  kGravityAcceleration = 1,
  kSourceTerm = 2,
  kUserClamp = 3,
};

[[nodiscard]] std::string_view timeStepCandidateSourceName(TimeStepCandidateSource source);

struct TimeStepCandidateSubmission {
  std::uint32_t element_index = 0;
  double dt_time_code = 0.0;
  TimeStepCandidateSource source = TimeStepCandidateSource::kUserClamp;
  std::string label;
};

struct TimeStepReconciliationResult {
  std::uint32_t submitted_candidates = 0;
  std::uint32_t elements_with_candidates = 0;
  std::uint32_t committed_transition_requests = 0;
  std::uint32_t clipped_to_min_dt = 0;
  std::uint32_t clipped_to_max_dt = 0;
};

struct PmSyncEvent {
  std::uint64_t gravity_kick_opportunity = 0;
  bool refresh_long_range_field = false;
  std::uint64_t field_version = 0;
  std::uint64_t last_refresh_opportunity = 0;
  std::uint64_t field_built_step_index = 0;
  double field_built_scale_factor = 1.0;
};

// Scheduler-owned PM cadence state. Gravity callbacks execute the solver, but cadence
// legality and refresh boundaries are represented here to avoid free-floating PM truth.
class PmSynchronizationState {
 public:
  void reset(std::uint64_t cadence_steps = 1);
  [[nodiscard]] PmSyncEvent registerKickOpportunity(
      std::uint64_t step_index,
      double scale_factor,
      bool has_long_range_field);
  void commitRefresh(const PmSyncEvent& event);

  [[nodiscard]] std::uint64_t cadenceSteps() const noexcept { return m_cadence_steps; }
  [[nodiscard]] std::uint64_t gravityKickOpportunity() const noexcept { return m_gravity_kick_opportunity; }
  [[nodiscard]] std::uint64_t fieldVersion() const noexcept { return m_field_version; }
  [[nodiscard]] std::uint64_t lastRefreshOpportunity() const noexcept { return m_last_refresh_opportunity; }
  [[nodiscard]] std::uint64_t lastRefreshStepIndex() const noexcept { return m_last_refresh_step_index; }
  [[nodiscard]] double lastRefreshScaleFactor() const noexcept { return m_last_refresh_scale_factor; }

 private:
  std::uint64_t m_cadence_steps = 1;
  std::uint64_t m_gravity_kick_opportunity = 0;
  std::uint64_t m_last_refresh_opportunity = 0;
  std::uint64_t m_field_version = 0;
  std::uint64_t m_last_refresh_step_index = 0;
  double m_last_refresh_scale_factor = 1.0;
  bool m_refresh_commit_pending = false;
  std::uint64_t m_pending_refresh_opportunity = 0;
  std::uint64_t m_pending_refresh_field_version = 0;
};

// Persisted scheduler state required for exact restart continuation.
struct TimeBinPersistentState {
  std::uint64_t current_tick = 0;
  std::uint8_t max_bin = 0;
  std::vector<std::uint8_t> bin_index;
  std::vector<std::uint64_t> next_activation_tick;
  std::vector<std::uint8_t> active_flag;
  std::vector<std::uint8_t> pending_bin_index;
};

// Typed limits that normalize physical timestep proposals into the discrete bin hierarchy.
struct TimeStepLimits {
  double min_dt_time_code = 0.0;
  double max_dt_time_code = 0.0;
  std::uint8_t max_bin = 0;
};

// Compact inputs for conservative CFL and gravity criteria hooks.
struct CflTimeStepInput {
  double cell_width_code = 0.0;
  double flow_speed_code = 0.0;
  double sound_speed_code = 0.0;
};

struct GravityTimeStepInput {
  double softening_length_code = 0.0;
  double acceleration_magnitude_code = 0.0;
};

using CriteriaHook = std::function<double(std::uint32_t)>;

struct TimeStepCriteriaHooks {
  CriteriaHook cfl_hook;
  CriteriaHook gravity_hook;
  CriteriaHook source_hook;
  CriteriaHook user_clamp_hook;
};

class TimeStepCriteriaRegistry {
 public:
  void registerCflHook(CriteriaHook hook);
  void registerGravityHook(CriteriaHook hook);
  void registerSourceHook(CriteriaHook hook);
  void registerUserClampHook(CriteriaHook hook);

  [[nodiscard]] const TimeStepCriteriaHooks& hooks() const noexcept;

 private:
  TimeStepCriteriaHooks m_hooks;
};

// Integer timeline scheduler with power-of-two bins and compact active set extraction.
class HierarchicalTimeBinScheduler {
 public:
  static constexpr std::uint8_t k_unset_pending_bin = 0xFF;

  explicit HierarchicalTimeBinScheduler(std::uint8_t max_bin = 0);

  void reset(std::uint32_t element_count, std::uint8_t initial_bin, std::uint64_t start_tick = 0);
  void setElementBin(std::uint32_t element_index, std::uint8_t bin_index, std::uint64_t current_tick);
  void submitCandidateTimeStep(
      std::uint32_t element_index,
      double dt_time_code,
      const TimeStepLimits& limits,
      TimeStepCandidateSource source,
      std::string_view label = {});
  void submitCandidateBin(
      std::uint32_t element_index,
      std::uint8_t target_bin,
      TimeStepCandidateSource source,
      std::string_view label = {});
  [[nodiscard]] TimeStepReconciliationResult reconcileCandidateTransitions();

  [[nodiscard]] std::span<const std::uint32_t> activeElements() const noexcept;

  [[nodiscard]] std::uint64_t currentTick() const noexcept;
  [[nodiscard]] std::uint8_t maxBin() const noexcept;
  [[nodiscard]] std::uint32_t elementCount() const noexcept;

  [[nodiscard]] bool isBinActiveAtTick(std::uint8_t bin_index, std::uint64_t tick) const;
  [[nodiscard]] std::uint64_t binPeriodTicks(std::uint8_t bin_index) const;

  std::span<const std::uint32_t> beginSubstep();
  void endSubstep();

  [[nodiscard]] const TimeBinHotMetadata& hotMetadata() const noexcept;
  [[nodiscard]] const TimeBinDiagnostics& diagnostics() const noexcept;

  [[nodiscard]] TimeBinPersistentState exportPersistentState() const;
  void importPersistentState(const TimeBinPersistentState& persistent_state);

 private:
  std::uint8_t clampBin(std::uint8_t requested) const noexcept;
  void eraseFromBin(std::uint32_t element_index, std::uint8_t bin_index);
  void insertIntoBin(std::uint32_t element_index, std::uint8_t bin_index);
  void rebuildActiveSet();
  void requestBinTransition(std::uint32_t element_index, std::uint8_t target_bin);
  void applyPendingTransitions();
  void validateInternalState(std::string_view source_label) const;
  void validateTransitionRequest(
      std::uint32_t element_index,
      std::uint8_t target_bin,
      std::string_view source_label) const;

  std::uint64_t m_current_tick = 0;
  std::uint8_t m_max_bin = 0;
  TimeBinHotMetadata m_hot;
  std::vector<std::vector<std::uint32_t>> m_elements_by_bin;
  std::vector<std::size_t> m_position_in_bin;
  std::vector<std::uint32_t> m_active_elements;
  TimeBinDiagnostics m_diagnostics;
  std::vector<std::uint8_t> m_candidate_bin_index;
  std::vector<std::string> m_candidate_label;
  TimeStepReconciliationResult m_last_reconciliation;
  bool m_substep_open = false;
};

[[nodiscard]] ActiveSetDescriptor makeSchedulerActiveSetDescriptor(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint32_t> active_cell_indices = {});

void debugAssertActiveSetDescriptorFresh(
    const ActiveSetDescriptor& active_set,
    const SimulationState& state);

void debugAssertActiveSetDescriptorFresh(
    const ActiveSetDescriptor& active_set,
    const SimulationState& state,
    std::uint64_t expected_scheduler_tick);

void debugAssertActiveSetDescriptorFresh(
    const ActiveSetDescriptor& active_set,
    const SimulationState& state,
    const HierarchicalTimeBinScheduler& scheduler);

struct TimeBinSchedulerIdentityRecord {
  // Full scheduler authority for one physical element identity. Migration packets
  // may carry time_bin as a diagnostic mirror, but exact continuation requires
  // these scheduler-owned lanes matched by stable particle_id/gas_cell_id.
  std::uint64_t element_id = 0;
  std::uint8_t bin_index = 0;
  std::uint64_t next_activation_tick = 0;
  std::uint8_t pending_bin_index = HierarchicalTimeBinScheduler::k_unset_pending_bin;
};

[[nodiscard]] std::vector<TimeBinSchedulerIdentityRecord> exportParticleSchedulerIdentityRecords(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    std::span<const std::uint32_t> particle_indices);

[[nodiscard]] TimeBinPersistentState remapSchedulerPersistentStateByParticleId(
    const TimeBinPersistentState& source_state,
    std::span<const std::uint64_t> source_particle_ids,
    std::span<const std::uint64_t> destination_particle_ids);

void remapSchedulerByParticleId(
    HierarchicalTimeBinScheduler& scheduler,
    std::span<const std::uint64_t> source_particle_ids,
    std::span<const std::uint64_t> destination_particle_ids);

void remapSchedulerByParticleReorderMap(
    HierarchicalTimeBinScheduler& scheduler,
    const ParticleReorderMap& reorder_map);

[[nodiscard]] TimeBinPersistentState rebuildSchedulerPersistentStateFromIdentityRecords(
    std::uint64_t current_tick,
    std::uint8_t max_bin,
    std::span<const TimeBinSchedulerIdentityRecord> records,
    std::span<const std::uint64_t> destination_element_ids);

void rebuildSchedulerFromParticleIdentityRecords(
    HierarchicalTimeBinScheduler& scheduler,
    std::span<const TimeBinSchedulerIdentityRecord> records,
    std::span<const std::uint64_t> destination_particle_ids);

void syncGasCellTimeBinMirrorsFromParticleScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    SimulationState& state);

[[nodiscard]] TimeBinMappingResult mapDtToTimeBin(double dt_time_code, const TimeStepLimits& limits);
[[nodiscard]] double binIndexToDt(std::uint8_t bin_index, const TimeStepLimits& limits);
enum class TimeBinMirrorDomain : std::uint8_t {
  kParticles = 0,
  kCells = 1,
  kParticlesAndCells = 2,
};

void syncTimeBinMirrorsFromScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    SimulationState& state,
    TimeBinMirrorDomain domain = TimeBinMirrorDomain::kParticles);
[[nodiscard]] bool timeBinMirrorsMatchScheduler(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    TimeBinMirrorDomain domain = TimeBinMirrorDomain::kParticles);
void debugAssertTimeBinMirrorAuthorityInvariant(
    const HierarchicalTimeBinScheduler& scheduler,
    const SimulationState& state,
    TimeBinMirrorDomain domain = TimeBinMirrorDomain::kParticles);

[[nodiscard]] double computeCflTimeStep(const CflTimeStepInput& input, double c_cfl);
[[nodiscard]] double computeGravityTimeStep(const GravityTimeStepInput& input, double eta);

[[nodiscard]] double combineTimeStepCriteria(
    std::uint32_t element_index,
    const TimeStepCriteriaHooks& hooks,
    double fallback_dt_time_code);

// da/dt = a H(a) for standard FLRW backgrounds.
[[nodiscard]] double computeScaleFactorRate(const LambdaCdmBackground& background, double scale_factor);

// Forward-Euler helper used by baseline tests and scheduler scaffolding.
[[nodiscard]] double advanceScaleFactorEuler(
    const LambdaCdmBackground& background,
    double scale_factor,
    double dt_time_code);

// dt estimate for an intended delta-a increment around the current scale factor.
[[nodiscard]] double estimateDeltaTimeFromScaleFactorStep(
    const LambdaCdmBackground& background,
    double scale_factor,
    double delta_scale_factor);

// Drift prefactor integral: integral_{a0}^{a1} da / (a^2 H(a)).
[[nodiscard]] double computeComovingDriftFactor(
    const LambdaCdmBackground& background,
    double scale_factor_begin,
    double scale_factor_end,
    std::uint32_t midpoint_samples = 16);

// Kick prefactor for comoving acceleration terms proportional to 1/a.
[[nodiscard]] double computeComovingKickFactor(
    const LambdaCdmBackground& background,
    double scale_factor_begin,
    double scale_factor_end,
    std::uint32_t midpoint_samples = 16);

// Hubble drag factor for dv/dt = -H(a) v over [a0, a1].
[[nodiscard]] double computeHubbleDragFactor(double scale_factor_begin, double scale_factor_end);

}  // namespace cosmosim::core
