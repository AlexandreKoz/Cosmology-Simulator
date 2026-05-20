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

class HierarchicalTimeBinScheduler;
struct IntegratorState;

// Explicit step stage contract shared by gravity, hydro, physics, analysis, and I/O modules.
enum class IntegrationStage : std::uint8_t {
  kGravityKickPre = 0,
  kDrift = 1,
  kForceRefresh = 2,
  kHydroUpdate = 3,
  kSourceTerms = 4,
  kGravityKickPost = 5,
  kAnalysisHooks = 6,
  kOutputCheck = 7,
};
enum class StageDataDomain : std::uint16_t {
  kNone = 0,
  kParticles = 1U << 0U,
  kGasCells = 1U << 1U,
  kMeshCells = 1U << 2U,
  kGhostCells = 1U << 3U,
  kPmField = 1U << 4U,
  kTreeState = 1U << 5U,
  kOutputState = 1U << 6U,
  kRestartState = 1U << 7U,
  kDiagnostics = 1U << 8U,
};
[[nodiscard]] constexpr StageDataDomain operator|(StageDataDomain lhs, StageDataDomain rhs) noexcept {
  return static_cast<StageDataDomain>(static_cast<std::uint16_t>(lhs) | static_cast<std::uint16_t>(rhs));
}
[[nodiscard]] constexpr bool hasStageDataDomain(StageDataDomain value, StageDataDomain mask) noexcept {
  return (static_cast<std::uint16_t>(value) & static_cast<std::uint16_t>(mask)) != 0U;
}
enum class StageSyncRequirement : std::uint8_t { kNone = 0, kLocalOnly = 1, kPmRefreshBoundary = 2, kGlobal = 3 };
enum class StageActiveSetFamily : std::uint8_t {
  kNone = 0,
  kAllParticles = 1,
  kActiveParticles = 2,
  kGasCells = 3,
  kMeshCells = 4,
  kGhostCells = 5,
  kOutputState = 6,
  kRestartState = 7,
};
enum class StageSafety : std::uint8_t { kUnsafe = 0, kSafe = 1 };
enum class StageSubsystem : std::uint8_t { kCore = 0, kGravity = 1, kHydro = 2, kSources = 3, kAnalysis = 4, kOutput = 5 };
struct StageContract {
  IntegrationStage stage = IntegrationStage::kGravityKickPre;
  StageDataDomain required_inputs = StageDataDomain::kNone;
  StageDataDomain mutated_state = StageDataDomain::kNone;
  StageDataDomain produced_outputs = StageDataDomain::kNone;
  StageDataDomain allowed_side_effects = StageDataDomain::kNone;
  StageSyncRequirement sync_requirements = StageSyncRequirement::kNone;
  StageActiveSetFamily active_set_family = StageActiveSetFamily::kNone;
  StageSafety restart_safety = StageSafety::kUnsafe;
  StageSafety output_safety = StageSafety::kUnsafe;
  StageSubsystem owner = StageSubsystem::kCore;
};

[[nodiscard]] std::string_view integrationStageName(IntegrationStage stage);
[[nodiscard]] constexpr std::size_t integrationStageCount() noexcept { return 8; }
[[nodiscard]] constexpr std::size_t integrationStageIndex(IntegrationStage stage) noexcept {
  return static_cast<std::size_t>(stage);
}

// Baseline stepping family; hierarchical bins can reuse the same stage contract.
enum class TimeStepScheme : std::uint8_t {
  kKickDriftKick = 0,
};

// Explicit integration-boundary classes distinguish local active-bin work from
// globally coherent restart/output points and legal PM-refresh surfaces.
enum class StepBoundaryKind : std::uint8_t {
  kLocalActiveBinStep = 0,
  kGlobalSynchronizationPoint = 1,
  kPmRefreshPoint = 2,
  kSnapshotPoint = 3,
  kCheckpointPoint = 4,
};

[[nodiscard]] std::string_view stepBoundaryKindName(StepBoundaryKind kind);
[[nodiscard]] bool isRestartSafeBoundary(StepBoundaryKind kind) noexcept;
[[nodiscard]] bool isOutputSafeBoundary(StepBoundaryKind kind) noexcept;

struct StepBoundaryState {
  StepBoundaryKind kind = StepBoundaryKind::kGlobalSynchronizationPoint;
  bool restart_safe = true;
  bool output_safe = true;
  bool pm_refresh_allowed = true;
  bool local_substep = false;
};

// Integrator-issued TreePM refresh directive.  Gravity callbacks may own solver
// buffers, but they must consume this placement contract instead of inventing
// legal refresh surfaces internally.  Local-bin force-refresh surfaces require a
// coherent predicted source view because inactive particles are not physically
// drifted in the persistent state during the local substep.
struct PmRefreshDirective {
  enum class Reason : std::uint8_t {
    kNone = 0,
    kInitialForceBootstrap = 1,
    kScheduledForceRefreshStage = 2,
  };
  bool force_refresh_surface = false;
  bool cadence_opportunity_allowed = false;
  bool initial_cache_bootstrap_allowed = false;
  bool requires_predicted_inactive_sources = false;
  bool has_sync_event = false;
  bool refresh_long_range_field = false;
  bool solver_executed = false;
  std::uint64_t gravity_kick_opportunity = 0;
  std::uint64_t field_version = 0;
  std::uint64_t last_refresh_opportunity = 0;
  std::uint64_t field_built_step_index = 0;
  double field_built_scale_factor = 1.0;
  double force_evaluation_scale_factor = 1.0;
  Reason reason = Reason::kNone;
};

struct CosmologicalStepFactors {
  bool cosmological = false;
  double time_begin_code = 0.0;
  double time_end_code = 0.0;
  double dt_time_code = 0.0;
  double time_si_per_code = 1.0;
  double dt_time_si = 0.0;
  double scale_factor_begin = 1.0;
  double scale_factor_midpoint = 1.0;
  double scale_factor_end = 1.0;
  double redshift_begin = 0.0;
  double redshift_end = 0.0;
  double hubble_begin_code = 0.0;
  double hubble_end_code = 0.0;
  double drift_factor_code = 0.0;
  // Particle velocity lanes named velocity_*_peculiar store physical peculiar
  // velocities u = a dx/dt.  The KDK kick sub-operators therefore apply the
  // exact homogeneous Hubble drag u <- (a_begin/a_end) u and a distinct
  // force kick for comoving acceleration kernels, not the drift integral.
  double first_kick_factor_code = 0.0;
  double second_kick_factor_code = 0.0;
  double first_hubble_drag_factor = 1.0;
  double second_hubble_drag_factor = 1.0;
  double hubble_drag_factor = 1.0;
};

class CosmologicalTimeline {
 public:
  explicit CosmologicalTimeline(const LambdaCdmBackground* background = nullptr, double time_si_per_code = 1.0);

  [[nodiscard]] CosmologicalStepFactors prepareStep(
      double current_time_code,
      double current_scale_factor,
      double dt_time_code) const;

  void commitStep(IntegratorState& integrator_state, const CosmologicalStepFactors& step) const;

 private:
  const LambdaCdmBackground* m_background = nullptr;
  double m_time_si_per_code = 1.0;
};

// Hierarchical stepping metadata kept outside particle arrays for auditable ownership.
struct TimeBinContext {
  bool hierarchical_enabled = false;
  std::uint8_t active_bin = 0;
  std::uint8_t max_bin = 0;
};

struct PmSyncEvent {
  std::uint64_t gravity_kick_opportunity = 0;
  bool refresh_long_range_field = false;
  std::uint64_t field_version = 0;
  std::uint64_t last_refresh_opportunity = 0;
  std::uint64_t field_built_step_index = 0;
  double field_built_scale_factor = 1.0;
};

// Restartable PM-cadence authority. The TreePM solver may own field buffers, but
// this state owns the legality of long-range refresh opportunities across restarts.
struct PmSynchronizationPersistentState {
  std::uint64_t cadence_steps = 1;
  std::uint64_t gravity_kick_opportunity = 0;
  std::uint64_t last_refresh_opportunity = 0;
  std::uint64_t field_version = 0;
  std::uint64_t last_refresh_step_index = 0;
  double last_refresh_scale_factor = 1.0;
  bool refresh_commit_pending = false;
  std::uint64_t pending_refresh_opportunity = 0;
  std::uint64_t pending_refresh_field_version = 0;
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
  [[nodiscard]] bool refreshCommitPending() const noexcept { return m_refresh_commit_pending; }

  [[nodiscard]] PmSynchronizationPersistentState exportPersistentState() const;
  void importPersistentState(const PmSynchronizationPersistentState& persistent_state);

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


// Persistent integrator state tracked by the orchestrator.
struct IntegratorState {
  double current_time_code = 0.0;
  double current_scale_factor = 1.0;
  double current_redshift = 0.0;
  double current_hubble_rate_code = 0.0;
  double time_si_per_code = 1.0;
  double dt_time_code = 0.0;
  double last_drift_factor_code = 0.0;
  double last_first_kick_factor_code = 0.0;
  double last_second_kick_factor_code = 0.0;
  double last_first_hubble_drag_factor = 1.0;
  double last_second_hubble_drag_factor = 1.0;
  StepBoundaryKind current_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint;
  StepBoundaryKind last_completed_boundary_kind = StepBoundaryKind::kCheckpointPoint;
  bool inside_kdk_step = false;
  bool last_completed_restart_safe = true;
  std::uint64_t step_index = 0;
  TimeStepScheme scheme = TimeStepScheme::kKickDriftKick;
  TimeBinContext time_bins;
  PmSynchronizationState pm_sync_state;
  bool pm_refresh_enabled = false;
  bool pm_long_range_field_valid = false;
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
  CosmologicalStepFactors timeline_step;
  StepBoundaryState boundary;
  PmRefreshDirective pm_refresh_directive;
  IntegrationStage stage = IntegrationStage::kGravityKickPre;
};

// Stage-bound handler interface implemented by gravity, hydro, source, analysis, and output modules.
// Ownership stays with the caller; StepOrchestrator stores non-owning pointers in
// deterministic registration order.  Each handler declares the exact stage set it
// may receive, and production dispatch only iterates the current stage bucket.
class IntegrationCallback {
 public:
  virtual ~IntegrationCallback() = default;

  [[nodiscard]] virtual std::string_view callbackName() const = 0;
  [[nodiscard]] virtual std::span<const IntegrationStage> integrationStages() const = 0;
  [[nodiscard]] virtual std::span<const StageContract> stageContracts() const { return {}; }
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
  [[nodiscard]] std::span<IntegrationCallback* const> handlersFor(IntegrationStage stage) const noexcept;
  [[nodiscard]] std::span<const StageContract> contractsFor(IntegrationStage stage) const noexcept;
  [[nodiscard]] std::optional<StageContract> contractForHandlerStage(
      const IntegrationCallback& callback,
      IntegrationStage stage) const noexcept;

  void executeOutputBoundary(
      SimulationState& state,
      IntegratorState& integrator_state,
      ProfilerSession* profiler_session = nullptr,
      StepBoundaryKind requested_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint) const;

  void executeSingleStep(
      SimulationState& state,
      IntegratorState& integrator_state,
      ActiveSetDescriptor active_set,
      const LambdaCdmBackground* cosmology_background,
      TransientStepWorkspace* workspace = nullptr,
      const ModePolicy* mode_policy = nullptr,
      ProfilerSession* profiler_session = nullptr,
      std::optional<std::uint64_t> expected_scheduler_tick = std::nullopt,
      StepBoundaryKind requested_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint) const;

  void executeSchedulerSubstep(
      SimulationState& state,
      IntegratorState& integrator_state,
      const HierarchicalTimeBinScheduler& scheduler,
      std::span<const std::uint32_t> active_particle_indices,
      std::span<const std::uint32_t> active_cell_indices,
      const LambdaCdmBackground* cosmology_background,
      TransientStepWorkspace* workspace = nullptr,
      const ModePolicy* mode_policy = nullptr,
      ProfilerSession* profiler_session = nullptr,
      StepBoundaryKind requested_boundary_kind = StepBoundaryKind::kGlobalSynchronizationPoint) const;

 private:
  StageScheduler m_scheduler;
  std::array<std::vector<IntegrationCallback*>, integrationStageCount()> m_handlers_by_stage;
  std::array<std::vector<StageContract>, integrationStageCount()> m_contracts_by_stage;
  std::size_t m_callback_count = 0;
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
  kCosmologyExpansion = 2,
  kSourceTerm = 3,
  kUserClamp = 4,
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
  std::array<std::uint32_t, 5> limiting_candidates_by_source{};
  TimeStepCandidateSource dominant_limiting_source = TimeStepCandidateSource::kUserClamp;
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
  void appendElements(std::uint32_t new_element_count, std::uint8_t initial_bin, std::uint64_t first_activation_tick);
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
  std::vector<TimeStepCandidateSource> m_candidate_source;
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

// Scheduler-backed reorder by rung. This is the only production-safe time-bin
// ordering path: it consumes scheduler authority directly instead of public
// particle/cell mirrors.
[[nodiscard]] ParticleReorderMap buildParticleReorderMapByScheduler(
    const SimulationState& state,
    const HierarchicalTimeBinScheduler& scheduler);


// Stable scheduler identity namespace. Current production scheduling is particle-row
// backed, with gas cells treated as particle-bound finite-volume carriers. Future AMR
// patches or mesh cells must enter through an explicit kind/stable_id pair rather
// than by overloading local row indices.
enum class ScheduledElementKind : std::uint8_t {
  kParticle = 0,
  kParticleBoundGasCell = 1,
  kAmrPatch = 2,
};

struct ScheduledElementKey {
  ScheduledElementKind kind = ScheduledElementKind::kParticle;
  std::uint64_t stable_id = 0;

  [[nodiscard]] friend bool operator==(const ScheduledElementKey&, const ScheduledElementKey&) = default;
};

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

[[nodiscard]] double computeCosmologyExpansionTimeStep(
    const LambdaCdmBackground& background,
    double scale_factor,
    double max_delta_ln_a,
    double max_hubble_time_fraction,
    double time_si_per_code = 1.0);

// da/dt = a H(a) for standard FLRW backgrounds.
[[nodiscard]] double computeScaleFactorRate(const LambdaCdmBackground& background, double scale_factor);

// Production scale-factor evolution: invert the FLRW time integral dt = integral da/(a H(a)).
[[nodiscard]] double advanceScaleFactorByCosmicTime(
    const LambdaCdmBackground& background,
    double scale_factor,
    double dt_time_code,
    std::uint32_t midpoint_samples = 64);

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

[[nodiscard]] StepBoundaryState classifyStepBoundary(
    const SimulationState& state,
    const ActiveSetDescriptor& active_set,
    bool scheduler_owned_substep,
    StepBoundaryKind requested_kind = StepBoundaryKind::kGlobalSynchronizationPoint);

void assertCanWriteSnapshotAtBoundary(const IntegratorState& integrator_state);
void assertCanWriteCheckpointAtBoundary(const IntegratorState& integrator_state);

}  // namespace cosmosim::core
