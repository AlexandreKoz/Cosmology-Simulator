#pragma once

#include <cstdint>
#include <string_view>

namespace cosmosim::core {
struct IntegratorState;
struct SimulationState;
struct StepContext;
class HierarchicalTimeBinScheduler;
}

namespace cosmosim::workflows {

class AnalysisRuntime;
class GravityRuntime;
class HydroAmrRuntime;
class SourceRuntime;
class TimeCoordinator;
namespace internal {
class DriftRuntime;
class OutputRestartRuntime;
}

// A transient stage view is derived from these runtime authorities. It is
// never restart truth: owners rebuild views after every invalidating event.
struct RuntimeResourceEpoch {
  std::uint64_t particle_index_generation = 0;
  std::uint64_t cell_index_generation = 0;
  std::uint64_t gas_identity_generation = 0;
  std::uint64_t scheduler_tick = 0;
  std::uint64_t step_index = 0;
};

enum class RuntimeEpochField : std::uint32_t {
  kNone = 0,
  kParticleIndex = 1U << 0U,
  kCellIndex = 1U << 1U,
  kGasIdentity = 1U << 2U,
  kSchedulerTick = 1U << 3U,
  kStepIndex = 1U << 4U,
};

[[nodiscard]] constexpr RuntimeEpochField operator|(
    RuntimeEpochField lhs,
    RuntimeEpochField rhs) noexcept {
  return static_cast<RuntimeEpochField>(
      static_cast<std::uint32_t>(lhs) | static_cast<std::uint32_t>(rhs));
}

[[nodiscard]] constexpr RuntimeEpochField operator&(
    RuntimeEpochField lhs,
    RuntimeEpochField rhs) noexcept {
  return static_cast<RuntimeEpochField>(
      static_cast<std::uint32_t>(lhs) & static_cast<std::uint32_t>(rhs));
}

class RuntimeResourceEpochSource {
 public:
  virtual ~RuntimeResourceEpochSource() = default;
  [[nodiscard]] virtual RuntimeResourceEpoch currentRuntimeEpoch() const noexcept = 0;
};

// Adapter over existing runtime authorities. It borrows them without exposing
// any of their mutation APIs through the stage-view surface.
class SimulationRuntimeEpochSource final : public RuntimeResourceEpochSource {
 public:
  SimulationRuntimeEpochSource(
      const core::SimulationState& state,
      const core::HierarchicalTimeBinScheduler& scheduler,
      const core::IntegratorState& integrator_state) noexcept;

  [[nodiscard]] RuntimeResourceEpoch currentRuntimeEpoch() const noexcept override;

 private:
  const core::SimulationState* m_state = nullptr;
  const core::HierarchicalTimeBinScheduler* m_scheduler = nullptr;
  const core::IntegratorState* m_integrator_state = nullptr;
};

// Capability carried by every production stage view. The capability exposes
// freshness validation only; it cannot mutate an epoch or recover the state
// object from which the view was derived.
class RuntimeResourceLease {
 public:
  RuntimeResourceLease(
      const RuntimeResourceEpochSource& source,
      RuntimeEpochField guarded_fields) noexcept;

  [[nodiscard]] const RuntimeResourceEpoch& capturedEpoch() const noexcept;
  [[nodiscard]] RuntimeEpochField guardedFields() const noexcept;
  [[nodiscard]] bool isFresh() const noexcept;
  void requireFresh(std::string_view view_name) const;

 private:
  const RuntimeResourceEpochSource* m_source = nullptr;
  RuntimeResourceEpoch m_captured_epoch{};
  RuntimeEpochField m_guarded_fields = RuntimeEpochField::kNone;
};

// These are deliberately distinct capabilities. A task accepting one view
// cannot be invoked with another view and cannot obtain SimulationState.
class DriftParticleStageView {
 public:
  explicit DriftParticleStageView(RuntimeResourceLease lease) noexcept;
  void requireFresh() const;

 private:
  friend class internal::DriftRuntime;
  friend class TimeCoordinator;
  DriftParticleStageView(
      RuntimeResourceLease lease,
      core::StepContext& context) noexcept;
  [[nodiscard]] core::StepContext& ownerContext() const;

  RuntimeResourceLease m_lease;
  core::StepContext* m_context = nullptr;
};

class GravityStageView {
 public:
  explicit GravityStageView(RuntimeResourceLease lease) noexcept;
  void requireFresh() const;

 private:
  friend class GravityRuntime;
  friend class TimeCoordinator;
  GravityStageView(
      RuntimeResourceLease lease,
      core::StepContext& context) noexcept;
  [[nodiscard]] core::StepContext& ownerContext() const;

  RuntimeResourceLease m_lease;
  core::StepContext* m_context = nullptr;
};

class HydroAmrStageView {
 public:
  explicit HydroAmrStageView(RuntimeResourceLease lease) noexcept;
  void requireFresh() const;

 private:
  friend class HydroAmrRuntime;
  friend class TimeCoordinator;
  HydroAmrStageView(
      RuntimeResourceLease lease,
      core::StepContext& context) noexcept;
  [[nodiscard]] core::StepContext& ownerContext() const;

  RuntimeResourceLease m_lease;
  core::StepContext* m_context = nullptr;
};

class SourceMutationStageView {
 public:
  explicit SourceMutationStageView(RuntimeResourceLease lease) noexcept;
  void requireFresh() const;

 private:
  friend class SourceRuntime;
  friend class TimeCoordinator;
  SourceMutationStageView(
      RuntimeResourceLease lease,
      core::StepContext& context) noexcept;
  [[nodiscard]] core::StepContext& ownerContext() const;

  RuntimeResourceLease m_lease;
  core::StepContext* m_context = nullptr;
};

class AnalysisStageView {
 public:
  explicit AnalysisStageView(RuntimeResourceLease lease) noexcept;
  void requireFresh() const;

 private:
  friend class AnalysisRuntime;
  friend class TimeCoordinator;
  AnalysisStageView(
      RuntimeResourceLease lease,
      core::StepContext& context) noexcept;
  [[nodiscard]] core::StepContext& ownerContext() const;

  RuntimeResourceLease m_lease;
  core::StepContext* m_context = nullptr;
};

class OutputRestartStageView {
 public:
  explicit OutputRestartStageView(RuntimeResourceLease lease) noexcept;
  void requireFresh() const;

 private:
  friend class TimeCoordinator;
  friend class internal::OutputRestartRuntime;
  OutputRestartStageView(
      RuntimeResourceLease lease,
      core::StepContext& context) noexcept;
  [[nodiscard]] core::StepContext& ownerContext() const;

  RuntimeResourceLease m_lease;
  core::StepContext* m_context = nullptr;
};

class MigrationOwnershipView {
 public:
  explicit MigrationOwnershipView(RuntimeResourceLease lease) noexcept;
  void requireFresh() const;

 private:
  RuntimeResourceLease m_lease;
};

class TimeCriteriaStageView {
 public:
  explicit TimeCriteriaStageView(RuntimeResourceLease lease) noexcept;
  void requireFresh() const;

 private:
  RuntimeResourceLease m_lease;
};

}  // namespace cosmosim::workflows
