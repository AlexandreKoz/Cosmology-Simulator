#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"

namespace cosmosim::workflows {

struct RuntimeServices;

enum class RuntimeStageViewKind : std::uint8_t {
  kStageAudit,
  kDriftParticles,
  kGravity,
  kHydroAmr,
  kSourceMutation,
  kAnalysis,
  kOutputRestart,
};

enum class RuntimeTaskPressureClass : std::uint8_t {
  kLow = 0,
  kModerate = 1,
  kHigh = 2,
};

enum class RuntimeTaskLifetimeBoundary : std::uint8_t {
  kTaskEnd = 0,
  kStageEnd = 1,
};

// Owner-managed tasks retain their existing phase/transaction reservations.
// Dispatcher-owned tasks are charged once by the execution plan. Neither mode
// grants permission to overlap tasks with unknown or incomplete peak models.
enum class RuntimeTaskMemoryOwnership : std::uint8_t {
  kOwnerManaged = 0,
  kDispatcherOwned = 1,
};

struct RuntimeTaskSchedulingProfile {
  std::uint64_t estimated_peak_bytes = 0U;
  // Zero is not an implicit proof of zero allocation. A complete incremental
  // peak model must opt in before this declaration can authorize overlap.
  bool peak_is_known = false;
  RuntimeTaskMemoryOwnership memory_ownership = RuntimeTaskMemoryOwnership::kOwnerManaged;
  core::MemoryClass memory_class = core::MemoryClass::kPhaseResident;
  RuntimeTaskLifetimeBoundary lifetime_boundary = RuntimeTaskLifetimeBoundary::kTaskEnd;
  RuntimeTaskPressureClass compute_pressure = RuntimeTaskPressureClass::kLow;
  RuntimeTaskPressureClass memory_bandwidth_pressure = RuntimeTaskPressureClass::kLow;
  RuntimeTaskPressureClass communication_pressure = RuntimeTaskPressureClass::kLow;
  bool optional = false;
};

struct RuntimeTaskDeclaration {
  std::string task_id;
  core::IntegrationStage stage = core::IntegrationStage::kGravityKickPre;
  std::int32_t ordinal = 0;
  RuntimeStageViewKind view_kind = RuntimeStageViewKind::kGravity;
  std::vector<RuntimeResourceAccess> resources;
  std::vector<std::string> dependencies;
  RuntimeTaskSchedulingProfile scheduling{};
};

using DriftStageTask = std::function<void(DriftParticleStageView&)>;
struct StageAuditTask {
  std::function<void(AnalysisStageView&)> execute;
  void operator()(AnalysisStageView& view) const { execute(view); }
};
using GravityStageTask = std::function<void(GravityStageView&)>;
using HydroAmrStageTask = std::function<void(HydroAmrStageView&)>;
using SourceMutationStageTask = std::function<void(SourceMutationStageView&)>;
using AnalysisStageTask = std::function<void(AnalysisStageView&)>;
using OutputRestartStageTask = std::function<void(OutputRestartStageView&)>;

using RuntimeStageTaskFunction = std::variant<
    StageAuditTask,
    DriftStageTask,
    GravityStageTask,
    HydroAmrStageTask,
    SourceMutationStageTask,
    AnalysisStageTask,
    OutputRestartStageTask>;

struct RuntimeStageTaskContribution {
  std::string task_id;
  // Optional bounded, state-dependent incremental peak. Evaluated once at
  // the task boundary; owner-managed tasks retain their existing reservations.
  std::function<std::uint64_t()> estimate_incremental_bytes;
  RuntimeStageTaskFunction task;
};

// The lifetime handle owns the concrete owner service captured by task
// callables. No caller can down-cast it through this interface.
struct RuntimeModuleInstance {
  std::shared_ptr<void> owner_lifetime;
  std::vector<RuntimeStageTaskContribution> stage_tasks;
};

struct RuntimeModuleFactoryContext {
  const RuntimeServices& services;
};

using RuntimeModuleFactory =
    std::function<RuntimeModuleInstance(const RuntimeModuleFactoryContext&)>;

struct RuntimeModuleDescriptor {
  std::string module_id;
  std::uint32_t schema_version = 1;
  std::int32_t construction_ordinal = 0;
  std::vector<std::string> prerequisites;
  std::vector<std::string> incompatibilities;
  std::vector<RuntimeTaskDeclaration> stage_tasks;
  RuntimeModuleFactory factory;
};

class RuntimeExecutionPlan {
 public:
  // Public only so the implementation can use generic typed dispatch without
  // granting tasks access to plan internals. Callers receive these only through
  // the read-only plan execution API.
  struct PlannedTask {
    std::string module_id;
    RuntimeTaskDeclaration declaration;
    RuntimeStageTaskFunction task;
    std::function<std::uint64_t()> estimate_incremental_bytes;
  };

  RuntimeExecutionPlan() = default;
  RuntimeExecutionPlan(RuntimeExecutionPlan&&) noexcept = default;
  RuntimeExecutionPlan& operator=(RuntimeExecutionPlan&&) noexcept = default;
  RuntimeExecutionPlan(const RuntimeExecutionPlan&) = delete;
  RuntimeExecutionPlan& operator=(const RuntimeExecutionPlan&) = delete;

  [[nodiscard]] std::size_t moduleCount() const noexcept;
  [[nodiscard]] std::size_t taskCount() const noexcept;
  [[nodiscard]] std::span<const std::string> orderedModuleIds() const noexcept;

  void executeAuditStage(core::IntegrationStage stage, AnalysisStageView& view) const;
  void executeStage(core::IntegrationStage stage, DriftParticleStageView& view) const;
  void executeStage(core::IntegrationStage stage, GravityStageView& view) const;
  void executeStage(core::IntegrationStage stage, HydroAmrStageView& view) const;
  void executeStage(core::IntegrationStage stage, SourceMutationStageView& view) const;
  void executeStage(core::IntegrationStage stage, AnalysisStageView& view) const;
  void executeStage(core::IntegrationStage stage, OutputRestartStageView& view) const;

 private:
  friend class RuntimeModuleRegistry;

  const RuntimeServices* m_services = nullptr;
  std::vector<std::string> m_ordered_module_ids;
  std::vector<RuntimeModuleInstance> m_module_instances;
  std::vector<PlannedTask> m_tasks;
};

class RuntimeModuleRegistry {
 public:
  void registerModule(RuntimeModuleDescriptor descriptor);
  [[nodiscard]] bool frozen() const noexcept;
  [[nodiscard]] std::size_t descriptorCount() const noexcept;
  [[nodiscard]] RuntimeExecutionPlan freezeAndInstantiate(
      const RuntimeModuleFactoryContext& context);

 private:
  std::vector<RuntimeModuleDescriptor> m_descriptors;
  bool m_frozen = false;
};

[[nodiscard]] RuntimeStageViewKind taskViewKind(
    const RuntimeStageTaskFunction& task) noexcept;
[[nodiscard]] std::string_view runtimeResourceKeyName(
    RuntimeResourceKey resource) noexcept;
[[nodiscard]] bool runtimeTasksMayOverlap(
    const RuntimeTaskDeclaration& lhs,
    const RuntimeTaskDeclaration& rhs,
    const core::MemoryGovernorSnapshot& memory_snapshot) noexcept;

}  // namespace cosmosim::workflows
