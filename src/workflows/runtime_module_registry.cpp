#include "cosmosim/workflows/runtime_module_registry.hpp"

#include <algorithm>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace cosmosim::workflows {

namespace internal {

class RuntimeTaskGrantBinder final {
 public:
  template <class View>
  static void bind(
      View& view,
      std::span<const RuntimeResourceAccess> resources) noexcept {
    view.bindDeclaredResources(resources);
  }

  template <class View>
  static void clear(View& view) noexcept {
    view.clearDeclaredResources();
  }
};

}  // namespace internal

namespace {

[[nodiscard]] bool viewSupportsStage(
    RuntimeStageViewKind view_kind,
    core::IntegrationStage stage) noexcept {
  switch (view_kind) {
    case RuntimeStageViewKind::kStageAudit:
      return true;
    case RuntimeStageViewKind::kDriftParticles:
      return stage == core::IntegrationStage::kDrift;
    case RuntimeStageViewKind::kGravity:
      return stage == core::IntegrationStage::kGravityKickPre ||
             stage == core::IntegrationStage::kForceRefresh ||
             stage == core::IntegrationStage::kGravityKickPost;
    case RuntimeStageViewKind::kHydroAmr:
      return stage == core::IntegrationStage::kHydroUpdate;
    case RuntimeStageViewKind::kSourceMutation:
      return stage == core::IntegrationStage::kSourceTerms;
    case RuntimeStageViewKind::kAnalysis:
      return stage == core::IntegrationStage::kAnalysisHooks;
    case RuntimeStageViewKind::kOutputRestart:
      return stage == core::IntegrationStage::kOutputCheck;
  }
  return false;
}

[[nodiscard]] std::optional<RuntimeResourceAccessMode> allowedModeForView(
    RuntimeStageViewKind view_kind,
    RuntimeResourceKey resource) noexcept {
  using Mode = RuntimeResourceAccessMode;
  using Key = RuntimeResourceKey;
  switch (view_kind) {
    case RuntimeStageViewKind::kStageAudit:
      return resource == Key::kDiagnostics
          ? std::optional<Mode>(Mode::kWrite)
          : std::nullopt;
    case RuntimeStageViewKind::kDriftParticles:
      switch (resource) {
        case Key::kParticlePosition: return Mode::kReadWrite;
        case Key::kParticleVelocity: return Mode::kRead;
        case Key::kMigrationOwnership: return Mode::kRead;
        default: return std::nullopt;
      }
    case RuntimeStageViewKind::kGravity:
      switch (resource) {
        case Key::kParticlePosition: return Mode::kRead;
        case Key::kParticleVelocity: return Mode::kReadWrite;
        case Key::kParticleGravitySource: return Mode::kRead;
        case Key::kHydroConservedState: return Mode::kRead;
        case Key::kMigrationOwnership: return Mode::kRead;
        case Key::kSchedulerTruth: return Mode::kRead;
        case Key::kGravityAcceleration: return Mode::kWrite;
        case Key::kIntegratorTruth: return Mode::kReadWrite;
        default: return std::nullopt;
      }
    case RuntimeStageViewKind::kHydroAmr:
      switch (resource) {
        case Key::kParticleVelocity: return Mode::kReadWrite;
        case Key::kHydroConservedState: return Mode::kReadWrite;
        case Key::kHydroPrimitiveState: return Mode::kReadWrite;
        case Key::kAmrPatchState: return Mode::kReadWrite;
        case Key::kGravityAcceleration: return Mode::kRead;
        case Key::kMigrationOwnership: return Mode::kRead;
        case Key::kIntegratorTruth: return Mode::kRead;
        default: return std::nullopt;
      }
    case RuntimeStageViewKind::kSourceMutation:
      switch (resource) {
        case Key::kSourceMutationState: return Mode::kReadWrite;
        case Key::kParticlePosition: return Mode::kReadWrite;
        case Key::kParticleVelocity: return Mode::kReadWrite;
        case Key::kHydroConservedState: return Mode::kReadWrite;
        case Key::kMigrationOwnership: return Mode::kReadWrite;
        case Key::kIntegratorTruth: return Mode::kRead;
        default: return std::nullopt;
      }
    case RuntimeStageViewKind::kAnalysis:
      switch (resource) {
        case Key::kParticlePosition:
        case Key::kParticleVelocity:
        case Key::kParticleGravitySource:
        case Key::kHydroPrimitiveState:
        case Key::kSourceMutationState:
        case Key::kMigrationOwnership:
        case Key::kIntegratorTruth:
          return Mode::kRead;
        case Key::kDiagnostics:
          return Mode::kWrite;
        default:
          return std::nullopt;
      }
    case RuntimeStageViewKind::kOutputRestart:
      switch (resource) {
        case Key::kParticlePosition:
        case Key::kParticleVelocity:
        case Key::kHydroConservedState:
        case Key::kSourceMutationState:
        case Key::kMigrationOwnership:
        case Key::kSchedulerTruth:
        case Key::kIntegratorTruth:
          return Mode::kRead;
        case Key::kOutputRestartState:
          return Mode::kReadWrite;
        case Key::kDiagnostics:
          return Mode::kWrite;
        default:
          return std::nullopt;
      }
  }
  return std::nullopt;
}

[[nodiscard]] bool viewAllowsAccess(
    RuntimeStageViewKind view_kind,
    RuntimeResourceAccess requested) noexcept {
  const auto allowed_mode = allowedModeForView(view_kind, requested.resource);
  return allowed_mode.has_value() && runtimeResourceAccessSatisfies(
      RuntimeResourceAccess{requested.resource, *allowed_mode}, requested);
}

void validateDescriptor(const RuntimeModuleDescriptor& descriptor) {
  if (descriptor.module_id.empty()) {
    throw std::invalid_argument("runtime module descriptor has an empty module_id");
  }
  if (!descriptor.factory) {
    throw std::invalid_argument(
        "runtime module '" + descriptor.module_id + "' has no factory");
  }
  std::set<std::string> task_ids;
  for (const RuntimeTaskDeclaration& task : descriptor.stage_tasks) {
    if (task.task_id.empty()) {
      throw std::invalid_argument(
          "runtime module '" + descriptor.module_id + "' has an empty task_id");
    }
    if (!task_ids.insert(task.task_id).second) {
      throw std::invalid_argument(
          "runtime module '" + descriptor.module_id +
          "' declares duplicate task '" + task.task_id + "'");
    }
    if (!viewSupportsStage(task.view_kind, task.stage)) {
      throw std::invalid_argument(
          "runtime task '" + task.task_id +
          "' declares a view that cannot execute its stage");
    }
    if (task.resources.empty()) {
      throw std::invalid_argument(
          "runtime task '" + task.task_id + "' declares no resource access");
    }
    std::set<RuntimeResourceKey> resources;
    for (const RuntimeResourceAccess access : task.resources) {
      if (!resources.insert(access.resource).second) {
        throw std::invalid_argument(
            "runtime task '" + task.task_id +
            "' declares duplicate access for resource '" +
            std::string(runtimeResourceKeyName(access.resource)) + "'");
      }
      if (!viewAllowsAccess(task.view_kind, access)) {
        throw std::invalid_argument(
            "runtime task '" + task.task_id +
            "' declares access outside its typed view for resource '" +
            std::string(runtimeResourceKeyName(access.resource)) + "'");
      }
    }
  }
}

[[nodiscard]] std::vector<std::size_t> resolveModuleOrder(
    const std::vector<RuntimeModuleDescriptor>& descriptors) {
  std::unordered_map<std::string, std::size_t> by_id;
  for (std::size_t index = 0; index < descriptors.size(); ++index) {
    if (!by_id.emplace(descriptors[index].module_id, index).second) {
      throw std::invalid_argument(
          "duplicate runtime module id '" + descriptors[index].module_id + "'");
    }
  }

  for (const RuntimeModuleDescriptor& descriptor : descriptors) {
    for (const std::string& incompatible : descriptor.incompatibilities) {
      if (by_id.contains(incompatible)) {
        throw std::invalid_argument(
            "runtime module '" + descriptor.module_id +
            "' is incompatible with registered module '" + incompatible + "'");
      }
    }
    for (const std::string& prerequisite : descriptor.prerequisites) {
      if (!by_id.contains(prerequisite)) {
        throw std::invalid_argument(
            "runtime module '" + descriptor.module_id +
            "' requires missing module '" + prerequisite + "'");
      }
      if (prerequisite == descriptor.module_id) {
        throw std::invalid_argument(
            "runtime module '" + descriptor.module_id + "' requires itself");
      }
    }
  }

  std::vector<std::size_t> ordered;
  std::set<std::size_t> emitted;
  while (ordered.size() != descriptors.size()) {
    std::vector<std::size_t> ready;
    for (std::size_t index = 0; index < descriptors.size(); ++index) {
      if (emitted.contains(index)) {
        continue;
      }
      const bool prerequisites_ready = std::all_of(
          descriptors[index].prerequisites.begin(),
          descriptors[index].prerequisites.end(),
          [&](const std::string& prerequisite) {
            return emitted.contains(by_id.at(prerequisite));
          });
      if (prerequisites_ready) {
        ready.push_back(index);
      }
    }
    if (ready.empty()) {
      throw std::invalid_argument("runtime module prerequisite graph contains a cycle");
    }
    std::sort(
        ready.begin(),
        ready.end(),
        [&](std::size_t lhs, std::size_t rhs) {
          return std::tuple{
                     descriptors[lhs].construction_ordinal,
                     descriptors[lhs].module_id} <
                 std::tuple{
                     descriptors[rhs].construction_ordinal,
                     descriptors[rhs].module_id};
        });
    for (const std::size_t index : ready) {
      emitted.insert(index);
      ordered.push_back(index);
    }
  }
  return ordered;
}

template <class View>
class TaskGrantScope final {
 public:
  TaskGrantScope(
      View& view,
      std::span<const RuntimeResourceAccess> resources) noexcept
      : m_view(view) {
    internal::RuntimeTaskGrantBinder::bind(m_view, resources);
  }

  ~TaskGrantScope() {
    internal::RuntimeTaskGrantBinder::clear(m_view);
  }

  TaskGrantScope(const TaskGrantScope&) = delete;
  TaskGrantScope& operator=(const TaskGrantScope&) = delete;

 private:
  View& m_view;
};

template <class View, class Task>
void executeTypedStage(
    std::span<const RuntimeExecutionPlan::PlannedTask> tasks,
    core::IntegrationStage stage,
    RuntimeStageViewKind expected_view,
    View& view) {
  view.requireFresh();
  for (const RuntimeExecutionPlan::PlannedTask& task : tasks) {
    if (task.declaration.stage != stage) {
      continue;
    }
    if (task.declaration.view_kind != expected_view) {
      continue;
    }
    if (!std::holds_alternative<Task>(task.task)) {
      throw std::logic_error(
          "frozen runtime task '" + task.declaration.task_id +
          "' has an inconsistent typed stage view");
    }
    TaskGrantScope<View> grant_scope(view, task.declaration.resources);
    std::get<Task>(task.task)(view);
  }
}

}  // namespace

RuntimeStageViewKind taskViewKind(const RuntimeStageTaskFunction& task) noexcept {
  return std::visit(
      [](const auto& typed_task) {
        using Task = std::decay_t<decltype(typed_task)>;
        if constexpr (std::is_same_v<Task, StageAuditTask>) {
          return RuntimeStageViewKind::kStageAudit;
        } else if constexpr (std::is_same_v<Task, DriftStageTask>) {
          return RuntimeStageViewKind::kDriftParticles;
        } else if constexpr (std::is_same_v<Task, GravityStageTask>) {
          return RuntimeStageViewKind::kGravity;
        } else if constexpr (std::is_same_v<Task, HydroAmrStageTask>) {
          return RuntimeStageViewKind::kHydroAmr;
        } else if constexpr (std::is_same_v<Task, SourceMutationStageTask>) {
          return RuntimeStageViewKind::kSourceMutation;
        } else if constexpr (std::is_same_v<Task, AnalysisStageTask>) {
          return RuntimeStageViewKind::kAnalysis;
        } else {
          return RuntimeStageViewKind::kOutputRestart;
        }
      },
      task);
}

std::string_view runtimeResourceKeyName(RuntimeResourceKey resource) noexcept {
  switch (resource) {
    case RuntimeResourceKey::kParticlePosition: return "particle_position";
    case RuntimeResourceKey::kParticleVelocity: return "particle_velocity";
    case RuntimeResourceKey::kParticleGravitySource: return "particle_gravity_source";
    case RuntimeResourceKey::kGravityAcceleration: return "gravity_acceleration";
    case RuntimeResourceKey::kHydroConservedState: return "hydro_conserved_state";
    case RuntimeResourceKey::kHydroPrimitiveState: return "hydro_primitive_state";
    case RuntimeResourceKey::kAmrPatchState: return "amr_patch_state";
    case RuntimeResourceKey::kSourceMutationState: return "source_mutation_state";
    case RuntimeResourceKey::kMigrationOwnership: return "migration_ownership";
    case RuntimeResourceKey::kSchedulerTruth: return "scheduler_truth";
    case RuntimeResourceKey::kIntegratorTruth: return "integrator_truth";
    case RuntimeResourceKey::kOutputRestartState: return "output_restart_state";
    case RuntimeResourceKey::kDiagnostics: return "diagnostics";
  }
  return "unknown";
}

std::size_t RuntimeExecutionPlan::moduleCount() const noexcept {
  return m_module_instances.size();
}

std::size_t RuntimeExecutionPlan::taskCount() const noexcept {
  return m_tasks.size();
}

std::span<const std::string> RuntimeExecutionPlan::orderedModuleIds() const noexcept {
  return m_ordered_module_ids;
}

void RuntimeExecutionPlan::executeAuditStage(
    core::IntegrationStage stage,
    AnalysisStageView& view) const {
  executeTypedStage<AnalysisStageView, StageAuditTask>(
      m_tasks, stage, RuntimeStageViewKind::kStageAudit, view);
}

void RuntimeExecutionPlan::executeStage(
    core::IntegrationStage stage,
    DriftParticleStageView& view) const {
  executeTypedStage<DriftParticleStageView, DriftStageTask>(
      m_tasks, stage, RuntimeStageViewKind::kDriftParticles, view);
}

void RuntimeExecutionPlan::executeStage(
    core::IntegrationStage stage,
    GravityStageView& view) const {
  executeTypedStage<GravityStageView, GravityStageTask>(
      m_tasks, stage, RuntimeStageViewKind::kGravity, view);
}

void RuntimeExecutionPlan::executeStage(
    core::IntegrationStage stage,
    HydroAmrStageView& view) const {
  executeTypedStage<HydroAmrStageView, HydroAmrStageTask>(
      m_tasks, stage, RuntimeStageViewKind::kHydroAmr, view);
}

void RuntimeExecutionPlan::executeStage(
    core::IntegrationStage stage,
    SourceMutationStageView& view) const {
  executeTypedStage<SourceMutationStageView, SourceMutationStageTask>(
      m_tasks, stage, RuntimeStageViewKind::kSourceMutation, view);
}

void RuntimeExecutionPlan::executeStage(
    core::IntegrationStage stage,
    AnalysisStageView& view) const {
  executeTypedStage<AnalysisStageView, AnalysisStageTask>(
      m_tasks, stage, RuntimeStageViewKind::kAnalysis, view);
}

void RuntimeExecutionPlan::executeStage(
    core::IntegrationStage stage,
    OutputRestartStageView& view) const {
  executeTypedStage<OutputRestartStageView, OutputRestartStageTask>(
      m_tasks, stage, RuntimeStageViewKind::kOutputRestart, view);
}

void RuntimeModuleRegistry::registerModule(RuntimeModuleDescriptor descriptor) {
  if (m_frozen) {
    throw std::logic_error("runtime module registry is already frozen");
  }
  validateDescriptor(descriptor);
  if (std::any_of(
          m_descriptors.begin(),
          m_descriptors.end(),
          [&](const RuntimeModuleDescriptor& existing) {
            return existing.module_id == descriptor.module_id;
          })) {
    throw std::invalid_argument(
        "duplicate runtime module id '" + descriptor.module_id + "'");
  }
  m_descriptors.push_back(std::move(descriptor));
}

bool RuntimeModuleRegistry::frozen() const noexcept { return m_frozen; }

std::size_t RuntimeModuleRegistry::descriptorCount() const noexcept {
  return m_descriptors.size();
}

RuntimeExecutionPlan RuntimeModuleRegistry::freezeAndInstantiate(
    const RuntimeModuleFactoryContext& context) {
  if (m_frozen) {
    throw std::logic_error("runtime module registry may be instantiated only once");
  }
  m_frozen = true;

  const std::vector<std::size_t> order = resolveModuleOrder(m_descriptors);
  RuntimeExecutionPlan plan;
  for (const std::size_t descriptor_index : order) {
    const RuntimeModuleDescriptor& descriptor = m_descriptors[descriptor_index];
    plan.m_ordered_module_ids.push_back(descriptor.module_id);
    RuntimeModuleInstance instance = descriptor.factory(context);

    std::unordered_map<std::string, RuntimeStageTaskContribution*> contribution_by_id;
    for (RuntimeStageTaskContribution& contribution : instance.stage_tasks) {
      if (contribution.task_id.empty() ||
          !contribution_by_id.emplace(contribution.task_id, &contribution).second) {
        throw std::invalid_argument(
            "runtime module factory '" + descriptor.module_id +
            "' returned an empty or duplicate task id");
      }
    }
    if (contribution_by_id.size() != descriptor.stage_tasks.size()) {
      throw std::invalid_argument(
          "runtime module factory '" + descriptor.module_id +
          "' did not return exactly its declared stage tasks");
    }
    for (const RuntimeTaskDeclaration& declaration : descriptor.stage_tasks) {
      const auto contribution_it = contribution_by_id.find(declaration.task_id);
      if (contribution_it == contribution_by_id.end()) {
        throw std::invalid_argument(
            "runtime module factory '" + descriptor.module_id +
            "' omitted declared task '" + declaration.task_id + "'");
      }
      if (taskViewKind(contribution_it->second->task) != declaration.view_kind) {
        throw std::invalid_argument(
            "runtime module factory task '" + declaration.task_id +
            "' does not match its declared typed view");
      }
      plan.m_tasks.push_back(RuntimeExecutionPlan::PlannedTask{
          .module_id = descriptor.module_id,
          .declaration = declaration,
          .task = std::move(contribution_it->second->task),
      });
    }
    plan.m_module_instances.push_back(std::move(instance));
  }

  std::stable_sort(
      plan.m_tasks.begin(),
      plan.m_tasks.end(),
      [](const RuntimeExecutionPlan::PlannedTask& lhs,
         const RuntimeExecutionPlan::PlannedTask& rhs) {
        return std::tuple{
                   static_cast<std::uint8_t>(lhs.declaration.stage),
                   lhs.declaration.ordinal,
                   lhs.module_id,
                   lhs.declaration.task_id} <
               std::tuple{
                   static_cast<std::uint8_t>(rhs.declaration.stage),
                   rhs.declaration.ordinal,
                   rhs.module_id,
                   rhs.declaration.task_id};
      });
  return plan;
}

}  // namespace cosmosim::workflows
