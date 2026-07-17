#include "cosmosim/workflows/runtime_module_registry.hpp"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/profiling.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/workflows/runtime_services.hpp"

namespace {

class StableEpochSource final
    : public cosmosim::workflows::RuntimeResourceEpochSource {
 public:
  [[nodiscard]] cosmosim::workflows::RuntimeResourceEpoch
  currentRuntimeEpoch() const noexcept override {
    return {};
  }
};

cosmosim::workflows::RuntimeModuleDescriptor makeAnalysisModule(
    std::string module_id,
    std::int32_t construction_ordinal,
    std::int32_t task_ordinal,
    std::vector<std::string> prerequisites,
    std::vector<std::string>* trace) {
  using namespace cosmosim::workflows;
  const std::string task_id = module_id + ".analysis";
  return RuntimeModuleDescriptor{
      .module_id = module_id,
      .schema_version = 1,
      .construction_ordinal = construction_ordinal,
      .prerequisites = std::move(prerequisites),
      .incompatibilities = {},
      .stage_tasks = {RuntimeTaskDeclaration{
          .task_id = task_id,
          .stage = cosmosim::core::IntegrationStage::kAnalysisHooks,
          .ordinal = task_ordinal,
          .view_kind = RuntimeStageViewKind::kAnalysis,
          .resources = {{
              .resource = RuntimeResourceKey::kDiagnostics,
              .mode = RuntimeResourceAccessMode::kWrite,
          }},
      }},
      .factory = [module_id = std::move(module_id), task_id, trace](
                     const RuntimeModuleFactoryContext&) {
        auto owner = std::make_shared<std::string>(module_id);
        RuntimeModuleInstance instance;
        instance.owner_lifetime = owner;
        instance.stage_tasks.push_back(RuntimeStageTaskContribution{
            .task_id = task_id,
            .task = AnalysisStageTask(
                [owner = std::move(owner), trace](AnalysisStageView& view) {
                  view.requireFresh();
                  trace->push_back(*owner);
                }),
        });
        return instance;
      },
  };
}

void expectInvalid(auto&& callback) {
  bool threw = false;
  try {
    callback();
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
}

}  // namespace

int main() {
  const cosmosim::parallel::MpiContext mpi_context(false, 1, 0);
  cosmosim::core::ProfilerSession profiler(true);
  cosmosim::workflows::RuntimeServices services{
      .mpi_context = mpi_context,
      .profiler = profiler,
      .deterministic_execution = true,
  };
  cosmosim::workflows::RuntimeModuleFactoryContext context{services};

  std::vector<std::string> trace;
  cosmosim::workflows::RuntimeModuleRegistry registry;
  registry.registerModule(makeAnalysisModule("base", 20, 20, {}, &trace));
  // The prerequisite fixes construction order while task ordinals independently
  // control deterministic execution order.
  registry.registerModule(
      makeAnalysisModule("probe", 10, 10, {"base"}, &trace));
  assert(registry.descriptorCount() == 2U);

  auto plan = registry.freezeAndInstantiate(context);
  assert(registry.frozen());
  assert(plan.moduleCount() == 2U);
  assert(plan.taskCount() == 2U);
  assert(plan.orderedModuleIds().size() == 2U);
  assert(plan.orderedModuleIds()[0] == "base");
  assert(plan.orderedModuleIds()[1] == "probe");

  StableEpochSource epoch_source;
  cosmosim::workflows::AnalysisStageView analysis_view(
      cosmosim::workflows::RuntimeResourceLease(
          epoch_source,
          cosmosim::workflows::RuntimeEpochField::kStepIndex));
  plan.executeStage(
      cosmosim::core::IntegrationStage::kAnalysisHooks,
      analysis_view);
  assert((trace == std::vector<std::string>{"probe", "base"}));

  bool rejected_after_freeze = false;
  try {
    registry.registerModule(makeAnalysisModule("late", 30, 30, {}, &trace));
  } catch (const std::logic_error&) {
    rejected_after_freeze = true;
  }
  assert(rejected_after_freeze);

  // Removing the probe descriptor removes its behavior; no execution-plan or
  // composition-root edit is necessary.
  trace.clear();
  cosmosim::workflows::RuntimeModuleRegistry without_probe;
  without_probe.registerModule(makeAnalysisModule("base", 20, 20, {}, &trace));
  auto plan_without_probe = without_probe.freezeAndInstantiate(context);
  plan_without_probe.executeStage(
      cosmosim::core::IntegrationStage::kAnalysisHooks,
      analysis_view);
  assert((trace == std::vector<std::string>{"base"}));

  cosmosim::workflows::RuntimeModuleRegistry missing_prerequisite;
  missing_prerequisite.registerModule(
      makeAnalysisModule("orphan", 0, 0, {"absent"}, &trace));
  expectInvalid([&]() {
    (void)missing_prerequisite.freezeAndInstantiate(context);
  });

  cosmosim::workflows::RuntimeModuleRegistry incompatible;
  auto first = makeAnalysisModule("first", 0, 0, {}, &trace);
  first.incompatibilities.push_back("second");
  incompatible.registerModule(std::move(first));
  incompatible.registerModule(makeAnalysisModule("second", 1, 1, {}, &trace));
  expectInvalid([&]() { (void)incompatible.freezeAndInstantiate(context); });

  cosmosim::workflows::RuntimeModuleRegistry bad_factory;
  auto mismatch = makeAnalysisModule("mismatch", 0, 0, {}, &trace);
  mismatch.factory = [](const cosmosim::workflows::RuntimeModuleFactoryContext&) {
    cosmosim::workflows::RuntimeModuleInstance instance;
    instance.stage_tasks.push_back({
        .task_id = "mismatch.analysis",
        .task = cosmosim::workflows::GravityStageTask(
            [](cosmosim::workflows::GravityStageView&) {}),
    });
    return instance;
  };
  bad_factory.registerModule(std::move(mismatch));
  expectInvalid([&]() { (void)bad_factory.freezeAndInstantiate(context); });
  return 0;
}
