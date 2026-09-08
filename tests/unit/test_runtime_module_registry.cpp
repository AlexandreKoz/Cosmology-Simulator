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
                  const auto resources = view.declaredResources();
                  assert(resources.size() == 1U);
                  assert(resources.front().resource == RuntimeResourceKey::kDiagnostics);
                  assert(resources.front().mode == RuntimeResourceAccessMode::kWrite);
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


void testTaskDependenciesAndConcurrencyGuards(
    const cosmosim::workflows::RuntimeModuleFactoryContext& context,
    std::vector<std::string>* trace) {
  using namespace cosmosim::workflows;

  RuntimeModuleRegistry dependency_registry;
  auto producer = makeAnalysisModule("producer", 0, 10, {}, trace);
  auto consumer = makeAnalysisModule("consumer", 1, 20, {}, trace);
  consumer.stage_tasks.front().dependencies = {"producer::producer.analysis"};
  dependency_registry.registerModule(std::move(producer));
  dependency_registry.registerModule(std::move(consumer));
  auto dependency_plan = dependency_registry.freezeAndInstantiate(context);
  assert(dependency_plan.taskCount() == 2U);

  RuntimeModuleRegistry invalid_order_registry;
  auto early = makeAnalysisModule("early", 0, 10, {}, trace);
  auto late = makeAnalysisModule("late_dependency", 1, 20, {}, trace);
  early.stage_tasks.front().dependencies = {"late_dependency::late_dependency.analysis"};
  invalid_order_registry.registerModule(std::move(early));
  invalid_order_registry.registerModule(std::move(late));
  expectInvalid([&]() { (void)invalid_order_registry.freezeAndInstantiate(context); });

  RuntimeTaskDeclaration bandwidth_a;
  bandwidth_a.task_id = "bandwidth_a";
  bandwidth_a.resources = {{RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kRead}};
  bandwidth_a.scheduling.estimated_peak_bytes = 256U;
  bandwidth_a.scheduling.peak_is_known = true;
  bandwidth_a.scheduling.memory_bandwidth_pressure = RuntimeTaskPressureClass::kHigh;

  RuntimeTaskDeclaration bandwidth_b;
  bandwidth_b.task_id = "bandwidth_b";
  bandwidth_b.resources = {{RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kRead}};
  bandwidth_b.scheduling.estimated_peak_bytes = 256U;
  bandwidth_b.scheduling.peak_is_known = true;
  bandwidth_b.scheduling.memory_bandwidth_pressure = RuntimeTaskPressureClass::kHigh;

  cosmosim::core::MemoryGovernorSnapshot roomy;
  roomy.headroom_bytes = 1024U;
  assert(!runtimeTasksMayOverlap(bandwidth_a, bandwidth_b, roomy));

  bandwidth_a.scheduling.memory_bandwidth_pressure = RuntimeTaskPressureClass::kModerate;
  bandwidth_b.scheduling.memory_bandwidth_pressure = RuntimeTaskPressureClass::kModerate;
  assert(runtimeTasksMayOverlap(bandwidth_a, bandwidth_b, roomy));

  cosmosim::core::MemoryGovernorSnapshot tight = roomy;
  tight.headroom_bytes = 400U;
  assert(!runtimeTasksMayOverlap(bandwidth_a, bandwidth_b, tight));

  bandwidth_b.resources = {{RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kWrite}};
  assert(!runtimeTasksMayOverlap(bandwidth_a, bandwidth_b, roomy));
  bandwidth_b.resources = {{RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kRead}};
  bandwidth_b.scheduling.peak_is_known = false;
  assert(!runtimeTasksMayOverlap(bandwidth_a, bandwidth_b, roomy));
  bandwidth_b.scheduling.peak_is_known = true;
  bandwidth_b.dependencies = {"producer::producer.analysis"};
  assert(!runtimeTasksMayOverlap(bandwidth_a, bandwidth_b, roomy));
}


void testDispatcherMemoryAdmission() {
  using namespace cosmosim::workflows;
  cosmosim::parallel::MpiContext mpi_context(false, 1, 0);
  cosmosim::core::ProfilerSession profiler(true);
  cosmosim::core::MemoryGovernor governor({.hard_limit_bytes = 256U});
  RuntimeServices services{.mpi_context = mpi_context, .profiler = profiler,
                           .memory_governor = &governor};
  RuntimeModuleFactoryContext context{services};
  std::vector<std::string> trace;
  auto makeModule = [&](std::string name, std::uint64_t bytes, bool owner_managed) {
    auto module = makeAnalysisModule(name, 0, 0, {}, &trace);
    module.stage_tasks.front().scheduling.peak_is_known = true;
    module.stage_tasks.front().scheduling.estimated_peak_bytes = bytes;
    module.stage_tasks.front().scheduling.memory_ownership = owner_managed
        ? RuntimeTaskMemoryOwnership::kOwnerManaged
        : RuntimeTaskMemoryOwnership::kDispatcherOwned;
    module.factory = [&, name, bytes, owner_managed](const RuntimeModuleFactoryContext& ctx) {
      RuntimeModuleInstance instance;
      instance.owner_lifetime = std::make_shared<std::string>(name);
      instance.stage_tasks.push_back(RuntimeStageTaskContribution{
          .task_id = name + ".analysis",
          .estimate_incremental_bytes = [bytes] { return bytes; },
          .task = AnalysisStageTask([&, bytes, owner_managed](AnalysisStageView& view) {
            view.requireFresh();
            const auto snapshot = governor.snapshot();
            assert(snapshot.committed_bytes == (owner_managed ? 0U : bytes));
            if (owner_managed) {
              auto physical = governor.reserve(cosmosim::core::MemoryClass::kDiagnostic,
                                               bytes, "test.owner");
              physical.commit();
              assert(governor.snapshot().committed_bytes == bytes);
            }
            trace.push_back("ran");
          }),
      });
      return instance;
    };
    return module;
  };
  StableEpochSource epoch_source;
  AnalysisStageView view(RuntimeResourceLease(epoch_source, RuntimeEpochField::kStepIndex));
  const auto execute = [&](std::uint64_t bytes, bool owner_managed) {
    RuntimeModuleRegistry registry;
    registry.registerModule(makeModule("admission", bytes, owner_managed));
    auto plan = registry.freezeAndInstantiate(context);
    plan.executeStage(cosmosim::core::IntegrationStage::kAnalysisHooks, view);
  };
  execute(128U, false);
  assert(governor.snapshot().committed_bytes == 0U);
  execute(128U, true);
  assert(governor.snapshot().committed_bytes == 0U);
  bool rejected = false;
  try { execute(300U, false); } catch (const std::exception&) { rejected = true; }
  assert(rejected && governor.snapshot().reserved_bytes == 0U);
  assert(governor.snapshot().committed_bytes == 0U);
}


void testDynamicOwnerMemoryContracts() {
  using namespace cosmosim::workflows;
  using cosmosim::core::MemoryClass;
  using cosmosim::core::MemoryGovernor;
  using cosmosim::core::MemoryGovernorSnapshot;
  assert(incrementalMemoryBeyondRetained(100U, 40U) == 60U);
  assert(incrementalMemoryBeyondRetained(40U, 100U) == 0U);
  assert(incrementalMemoryBeyondRetained(UINT64_MAX, 0U) == UINT64_MAX);

  cosmosim::parallel::MpiContext mpi_context(false, 1, 0);
  cosmosim::core::ProfilerSession profiler(true);
  MemoryGovernor governor({.hard_limit_bytes = 256U});
  RuntimeServices services{.mpi_context = mpi_context, .profiler = profiler,
                           .memory_governor = &governor};
  RuntimeModuleFactoryContext context{services};
  StableEpochSource epoch_source;
  AnalysisStageView view(RuntimeResourceLease(epoch_source, RuntimeEpochField::kStepIndex));
  std::vector<std::string> trace;
  std::uint64_t requested_bytes = 128U;
  bool complete = true;
  auto makeDynamicModule = [&](bool dispatcher_owned) {
    auto module = makeAnalysisModule("dynamic", 0, 0, {}, &trace);
    module.stage_tasks.front().scheduling.memory_ownership = dispatcher_owned
        ? RuntimeTaskMemoryOwnership::kDispatcherOwned
        : RuntimeTaskMemoryOwnership::kOwnerManaged;
    module.stage_tasks.front().scheduling.peak_is_known = dispatcher_owned;
    module.stage_tasks.front().scheduling.estimated_peak_bytes = dispatcher_owned ? 128U : 0U;
    module.factory = [&, dispatcher_owned](const RuntimeModuleFactoryContext&) {
      RuntimeModuleInstance instance;
      instance.owner_lifetime = std::make_shared<int>(1);
      instance.stage_tasks.push_back(RuntimeStageTaskContribution{
          .task_id = "dynamic.analysis",
          .estimate_memory = [&]() {
            return RuntimeTaskMemoryEstimate{requested_bytes, complete,
                complete ? std::string_view{} : std::string_view{"partial"}};
          },
          .task = AnalysisStageTask([&, dispatcher_owned](AnalysisStageView& stage_view) {
            stage_view.requireFresh();
            if (dispatcher_owned) {
              assert(governor.snapshot().committed_bytes == requested_bytes);
            } else {
              assert(governor.snapshot().committed_bytes == 0U);
              auto physical = governor.reserve(MemoryClass::kDiagnostic, 16U, "test.physical");
              physical.commit();
            }
            trace.push_back("dynamic");
          }),
      });
      return instance;
    };
    return module;
  };
  auto makePlan = [&](bool dispatcher_owned) {
    RuntimeModuleRegistry registry;
    registry.registerModule(makeDynamicModule(dispatcher_owned));
    return registry.freezeAndInstantiate(context);
  };
  auto owner_plan = makePlan(false);
  assert(owner_plan.taskMemoryEstimate("dynamic::dynamic.analysis").complete);
  requested_bytes = 300U;
  assert(owner_plan.taskMemoryEstimate("dynamic::dynamic.analysis").incremental_bytes == 300U);
  bool rejected = false;
  try { owner_plan.executeStage(cosmosim::core::IntegrationStage::kAnalysisHooks, view); }
  catch (const cosmosim::core::MemoryAdmissionError&) { rejected = true; }
  assert(rejected && trace.empty());
  complete = false;
  // Partial estimates cannot reject a legal smaller physical allocation.
  owner_plan.executeStage(cosmosim::core::IntegrationStage::kAnalysisHooks, view);
  assert(trace.size() == 1U);
  assert(governor.snapshot().committed_bytes == 0U);
  assert(governor.snapshot().reserved_bytes == 0U);
  auto dispatcher_plan = makePlan(true);
  rejected = false;
  try { dispatcher_plan.executeStage(cosmosim::core::IntegrationStage::kAnalysisHooks, view); }
  catch (const std::logic_error&) { rejected = true; }
  assert(rejected && trace.size() == 1U);
  complete = true;
  requested_bytes = 128U;
  dispatcher_plan.executeStage(cosmosim::core::IntegrationStage::kAnalysisHooks, view);
  assert(trace.size() == 2U);
  assert(governor.snapshot().committed_bytes == 0U);
  assert(governor.snapshot().reserved_bytes == 0U);
  rejected = false;
  try { (void)owner_plan.taskMemoryEstimate("dynamic.analysis"); }
  catch (const std::invalid_argument&) { rejected = true; }
  assert(rejected);
  rejected = false;
  try { (void)owner_plan.taskMemoryEstimate("missing::dynamic.analysis"); }
  catch (const std::out_of_range&) { rejected = true; }
  assert(rejected);
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
  testTaskDependenciesAndConcurrencyGuards(context, &trace);
  testDispatcherMemoryAdmission();
  testDynamicOwnerMemoryContracts();
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
  assert(analysis_view.declaredResources().empty());

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


  cosmosim::workflows::RuntimeModuleRegistry forbidden_analysis_write;
  auto forbidden = makeAnalysisModule("forbidden", 0, 0, {}, &trace);
  forbidden.stage_tasks.front().resources = {{
      .resource = cosmosim::workflows::RuntimeResourceKey::kParticlePosition,
      .mode = cosmosim::workflows::RuntimeResourceAccessMode::kWrite,
  }};
  expectInvalid([&]() {
    forbidden_analysis_write.registerModule(std::move(forbidden));
  });

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
