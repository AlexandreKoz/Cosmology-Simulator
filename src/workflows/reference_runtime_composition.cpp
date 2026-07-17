#include "workflows/internal/reference_runtime_composition.hpp"

#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/workflows/analysis_runtime.hpp"
#include "cosmosim/workflows/gravity_runtime.hpp"
#include "cosmosim/workflows/hydro_amr_runtime.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "cosmosim/workflows/source_runtime.hpp"
#include "cosmosim/workflows/output_restart_runtime.hpp"

namespace cosmosim::workflows::internal {

struct CompositionAssembly {
  std::shared_ptr<GravityRuntime> gravity;
  std::shared_ptr<HydroAmrRuntime> hydro_amr;
};

class DriftRuntime final {
 public:
  explicit DriftRuntime(const RuntimeServices& services) noexcept
      : m_services(services) {}

  void execute(DriftParticleStageView& view) const {
    view.requireFresh();
    core::StepContext& context = view.ownerContext();
    if (context.stage != core::IntegrationStage::kDrift) {
      throw std::logic_error("drift task received a non-drift stage");
    }
    const double drift_factor = context.timeline_step.drift_factor_code;
    const std::uint32_t world_rank =
        static_cast<std::uint32_t>(m_services.mpi_context.worldRank());
    for (const std::uint32_t particle_index : context.active_set.particle_indices) {
      if (particle_index >= context.state.particles.size() ||
          particle_index >= context.state.particle_sidecar.owning_rank.size()) {
        throw std::out_of_range("drift task active particle index out of range");
      }
      if (context.state.particle_sidecar.owning_rank[particle_index] != world_rank) {
        continue;
      }
      context.state.particles.position_x_comoving[particle_index] +=
          context.state.particles.velocity_x_peculiar[particle_index] * drift_factor;
      context.state.particles.position_y_comoving[particle_index] +=
          context.state.particles.velocity_y_peculiar[particle_index] * drift_factor;
      context.state.particles.position_z_comoving[particle_index] +=
          context.state.particles.velocity_z_peculiar[particle_index] * drift_factor;
    }
    context.state.requireGasCellIdentityMapCoversDenseRows("drift task");
  }

 private:
  const RuntimeServices& m_services;
};

namespace {

[[nodiscard]] RuntimeResourceAccess read(RuntimeResourceKey resource) {
  return {.resource = resource, .mode = RuntimeResourceAccessMode::kRead};
}

[[nodiscard]] RuntimeResourceAccess write(RuntimeResourceKey resource) {
  return {.resource = resource, .mode = RuntimeResourceAccessMode::kWrite};
}

[[nodiscard]] RuntimeResourceAccess readWrite(RuntimeResourceKey resource) {
  return {.resource = resource, .mode = RuntimeResourceAccessMode::kReadWrite};
}

[[nodiscard]] RuntimeModuleDescriptor makeAnalysisDescriptor(
    const ReferenceRuntimeCompositionInputs& inputs) {
  static constexpr std::array stages{
      core::IntegrationStage::kGravityKickPre,
      core::IntegrationStage::kDrift,
      core::IntegrationStage::kForceRefresh,
      core::IntegrationStage::kHydroUpdate,
      core::IntegrationStage::kSourceTerms,
      core::IntegrationStage::kGravityKickPost,
      core::IntegrationStage::kAnalysisHooks,
      core::IntegrationStage::kOutputCheck,
  };
  std::vector<RuntimeTaskDeclaration> declarations;
  declarations.reserve(stages.size() + 1U);
  for (const core::IntegrationStage stage : stages) {
    declarations.push_back(RuntimeTaskDeclaration{
        .task_id = "analysis.stage_audit." + std::string(core::integrationStageName(stage)),
        .stage = stage,
        .ordinal = -100,
        .view_kind = RuntimeStageViewKind::kStageAudit,
        .resources = {write(RuntimeResourceKey::kDiagnostics)},
    });
  }
  declarations.push_back(RuntimeTaskDeclaration{
      .task_id = "analysis.diagnostics",
      .stage = core::IntegrationStage::kAnalysisHooks,
      .ordinal = 100,
      .view_kind = RuntimeStageViewKind::kAnalysis,
      .resources = {read(RuntimeResourceKey::kParticlePosition),
                    read(RuntimeResourceKey::kHydroPrimitiveState),
                    write(RuntimeResourceKey::kDiagnostics)},
  });

  return RuntimeModuleDescriptor{
      .module_id = "analysis",
      .schema_version = 1,
      .construction_ordinal = 10,
      .prerequisites = {},
      .incompatibilities = {},
      .stage_tasks = std::move(declarations),
      .factory = [&config = inputs.config, &report = inputs.report](
                     const RuntimeModuleFactoryContext&) {
        std::shared_ptr<AnalysisRuntime> owner(
            makeAnalysisRuntime(config, report.stage_sequence));
        RuntimeModuleInstance instance;
        instance.owner_lifetime = owner;
        static constexpr std::array factory_stages{
            core::IntegrationStage::kGravityKickPre,
            core::IntegrationStage::kDrift,
            core::IntegrationStage::kForceRefresh,
            core::IntegrationStage::kHydroUpdate,
            core::IntegrationStage::kSourceTerms,
            core::IntegrationStage::kGravityKickPost,
            core::IntegrationStage::kAnalysisHooks,
            core::IntegrationStage::kOutputCheck,
        };
        for (const core::IntegrationStage stage : factory_stages) {
          instance.stage_tasks.push_back(RuntimeStageTaskContribution{
              .task_id = "analysis.stage_audit." +
                  std::string(core::integrationStageName(stage)),
              .task = StageAuditTask{[owner](AnalysisStageView& view) {
                owner->audit(view);
              }},
          });
        }
        instance.stage_tasks.push_back(RuntimeStageTaskContribution{
            .task_id = "analysis.diagnostics",
            .task = AnalysisStageTask([owner](AnalysisStageView& view) {
              owner->executeDiagnostics(view);
            }),
        });
        return instance;
      },
  };
}

[[nodiscard]] RuntimeModuleDescriptor makeDriftDescriptor() {
  return RuntimeModuleDescriptor{
      .module_id = "drift",
      .schema_version = 1,
      .construction_ordinal = 20,
      .prerequisites = {},
      .incompatibilities = {},
      .stage_tasks = {RuntimeTaskDeclaration{
          .task_id = "time.drift_particles",
          .stage = core::IntegrationStage::kDrift,
          .ordinal = 0,
          .view_kind = RuntimeStageViewKind::kDriftParticles,
          .resources = {readWrite(RuntimeResourceKey::kParticlePosition),
                        read(RuntimeResourceKey::kParticleVelocity),
                        read(RuntimeResourceKey::kMigrationOwnership)},
      }},
      .factory = [](const RuntimeModuleFactoryContext& context) {
        auto owner = std::make_shared<DriftRuntime>(context.services);
        RuntimeModuleInstance instance;
        instance.owner_lifetime = owner;
        instance.stage_tasks.push_back(RuntimeStageTaskContribution{
            .task_id = "time.drift_particles",
            .task = DriftStageTask([owner](DriftParticleStageView& view) {
              owner->execute(view);
            }),
        });
        return instance;
      },
  };
}

[[nodiscard]] RuntimeModuleDescriptor makeGravityDescriptor(
    const ReferenceRuntimeCompositionInputs& inputs,
    std::shared_ptr<CompositionAssembly> assembly) {
  std::vector<RuntimeTaskDeclaration> declarations;
  for (const core::IntegrationStage stage : {
           core::IntegrationStage::kGravityKickPre,
           core::IntegrationStage::kForceRefresh,
           core::IntegrationStage::kGravityKickPost}) {
    declarations.push_back(RuntimeTaskDeclaration{
        .task_id = "gravity." + std::string(core::integrationStageName(stage)),
        .stage = stage,
        .ordinal = 0,
        .view_kind = RuntimeStageViewKind::kGravity,
        .resources = {read(RuntimeResourceKey::kParticleGravitySource),
                      write(RuntimeResourceKey::kGravityAcceleration),
                      readWrite(RuntimeResourceKey::kIntegratorTruth)},
    });
  }
  return RuntimeModuleDescriptor{
      .module_id = "gravity",
      .schema_version = 1,
      .construction_ordinal = 30,
      .prerequisites = {},
      .incompatibilities = {},
      .stage_tasks = std::move(declarations),
      .factory = [inputs, assembly](const RuntimeModuleFactoryContext&) {
        assembly->gravity = std::shared_ptr<GravityRuntime>(makeGravityRuntime(
            inputs.config,
            inputs.mode_policy,
            inputs.services,
            inputs.zoom_region_path));
        RuntimeModuleInstance instance;
        instance.owner_lifetime = assembly->gravity;
        for (const core::IntegrationStage stage : {
                 core::IntegrationStage::kGravityKickPre,
                 core::IntegrationStage::kForceRefresh,
                 core::IntegrationStage::kGravityKickPost}) {
          instance.stage_tasks.push_back(RuntimeStageTaskContribution{
              .task_id = "gravity." + std::string(core::integrationStageName(stage)),
              .task = GravityStageTask([owner = assembly->gravity](GravityStageView& view) {
                owner->execute(view);
              }),
          });
        }
        return instance;
      },
  };
}

[[nodiscard]] RuntimeModuleDescriptor makeHydroAmrDescriptor(
    const ReferenceRuntimeCompositionInputs& inputs,
    std::shared_ptr<CompositionAssembly> assembly) {
  return RuntimeModuleDescriptor{
      .module_id = "hydro_amr",
      .schema_version = 1,
      .construction_ordinal = 40,
      .prerequisites = {"gravity"},
      .incompatibilities = {},
      .stage_tasks = {RuntimeTaskDeclaration{
          .task_id = "hydro_amr.update",
          .stage = core::IntegrationStage::kHydroUpdate,
          .ordinal = 0,
          .view_kind = RuntimeStageViewKind::kHydroAmr,
          .resources = {readWrite(RuntimeResourceKey::kHydroConservedState),
                        write(RuntimeResourceKey::kHydroPrimitiveState),
                        readWrite(RuntimeResourceKey::kAmrPatchState),
                        read(RuntimeResourceKey::kGravityAcceleration)},
      }},
      .factory = [inputs, assembly](const RuntimeModuleFactoryContext&) {
        if (!assembly->gravity) {
          throw std::logic_error("hydro/AMR factory requires the gravity owner");
        }
        assembly->hydro_amr = std::shared_ptr<HydroAmrRuntime>(makeHydroAmrRuntime(
            inputs.config,
            inputs.mode_policy,
            *assembly->gravity,
            inputs.services));
        RuntimeModuleInstance instance;
        instance.owner_lifetime = assembly->hydro_amr;
        instance.stage_tasks.push_back(RuntimeStageTaskContribution{
            .task_id = "hydro_amr.update",
            .task = HydroAmrStageTask([owner = assembly->hydro_amr](HydroAmrStageView& view) {
              owner->execute(view);
            }),
        });
        return instance;
      },
  };
}

[[nodiscard]] RuntimeModuleDescriptor makeSourceDescriptor(
    const ReferenceRuntimeCompositionInputs& inputs) {
  return RuntimeModuleDescriptor{
      .module_id = "sources",
      .schema_version = 1,
      .construction_ordinal = 50,
      .prerequisites = {},
      .incompatibilities = {},
      .stage_tasks = {RuntimeTaskDeclaration{
          .task_id = "sources.update",
          .stage = core::IntegrationStage::kSourceTerms,
          .ordinal = 0,
          .view_kind = RuntimeStageViewKind::kSourceMutation,
          .resources = {readWrite(RuntimeResourceKey::kSourceMutationState),
                        readWrite(RuntimeResourceKey::kParticlePosition),
                        readWrite(RuntimeResourceKey::kHydroConservedState)},
      }},
      .factory = [&config = inputs.config, &units = inputs.units,
                  world_rank = inputs.world_rank](const RuntimeModuleFactoryContext&) {
        std::shared_ptr<SourceRuntime> owner(
            makeSourceRuntime(config, units, world_rank));
        RuntimeModuleInstance instance;
        instance.owner_lifetime = owner;
        instance.stage_tasks.push_back(RuntimeStageTaskContribution{
            .task_id = "sources.update",
            .task = SourceMutationStageTask([owner](SourceMutationStageView& view) {
              owner->execute(view);
            }),
        });
        return instance;
      },
  };
}

[[nodiscard]] RuntimeModuleDescriptor makeOutputRestartDescriptor(
    const ReferenceRuntimeCompositionInputs& inputs,
    std::shared_ptr<CompositionAssembly> assembly) {
  return RuntimeModuleDescriptor{
      .module_id = "output_restart",
      .schema_version = 1,
      .construction_ordinal = 60,
      .prerequisites = {"gravity"},
      .incompatibilities = {},
      .stage_tasks = {RuntimeTaskDeclaration{
          .task_id = "output_restart.boundary",
          .stage = core::IntegrationStage::kOutputCheck,
          .ordinal = 0,
          .view_kind = RuntimeStageViewKind::kOutputRestart,
          .resources = {read(RuntimeResourceKey::kParticlePosition),
                        read(RuntimeResourceKey::kHydroConservedState),
                        read(RuntimeResourceKey::kMigrationOwnership),
                        read(RuntimeResourceKey::kSchedulerTruth),
                        read(RuntimeResourceKey::kIntegratorTruth),
                        readWrite(RuntimeResourceKey::kOutputRestartState),
                        write(RuntimeResourceKey::kDiagnostics)},
      }},
      .factory = [inputs, assembly](const RuntimeModuleFactoryContext&) {
        if (!assembly->gravity) {
          throw std::logic_error("output/restart factory requires the gravity owner");
        }
        auto owner = std::make_shared<OutputRestartRuntime>(
            inputs.frozen_config,
            inputs.config,
            inputs.particle_scheduler,
            inputs.gas_cell_scheduler,
            *assembly->gravity,
            inputs.report,
            inputs.profiler,
            inputs.pending_output,
            inputs.options.write_outputs);
        RuntimeModuleInstance instance;
        instance.owner_lifetime = owner;
        instance.stage_tasks.push_back(RuntimeStageTaskContribution{
            .task_id = "output_restart.boundary",
            .task = OutputRestartStageTask([owner](OutputRestartStageView& view) {
              owner->execute(view);
            }),
        });
        return instance;
      },
  };
}

}  // namespace

ReferenceRuntimeComposition buildReferenceRuntimeComposition(
    ReferenceRuntimeCompositionInputs inputs) {
  auto assembly = std::make_shared<CompositionAssembly>();
  RuntimeModuleRegistry registry;
  registry.registerModule(makeAnalysisDescriptor(inputs));
  registry.registerModule(makeDriftDescriptor());
  registry.registerModule(makeGravityDescriptor(inputs, assembly));
  registry.registerModule(makeHydroAmrDescriptor(inputs, assembly));
  registry.registerModule(makeSourceDescriptor(inputs));
  registry.registerModule(makeOutputRestartDescriptor(inputs, assembly));
  if (inputs.options.register_runtime_modules) {
    inputs.options.register_runtime_modules(registry);
  }
  RuntimeExecutionPlan execution_plan = registry.freezeAndInstantiate(
      RuntimeModuleFactoryContext{inputs.services});
  if (!assembly->gravity || !assembly->hydro_amr) {
    throw std::logic_error("reference runtime composition omitted a required owner");
  }
  return ReferenceRuntimeComposition{
      .gravity = std::move(assembly->gravity),
      .hydro_amr = std::move(assembly->hydro_amr),
      .execution_plan = std::move(execution_plan),
  };
}

}  // namespace cosmosim::workflows::internal
