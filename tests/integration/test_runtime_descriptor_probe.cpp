#include <cassert>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <sstream>
#include <string>

#include "cosmosim/core/config.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"
#include "cosmosim/workflows/runtime_module_registry.hpp"

int main() {
  std::stringstream config_text;
  config_text << "schema_version = 1\n\n";
  config_text << "[mode]\n";
  config_text << "mode = zoom_in\n";
  config_text << "ic_file = generated\n";
  config_text << "zoom_high_res_region = false\n\n";
  config_text << "[numerics]\n";
  config_text << "time_begin_code = 0.01\n";
  config_text << "time_end_code = 0.0101\n";
  config_text << "max_global_steps = 1\n";
  config_text << "hierarchical_max_rung = 0\n";
  config_text << "treepm_pm_grid = 9\n";
  config_text << "treepm_rcut_cells = 3.9\n\n";
  config_text << "[output]\n";
  config_text << "run_name = runtime_descriptor_probe\n";
  config_text << "output_directory = integration_outputs\n";

  const cosmosim::core::FrozenConfig frozen =
      cosmosim::core::loadFrozenConfigFromString(
          config_text.str(), "test_runtime_descriptor_probe");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);

  std::uint64_t probe_invocations = 0;
  cosmosim::workflows::ReferenceWorkflowOptions options;
  options.write_outputs = false;
  options.max_steps_override = 1;
  options.register_runtime_modules = [&](cosmosim::workflows::RuntimeModuleRegistry& registry) {
    using namespace cosmosim::workflows;
    registry.registerModule(RuntimeModuleDescriptor{
        .module_id = "test_probe",
        .schema_version = 1,
        .construction_ordinal = 1000,
        .prerequisites = {"analysis"},
        .incompatibilities = {},
        .stage_tasks = {RuntimeTaskDeclaration{
            .task_id = "test_probe.analysis_diagnostic",
            .stage = cosmosim::core::IntegrationStage::kAnalysisHooks,
            .ordinal = 1000,
            .view_kind = RuntimeStageViewKind::kAnalysis,
            .resources = {RuntimeResourceAccess{
                .resource = RuntimeResourceKey::kDiagnostics,
                .mode = RuntimeResourceAccessMode::kWrite,
            }},
        }},
        .factory = [&probe_invocations](const RuntimeModuleFactoryContext&) {
          auto owner = std::make_shared<std::string>("test_probe");
          RuntimeModuleInstance instance;
          instance.owner_lifetime = owner;
          instance.stage_tasks.push_back(RuntimeStageTaskContribution{
              .task_id = "test_probe.analysis_diagnostic",
              .task = AnalysisStageTask(
                  [owner, &probe_invocations](AnalysisStageView& view) {
                    view.requireFresh();
                    assert(*owner == "test_probe");
                    ++probe_invocations;
                  }),
          });
          return instance;
        },
    });
  };

  const std::filesystem::path output_root =
      std::filesystem::temp_directory_path() / "cosmosim_runtime_descriptor_probe";
  const cosmosim::workflows::ReferenceWorkflowReport report =
      runner.run(output_root, options);
  assert(report.completed_steps == 1U);
  assert(report.canonical_stage_order);
  assert(probe_invocations == 1U);
  return 0;
}
