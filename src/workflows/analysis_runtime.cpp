#include "cosmosim/workflows/analysis_runtime.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"

namespace cosmosim::workflows {
namespace {

class AnalysisRuntimeImpl final : public AnalysisRuntime {
 public:
  AnalysisRuntimeImpl(
      const core::SimulationConfig& config,
      std::vector<std::string>& stage_sequence)
      : m_config(config),
        m_stage_sequence(&stage_sequence),
        m_diagnostics(config) {}

  void audit(AnalysisStageView& view) override {
    view.requireFresh();
    const core::StepContext& context = internal::RuntimeStageAccess::analysisContext(
        view,
        {{RuntimeResourceKey::kDiagnostics, RuntimeResourceAccessMode::kWrite}});
    m_stage_sequence->push_back(
        std::string(core::integrationStageName(context.stage)));
  }

  void executeDiagnostics(AnalysisStageView& view) override {
    view.requireFresh();
    core::StepContext& context = internal::RuntimeStageAccess::analysisContext(
        view,
        {{RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kParticleGravitySource, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kHydroPrimitiveState, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kSourceMutationState, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kIntegratorTruth, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kDiagnostics, RuntimeResourceAccessMode::kWrite}});
    if (context.stage != core::IntegrationStage::kAnalysisHooks) {
      throw std::logic_error("analysis diagnostics handler received an unregistered stage");
    }
    if (!m_config.analysis.enable_diagnostics) {
      return;
    }
    const std::uint64_t step = context.integrator_state.step_index;
    const auto run = [&](analysis::DiagnosticClass diagnostic_class) {
      const analysis::DiagnosticsBundle bundle = m_diagnostics.generateBundle(
          context.state,
          step,
          context.integrator_state.current_scale_factor,
          diagnostic_class,
          context.workspace);
      m_diagnostics.writeBundle(bundle);
    };
    if (step % static_cast<std::uint64_t>(
                   m_config.analysis.run_health_interval_steps) == 0) {
      run(analysis::DiagnosticClass::kRunHealth);
    }
    if (step % static_cast<std::uint64_t>(
                   m_config.analysis.science_light_interval_steps) == 0 &&
        m_config.analysis.diagnostics_execution_policy !=
            core::AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthOnly) {
      run(analysis::DiagnosticClass::kScienceLight);
    }
    if (step % static_cast<std::uint64_t>(
                   m_config.analysis.science_heavy_interval_steps) == 0 &&
        m_config.analysis.diagnostics_execution_policy ==
            core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional) {
      run(analysis::DiagnosticClass::kScienceHeavy);
    }
    m_diagnostics.enforceRetentionPolicy();
  }

 private:
  core::SimulationConfig m_config;
  std::vector<std::string>* m_stage_sequence = nullptr;
  analysis::DiagnosticsEngine m_diagnostics;
};

}  // namespace

std::unique_ptr<AnalysisRuntime> makeAnalysisRuntime(
    const core::SimulationConfig& config,
    std::vector<std::string>& stage_sequence) {
  return std::make_unique<AnalysisRuntimeImpl>(config, stage_sequence);
}

}  // namespace cosmosim::workflows
