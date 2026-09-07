#include "cosmosim/workflows/analysis_runtime.hpp"

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/optional_diagnostic_cadence.hpp"

namespace cosmosim::workflows {
namespace {

class AnalysisRuntimeImpl final : public AnalysisRuntime {
 public:
  AnalysisRuntimeImpl(
      const core::SimulationConfig& config,
      std::vector<std::string>& stage_sequence,
      const RuntimeServices& services)
      : m_config(config),
        m_stage_sequence(&stage_sequence),
        m_services(services),
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
    if (context.integrator_state.step_index ==
        std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error("analysis completed-step index overflows uint64");
    }
    // AnalysisHooks observes the fully kicked end-of-step state before the
    // integrator commits that epoch. Label cadence and science products by the
    // state being measured, not by the still-committed step-begin metadata.
    const std::uint64_t step = context.integrator_state.step_index + 1U;
    const double scale_factor = context.timeline_step.scale_factor_end;
    if (!std::isfinite(scale_factor) || scale_factor <= 0.0) {
      throw std::runtime_error(
          "analysis completed-state scale factor must be finite and positive");
    }
    const auto run = [&](analysis::DiagnosticClass diagnostic_class) {
      const analysis::DiagnosticsBundle bundle = m_diagnostics.generateBundle(
          context.state,
          step,
          scale_factor,
          diagnostic_class,
          context.workspace);
      m_diagnostics.writeBundle(bundle);
    };
    if (step % static_cast<std::uint64_t>(
                   m_config.analysis.run_health_interval_steps) == 0) {
      run(analysis::DiagnosticClass::kRunHealth);
    }

    const core::MemoryPressure memory_pressure =
        m_services.memory_governor != nullptr
        ? m_services.memory_governor->snapshot().pressure
        : core::MemoryPressure::kGreen;
    const bool defer_optional_analysis =
        memory_pressure == core::MemoryPressure::kRed ||
        memory_pressure == core::MemoryPressure::kTrip;
    const bool science_light_due_now =
        step % static_cast<std::uint64_t>(
                   m_config.analysis.science_light_interval_steps) == 0 &&
        m_config.analysis.diagnostics_execution_policy !=
            core::AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthOnly;
    const bool science_heavy_due_now =
        step % static_cast<std::uint64_t>(
                   m_config.analysis.science_heavy_interval_steps) == 0 &&
        m_config.analysis.diagnostics_execution_policy ==
            core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional;
    const bool science_light_requested = science_light_due_now || m_science_light_pending.pending;
    const bool science_heavy_requested = science_heavy_due_now || m_science_heavy_pending.pending;
    if (defer_optional_analysis) {
      if (science_light_due_now) { m_science_light_pending.recordDue(step); }
      if (science_heavy_due_now) { m_science_heavy_pending.recordDue(step); }
      if (science_light_requested || science_heavy_requested) {
        m_services.profiler.recordEvent(core::RuntimeEvent{
            .event_kind = "analysis.memory_pressure_deferral",
            .severity = core::RuntimeEventSeverity::kInfo,
            .subsystem = "analysis.diagnostics",
            .step_index = step,
            .simulation_time_code = context.timeline_step.time_end_code,
            .scale_factor = scale_factor,
            .message = "optional science diagnostics deferred under process memory pressure",
            .payload = {
                {"pressure", std::string(core::memoryPressureLabel(memory_pressure))},
                {"science_light_due", science_light_due_now ? "true" : "false"},
                {"science_heavy_due", science_heavy_due_now ? "true" : "false"},
                {"science_light_pending", m_science_light_pending.pending ? "true" : "false"},
                {"science_heavy_pending", m_science_heavy_pending.pending ? "true" : "false"},
                {"science_light_missed_count", std::to_string(m_science_light_pending.missed_count)},
                {"science_heavy_missed_count", std::to_string(m_science_heavy_pending.missed_count)},
            },
        });
      }
    } else {
      const auto executeOptional = [&](analysis::DiagnosticClass diagnostic_class,
                                       internal::OptionalDiagnosticCadence& pending,
                                       bool due_now, const char* class_name) {
        if (!due_now && !pending.pending) { return; }
        const auto previous = pending;
        // Do not clear pending state until the physical product was written.
        // The engine's existing memory admission remains authoritative.
        run(diagnostic_class);
        pending.clear();
        if (previous.pending) {
          m_services.profiler.recordEvent(core::RuntimeEvent{
              .event_kind = "analysis.memory_pressure_catchup",
              .severity = core::RuntimeEventSeverity::kInfo,
              .subsystem = "analysis.diagnostics",
              .step_index = step,
              .simulation_time_code = context.timeline_step.time_end_code,
              .scale_factor = scale_factor,
              .message = "coalesced optional science diagnostic completed at the current physical epoch",
              .payload = {
                  {"diagnostic_class", class_name},
                  {"first_due_step", std::to_string(previous.first_due_step)},
                  {"latest_due_step", std::to_string(previous.latest_due_step)},
                  {"missed_count", std::to_string(previous.missed_count)},
                  {"coalesced_count", std::to_string(previous.missed_count - 1U)},
                  {"actual_execution_step", std::to_string(step)},
                  {"science_light_catchup", diagnostic_class == analysis::DiagnosticClass::kScienceLight ? "true" : "false"},
                  {"science_heavy_catchup", diagnostic_class == analysis::DiagnosticClass::kScienceHeavy ? "true" : "false"},
                  {"historical_state_replayed", "false"},
              },
          });
        }
      };
      executeOptional(analysis::DiagnosticClass::kScienceLight,
                      m_science_light_pending, science_light_due_now, "science_light");
      executeOptional(analysis::DiagnosticClass::kScienceHeavy,
                      m_science_heavy_pending, science_heavy_due_now, "science_heavy");
    }
    m_diagnostics.enforceRetentionPolicy();
  }

  [[nodiscard]] std::string optionalCadenceProvenance() const override {
    const auto describe = [](const internal::OptionalDiagnosticCadence& state) {
      return std::string("pending=") + (state.pending ? "true" : "false") +
          ",first_due_step=" + std::to_string(state.first_due_step) +
          ",latest_due_step=" + std::to_string(state.latest_due_step) +
          ",missed_count=" + std::to_string(state.missed_count);
    };
    return "optional_diagnostic_cadence_policy=coalesced_nonpersistent\n"
           "optional_diagnostic_science_light=" + describe(m_science_light_pending) + "\n" +
           "optional_diagnostic_science_heavy=" + describe(m_science_heavy_pending) + "\n";
  }

  void finalizePending(std::uint64_t completed_step, double time_code,
                       double scale_factor) override {
    const auto drop = [&](internal::OptionalDiagnosticCadence& pending,
                          const char* class_name) {
      if (!pending.pending) { return; }
      m_services.profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "analysis.optional_cadence_dropped",
          .severity = core::RuntimeEventSeverity::kWarning,
          .subsystem = "analysis.diagnostics",
          .step_index = completed_step,
          .simulation_time_code = time_code,
          .scale_factor = scale_factor,
          .message = "pending optional science cadence discarded at run termination; no historical state replay",
          .payload = {
              {"diagnostic_class", class_name},
              {"first_due_step", std::to_string(pending.first_due_step)},
              {"latest_due_step", std::to_string(pending.latest_due_step)},
              {"dropped_count", std::to_string(pending.missed_count)},
              {"reason", "run_termination"},
          },
      });
      pending.clear();
    };
    drop(m_science_light_pending, "science_light");
    drop(m_science_heavy_pending, "science_heavy");
  }

 private:
  core::SimulationConfig m_config;
  std::vector<std::string>* m_stage_sequence = nullptr;
  const RuntimeServices& m_services;
  analysis::DiagnosticsEngine m_diagnostics;
  internal::OptionalDiagnosticCadence m_science_light_pending;
  internal::OptionalDiagnosticCadence m_science_heavy_pending;
};

}  // namespace

std::unique_ptr<AnalysisRuntime> makeAnalysisRuntime(
    const core::SimulationConfig& config,
    std::vector<std::string>& stage_sequence,
    const RuntimeServices& services) {
  return std::make_unique<AnalysisRuntimeImpl>(config, stage_sequence, services);
}

}  // namespace cosmosim::workflows
