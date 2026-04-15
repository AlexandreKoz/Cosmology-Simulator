#include "cosmosim/workflows/reference_workflow.hpp"

#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"
#include "cosmosim/physics/black_hole_agn.hpp"
#include "cosmosim/physics/star_formation.hpp"

namespace cosmosim::workflows {

ReferenceWorkflowRunner::ReferenceWorkflowRunner(core::FrozenConfig frozen_config)
    : m_frozen_config(std::move(frozen_config)) {}

const core::FrozenConfig& ReferenceWorkflowRunner::frozenConfig() const noexcept {
  return m_frozen_config;
}

ReferenceWorkflowReport ReferenceWorkflowRunner::run(
    const std::filesystem::path& run_directory,
    const ReferenceWorkflowOptions& options) const {
  class StageAuditCallback final : public core::IntegrationCallback {
   public:
    explicit StageAuditCallback(std::vector<std::string>* stage_sequence)
        : m_stage_sequence(stage_sequence) {}
    [[nodiscard]] std::string_view callbackName() const override { return "stage_audit"; }
    void onStage(core::StepContext& context) override {
      m_stage_sequence->push_back(std::string(core::integrationStageName(context.stage)));
    }

   private:
    std::vector<std::string>* m_stage_sequence = nullptr;
  };

  const core::SimulationConfig& config = m_frozen_config.config;

  ReferenceWorkflowReport report;
  report.config_compatible =
      config.numerics.gravity_solver == core::GravitySolver::kTreePm &&
      config.numerics.hydro_solver == core::HydroSolver::kGodunovFv &&
      config.units.coordinate_frame == core::CoordinateFrame::kComoving;
  report.schema_compatible =
      config.schema_version == 1 && io::gadgetArepoSchemaMap().schema_version == 1 &&
      io::isRestartSchemaCompatible(io::restartSchema().version);

  core::ProfilerSession profiler(true);
  profiler.recordEvent(core::RuntimeEvent{
      .event_kind = "config.freeze",
      .severity = core::RuntimeEventSeverity::kInfo,
      .subsystem = "core.config",
      .step_index = options.step_index,
      .simulation_time_code = config.numerics.time_begin_code,
      .scale_factor = 1.0,
      .message = "frozen configuration accepted for reference workflow",
      .payload =
          {
              {"schema_version", std::to_string(config.schema_version)},
              {"config_hash_hex", m_frozen_config.provenance.config_hash_hex},
          },
  });

  if (!report.config_compatible || !report.schema_compatible) {
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "validation.fatal",
        .severity = core::RuntimeEventSeverity::kFatal,
        .subsystem = "workflows.reference",
        .step_index = options.step_index,
        .simulation_time_code = config.numerics.time_begin_code,
        .scale_factor = 1.0,
        .message = "reference workflow compatibility validation failed",
        .payload =
            {
                {"config_compatible", report.config_compatible ? "true" : "false"},
                {"schema_compatible", report.schema_compatible ? "true" : "false"},
            },
    });
    report.operational_report_json_path = run_directory / "reference_operational_events.json";
    core::writeOperationalReportJson(
        profiler,
        report.operational_report_json_path,
        config.output.run_name,
        m_frozen_config.provenance.config_hash_hex);
    throw std::runtime_error("reference workflow compatibility validation failed");
  }

  std::filesystem::create_directories(run_directory);

  core::SimulationConfig generated_ic_config = config;
  generated_ic_config.output.output_directory = run_directory.string();
  io::IcReadResult ic_result = io::buildGeneratedIsolatedIc(generated_ic_config, 16, 16, 1000);
  core::SimulationState state = std::move(ic_result.state);
  state.rebuildSpeciesIndex();

  core::IntegratorState integrator_state;
  integrator_state.step_index = options.step_index;
  integrator_state.current_time_code = config.numerics.time_begin_code;
  integrator_state.dt_time_code = options.dt_time_code;
  integrator_state.current_scale_factor = 1.0;
  integrator_state.time_bins.hierarchical_enabled = true;
  integrator_state.time_bins.max_bin =
      static_cast<std::uint8_t>(std::max(0, std::min(config.numerics.hierarchical_max_rung, 255)));

  const core::ModePolicy mode_policy = core::buildModePolicy(config.mode);
  core::CosmologyBackgroundConfig background_config;
  background_config.hubble_param = config.cosmology.hubble_param;
  background_config.omega_matter = config.cosmology.omega_matter;
  background_config.omega_lambda = config.cosmology.omega_lambda;
  const core::LambdaCdmBackground background(background_config);

  std::vector<std::uint32_t> particle_indices(state.particles.size());
  std::vector<std::uint32_t> cell_indices(state.cells.size());
  for (std::size_t i = 0; i < particle_indices.size(); ++i) {
    particle_indices[i] = static_cast<std::uint32_t>(i);
  }
  for (std::size_t i = 0; i < cell_indices.size(); ++i) {
    cell_indices[i] = static_cast<std::uint32_t>(i);
  }

  core::ActiveSetDescriptor active_set{.particle_indices = particle_indices,
                                        .cell_indices = cell_indices,
                                        .particles_are_subset = false,
                                        .cells_are_subset = false};

  core::StepOrchestrator orchestrator;
  StageAuditCallback stage_audit(&report.stage_sequence);
  orchestrator.registerCallback(stage_audit);

  analysis::DiagnosticsCallback diagnostics_callback(config);
  orchestrator.registerCallback(diagnostics_callback);

  physics::StarFormationCallback star_formation_callback(
      physics::StarFormationModel(physics::makeStarFormationConfig(config.physics)));
  orchestrator.registerCallback(star_formation_callback);

  physics::BlackHoleAgnCallback bh_callback(
      physics::BlackHoleAgnModel(physics::makeBlackHoleAgnConfig(config.physics)));
  orchestrator.registerCallback(bh_callback);

  core::TransientStepWorkspace workspace;
  orchestrator.executeSingleStep(
      state,
      integrator_state,
      active_set,
      &background,
      &workspace,
      &mode_policy,
      &profiler);

  std::vector<core::IntegrationStage> observed_stages;
  observed_stages.reserve(report.stage_sequence.size());
  for (const std::string& stage_name : report.stage_sequence) {
    for (const core::IntegrationStage expected_stage : core::StageScheduler::kickDriftKickOrder()) {
      if (stage_name == core::integrationStageName(expected_stage)) {
        observed_stages.push_back(expected_stage);
        break;
      }
    }
  }
  report.canonical_stage_order = core::isCanonicalIntegrationStageOrder(observed_stages);

  report.profiler_json_path = run_directory / "reference_profile.json";
  report.profiler_csv_path = run_directory / "reference_profile.csv";
  core::writeProfilerReportJson(profiler, report.profiler_json_path);
  core::writeProfilerReportCsv(profiler, report.profiler_csv_path);

#if COSMOSIM_ENABLE_HDF5
  if (options.write_outputs) {
    report.restart_roundtrip_executed = true;
    report.snapshot_roundtrip_executed = true;

    report.restart_path = run_directory / "reference_restart_0000.hdf5";
    io::RestartWritePayload restart_payload;
    restart_payload.state = &state;
    restart_payload.integrator_state = &integrator_state;
    restart_payload.provenance =
        core::makeProvenanceRecord(m_frozen_config.provenance.config_hash_hex, "reference_workflow");
    restart_payload.normalized_config_text = m_frozen_config.normalized_text;
    restart_payload.normalized_config_hash_hex = m_frozen_config.provenance.config_hash_hex;

    core::HierarchicalTimeBinScheduler scheduler_state(2);
    scheduler_state.reset(static_cast<std::uint32_t>(state.particles.size()), 0);
    restart_payload.scheduler = &scheduler_state;

    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "restart.write.begin",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "io.restart",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "writing restart checkpoint",
        .payload = {{"path", report.restart_path.string()}},
    });
    try {
      io::writeRestartCheckpointHdf5(report.restart_path, restart_payload);
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "restart.write.complete",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "io.restart",
          .step_index = integrator_state.step_index,
          .simulation_time_code = integrator_state.current_time_code,
          .scale_factor = integrator_state.current_scale_factor,
          .message = "restart checkpoint written",
          .payload = {{"path", report.restart_path.string()}},
      });
    } catch (const std::exception& ex) {
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "restart.write.failure",
          .severity = core::RuntimeEventSeverity::kError,
          .subsystem = "io.restart",
          .step_index = integrator_state.step_index,
          .simulation_time_code = integrator_state.current_time_code,
          .scale_factor = integrator_state.current_scale_factor,
          .message = "restart checkpoint write failed",
          .payload = {{"path", report.restart_path.string()}, {"error", ex.what()}},
      });
      throw;
    }
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "restart.read.begin",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "io.restart",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "reading restart checkpoint",
        .payload = {{"path", report.restart_path.string()}},
    });
    try {
      const io::RestartReadResult restart_read = io::readRestartCheckpointHdf5(report.restart_path);
      report.restart_roundtrip_ok = restart_read.state.particles.size() == state.particles.size();
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "restart.read.complete",
          .severity = report.restart_roundtrip_ok ? core::RuntimeEventSeverity::kInfo
                                                  : core::RuntimeEventSeverity::kWarning,
          .subsystem = "io.restart",
          .step_index = integrator_state.step_index,
          .simulation_time_code = integrator_state.current_time_code,
          .scale_factor = integrator_state.current_scale_factor,
          .message = "restart checkpoint read complete",
          .payload =
              {
                  {"path", report.restart_path.string()},
                  {"roundtrip_ok", report.restart_roundtrip_ok ? "true" : "false"},
              },
      });
    } catch (const std::exception& ex) {
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "restart.read.failure",
          .severity = core::RuntimeEventSeverity::kError,
          .subsystem = "io.restart",
          .step_index = integrator_state.step_index,
          .simulation_time_code = integrator_state.current_time_code,
          .scale_factor = integrator_state.current_scale_factor,
          .message = "restart checkpoint read failed",
          .payload = {{"path", report.restart_path.string()}, {"error", ex.what()}},
      });
      throw;
    }

    report.snapshot_path = run_directory / "reference_snapshot_0000.hdf5";
    io::SnapshotWritePayload snapshot_payload;
    snapshot_payload.state = &state;
    snapshot_payload.config = &config;
    snapshot_payload.normalized_config_text = m_frozen_config.normalized_text;
    snapshot_payload.provenance =
        core::makeProvenanceRecord(m_frozen_config.provenance.config_hash_hex, "reference_workflow");
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "snapshot.write.begin",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "io.snapshot",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "writing snapshot output",
        .payload = {{"path", report.snapshot_path.string()}},
    });
    try {
      io::writeGadgetArepoSnapshotHdf5(report.snapshot_path, snapshot_payload);
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "snapshot.write.complete",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "io.snapshot",
          .step_index = integrator_state.step_index,
          .simulation_time_code = integrator_state.current_time_code,
          .scale_factor = integrator_state.current_scale_factor,
          .message = "snapshot output written",
          .payload = {{"path", report.snapshot_path.string()}},
      });
    } catch (const std::exception& ex) {
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "snapshot.write.failure",
          .severity = core::RuntimeEventSeverity::kError,
          .subsystem = "io.snapshot",
          .step_index = integrator_state.step_index,
          .simulation_time_code = integrator_state.current_time_code,
          .scale_factor = integrator_state.current_scale_factor,
          .message = "snapshot write failed",
          .payload = {{"path", report.snapshot_path.string()}, {"error", ex.what()}},
      });
      throw;
    }

    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "snapshot.read.begin",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "io.snapshot",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "reading snapshot output",
        .payload = {{"path", report.snapshot_path.string()}},
    });
    try {
      const io::SnapshotReadResult snapshot_read =
          io::readGadgetArepoSnapshotHdf5(report.snapshot_path, config);
      report.snapshot_roundtrip_ok = snapshot_read.state.particles.size() == state.particles.size();
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "snapshot.read.complete",
          .severity = report.snapshot_roundtrip_ok ? core::RuntimeEventSeverity::kInfo
                                                   : core::RuntimeEventSeverity::kWarning,
          .subsystem = "io.snapshot",
          .step_index = integrator_state.step_index,
          .simulation_time_code = integrator_state.current_time_code,
          .scale_factor = integrator_state.current_scale_factor,
          .message = "snapshot read complete",
          .payload =
              {
                  {"path", report.snapshot_path.string()},
                  {"roundtrip_ok", report.snapshot_roundtrip_ok ? "true" : "false"},
              },
      });
    } catch (const std::exception& ex) {
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "snapshot.read.failure",
          .severity = core::RuntimeEventSeverity::kError,
          .subsystem = "io.snapshot",
          .step_index = integrator_state.step_index,
          .simulation_time_code = integrator_state.current_time_code,
          .scale_factor = integrator_state.current_scale_factor,
          .message = "snapshot read failed",
          .payload = {{"path", report.snapshot_path.string()}, {"error", ex.what()}},
      });
      throw;
    }
  }
#endif
  report.operational_report_json_path = run_directory / "reference_operational_events.json";
  core::writeOperationalReportJson(
      profiler,
      report.operational_report_json_path,
      config.output.run_name,
      m_frozen_config.provenance.config_hash_hex);

  return report;
}

}  // namespace cosmosim::workflows
