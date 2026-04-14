#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>

#include "cosmosim/analysis/diagnostics.hpp"

namespace {

cosmosim::core::SimulationConfig makeConfig() {
  cosmosim::core::SimulationConfig config;
  config.output.run_name = "unit_analysis";
  config.output.output_directory = "test_outputs";
  config.cosmology.box_size_mpc_comoving = 1.0;
  config.analysis.power_spectrum_mesh_n = 4;
  config.analysis.power_spectrum_bin_count = 4;
  config.analysis.quicklook_grid_n = 4;
  config.analysis.sf_history_bin_count = 4;
  config.analysis.diagnostics_stem = "diag";
  return config;
}

cosmosim::core::SimulationState makeSingleModeState() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(4);
  state.resizeCells(4);
  state.species.count_by_species = {4, 0, 0, 0, 0};
  for (std::size_t i = 0; i < 4; ++i) {
    state.particle_sidecar.particle_id[i] = 100 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particles.mass_code[i] = 1.0;
    state.particles.position_y_comoving[i] = 0.125;
    state.particles.position_z_comoving[i] = 0.125;
  }

  state.particles.position_x_comoving[0] = 0.125;
  state.particles.position_x_comoving[1] = 0.125;
  state.particles.position_x_comoving[2] = 0.625;
  state.particles.position_x_comoving[3] = 0.625;

  for (std::size_t i = 0; i < 4; ++i) {
    state.cells.center_x_comoving[i] = (static_cast<double>(i) + 0.5) / 4.0;
    state.cells.center_y_comoving[i] = 0.125;
    state.cells.center_z_comoving[i] = 0.5;
    state.cells.mass_code[i] = 1.0;
    state.gas_cells.density_code[i] = 2.0 + static_cast<double>(i);
  }

  state.rebuildSpeciesIndex();
  return state;
}

void testPowerSpectrumHasSignalAndFiniteModes() {
  cosmosim::analysis::DiagnosticsEngine engine(makeConfig());
  const cosmosim::core::SimulationState state = makeSingleModeState();

  const auto bins = engine.computePowerSpectrum(state, 4, 4);
  assert(!bins.empty());

  double summed_power = 0.0;
  std::uint64_t total_modes = 0;
  for (const auto& bin : bins) {
    assert(std::isfinite(bin.k_center_code));
    assert(std::isfinite(bin.power_code_volume));
    assert(bin.mode_count > 0);
    summed_power += bin.power_code_volume;
    total_modes += bin.mode_count;
  }

  assert(summed_power > 0.0);
  assert(total_modes > 0);
}

void testDerivedDiagnosticsSanity() {
  cosmosim::core::SimulationConfig config = makeConfig();
  cosmosim::analysis::DiagnosticsEngine engine(config);
  cosmosim::core::SimulationState state = makeSingleModeState();

  state.star_particles.resize(2);
  state.star_particles.formation_scale_factor[0] = 0.2;
  state.star_particles.formation_scale_factor[1] = 0.8;
  state.star_particles.birth_mass_code[0] = 2.0;
  state.star_particles.birth_mass_code[1] = 3.0;

  state.particles.velocity_y_peculiar[0] = 1.0;
  state.particles.position_x_comoving[0] = 0.25;

  const auto sfh = engine.computeStarFormationHistory(state, 4);
  double formed = 0.0;
  for (const auto& bin : sfh) {
    formed += bin.formed_mass_code;
  }
  assert(std::abs(formed - 5.0) < 1.0e-12);

  const auto angular = engine.computeAngularMomentumBudget(state);
  const double l_norm =
      std::sqrt(angular.total_l_code[0] * angular.total_l_code[0] + angular.total_l_code[1] * angular.total_l_code[1] +
                angular.total_l_code[2] * angular.total_l_code[2]);
  assert(std::isfinite(l_norm));
  assert(l_norm > 0.0);

  const auto bundle = engine.generateBundle(state, 4, 0.7, cosmosim::analysis::DiagnosticClass::kScienceLight);
  assert(bundle.xy_projection_density_code.size() == 16);
  assert(bundle.xy_slice_density_code.size() == 16);
  assert(!bundle.records.empty());
  assert(bundle.records.front().name == "run_health_counters");
}

void testProvisionalHeavyDiagnosticsRequireExplicitPolicy() {
  cosmosim::core::SimulationConfig config = makeConfig();
  config.analysis.diagnostics_execution_policy =
      cosmosim::core::AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthAndLightScience;
  config.analysis.enable_diagnostics = true;
  config.analysis.run_health_interval_steps = 1;
  config.analysis.science_light_interval_steps = 1;
  config.analysis.science_heavy_interval_steps = 1;
  config.output.run_name = "unit_analysis_policy_default";

  const std::filesystem::path output_root =
      std::filesystem::path(config.output.output_directory) / config.output.run_name;
  std::filesystem::remove_all(output_root);

  cosmosim::analysis::DiagnosticsCallback callback(config);
  cosmosim::core::SimulationState state = makeSingleModeState();
  cosmosim::core::IntegratorState integrator_state;
  integrator_state.current_scale_factor = 0.5;
  integrator_state.step_index = 0;

  std::vector<std::uint32_t> particle_indices(state.particles.size());
  std::vector<std::uint32_t> cell_indices(state.cells.size());
  for (std::size_t i = 0; i < particle_indices.size(); ++i) {
    particle_indices[i] = static_cast<std::uint32_t>(i);
  }
  for (std::size_t i = 0; i < cell_indices.size(); ++i) {
    cell_indices[i] = static_cast<std::uint32_t>(i);
  }

  const cosmosim::core::ActiveSetDescriptor active{
      .particle_indices = particle_indices,
      .cell_indices = cell_indices,
      .particles_are_subset = false,
      .cells_are_subset = false,
  };
  cosmosim::core::StepContext context{
      .state = state,
      .integrator_state = integrator_state,
      .active_set = active,
      .workspace = nullptr,
      .cosmology_background = nullptr,
      .mode_policy = nullptr,
      .stage = cosmosim::core::IntegrationStage::kAnalysisHooks,
  };

  callback.onStage(context);
  const auto timing = callback.timing();
  assert(timing.heavy_calls == 0);
  assert(timing.light_calls == 1);
  assert(timing.run_health_calls == 1);

  bool found_heavy_bundle = false;
  for (const auto& entry : std::filesystem::directory_iterator(output_root / "diagnostics")) {
    if (entry.path().filename().string().find("science_heavy") != std::string::npos) {
      found_heavy_bundle = true;
    }
  }
  assert(!found_heavy_bundle);

  std::filesystem::remove_all(output_root);
}

void testEngineLevelHeavyDiagnosticsQuarantineAndTruthfulRecords() {
  cosmosim::core::SimulationConfig blocked_config = makeConfig();
  blocked_config.analysis.diagnostics_execution_policy =
      cosmosim::core::AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthAndLightScience;

  cosmosim::analysis::DiagnosticsEngine blocked_engine(blocked_config);
  const cosmosim::core::SimulationState state = makeSingleModeState();

  const auto blocked_bundle =
      blocked_engine.generateBundle(state, 7, 0.6, cosmosim::analysis::DiagnosticClass::kScienceHeavy);
  assert(blocked_bundle.power_spectrum.empty());

  bool found_power_record = false;
  for (const auto& record : blocked_bundle.records) {
    if (record.name == "power_spectrum") {
      found_power_record = true;
      assert(record.tier == cosmosim::analysis::DiagnosticTier::kReferenceScience);
      assert(record.maturity == cosmosim::analysis::DiagnosticMaturity::kProvisional);
      assert(record.scalability == cosmosim::analysis::DiagnosticScalability::kHeavyReference);
      assert(!record.executed);
      assert(record.policy_note == "blocked_by_execution_policy");
    }
  }
  assert(found_power_record);

  cosmosim::core::SimulationConfig allowed_config = makeConfig();
  allowed_config.analysis.diagnostics_execution_policy =
      cosmosim::core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional;
  cosmosim::analysis::DiagnosticsEngine allowed_engine(allowed_config);
  const auto allowed_bundle =
      allowed_engine.generateBundle(state, 7, 0.6, cosmosim::analysis::DiagnosticClass::kScienceHeavy);
  assert(!allowed_bundle.power_spectrum.empty());
  found_power_record = false;
  for (const auto& record : allowed_bundle.records) {
    if (record.name == "power_spectrum") {
      found_power_record = true;
      assert(record.executed);
      assert(record.policy_note == "reference_only_non_default");
    }
  }
  assert(found_power_record);
}

}  // namespace

int main() {
  testPowerSpectrumHasSignalAndFiniteModes();
  testDerivedDiagnosticsSanity();
  testProvisionalHeavyDiagnosticsRequireExplicitPolicy();
  testEngineLevelHeavyDiagnosticsQuarantineAndTruthfulRecords();
  return 0;
}
