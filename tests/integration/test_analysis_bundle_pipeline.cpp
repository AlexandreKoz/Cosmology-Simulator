#include <cassert>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/core/time_integration.hpp"

namespace {

class NoopCallback final : public cosmosim::core::IntegrationCallback {
 public:
  std::string_view callbackName() const override { return "noop"; }
  void onStage(cosmosim::core::StepContext& /*context*/) override {}
};

}  // namespace

int main() {
  cosmosim::core::SimulationConfig config;
  config.output.output_directory = "test_outputs";
  config.output.run_name = "integration_analysis";
  config.cosmology.box_size_mpc_comoving = 1.0;
  config.analysis.enable_diagnostics = true;
  config.analysis.diagnostics_execution_policy =
      cosmosim::core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional;
  config.analysis.run_health_interval_steps = 1;
  config.analysis.science_light_interval_steps = 1;
  config.analysis.science_heavy_interval_steps = 2;
  config.analysis.retention_bundle_count = 16;
  config.analysis.quicklook_grid_n = 4;
  config.analysis.power_spectrum_mesh_n = 4;
  config.analysis.power_spectrum_bin_count = 4;
  config.analysis.sf_history_bin_count = 4;
  config.analysis.diagnostics_stem = "diag";

  const std::filesystem::path output_root =
      std::filesystem::path(config.output.output_directory) / config.output.run_name;
  std::filesystem::remove_all(output_root);

  cosmosim::core::SimulationState state;
  state.resizeParticles(8);
  state.resizeCells(8);
  state.species.count_by_species = {4, 2, 2, 0, 0};

  for (std::size_t i = 0; i < 8; ++i) {
    state.particle_sidecar.particle_id[i] = 1000 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(
        (i < 4) ? cosmosim::core::ParticleSpecies::kDarkMatter
                : (i < 6) ? cosmosim::core::ParticleSpecies::kGas : cosmosim::core::ParticleSpecies::kStar);
    state.particles.mass_code[i] = 1.0;
    state.particles.position_x_comoving[i] = (static_cast<double>(i) + 0.5) / 8.0;
    state.particles.position_y_comoving[i] = 0.2;
    state.particles.position_z_comoving[i] = 0.7;
    state.particles.velocity_x_peculiar[i] = 0.1;
    state.particles.velocity_y_peculiar[i] = 0.2;
    state.particles.velocity_z_peculiar[i] = 0.3;
  }

  state.star_particles.resize(2);
  state.star_particles.formation_scale_factor[0] = 0.4;
  state.star_particles.formation_scale_factor[1] = 0.6;
  state.star_particles.birth_mass_code[0] = 1.0;
  state.star_particles.birth_mass_code[1] = 1.0;

  for (std::size_t i = 0; i < 8; ++i) {
    state.cells.center_x_comoving[i] = (static_cast<double>(i) + 0.5) / 8.0;
    state.cells.center_y_comoving[i] = 0.5;
    state.cells.center_z_comoving[i] = 0.5;
    state.cells.mass_code[i] = 1.0;
    state.gas_cells.density_code[i] = 2.0;
  }

  state.rebuildSpeciesIndex();

  cosmosim::analysis::DiagnosticsCallback diagnostics_callback(config);
  NoopCallback noop;

  cosmosim::core::StepOrchestrator orchestrator;
  orchestrator.registerCallback(noop);
  orchestrator.registerCallback(diagnostics_callback);

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.current_scale_factor = 0.5;
  integrator_state.dt_time_code = 0.01;

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

  orchestrator.executeSingleStep(state, integrator_state, active, nullptr, nullptr, nullptr);
  orchestrator.executeSingleStep(state, integrator_state, active, nullptr, nullptr, nullptr);

  const std::filesystem::path diagnostics_dir = output_root / "diagnostics";
  assert(std::filesystem::exists(diagnostics_dir));

  std::size_t json_count = 0;
  bool found_heavy = false;
  bool found_heavy_metadata = false;
  for (const auto& entry : std::filesystem::directory_iterator(diagnostics_dir)) {
    if (entry.path().extension() == ".json") {
      ++json_count;
      const std::string filename = entry.path().filename().string();
      if (filename.find("science_heavy") != std::string::npos) {
        found_heavy = true;
        std::ifstream in(entry.path());
        const std::string body((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        if (body.find("\"diagnostics_execution_policy\": \"all_including_provisional\"") !=
                std::string::npos &&
            body.find("\"tier\": \"reference_science\"") != std::string::npos &&
            body.find("\"maturity\": \"provisional\"") != std::string::npos &&
            body.find("\"name\": \"power_spectrum\"") != std::string::npos) {
          found_heavy_metadata = true;
        }
      }
    }
  }
  assert(json_count >= 3);
  assert(found_heavy);
  assert(found_heavy_metadata);

  const auto timing = diagnostics_callback.timing();
  assert(timing.run_health_calls >= 2);
  assert(timing.light_calls >= 2);
  assert(timing.heavy_calls >= 1);

  std::filesystem::remove_all(output_root);
  return 0;
}
