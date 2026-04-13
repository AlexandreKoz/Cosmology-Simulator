#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>

#include "cosmosim/analysis/halo_workflow.hpp"

namespace {

cosmosim::core::SimulationState makePipelineState() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(10);
  state.species.count_by_species = {10, 0, 0, 0, 0};
  state.metadata.normalized_config_hash = 1234567;

  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 1000 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particles.mass_code[i] = 1.0;
    state.particles.velocity_x_peculiar[i] = 0.01;
    state.particles.velocity_y_peculiar[i] = 0.02;
    state.particles.velocity_z_peculiar[i] = 0.03;
    state.particles.position_x_comoving[i] = (i < 5) ? (0.1 + 0.01 * static_cast<double>(i))
                                                     : (0.6 + 0.01 * static_cast<double>(i - 5));
    state.particles.position_y_comoving[i] = (i < 5) ? 0.10 : 0.60;
    state.particles.position_z_comoving[i] = (i < 5) ? 0.10 : 0.60;
  }

  return state;
}

}  // namespace

int main() {
  cosmosim::core::SimulationConfig config;
  config.output.output_directory = "test_outputs";
  config.output.run_name = "integration_halo";
  config.analysis.enable_halo_workflow = true;
  config.analysis.halo_on_the_fly = false;
  config.analysis.halo_fof_linking_length_factor = 0.3;
  config.analysis.halo_fof_min_group_size = 3;
  config.cosmology.box_size_mpc_comoving = 1.0;

  const auto output_root = std::filesystem::path(config.output.output_directory) / config.output.run_name;
  std::filesystem::remove_all(output_root);

  cosmosim::analysis::HaloWorkflowPlanner planner(config);
  const auto report = planner.runSnapshotWorkflow(makePipelineState(), 11, 0.25, nullptr);

  assert(std::filesystem::exists(report.halo_catalog_path));
  assert(std::filesystem::exists(report.merger_tree_plan_path));
  assert(report.catalog.halos.size() == 2);

  std::ifstream in(report.halo_catalog_path);
  std::string json_text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  assert(json_text.find("cosmosim_halo_catalog_v1") != std::string::npos);
  assert(json_text.find("\"schema_version\": 1") != std::string::npos);

  std::filesystem::remove_all(output_root);
  return 0;
}
