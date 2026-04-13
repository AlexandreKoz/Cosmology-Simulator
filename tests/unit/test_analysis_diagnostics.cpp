#include <cassert>
#include <cmath>
#include <filesystem>
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
}

}  // namespace

int main() {
  testPowerSpectrumHasSignalAndFiniteModes();
  testDerivedDiagnosticsSanity();
  return 0;
}
