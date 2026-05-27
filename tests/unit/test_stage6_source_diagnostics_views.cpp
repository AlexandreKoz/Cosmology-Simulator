#include <cassert>
#include <cstdint>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/physics/black_hole_agn.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"
#include "cosmosim/physics/tracer_support.hpp"

namespace {

cosmosim::core::SimulationState makeSourceState() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(2);
  state.resizeCells(2);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 100 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
    state.particles.mass_code[i] = 1.0;
    state.particles.position_x_comoving[i] = static_cast<double>(i);
    state.particles.velocity_y_peculiar[i] = 1.0;
  }
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.center_y_comoving[i] = 0.0;
    state.cells.center_z_comoving[i] = 0.0;
    state.cells.mass_code[i] = 10.0;
    state.gas_cells.density_code[i] = 100.0;
    state.gas_cells.temperature_code[i] = 100.0;
    state.gas_cells.internal_energy_code[i] = 1.0;
  }
  state.rebuildSpeciesIndex();
  return state;
}

void testStarFormationUsesNarrowRuntimeViewForCellHotFields() {
  cosmosim::core::SimulationState state = makeSourceState();
  cosmosim::physics::StarFormationConfig config;
  config.density_threshold_code = 1.0;
  config.temperature_threshold_k = 1.0e6;
  config.min_converging_flow_rate_code = 0.0;
  config.epsilon_ff = 1.0;
  config.min_star_particle_mass_code = 0.5;
  config.stochastic_spawning = false;
  config.newton_g_code = 1.0;
  cosmosim::physics::StarFormationModel model(config);

  std::vector<std::uint32_t> active_cells{1};
  std::vector<double> divergence(2, -1.0);
  std::vector<double> metallicity(2, 0.02);
  auto mass_before = state.cells.mass_code[0];

  cosmosim::physics::StarFormationRuntimeView view{
      .active_cell_indices = active_cells,
      .center_x_comoving = state.cells.center_x_comoving,
      .center_y_comoving = state.cells.center_y_comoving,
      .center_z_comoving = state.cells.center_z_comoving,
      .gas_mass_code = state.cells.mass_code,
      .gas_density_code = state.gas_cells.density_code,
      .gas_temperature_k = state.gas_cells.temperature_code,
      .velocity_divergence_code = divergence,
      .metallicity_mass_fraction = metallicity,
  };

  const auto report = model.applyFromView(state, view, 0.05, 0.5, 7, 0);
  assert(report.counters.scanned_cells == 1);
  assert(report.counters.spawned_particles == 1);
  assert(state.cells.mass_code[0] == mass_before);
  assert(state.cells.mass_code[1] < 10.0);
  assert(state.star_particles.size() == 1);
}

void testStellarFeedbackGeometryAndDepositionViewsAreSufficientForHotLoop() {
  cosmosim::core::SimulationState state = makeSourceState();
  const std::size_t star_particle = state.particles.size();
  state.resizeParticles(star_particle + 1);
  state.particles.position_x_comoving[star_particle] = 0.0;
  state.particles.mass_code[star_particle] = 1.0;
  state.particle_sidecar.particle_id[star_particle] = 900;
  state.particle_sidecar.species_tag[star_particle] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar);
  state.star_particles.resize(1);
  state.star_particles.particle_index[0] = static_cast<std::uint32_t>(star_particle);
  state.star_particles.birth_mass_code[0] = 1.0;

  cosmosim::physics::StellarFeedbackConfig config;
  config.neighbor_count = 1;
  config.variant = cosmosim::physics::StellarFeedbackVariant::kNone;
  config.mode = cosmosim::physics::StellarFeedbackMode::kThermal;
  config.sn_energy_erg_per_mass_code = 10.0;
  cosmosim::physics::StellarFeedbackModel model(config);
  cosmosim::physics::StellarFeedbackModuleState module_state;

  cosmosim::physics::StellarFeedbackGeometryView geometry{
      .particle_position_x_comoving = state.particles.position_x_comoving,
      .particle_position_y_comoving = state.particles.position_y_comoving,
      .particle_position_z_comoving = state.particles.position_z_comoving,
      .cell_center_x_comoving = state.cells.center_x_comoving,
      .cell_center_y_comoving = state.cells.center_y_comoving,
      .cell_center_z_comoving = state.cells.center_z_comoving,
  };
  cosmosim::physics::StellarFeedbackDepositionView deposition{
      .cell_mass_code = state.cells.mass_code,
      .gas_density_code = state.gas_cells.density_code,
      .gas_internal_energy_code = state.gas_cells.internal_energy_code,
  };
  const std::uint32_t active_star = 0;
  const double returned_mass = 0.5;
  const double returned_metals = 0.1;
  const auto report = model.applyWithViews(
      state, module_state, geometry, deposition, std::span<const std::uint32_t>(&active_star, 1),
      std::span<const double>(&returned_mass, 1), std::span<const double>(&returned_metals, 1), 1.0);
  assert(report.counters.feedback_stars == 1);
  assert(report.counters.target_cells_visited == 1);
  assert(state.cells.mass_code[0] > 10.0 || state.cells.mass_code[1] > 10.0);
}


void testBlackHoleAccretionUsesNarrowRuntimeView() {
  cosmosim::core::SimulationState state = makeSourceState();
  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 0;
  state.black_holes.host_cell_index[0] = 0;
  state.black_holes.subgrid_mass_code[0] = 10.0;
  state.black_holes.duty_cycle_total_time_code[0] = 0.0;
  state.gas_cells.density_code[0] = 1.0e6;
  state.gas_cells.sound_speed_code[0] = 1.0e-3;
  state.gas_cells.internal_energy_code[0] = 1.0;

  cosmosim::physics::BlackHoleAgnConfig config;
  config.enabled = true;
  config.seed_halo_mass_threshold_code = 1.0;
  config.seed_mass_code = 1.0;
  config.alpha_bondi = 1.0;
  config.use_eddington_cap = false;
  config.newton_g_code = 1.0;
  config.speed_of_light_code = 1.0;
  config.proton_mass_code = 1.0;
  config.thomson_cross_section_code = 1.0;
  cosmosim::physics::BlackHoleAgnModel model(config);

  const std::uint32_t active_bh = 0;
  cosmosim::physics::BlackHoleAgnAccretionView view{
      .active_black_hole_indices = std::span<const std::uint32_t>(&active_bh, 1),
      .particle_index = state.black_holes.particle_index,
      .host_cell_index = state.black_holes.host_cell_index,
      .subgrid_mass_code = state.black_holes.subgrid_mass_code,
      .accretion_rate_code = state.black_holes.accretion_rate_code,
      .feedback_energy_code = state.black_holes.feedback_energy_code,
      .eddington_ratio = state.black_holes.eddington_ratio,
      .cumulative_accreted_mass_code = state.black_holes.cumulative_accreted_mass_code,
      .cumulative_feedback_energy_code = state.black_holes.cumulative_feedback_energy_code,
      .duty_cycle_active_time_code = state.black_holes.duty_cycle_active_time_code,
      .duty_cycle_total_time_code = state.black_holes.duty_cycle_total_time_code,
      .gas_density_code = state.gas_cells.density_code,
      .gas_sound_speed_code = state.gas_cells.sound_speed_code,
      .gas_internal_energy_code = state.gas_cells.internal_energy_code,
      .particle_mass_code = state.particles.mass_code,
  };
  const auto counters = model.applyAccretionFromView(view, 0.01);
  assert(counters.scanned_bh == 1);
  assert(counters.active_bh == 1);
  assert(state.black_holes.subgrid_mass_code[0] > 10.0);
  assert(state.particles.mass_code[0] == state.black_holes.subgrid_mass_code[0]);
}

void testTracerMassUpdateUsesNarrowHostView() {
  cosmosim::core::SimulationState state = makeSourceState();
  state.resizeParticles(3);
  state.particles.mass_code[2] = 0.1;
  state.particle_sidecar.species_tag[2] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kTracer);
  state.tracers.resize(1);
  state.tracers.particle_index[0] = 2;
  state.tracers.host_cell_index[0] = 1;
  state.tracers.mass_fraction_of_host[0] = 0.25;

  cosmosim::physics::TracerConfig config;
  config.enabled = true;
  config.track_mass = true;
  cosmosim::physics::TracerModel model(config);
  const std::uint32_t active_cell = 1;
  cosmosim::physics::TracerHostMassView view{
      .active_cell_indices = std::span<const std::uint32_t>(&active_cell, 1),
      .tracer_particle_index = state.tracers.particle_index,
      .host_cell_index = state.tracers.host_cell_index,
      .mass_fraction_of_host = state.tracers.mass_fraction_of_host,
      .last_host_mass_code = state.tracers.last_host_mass_code,
      .cumulative_exchanged_mass_code = state.tracers.cumulative_exchanged_mass_code,
      .host_mass_code = state.cells.mass_code,
      .particle_mass_code = state.particles.mass_code,
  };
  const auto counters = model.updateMassFromHostCellsView(view);
  assert(counters.updated_tracers == 1);
  assert(state.particles.mass_code[2] == 2.5);
  assert(state.tracers.last_host_mass_code[0] == 10.0);
}


void testStellarEvolutionUsesNarrowRuntimeView() {
  cosmosim::core::SimulationState state = makeSourceState();
  state.resizeParticles(3);
  state.particles.mass_code[2] = 5.0;
  state.particle_sidecar.species_tag[2] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar);
  state.star_particles.resize(1);
  state.star_particles.particle_index[0] = 2;
  state.star_particles.birth_mass_code[0] = 5.0;
  state.star_particles.formation_scale_factor[0] = 0.5;

  cosmosim::physics::StellarEvolutionConfig config;
  config.enabled = true;
  config.hubble_time_years = 1.0e10;
  cosmosim::physics::StellarEvolutionBookkeeper bookkeeper(
      config, cosmosim::physics::StellarEvolutionTable::makeBuiltinReference());

  const std::uint32_t active_star = 0;
  cosmosim::physics::StellarEvolutionRuntimeView view{
      .active_star_indices = std::span<const std::uint32_t>(&active_star, 1),
      .particle_index = state.star_particles.particle_index,
      .birth_mass_code = state.star_particles.birth_mass_code,
      .formation_scale_factor = state.star_particles.formation_scale_factor,
      .stellar_age_years_last = state.star_particles.stellar_age_years_last,
      .stellar_returned_mass_cumulative_code = state.star_particles.stellar_returned_mass_cumulative_code,
      .stellar_returned_metals_cumulative_code = state.star_particles.stellar_returned_metals_cumulative_code,
      .stellar_feedback_energy_cumulative_erg = state.star_particles.stellar_feedback_energy_cumulative_erg,
      .returned_mass_channel_cumulative_code = {
          state.star_particles.stellar_returned_mass_channel_cumulative_code[0],
          state.star_particles.stellar_returned_mass_channel_cumulative_code[1],
          state.star_particles.stellar_returned_mass_channel_cumulative_code[2]},
      .returned_metals_channel_cumulative_code = {
          state.star_particles.stellar_returned_metals_channel_cumulative_code[0],
          state.star_particles.stellar_returned_metals_channel_cumulative_code[1],
          state.star_particles.stellar_returned_metals_channel_cumulative_code[2]},
      .feedback_energy_channel_cumulative_erg = {
          state.star_particles.stellar_feedback_energy_channel_cumulative_erg[0],
          state.star_particles.stellar_feedback_energy_channel_cumulative_erg[1],
          state.star_particles.stellar_feedback_energy_channel_cumulative_erg[2]},
      .particle_mass_code = state.particles.mass_code,
  };
  const auto report = bookkeeper.applyFromView(view, 1.0, 0.1);
  assert(report.counters.scanned_stars == 1);
  assert(report.counters.evolved_stars == 1);
  assert(state.particles.mass_code[2] < 5.0);
  assert(state.star_particles.stellar_returned_mass_cumulative_code[0] > 0.0);
}

void testDiagnosticsNarrowViewsProduceSameResultsAsStateWrappers() {
  cosmosim::core::SimulationConfig config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.cosmology.box_size_mpc_comoving = 4.0;
  config.analysis.quicklook_grid_n = 4;
  cosmosim::analysis::DiagnosticsEngine engine(config);
  const cosmosim::core::SimulationState state = makeSourceState();
  const auto view = cosmosim::analysis::buildDiagnosticsStateView(state);

  const auto health_from_view = engine.computeRunHealth(view);
  const auto health_from_state = engine.computeRunHealth(state);
  assert(health_from_view.particle_count == health_from_state.particle_count);
  assert(health_from_view.cell_count == health_from_state.cell_count);

  const auto angular_from_view = engine.computeAngularMomentumBudget(view.particles);
  const auto angular_from_state = engine.computeAngularMomentumBudget(state);
  assert(angular_from_view.total_l_code[2] == angular_from_state.total_l_code[2]);

  const auto projection_from_view = engine.computeGasXyProjectionDensity(view.gas_cells, 4);
  const auto projection_from_state = engine.computeGasXyProjectionDensity(state, 4);
  assert(projection_from_view == projection_from_state);
}

}  // namespace

int main() {
  testStarFormationUsesNarrowRuntimeViewForCellHotFields();
  testStellarFeedbackGeometryAndDepositionViewsAreSufficientForHotLoop();
  testBlackHoleAccretionUsesNarrowRuntimeView();
  testTracerMassUpdateUsesNarrowHostView();
  testStellarEvolutionUsesNarrowRuntimeView();
  testDiagnosticsNarrowViewsProduceSameResultsAsStateWrappers();
  return 0;
}
