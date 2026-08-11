#include <array>
#include <cassert>
#include <cmath>
#include <string>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/black_hole_agn.hpp"

namespace {

void testAccretionFormulaAndEddingtonCap() {
  cosmosim::physics::BlackHoleAgnConfig config;
  config.enabled = true;
  config.alpha_bondi = 2.0;
  config.use_eddington_cap = true;
  cosmosim::physics::BlackHoleAgnModel model(config);

  const auto rates = model.computeAccretionRates(5.0, 20.0, 3.0, 4.0);
  const double denom = std::pow((3.0 * 3.0) + (4.0 * 4.0), 1.5);
  const double expected_bondi =
      config.alpha_bondi * 4.0 * cosmosim::core::constants::k_pi * config.newton_g_si * config.newton_g_si *
      25.0 * 20.0 / denom;
  assert(std::abs(rates.mdot_bondi_code - expected_bondi) / expected_bondi < 1.0e-12);
  assert(rates.mdot_acc_code <= rates.mdot_edd_code + 1.0e-30);
}

void testSeedEligibilityRespectsThresholdAndMultiplicity() {
  cosmosim::core::SimulationState state;
  state.resizeCells(2);
  state.cells.center_x_comoving[0] = 0.0;
  state.cells.center_x_comoving[1] = 1.0;
  state.cells.mass_code[0] = 10.0;
  state.cells.mass_code[1] = 10.0;

  cosmosim::physics::BlackHoleAgnConfig config;
  config.enabled = true;
  config.seed_halo_mass_threshold_code = 100.0;
  config.seed_max_per_cell = 1;
  cosmosim::physics::BlackHoleAgnModel model(config);

  cosmosim::physics::BlackHoleSeedCandidate under;
  under.cell_index = 0;
  under.host_halo_mass_code = 50.0;
  assert(!model.isSeedEligible(state, under));

  cosmosim::physics::BlackHoleSeedCandidate ok;
  ok.cell_index = 0;
  ok.host_halo_mass_code = 100.0;
  assert(model.isSeedEligible(state, ok));

  state.resizeParticles(1);
  state.particle_sidecar.species_tag[0] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kBlackHole);
  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 0;
  state.black_holes.host_cell_index[0] = 0;
  assert(!model.isSeedEligible(state, ok));
}

void testApplyMassGrowthFeedbackAndMetadata() {
  cosmosim::core::SimulationState state;
  state.resizeCells(1);
  state.cells.center_x_comoving[0] = 1.0;
  state.cells.center_y_comoving[0] = 2.0;
  state.cells.center_z_comoving[0] = 3.0;
  state.cells.time_bin[0] = 2;
  state.cells.mass_code[0] = 100.0;
  state.gas_cells.gas_cell_id[0] = 1234U;
  state.gas_cells.density_code[0] = 30.0;
  state.gas_cells.metal_mass_code[0] = 2.0;
  state.gas_cells.velocity_x_peculiar[0] = 1.5;
  state.gas_cells.sound_speed_code[0] = 5.0;
  state.gas_cells.internal_energy_code[0] = 10.0;

  cosmosim::physics::BlackHoleAgnConfig config;
  config.enabled = true;
  config.seed_halo_mass_threshold_code = 200.0;
  config.seed_mass_code = 4.0;
  cosmosim::physics::BlackHoleAgnModel model(config);

  const double total_mass_before_seed = state.cells.mass_code[0];
  const std::array<cosmosim::physics::BlackHoleSeedCandidate, 1> seeds{{{0, 250.0, 0}}};
  const auto seed_report = model.apply(state, seeds, 1.0, 0);
  assert(seed_report.counters.seeded_bh == 1);
  assert(state.black_holes.size() == 1);
  assert(std::abs((state.cells.mass_code[0] + state.particles.mass_code.back()) -
                  total_mass_before_seed) < 1.0e-12);
  assert(state.particles.velocity_x_peculiar.back() == 1.5);

  const double total_mass_before_growth = state.cells.mass_code[0] + state.particles.mass_code.back();
  const double energy_before = state.gas_cells.internal_energy_code[0];
  const auto growth_report = model.apply(state, std::array<cosmosim::physics::BlackHoleSeedCandidate, 0>{}, 2.0, 1);
  assert(growth_report.counters.active_bh == 1);
  assert(state.black_holes.cumulative_accreted_mass_code[0] > 0.0);
  assert(std::abs((state.cells.mass_code[0] + state.particles.mass_code.back()) -
                  total_mass_before_growth) < 1.0e-10);
  assert(growth_report.counters.gas_mass_removed_code ==
         growth_report.counters.integrated_accreted_mass_code);
  assert(state.gas_cells.internal_energy_code[0] > energy_before);
  assert(state.black_holes.duty_cycle_total_time_code[0] >= 2.0);

  const auto* sidecar = state.sidecars.find("black_hole_agn");
  assert(sidecar != nullptr);
  const std::string payload(reinterpret_cast<const char*>(sidecar->payload.data()), sidecar->payload.size());
  assert(payload.find("seeded_bh=") != std::string::npos);
}


void testCosmologicalDensityConversionAndGasLimitedAccretion() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(1);
  state.particles.mass_code[0] = 2.0;
  state.particles.velocity_x_peculiar[0] = 3.0;
  state.cells.mass_code[0] = 5.0;
  state.gas_cells.density_code[0] = 1.0;
  state.gas_cells.sound_speed_code[0] = 1.0;
  state.gas_cells.velocity_x_peculiar[0] = 1.0;
  state.gas_cells.metal_mass_code[0] = 0.5;
  state.gas_cells.internal_energy_code[0] = 1.0;
  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 0;
  state.black_holes.host_cell_index[0] = 0;
  state.black_holes.subgrid_mass_code[0] = 2.0;

  cosmosim::physics::BlackHoleAgnConfig config;
  config.enabled = true;
  config.alpha_bondi = 1.0;
  config.use_eddington_cap = false;
  config.newton_g_code = 1.0;
  config.speed_of_light_code = 1.0;
  config.epsilon_r = 0.1;
  config.epsilon_f = 0.0;
  cosmosim::physics::BlackHoleAgnModel model(config);

  const double expected_physical_density = 8.0;  // rho_comov / 0.5^3
  const double relative_velocity = 2.0;
  const auto expected = model.computeAccretionRates(
      2.0, expected_physical_density, 1.0, relative_velocity);
  const std::array<std::uint32_t, 1> active{0U};
  cosmosim::physics::BlackHoleAgnAccretionView view{
      .active_black_hole_indices = active,
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
      .gas_mass_code = state.cells.mass_code,
      .gas_density_code = state.gas_cells.density_code,
      .gas_metal_mass_code = state.gas_cells.metal_mass_code,
      .gas_sound_speed_code = state.gas_cells.sound_speed_code,
      .gas_velocity_x_peculiar = state.gas_cells.velocity_x_peculiar,
      .gas_velocity_y_peculiar = state.gas_cells.velocity_y_peculiar,
      .gas_velocity_z_peculiar = state.gas_cells.velocity_z_peculiar,
      .gas_internal_energy_code = state.gas_cells.internal_energy_code,
      .particle_mass_code = state.particles.mass_code,
      .particle_velocity_x_peculiar = state.particles.velocity_x_peculiar,
      .particle_velocity_y_peculiar = state.particles.velocity_y_peculiar,
      .particle_velocity_z_peculiar = state.particles.velocity_z_peculiar,
  };
  const double dt = 1.0e-4;
  const double total_mass_before = state.cells.mass_code[0] + state.particles.mass_code[0];
  const auto counters = model.applyAccretionFromView(view, dt, 0.5, true);
  assert(std::abs(state.black_holes.accretion_rate_code[0] - expected.mdot_acc_code) /
             std::max(expected.mdot_acc_code, 1.0e-30) < 1.0e-12);
  assert(std::abs((state.cells.mass_code[0] + state.particles.mass_code[0]) - total_mass_before) < 1.0e-12);
  assert(counters.gas_mass_removed_code == counters.integrated_accreted_mass_code);
}

}  // namespace

int main() {
  testAccretionFormulaAndEddingtonCap();
  testSeedEligibilityRespectsThresholdAndMultiplicity();
  testApplyMassGrowthFeedbackAndMetadata();
  testCosmologicalDensityConversionAndGasLimitedAccretion();
  return 0;
}
