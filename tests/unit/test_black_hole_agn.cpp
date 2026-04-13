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
  state.gas_cells.density_code[0] = 30.0;
  state.gas_cells.sound_speed_code[0] = 5.0;
  state.gas_cells.internal_energy_code[0] = 10.0;

  cosmosim::physics::BlackHoleAgnConfig config;
  config.enabled = true;
  config.seed_halo_mass_threshold_code = 200.0;
  config.seed_mass_code = 4.0;
  cosmosim::physics::BlackHoleAgnModel model(config);

  const std::array<cosmosim::physics::BlackHoleSeedCandidate, 1> seeds{{{0, 250.0, 0}}};
  const auto seed_report = model.apply(state, seeds, 1.0, 0);
  assert(seed_report.counters.seeded_bh == 1);
  assert(state.black_holes.size() == 1);

  const double energy_before = state.gas_cells.internal_energy_code[0];
  const auto growth_report = model.apply(state, std::array<cosmosim::physics::BlackHoleSeedCandidate, 0>{}, 2.0, 1);
  assert(growth_report.counters.active_bh == 1);
  assert(state.black_holes.cumulative_accreted_mass_code[0] > 0.0);
  assert(state.gas_cells.internal_energy_code[0] > energy_before);
  assert(state.black_holes.duty_cycle_total_time_code[0] >= 2.0);

  const auto* sidecar = state.sidecars.find("black_hole_agn");
  assert(sidecar != nullptr);
  const std::string payload(reinterpret_cast<const char*>(sidecar->payload.data()), sidecar->payload.size());
  assert(payload.find("seeded_bh=") != std::string::npos);
}

}  // namespace

int main() {
  testAccretionFormulaAndEddingtonCap();
  testSeedEligibilityRespectsThresholdAndMultiplicity();
  testApplyMassGrowthFeedbackAndMetadata();
  return 0;
}
