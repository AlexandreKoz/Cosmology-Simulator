#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/physics/star_formation.hpp"

namespace {

[[nodiscard]] bool nearlyEqual(double lhs, double rhs, double relative = 1.0e-12) {
  return std::abs(lhs - rhs) <= relative * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

[[nodiscard]] cosmosim::physics::StarFormationConfig makeAdaptiveConfig() {
  cosmosim::physics::StarFormationConfig config;
  config.model = cosmosim::core::StarFormationModelKind::kAdaptiveBoundJeans;
  config.epsilon_ff = 1.0;
  config.bound_alpha_vir_max = 1.0;
  config.require_converging_flow = true;
  config.collapse_timescale =
      cosmosim::core::StarFormationCollapseTimescale::kMinimumFreeFallOrCompression;
  config.jeans_mass_floor_code = 0.1;
  config.star_particle_mass_policy = cosmosim::core::StarParticleMassPolicy::kFixed;
  config.target_star_particle_mass_code = 0.5;
  config.min_star_particle_mass_code = 0.1;
  config.max_star_particle_mass_code = 10.0;
  config.max_spawn_particles_per_cell_step = 8;
  config.max_fractional_mass_conversion = 0.25;
  config.min_remaining_gas_fraction = 0.01;
  config.min_remaining_gas_mass_code = 0.0;
  config.stochastic_spawning = true;
  config.random_seed = 123456789ULL;
  config.newton_g_code = 1.0;
  config.density_is_comoving = false;
  config.geometry_is_comoving = false;
  return config;
}

[[nodiscard]] cosmosim::physics::StarFormationCellInput makeEligibleAdaptiveCell(
    std::uint32_t row = 0,
    std::uint64_t gas_cell_id = 42) {
  cosmosim::physics::StarFormationCellInput cell;
  cell.cell_index = row;
  cell.gas_cell_id = gas_cell_id;
  cell.owning_rank = 0;
  cell.is_active = true;
  cell.is_owned = true;
  cell.is_leaf = true;
  cell.is_ghost = false;
  cell.gas_mass_code = 10.0;
  cell.gas_density_code = 100.0;
  cell.cell_volume_code = 0.1;
  cell.gas_temperature_k = 100.0;
  cell.gas_sound_speed_code = 0.01;
  cell.velocity_x_peculiar = 2.0;
  cell.velocity_y_peculiar = -3.0;
  cell.velocity_z_peculiar = 4.0;
  cell.velocity_divergence_code = -20.0;
  cell.velocity_gradient_frobenius_sq_code = 0.0;
  cell.gas_metal_mass_code = 0.2;
  cell.center_x_comoving = 1.0 + row;
  cell.center_y_comoving = 2.0 + row;
  cell.center_z_comoving = 3.0 + row;
  return cell;
}

void initializeAdaptiveState(
    cosmosim::core::SimulationState& state,
    std::span<const cosmosim::physics::StarFormationCellInput> cells) {
  state.resizeCells(cells.size());
  for (const auto& cell : cells) {
    const std::size_t row = cell.cell_index;
    state.cells.center_x_comoving[row] = cell.center_x_comoving;
    state.cells.center_y_comoving[row] = cell.center_y_comoving;
    state.cells.center_z_comoving[row] = cell.center_z_comoving;
    state.cells.mass_code[row] = cell.gas_mass_code;
    state.gas_cells.gas_cell_id[row] = cell.gas_cell_id;
    state.gas_cells.density_code[row] = cell.gas_density_code;
    state.gas_cells.temperature_code[row] = cell.gas_temperature_k;
    state.gas_cells.sound_speed_code[row] = cell.gas_sound_speed_code;
    state.gas_cells.pressure_code[row] = 5.0;
    state.gas_cells.internal_energy_code[row] = 2.0;
    state.gas_cells.velocity_x_peculiar[row] = cell.velocity_x_peculiar;
    state.gas_cells.velocity_y_peculiar[row] = cell.velocity_y_peculiar;
    state.gas_cells.velocity_z_peculiar[row] = cell.velocity_z_peculiar;
    state.gas_cells.metal_mass_code[row] = cell.gas_metal_mass_code;
  }
}

void testLegacyThresholdEligibility() {
  cosmosim::physics::StarFormationConfig config;
  config.density_threshold_code = 5.0;
  config.temperature_threshold_k = 1.0e4;
  config.min_converging_flow_rate_code = -0.2;
  cosmosim::physics::StarFormationModel model(config);

  cosmosim::physics::StarFormationCellInput cell;
  cell.gas_mass_code = 1.0;
  cell.gas_density_code = 5.2;
  cell.gas_temperature_k = 9.0e3;
  cell.velocity_divergence_code = -0.3;
  assert(model.isEligible(cell));

  cell.gas_temperature_k = 2.0e4;
  assert(!model.isEligible(cell));
}

void testPhysicalFrameConversions() {
  auto config = makeAdaptiveConfig();
  config.density_is_comoving = true;
  config.geometry_is_comoving = true;
  cosmosim::physics::StarFormationModel model(config);
  assert(nearlyEqual(model.physicalDensityCode(8.0, 0.5), 64.0));
  assert(nearlyEqual(model.physicalCellScaleCode(8.0, 0.5), 1.0));
  assert(model.physicalDensityCode(8.0, 0.0) == 0.0);
}

void testVirialJeansAndFlowEligibility() {
  auto config = makeAdaptiveConfig();
  cosmosim::physics::StarFormationModel model(config);
  auto cell = makeEligibleAdaptiveCell();

  const double dx = std::cbrt(cell.cell_volume_code);
  const double expected_alpha =
      (2.0 * std::pow(cell.gas_sound_speed_code / dx, 2)) /
      (8.0 * std::numbers::pi * config.newton_g_code * cell.gas_density_code);
  assert(nearlyEqual(
      model.virialParameter(
          cell.gas_density_code,
          dx,
          cell.gas_sound_speed_code,
          cell.velocity_gradient_frobenius_sq_code),
      expected_alpha));
  assert(model.jeansMassCode(cell.gas_density_code, cell.gas_sound_speed_code) > 0.0);
  assert(model.isEligible(cell));

  cell.velocity_divergence_code = 0.0;
  assert(!model.isEligible(cell));
  cell = makeEligibleAdaptiveCell();
  cell.velocity_gradient_frobenius_sq_code = 1.0e6;
  assert(!model.isEligible(cell));
  cell = makeEligibleAdaptiveCell();
  cell.gas_sound_speed_code = 100.0;
  assert(!model.isEligible(cell));
  cell = makeEligibleAdaptiveCell();
  cell.gas_cell_id = 0;
  assert(!model.isEligible(cell));
}

void testRateExactConversionAndTimeStepLimit() {
  auto config = makeAdaptiveConfig();
  config.epsilon_ff = 0.2;
  config.max_fractional_mass_conversion = 0.9;
  config.min_remaining_gas_fraction = 0.0;
  config.stochastic_spawning = false;
  cosmosim::physics::StarFormationModel model(config);
  auto cell = makeEligibleAdaptiveCell();
  cell.velocity_divergence_code = -1.0e-6;

  const double dt = 0.01;
  const double t_ff = model.freeFallTimeCode(cell.gas_density_code);
  const double expected_mass = cell.gas_mass_code *
      (-std::expm1(-config.epsilon_ff * dt / t_ff));
  const auto outcome = model.sampleCellOutcome(cell, dt, 99, 0, 1.0);
  assert(nearlyEqual(outcome.expected_spawn_mass_code, expected_mass));

  const double expected_limit =
      -t_ff * std::log1p(-config.max_fractional_mass_conversion) / config.epsilon_ff;
  assert(nearlyEqual(model.sourceTimeStepLimitCode(t_ff), expected_limit));

  const double tiny_dt = 1.0e-14;
  const auto tiny = model.sampleCellOutcome(cell, tiny_dt, 99, 0, 1.0);
  const double tiny_expected = cell.gas_mass_code *
      (-std::expm1(-config.epsilon_ff * tiny_dt / t_ff));
  assert(nearlyEqual(tiny.expected_spawn_mass_code, tiny_expected, 1.0e-15));
}

void testCounterRngGoldenVectorsAndRankInvariance() {
  const double u0 = cosmosim::physics::starFormationUniform01(123456789ULL, 42ULL, 99ULL, 0U);
  const double u1 = cosmosim::physics::starFormationUniform01(123456789ULL, 42ULL, 99ULL, 1U);
  assert(nearlyEqual(u0, 0.8626457080226616, 1.0e-15));
  assert(nearlyEqual(u1, 0.19965097066605675, 1.0e-15));
  assert(cosmosim::physics::starFormationBirthKey(42ULL, 99ULL, 0U) ==
         0x7572a540acdbf584ULL);
  assert(cosmosim::physics::starFormationBirthKey(42ULL, 99ULL, 1U) ==
         0x4d4f92f0b00c40d8ULL);

  const auto config = makeAdaptiveConfig();
  cosmosim::physics::StarFormationModel model(config);
  const auto cell = makeEligibleAdaptiveCell();
  const auto rank0 = model.sampleCellOutcome(cell, 0.01, 99, 0, 1.0);
  const auto rank17 = model.sampleCellOutcome(cell, 0.01, 99, 17, 1.0);
  assert(rank0.random_u01 == rank17.random_u01);
  assert(rank0.spawned_particle_count == rank17.spawned_particle_count);
  assert(rank0.spawned_mass_code == rank17.spawned_mass_code);
}


void testUnbiasedSamplingNearMassCap() {
  auto config = makeAdaptiveConfig();
  config.target_star_particle_mass_code = 1.0;
  config.min_star_particle_mass_code = 0.1;
  config.max_fractional_mass_conversion = 0.25;
  config.min_remaining_gas_fraction = 0.0;
  cosmosim::physics::StarFormationModel model(config);
  auto cell = makeEligibleAdaptiveCell();

  // A very long step saturates the 2.5-code-mass conversion cap. The old
  // post-draw clipping law produced a biased mean for target mass 1.0. The
  // feasible-target policy represents the cap exactly without exceeding it.
  const auto capped = model.sampleCellOutcome(cell, 100.0, 17, 0, 1.0);
  assert(nearlyEqual(capped.expected_spawn_mass_code, 2.5));
  assert(nearlyEqual(capped.spawned_mass_code, 2.5));
  assert(capped.spawned_particle_count == 3);
  assert(nearlyEqual(capped.target_particle_mass_code, 2.5 / 3.0));

  // In the low-expectation Bernoulli regime, the ensemble mean must recover the
  // exact continuous expectation without a seed/rank offset bias.
  const double dt = 1.0e-4;
  const auto reference = model.sampleCellOutcome(cell, dt, 0, 0, 1.0);
  double realized_sum = 0.0;
  constexpr std::uint64_t sample_count = 200000;
  for (std::uint64_t tick = 0; tick < sample_count; ++tick) {
    realized_sum += model.sampleCellOutcome(cell, dt, tick, 99, 1.0).spawned_mass_code;
  }
  const double realized_mean = realized_sum / static_cast<double>(sample_count);
  const double standard_error = std::sqrt(
      reference.spawn_probability * (1.0 - reference.spawn_probability) /
      static_cast<double>(sample_count)) * reference.target_particle_mass_code;
  assert(std::abs(realized_mean - reference.expected_spawn_mass_code) < 6.0 * standard_error);
}

void testConservativeMultiParticleBirth() {
  auto config = makeAdaptiveConfig();
  cosmosim::physics::StarFormationModel model(config);
  const auto cell = makeEligibleAdaptiveCell();
  cosmosim::core::SimulationState state;
  initializeAdaptiveState(state, std::span<const cosmosim::physics::StarFormationCellInput>(&cell, 1));

  const auto report = model.applyFromInputs(state, {&cell, 1}, 1.0, 1.0, 77);
  assert(report.counters.spawn_events == 1);
  assert(report.counters.spawned_particles == 5);
  assert(state.particles.size() == 5);
  assert(state.star_particles.size() == 5);
  assert(nearlyEqual(state.cells.mass_code[0], 7.5));
  assert(nearlyEqual(state.gas_cells.density_code[0], 75.0));
  assert(nearlyEqual(state.gas_cells.metal_mass_code[0], 0.15));
  assert(nearlyEqual(report.counters.mass_residual_code, 0.0));
  assert(nearlyEqual(report.counters.momentum_residual_norm_code, 0.0));
  assert(nearlyEqual(report.counters.metal_mass_residual_code, 0.0));
  assert(nearlyEqual(report.counters.gas_internal_energy_removed_code, 5.0));
  assert(nearlyEqual(report.counters.star_formation_internal_energy_sink_code, 5.0));
  assert(nearlyEqual(report.counters.star_kinetic_energy_created_code, 36.25));

  double star_mass = 0.0;
  double stellar_metal_mass = 0.0;
  for (std::size_t star_row = 0; star_row < state.star_particles.size(); ++star_row) {
    const std::size_t particle_index = state.star_particles.particle_index[star_row];
    star_mass += state.particles.mass_code[particle_index];
    stellar_metal_mass += state.particles.mass_code[particle_index] *
        state.star_particles.metallicity_mass_fraction[star_row];
    assert(state.particles.velocity_x_peculiar[particle_index] == 2.0);
    assert(state.particles.velocity_y_peculiar[particle_index] == -3.0);
    assert(state.particles.velocity_z_peculiar[particle_index] == 4.0);
    assert(state.particles.position_x_comoving[particle_index] == 1.0);
    assert(state.star_particles.parent_gas_cell_id[star_row] == 42);
    assert(state.star_particles.birth_tick[star_row] == 77);
    assert(state.star_particles.birth_ordinal[star_row] == star_row);
  }
  assert(nearlyEqual(star_mass, 2.5));
  assert(nearlyEqual(stellar_metal_mass, 0.05));
  assert(nearlyEqual(state.cells.mass_code[0] + star_mass, 10.0));
  assert(state.validateUniqueParticleIds());

  const auto* sidecar = state.sidecars.find("star_formation");
  assert(sidecar != nullptr);
  const std::string payload(reinterpret_cast<const char*>(sidecar->payload.data()), sidecar->payload.size());
  assert(payload.find("model_name=adaptive_bound_jeans") != std::string::npos);
  assert(payload.find("spawned_particles=5") != std::string::npos);
}

void testDenseRowReorderInvariance() {
  auto config = makeAdaptiveConfig();
  config.target_star_particle_mass_code = 2.5;
  cosmosim::physics::StarFormationModel model(config);

  std::array<cosmosim::physics::StarFormationCellInput, 2> inputs_a{
      makeEligibleAdaptiveCell(0, 1001), makeEligibleAdaptiveCell(1, 1002)};
  inputs_a[1].center_x_comoving = 9.0;
  std::array<cosmosim::physics::StarFormationCellInput, 2> inputs_b{
      makeEligibleAdaptiveCell(1, 1001), makeEligibleAdaptiveCell(0, 1002)};
  inputs_b[0].center_x_comoving = inputs_a[0].center_x_comoving;
  inputs_b[1].center_x_comoving = inputs_a[1].center_x_comoving;

  cosmosim::core::SimulationState state_a;
  cosmosim::core::SimulationState state_b;
  initializeAdaptiveState(state_a, inputs_a);
  initializeAdaptiveState(state_b, inputs_b);
  const auto report_a = model.applyFromInputs(state_a, inputs_a, 1.0, 1.0, 1234);
  const auto report_b = model.applyFromInputs(state_b, inputs_b, 1.0, 1.0, 1234);
  assert(report_a.birth_keys == report_b.birth_keys);
  assert(report_a.counters.spawned_particles == report_b.counters.spawned_particles);
  assert(report_a.counters.spawned_mass_code == report_b.counters.spawned_mass_code);
  assert(state_a.particle_sidecar.particle_id == state_b.particle_sidecar.particle_id);
}

void testInvalidAndDisabledNoOp() {
  auto config = makeAdaptiveConfig();
  cosmosim::physics::StarFormationModel model(config);
  auto invalid = makeEligibleAdaptiveCell();
  invalid.gas_density_code = std::numeric_limits<double>::quiet_NaN();
  assert(!model.isEligible(invalid));
  invalid = makeEligibleAdaptiveCell();
  invalid.cell_volume_code = -1.0;
  assert(!model.isEligible(invalid));

  config.enabled = false;
  cosmosim::physics::StarFormationModel disabled(config);
  const auto cell = makeEligibleAdaptiveCell();
  cosmosim::core::SimulationState state;
  initializeAdaptiveState(state, {&cell, 1});
  const double gas_mass_before = state.cells.mass_code[0];
  const auto report = disabled.applyFromInputs(state, {&cell, 1}, 1.0, 1.0, 1);
  assert(report.counters.scanned_cells == 0);
  assert(state.cells.mass_code[0] == gas_mass_before);
  assert(state.particles.size() == 0);
  assert(state.star_particles.size() == 0);
  assert(state.sidecars.find("star_formation") == nullptr);
}

void testLegacyTimeIntegrationCallbackHook() {
  cosmosim::core::SimulationState state;
  state.resizeCells(1);
  state.cells.mass_code[0] = 2.0;
  state.gas_cells.density_code[0] = 20.0;
  state.gas_cells.temperature_code[0] = 5.0e3;

  cosmosim::physics::StarFormationConfig config;
  config.stochastic_spawning = false;
  config.epsilon_ff = 0.1;
  cosmosim::physics::StarFormationCallback callback{cosmosim::physics::StarFormationModel(config)};
  callback.setVelocityDivergenceCode(std::array<double, 1>{-0.5});
  callback.setMetallicityMassFraction(std::array<double, 1>{0.02});

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 1.0e8;
  const std::array<std::uint32_t, 1> active_cells{0};
  cosmosim::core::ActiveSetDescriptor active_set{
      .cell_indices = active_cells,
      .cells_are_subset = true,
  };
  cosmosim::core::StepContext context{
      .state = state,
      .integrator_state = integrator_state,
      .active_set = active_set,
      .stage = cosmosim::core::IntegrationStage::kSourceTerms,
  };

  callback.onStage(context);
  assert(callback.lastStepReport().counters.spawn_events == 1);
  assert(state.particles.size() == 1);
}

}  // namespace

int main() {
  testLegacyThresholdEligibility();
  testPhysicalFrameConversions();
  testVirialJeansAndFlowEligibility();
  testRateExactConversionAndTimeStepLimit();
  testCounterRngGoldenVectorsAndRankInvariance();
  testUnbiasedSamplingNearMassCap();
  testConservativeMultiParticleBirth();
  testDenseRowReorderInvariance();
  testInvalidAndDisabledNoOp();
  testLegacyTimeIntegrationCallbackHook();
  return 0;
}
