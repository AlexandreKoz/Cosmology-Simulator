#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"

namespace {

void seedSingleStarCase(cosmosim::core::SimulationState& state, std::size_t cell_count = 4U) {
  state.resizeCells(cell_count);
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.center_y_comoving[i] = 0.0;
    state.cells.center_z_comoving[i] = 0.0;
    state.cells.mass_code[i] = 10.0;
    state.gas_cells.gas_cell_id[i] = 100U + i;
    state.gas_cells.density_code[i] = 5.0;  // volume = 2
    state.gas_cells.internal_energy_code[i] = 2.0;
    state.gas_cells.metal_mass_code[i] = 0.1;
  }

  state.resizeParticles(1);
  state.star_particles.resize(1);
  state.particles.position_x_comoving[0] = 0.1;
  state.particles.position_y_comoving[0] = 0.0;
  state.particles.position_z_comoving[0] = 0.0;
  state.star_particles.particle_index[0] = 0;
  state.star_particles.birth_mass_code[0] = 1.0;
}

cosmosim::physics::StellarFeedbackConfig baseConfig() {
  cosmosim::physics::StellarFeedbackConfig config;
  config.variant = cosmosim::physics::StellarFeedbackVariant::kNone;
  config.mode = cosmosim::physics::StellarFeedbackMode::kThermal;
  config.epsilon_thermal = 1.0;
  config.epsilon_kinetic = 0.0;
  config.epsilon_momentum = 0.0;
  config.neighbor_count = 1;
  return config;
}

cosmosim::physics::StellarFeedbackGeometryView geometryView(
    const cosmosim::core::SimulationState& state) {
  return cosmosim::physics::StellarFeedbackGeometryView{
      .particle_position_x_comoving = state.particles.position_x_comoving,
      .particle_position_y_comoving = state.particles.position_y_comoving,
      .particle_position_z_comoving = state.particles.position_z_comoving,
      .cell_center_x_comoving = state.cells.center_x_comoving,
      .cell_center_y_comoving = state.cells.center_y_comoving,
      .cell_center_z_comoving = state.cells.center_z_comoving,
      .gas_cell_id = state.gas_cells.gas_cell_id,
  };
}

cosmosim::physics::StellarFeedbackDepositionView depositionView(
    cosmosim::core::SimulationState& state) {
  return cosmosim::physics::StellarFeedbackDepositionView{
      .cell_mass_code = state.cells.mass_code,
      .gas_density_code = state.gas_cells.density_code,
      .gas_internal_energy_code = state.gas_cells.internal_energy_code,
      .gas_metal_mass_code = state.gas_cells.metal_mass_code,
  };
}

void testMassDensityEnergyAndMetalConservation() {
  cosmosim::core::SimulationState state;
  seedSingleStarCase(state);
  auto config = baseConfig();
  config.total_energy_code_per_erg = 2.0;
  cosmosim::physics::StellarFeedbackModel model(config);
  cosmosim::physics::StellarFeedbackModuleState module_state;

  const std::vector<std::uint32_t> active = {0};
  const std::vector<double> returned_mass = {0.5};
  const std::vector<double> returned_metals = {0.05};
  const std::vector<double> energy = {4.0};
  const auto report = model.apply(
      state, module_state, active, returned_mass, returned_metals, 1.0, energy);

  assert(report.counters.feedback_stars == 1);
  assert(std::abs(report.counters.deposited_mass_code - 0.5) < 1.0e-12);
  assert(std::abs(report.counters.deposited_metals_code - 0.05) < 1.0e-12);
  assert(std::abs(state.cells.mass_code[0] - 10.5) < 1.0e-12);
  assert(std::abs(state.gas_cells.density_code[0] - 5.25) < 1.0e-12);
  assert(std::abs(state.gas_cells.metal_mass_code[0] - 0.15) < 1.0e-12);
  const double expected_u = (10.0 * 2.0 + 4.0 * 2.0) / 10.5;
  assert(std::abs(state.gas_cells.internal_energy_code[0] - expected_u) < 1.0e-12);
  assert(state.star_particles.enrichment_carry_mass_code[0] == 0.0);
  assert(state.star_particles.stellar_deposited_mass_cumulative_code[0] == 0.5);
}

void testNoNeighborBudgetIsDurablyCarried() {
  cosmosim::core::SimulationState state;
  seedSingleStarCase(state, 0U);
  cosmosim::physics::StellarFeedbackModel model(baseConfig());
  cosmosim::physics::StellarFeedbackModuleState module_state;
  const std::vector<std::uint32_t> active = {0};
  const std::vector<double> returned_mass = {0.25};
  const std::vector<double> returned_metals = {0.02};
  const std::vector<double> energy = {3.0};
  const auto report = model.apply(
      state, module_state, active, returned_mass, returned_metals, 1.0, energy);
  assert(report.counters.deposited_mass_code == 0.0);
  assert(std::abs(state.star_particles.enrichment_carry_mass_code[0] - 0.25) < 1.0e-14);
  assert(std::abs(state.star_particles.enrichment_carry_metals_code[0] - 0.02) < 1.0e-14);
  assert(std::abs(state.star_particles.enrichment_carry_feedback_energy_erg[0] - 3.0) < 1.0e-14);
}

void testModeSelectionMomentumOnly() {
  auto config = baseConfig();
  config.mode = cosmosim::physics::StellarFeedbackMode::kMomentum;
  config.epsilon_thermal = 0.0;
  config.epsilon_momentum = 1.0;
  cosmosim::physics::StellarFeedbackModel model(config);
  const auto budget = model.computeBudget(1.0, 0.1, 0.01);
  assert(budget.thermal_energy_erg == 0.0);
  assert(budget.kinetic_energy_erg == 0.0);
  assert(budget.momentum_budget_code > 0.0);
}

void testRuntimeCarryMirrorIsEliminated() {
  cosmosim::physics::StellarFeedbackModuleState module_state;
  module_state.ensureStarCapacity(100000000U);
  assert(module_state.ownedCapacityBytes() == 0U);
}

void testSpatialIndexMatchesExactNearestSelection() {
  cosmosim::core::SimulationState state;
  seedSingleStarCase(state, 257U);
  state.particles.position_x_comoving[0] = 4.375;
  state.particles.position_y_comoving[0] = -1.25;
  state.particles.position_z_comoving[0] = 0.75;
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    const std::uint32_t ii = static_cast<std::uint32_t>(i);
    state.cells.center_x_comoving[i] =
        0.125 * static_cast<double>((ii * 73U) % 257U);
    state.cells.center_y_comoving[i] =
        0.075 * static_cast<double>((ii * 29U) % 43U) - 1.6;
    state.cells.center_z_comoving[i] =
        0.05 * static_cast<double>((ii * 17U) % 31U) - 0.7;
    state.gas_cells.gas_cell_id[i] = 100000U + (257U - ii);
  }
  auto config = baseConfig();
  config.neighbor_count = 11U;
  cosmosim::physics::StellarFeedbackModel model(config);
  const auto geometry = geometryView(state);
  cosmosim::physics::StellarFeedbackSpatialIndex spatial_index;
  spatial_index.rebuild(geometry);

  const auto direct = model.selectTargets(geometry, 0U);
  const auto indexed = model.selectTargets(geometry, spatial_index, 0U);
  assert(direct.size() == indexed.size());
  for (std::size_t i = 0; i < direct.size(); ++i) {
    assert(direct[i].cell_index == indexed[i].cell_index);
    assert(std::abs(direct[i].weight - indexed[i].weight) < 1.0e-15);
    assert(std::abs(direct[i].radial_dx_comoving - indexed[i].radial_dx_comoving) < 1.0e-15);
    assert(std::abs(direct[i].radial_dy_comoving - indexed[i].radial_dy_comoving) < 1.0e-15);
    assert(std::abs(direct[i].radial_dz_comoving - indexed[i].radial_dz_comoving) < 1.0e-15);
  }
  assert(spatial_index.ownedCapacityBytes() <=
         state.cells.size() * sizeof(std::uint32_t));
  assert(spatial_index.highWaterBytes() == spatial_index.ownedCapacityBytes());
}

void testEventBatchingMatchesSequentialApplicationWithoutHistoryGrowth() {
  cosmosim::core::SimulationState batched;
  cosmosim::core::SimulationState sequential;
  seedSingleStarCase(batched, 16U);
  seedSingleStarCase(sequential, 16U);
  for (std::size_t i = 0; i < batched.cells.size(); ++i) {
    const double y = 0.1 * static_cast<double>(i % 3U);
    const double z = -0.05 * static_cast<double>(i % 5U);
    batched.cells.center_y_comoving[i] = y;
    batched.cells.center_z_comoving[i] = z;
    sequential.cells.center_y_comoving[i] = y;
    sequential.cells.center_z_comoving[i] = z;
  }

  auto config = baseConfig();
  config.neighbor_count = 4U;
  config.total_energy_code_per_erg = 1.5;
  cosmosim::physics::StellarFeedbackModel model(config);
  cosmosim::physics::StellarFeedbackModuleState batched_state;
  cosmosim::physics::StellarFeedbackModuleState sequential_state;
  auto batched_geometry = geometryView(batched);
  auto sequential_geometry = geometryView(sequential);
  cosmosim::physics::StellarFeedbackSpatialIndex batched_index;
  cosmosim::physics::StellarFeedbackSpatialIndex sequential_index;
  batched_index.rebuild(batched_geometry);
  sequential_index.rebuild(sequential_geometry);
  const auto initial_index_bytes = batched_index.ownedCapacityBytes();

  const std::vector<cosmosim::physics::StellarFeedbackEvent> events = {
      {.star_index = 0U,
       .returned_mass_code = 0.2,
       .returned_metals_code = 0.02,
       .feedback_energy_erg = 1.0},
      {.star_index = 0U,
       .returned_mass_code = 0.3,
       .returned_metals_code = 0.03,
       .feedback_energy_erg = 2.0},
  };
  const auto batched_report = model.applyEventsWithViews(
      batched,
      batched_state,
      batched_geometry,
      &batched_index,
      depositionView(batched),
      events,
      1.0);
  cosmosim::physics::StellarFeedbackStepCounters sequential_counters;
  for (const auto& event : events) {
    const std::array<cosmosim::physics::StellarFeedbackEvent, 1> one_event{event};
    const auto report = model.applyEventsWithViews(
        sequential,
        sequential_state,
        sequential_geometry,
        &sequential_index,
        depositionView(sequential),
        one_event,
        1.0);
    sequential_counters.feedback_stars += report.counters.feedback_stars;
    sequential_counters.target_cells_visited += report.counters.target_cells_visited;
    sequential_counters.deposited_mass_code += report.counters.deposited_mass_code;
    sequential_counters.deposited_metals_code += report.counters.deposited_metals_code;
    sequential_counters.deposited_thermal_energy_erg +=
        report.counters.deposited_thermal_energy_erg;
  }

  for (std::size_t i = 0; i < batched.cells.size(); ++i) {
    assert(std::abs(batched.cells.mass_code[i] - sequential.cells.mass_code[i]) < 1.0e-14);
    assert(std::abs(batched.gas_cells.density_code[i] -
                    sequential.gas_cells.density_code[i]) < 1.0e-14);
    assert(std::abs(batched.gas_cells.internal_energy_code[i] -
                    sequential.gas_cells.internal_energy_code[i]) < 1.0e-14);
    assert(std::abs(batched.gas_cells.metal_mass_code[i] -
                    sequential.gas_cells.metal_mass_code[i]) < 1.0e-14);
  }
  assert(std::abs(
             batched.star_particles.stellar_deposited_mass_cumulative_code[0] -
             sequential.star_particles.stellar_deposited_mass_cumulative_code[0]) < 1.0e-14);
  assert(std::abs(batched_report.counters.deposited_mass_code -
                  sequential_counters.deposited_mass_code) < 1.0e-14);
  assert(std::abs(batched_report.counters.deposited_metals_code -
                  sequential_counters.deposited_metals_code) < 1.0e-14);
  assert(std::abs(batched_report.counters.deposited_thermal_energy_erg -
                  sequential_counters.deposited_thermal_energy_erg) < 1.0e-14);
  assert(batched_state.ownedCapacityBytes() == 0U);
  assert(batched_index.ownedCapacityBytes() == initial_index_bytes);
  assert(batched_index.highWaterBytes() == initial_index_bytes);
}

}  // namespace

int main() {
  testMassDensityEnergyAndMetalConservation();
  testNoNeighborBudgetIsDurablyCarried();
  testModeSelectionMomentumOnly();
  testRuntimeCarryMirrorIsEliminated();
  testSpatialIndexMatchesExactNearestSelection();
  testEventBatchingMatchesSequentialApplicationWithoutHistoryGrowth();
  return 0;
}
