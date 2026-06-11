#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace {

bool throwsWithContext(const std::function<void()>& action, const std::string& required_context) {
  try {
    action();
  } catch (const std::exception& ex) {
    return std::string(ex.what()).find(required_context) != std::string::npos;
  }
  return false;
}

void testMappingRegressionReference() {
  const cosmosim::core::TimeStepLimits limits{
      .min_dt_time_code = 0.0625,
      .max_dt_time_code = 1.0,
      .max_bin = 4,
  };

  const auto coarse = cosmosim::core::mapDtToTimeBin(0.999, limits);
  const auto medium = cosmosim::core::mapDtToTimeBin(0.2, limits);
  const auto fine = cosmosim::core::mapDtToTimeBin(0.0625, limits);

  assert(coarse.bin_index == 3);
  assert(medium.bin_index == 1);
  assert(fine.bin_index == 0);

  assert(cosmosim::core::binIndexToDt(coarse.bin_index, limits) == 0.5);
  assert(cosmosim::core::binIndexToDt(medium.bin_index, limits) == 0.125);
}

void testHydroCflRejectsTooLargeAcceptedDtBeforeMutation() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(1);
  state.particle_sidecar.particle_id[0] = 1001;
  state.particle_sidecar.species_tag[0] =
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
  state.species.count_by_species = {};
  state.species.count_by_species[
      static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)] = 1;
  state.rebuildSpeciesIndex();
  state.refreshGasCellIdentityFromParticleOrder();
  state.gas_cells.density_code[0] = 1.0;
  state.gas_cells.pressure_code[0] = 1.0;
  state.gas_cells.internal_energy_code[0] = 1.5;
  state.gas_cells.sound_speed_code[0] = 1.0;
  state.cells.mass_code[0] = 1.0;
  state.particles.mass_code[0] = 1.0;
  state.particles.velocity_x_peculiar[0] = 100.0;
  state.particles.velocity_y_peculiar[0] = 0.0;
  state.particles.velocity_z_peculiar[0] = 0.0;

  const double density_before = state.gas_cells.density_code[0];
  const double pressure_before = state.gas_cells.pressure_code[0];
  const double mass_before = state.cells.mass_code[0];
  const cosmosim::core::DirectionalCflTimeStepInput cfl_input{
      .cell_width_axis_code = {0.25, 0.25, 0.25},
      .velocity_axis_code = {
          state.particles.velocity_x_peculiar[0],
          state.particles.velocity_y_peculiar[0],
          state.particles.velocity_z_peculiar[0]},
      .sound_speed_code = state.gas_cells.sound_speed_code[0],
  };
  const double accepted_dt_time_code = 0.1;
  const auto diagnostics = cosmosim::core::makeHydroCflDiagnostics(
      0,
      cfl_input,
      0.4,
      accepted_dt_time_code,
      state.gas_cells.gas_cell_id[0],
      std::nullopt,
      std::nullopt);
  assert(diagnostics.proposed_dt_time_code < accepted_dt_time_code);
  assert(throwsWithContext(
      [&]() { cosmosim::core::assertHydroCflStable(diagnostics); },
      "hydro CFL violation"));

  assert(state.gas_cells.density_code[0] == density_before);
  assert(state.gas_cells.pressure_code[0] == pressure_before);
  assert(state.cells.mass_code[0] == mass_before);
}

}  // namespace

int main() {
  testMappingRegressionReference();
  testHydroCflRejectsTooLargeAcceptedDtBeforeMutation();
  return 0;
}
