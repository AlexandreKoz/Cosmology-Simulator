#include <cassert>
#include <cmath>
#include <filesystem>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "restart_equivalence_harness.hpp"
#include "restart_equivalence_scenarios.hpp"

namespace {

cosmosim::hydro::HydroPatchGeometry makePeriodic1dGeometry(std::size_t cell_count) {
  cosmosim::hydro::HydroPatchGeometry geometry;
  geometry.cell_volume_comoving = 1.0;
  geometry.faces.reserve(cell_count);
  for (std::size_t i = 0; i < cell_count; ++i) {
    geometry.faces.push_back(cosmosim::hydro::HydroFace{
        .owner_cell = i,
        .neighbor_cell = (i + 1U) % cell_count,
        .area_comoving = 1.0,
        .normal_x = 1.0,
        .normal_y = 0.0,
        .normal_z = 0.0});
  }
  return geometry;
}

void applyProductionHydroStep(cosmosim::tests::RestartEquivalenceStepContext& context) {
  constexpr double gamma = 1.4;
  cosmosim::hydro::HydroConservedStateSoa conserved(context.state.cells.size());
  for (std::size_t cidx = 0; cidx < context.state.cells.size(); ++cidx) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = context.state.gas_cells.density_code[cidx];
    primitive.vel_x_peculiar = context.state.particles.velocity_x_peculiar[cidx];
    primitive.vel_y_peculiar = context.state.particles.velocity_y_peculiar[cidx];
    primitive.vel_z_peculiar = context.state.particles.velocity_z_peculiar[cidx];
    primitive.pressure_comoving = context.state.gas_cells.pressure_code[cidx];
    conserved.storeCell(cidx, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }

  const auto geometry = makePeriodic1dGeometry(context.state.cells.size());
  const cosmosim::hydro::HydroUpdateContext update{
      .dt_code = 1.0e-3,
      .scale_factor = context.integrator_state.current_scale_factor,
      .hubble_rate_code = 0.0,
  };
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::HydroCoreSolver solver(gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code,
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = true});
  cosmosim::hydro::HllcRiemannSolver riemann_solver;
  solver.advancePatch(conserved, geometry, update, reconstruction, riemann_solver, {}, source_context, nullptr);

  for (std::size_t cidx = 0; cidx < context.state.cells.size(); ++cidx) {
    const auto primitive = cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(cidx), gamma);
    context.state.gas_cells.density_code[cidx] = primitive.rho_comoving;
    context.state.gas_cells.pressure_code[cidx] = primitive.pressure_comoving;
    context.state.cells.mass_code[cidx] = primitive.rho_comoving * geometry.cell_volume_comoving;
    context.state.particles.mass_code[cidx] = context.state.cells.mass_code[cidx];
    context.state.particles.velocity_x_peculiar[cidx] = primitive.vel_x_peculiar;
    context.state.particles.velocity_y_peculiar[cidx] = primitive.vel_y_peculiar;
    context.state.particles.velocity_z_peculiar[cidx] = primitive.vel_z_peculiar;
    context.state.gas_cells.internal_energy_code[cidx] =
        primitive.pressure_comoving / std::max((gamma - 1.0) * primitive.rho_comoving, 1.0e-30);
    context.state.gas_cells.temperature_code[cidx] =
        primitive.pressure_comoving / std::max(primitive.rho_comoving, 1.0e-30);
    context.state.gas_cells.sound_speed_code[cidx] =
        std::sqrt(gamma * primitive.pressure_comoving / std::max(primitive.rho_comoving, 1.0e-30));
  }
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_HDF5
  const auto restart_path = cosmosim::tests::stage8RestartPath("restart_equivalence_hydro_toy");
  auto state = cosmosim::tests::makeStage8HydroToyState(24, "restart_equivalence_hydro_toy");
  for (std::size_t cidx = 0; cidx < state.cells.size(); ++cidx) {
    const bool left = cidx < state.cells.size() / 2U;
    state.gas_cells.density_code[cidx] = left ? 1.0 : 0.125;
    state.gas_cells.pressure_code[cidx] = left ? 1.0 : 0.1;
    state.gas_cells.temperature_code[cidx] = state.gas_cells.pressure_code[cidx] / state.gas_cells.density_code[cidx];
    state.gas_cells.internal_energy_code[cidx] = state.gas_cells.pressure_code[cidx] /
        ((1.4 - 1.0) * state.gas_cells.density_code[cidx]);
    state.gas_cells.sound_speed_code[cidx] = std::sqrt(1.4 * state.gas_cells.pressure_code[cidx] / state.gas_cells.density_code[cidx]);
    state.cells.mass_code[cidx] = state.gas_cells.density_code[cidx];
    state.particles.mass_code[cidx] = state.cells.mass_code[cidx];
  }
  auto scheduler = cosmosim::tests::makeStage8Scheduler(static_cast<std::uint32_t>(state.particles.size()), 2);
  auto integrator_state = cosmosim::tests::makeStage8IntegratorState(2, 2);
  auto output_state = cosmosim::tests::makeStage8OutputCadenceState(false);
  auto scenario = cosmosim::tests::makeStage8Scenario(
      std::move(state), integrator_state, std::move(scheduler), std::move(output_state), restart_path, 30, 12);
  scenario.step_kernel = applyProductionHydroStep;
  scenario.tolerances.position_abs = 1.0e-12;
  scenario.tolerances.velocity_abs = 1.0e-12;
  scenario.tolerances.scalar_abs = 1.0e-12;
  const auto result = cosmosim::tests::runRestartEquivalenceScenario(std::move(scenario));
  assert(result.direct_state.gas_cells.density_code[11U] > 0.2);
  assert(result.direct_state.gas_cells.density_code[12U] < 0.9);
  std::filesystem::remove(restart_path);
#endif
  return 0;
}
