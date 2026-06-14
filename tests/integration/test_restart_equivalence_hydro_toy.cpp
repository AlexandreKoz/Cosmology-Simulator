#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "restart_equivalence_harness.hpp"
#include "restart_equivalence_scenarios.hpp"

namespace {

constexpr std::size_t k_nx = 4;
constexpr std::size_t k_ny = 3;
constexpr std::size_t k_nz = 2;
constexpr double k_dx_code = 0.25;
constexpr double k_dy_code = 0.5;
constexpr double k_dz_code = 0.75;

cosmosim::hydro::HydroCartesianPatchSpec makeRestartHydroPatchSpec() {
  return cosmosim::hydro::HydroCartesianPatchSpec{
      .nx = k_nx,
      .ny = k_ny,
      .nz = k_nz,
      .origin_x_comoving = 0.125,
      .origin_y_comoving = 0.25,
      .origin_z_comoving = 0.375,
      .cell_width_x_comoving = k_dx_code,
      .cell_width_y_comoving = k_dy_code,
      .cell_width_z_comoving = k_dz_code};
}

cosmosim::core::SimulationState makeH1CartesianHydroState(const std::string& run_name) {
  const auto spec = makeRestartHydroPatchSpec();
  const auto geometry = cosmosim::hydro::makeCartesianPatchGeometry(spec);
  auto state = cosmosim::tests::makeStage8HydroToyState(geometry.cellCount(), run_name);
  state.resizePatches(1);
  state.patches.patch_id[0] = 4242;
  state.patches.level[0] = 1;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = static_cast<std::uint32_t>(geometry.cellCount());
  state.patches.owning_rank[0] = 0;

  for (std::size_t cidx = 0; cidx < geometry.cellCount(); ++cidx) {
    const auto ijk = geometry.cellIjk(cidx);
    const double x = spec.origin_x_comoving + (static_cast<double>(ijk[0]) + 0.5) * spec.cell_width_x_comoving;
    const double y = spec.origin_y_comoving + (static_cast<double>(ijk[1]) + 0.5) * spec.cell_width_y_comoving;
    const double z = spec.origin_z_comoving + (static_cast<double>(ijk[2]) + 0.5) * spec.cell_width_z_comoving;
    const double pattern =
        0.03 * static_cast<double>(ijk[0]) - 0.02 * static_cast<double>(ijk[1]) +
        0.015 * static_cast<double>(ijk[2]);
    state.cells.center_x_comoving[cidx] = x;
    state.cells.center_y_comoving[cidx] = y;
    state.cells.center_z_comoving[cidx] = z;
    state.cells.patch_index[cidx] = 0;
    state.particles.position_x_comoving[cidx] = x;
    state.particles.position_y_comoving[cidx] = y;
    state.particles.position_z_comoving[cidx] = z;
    state.particles.velocity_x_peculiar[cidx] = 0.02 + 0.004 * static_cast<double>(ijk[1]);
    state.particles.velocity_y_peculiar[cidx] = -0.015 + 0.003 * static_cast<double>(ijk[2]);
    state.particles.velocity_z_peculiar[cidx] = 0.01 - 0.002 * static_cast<double>(ijk[0]);
    state.gas_cells.velocity_x_peculiar[cidx] = state.particles.velocity_x_peculiar[cidx];
    state.gas_cells.velocity_y_peculiar[cidx] = state.particles.velocity_y_peculiar[cidx];
    state.gas_cells.velocity_z_peculiar[cidx] = state.particles.velocity_z_peculiar[cidx];
    state.gas_cells.density_code[cidx] = 1.0 + pattern;
    state.gas_cells.pressure_code[cidx] = 0.8 + 0.5 * pattern;
    state.gas_cells.temperature_code[cidx] =
        state.gas_cells.pressure_code[cidx] / state.gas_cells.density_code[cidx];
    state.gas_cells.internal_energy_code[cidx] =
        state.gas_cells.pressure_code[cidx] / ((1.4 - 1.0) * state.gas_cells.density_code[cidx]);
    state.gas_cells.sound_speed_code[cidx] =
        std::sqrt(1.4 * state.gas_cells.pressure_code[cidx] / state.gas_cells.density_code[cidx]);
    state.cells.mass_code[cidx] = state.gas_cells.density_code[cidx] * geometry.cell_volume_comoving;
    state.particles.mass_code[cidx] = state.cells.mass_code[cidx];
  }
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(state.cells.size());
  for (std::uint32_t cidx = 0; cidx < state.cells.size(); ++cidx) {
    std::optional<std::uint64_t> parent_particle_id = state.gas_cells.parent_particle_id[cidx];
    if (cidx == 0U) {
      parent_particle_id = std::nullopt;
    } else if (cidx == 2U) {
      parent_particle_id = state.gas_cells.parent_particle_id[1U];
    }
    records.push_back(cosmosim::core::GasCellIdentityRecord{
        .gas_cell_id = state.gas_cells.gas_cell_id[cidx],
        .parent_particle_id = parent_particle_id,
        .owning_patch_id = state.patches.patch_id[0],
        .local_cell_row = cidx,
    });
    state.gas_cells.parent_particle_id[cidx] = parent_particle_id.value_or(0U);
  }
  state.gas_cell_identity.assign(std::move(records));
  assert(state.validateOwnershipInvariants());
  return state;
}

void applyProductionHydroStep(cosmosim::tests::RestartEquivalenceStepContext& context) {
  constexpr double gamma = 1.4;
  cosmosim::hydro::HydroConservedStateSoa conserved(context.state.cells.size());
  for (std::size_t cidx = 0; cidx < context.state.cells.size(); ++cidx) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = context.state.gas_cells.density_code[cidx];
    primitive.vel_x_peculiar = context.state.gas_cells.velocity_x_peculiar[cidx];
    primitive.vel_y_peculiar = context.state.gas_cells.velocity_y_peculiar[cidx];
    primitive.vel_z_peculiar = context.state.gas_cells.velocity_z_peculiar[cidx];
    primitive.pressure_comoving = context.state.gas_cells.pressure_code[cidx];
    conserved.storeCell(cidx, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }

  const auto geometry = cosmosim::hydro::makeCartesianPatchGeometry(makeRestartHydroPatchSpec());
  assert(geometry.cellCount() == context.state.cells.size());
  const cosmosim::hydro::HydroUpdateContext update{
      .dt_code = 1.0e-3,
      .scale_factor = context.integrator_state.current_scale_factor,
      .hubble_rate_code = 0.0,
  };
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::HydroCoreSolver solver(gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code / geometry.cell_width_x_comoving,
      .dt_over_cell_width_code = {
          update.dt_code / geometry.cell_width_x_comoving,
          update.dt_code / geometry.cell_width_y_comoving,
          update.dt_code / geometry.cell_width_z_comoving},
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
    context.state.gas_cells.velocity_x_peculiar[cidx] = primitive.vel_x_peculiar;
    context.state.gas_cells.velocity_y_peculiar[cidx] = primitive.vel_y_peculiar;
    context.state.gas_cells.velocity_z_peculiar[cidx] = primitive.vel_z_peculiar;
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
  auto state = makeH1CartesianHydroState("restart_equivalence_hydro_toy");
  auto scheduler = cosmosim::tests::makeStage8Scheduler(static_cast<std::uint32_t>(state.particles.size()), 2);
  scheduler.setElementBin(2, 1, scheduler.currentTick());
  auto integrator_state = cosmosim::tests::makeStage8IntegratorState(2, 2);
  auto output_state = cosmosim::tests::makeStage8OutputCadenceState(false);
  auto scenario = cosmosim::tests::makeStage8Scenario(
      std::move(state), integrator_state, std::move(scheduler), std::move(output_state), restart_path, 30, 12);
  scenario.step_kernel = applyProductionHydroStep;
  scenario.tolerances.position_abs = 1.0e-12;
  scenario.tolerances.velocity_abs = 1.0e-12;
  scenario.tolerances.scalar_abs = 1.0e-12;
  const auto result = cosmosim::tests::runRestartEquivalenceScenario(std::move(scenario));
  assert(result.direct_state.gas_cells.density_code[5U] != result.direct_state.gas_cells.density_code[6U]);
  assert(result.direct_state.patches.patch_id[0] == 4242);
  assert(result.direct_state.patches.cell_count[0] == k_nx * k_ny * k_nz);
  std::filesystem::remove(restart_path);
#endif
  return 0;
}
