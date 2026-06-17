#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "amr_hydro_validation_helpers.hpp"

namespace {

namespace ahv = cosmosim::tests::amr_hydro_validation;

[[nodiscard]] cosmosim::amr::PatchDescriptor sedovParentPatch() {
  return cosmosim::amr::PatchDescriptor{
      .patch_id = 301,
      .parent_patch_id = 0,
      .level = 0,
      .morton_key = 301,
      .origin_comov = {0.0, 0.0, 0.0},
      .extent_comov = {1.0, 1.0, 1.0},
      .cell_dims = {2, 2, 2}};
}

[[nodiscard]] double radiusFromCenter(double x, double y, double z) {
  const double dx = x - 0.5;
  const double dy = y - 0.5;
  const double dz = z - 0.5;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void imposeSedovDeposit(cosmosim::core::SimulationState& state) {
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const double radius = radiusFromCenter(
        state.cells.center_x_comoving[row],
        state.cells.center_y_comoving[row],
        state.cells.center_z_comoving[row]);
    const double pressure = radius < 0.24 ? 6.0 : 1.0e-3;
    state.gas_cells.density_code[row] = 1.0;
    state.gas_cells.pressure_code[row] = pressure;
    state.gas_cells.internal_energy_code[row] = pressure / ((ahv::k_gamma - 1.0) * state.gas_cells.density_code[row]);
    state.gas_cells.velocity_x_peculiar[row] = 0.0;
    state.gas_cells.velocity_y_peculiar[row] = 0.0;
    state.gas_cells.velocity_z_peculiar[row] = 0.0;
    state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
    state.gas_cells.sound_speed_code[row] =
        std::sqrt(ahv::k_gamma * pressure / state.gas_cells.density_code[row]);
    state.cells.mass_code[row] = state.gas_cells.density_code[row] * (1.0 / 64.0);
  }
}

void testAmrSedovRefinement() {
  const auto parent = sedovParentPatch();
  cosmosim::core::SimulationState state = ahv::makePatchState(
      {parent},
      [](const cosmosim::amr::PatchDescriptor&, std::size_t, const ahv::CellCenter&) {
        return cosmosim::hydro::HydroPrimitiveState{
            .rho_comoving = 1.0,
            .pressure_comoving = 1.0e-3};
      });

  const auto refine = cosmosim::amr::refineProductionPatchInSimulationState(state, parent, 400, 30000);
  ahv::requireOrThrow(refine.refined_patch_count == 1U, "Sedov: parent patch was not refined");
  ahv::requireOrThrow(refine.created_gas_cell_count == 64U, "Sedov: unexpected refined cell count");
  ahv::requireFinitePositiveState(state, "Sedov after refine");
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    ahv::requireOrThrow(state.cells.mass_code[row] > 0.0, "Sedov: zero child mass after refine");
    ahv::requireOrThrow(state.gas_cells.pressure_code[row] > 0.0, "Sedov: zero child pressure after refine");
  }

  imposeSedovDeposit(state);
  const auto before = ahv::totalState(state, cosmosim::amr::buildProductionAmrPatchDescriptors(state));
  const auto diagnostics = ahv::advanceProductionHydroSteps(state, 12U, 1.5e-5);
  ahv::requireOrThrow(diagnostics.advanced_patch_count == 8U, "Sedov: refined child patches were not advanced");
  ahv::requireOrThrow(diagnostics.ghost_fill.same_level_ghosts_filled > 0U, "Sedov: same-level refined ghosts were not filled");
  ahv::requireFinitePositiveState(state, "Sedov after evolution");

  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  const auto after = ahv::totalState(state, descriptors);
  ahv::requireOrThrow(
      ahv::relativeDifference(after.total_energy, before.total_energy) < 2.0e-3,
      "Sedov: energy left the CI bound");

  double central_pressure = 0.0;
  double outer_pressure = 0.0;
  double shell_density = 0.0;
  double shell_density_sq = 0.0;
  double radial_velocity_shell = 0.0;
  std::size_t central_count = 0;
  std::size_t outer_count = 0;
  std::size_t shell_count = 0;
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const double radius = radiusFromCenter(
        state.cells.center_x_comoving[row],
        state.cells.center_y_comoving[row],
        state.cells.center_z_comoving[row]);
    const auto primitive = ahv::primitiveForRow(state, row);
    if (radius < 0.24) {
      central_pressure += primitive.pressure_comoving;
      ++central_count;
    }
    if (radius > 0.48) {
      outer_pressure += primitive.pressure_comoving;
      ++outer_count;
    }
    if (radius > 0.24 && radius < 0.48) {
      const double dx = state.cells.center_x_comoving[row] - 0.5;
      const double dy = state.cells.center_y_comoving[row] - 0.5;
      const double dz = state.cells.center_z_comoving[row] - 0.5;
      const double inv_radius = 1.0 / std::max(radius, 1.0e-14);
      shell_density += primitive.rho_comoving;
      shell_density_sq += primitive.rho_comoving * primitive.rho_comoving;
      radial_velocity_shell +=
          (primitive.vel_x_peculiar * dx + primitive.vel_y_peculiar * dy + primitive.vel_z_peculiar * dz) *
          inv_radius;
      ++shell_count;
    }
  }
  central_pressure /= static_cast<double>(central_count);
  outer_pressure /= static_cast<double>(outer_count);
  shell_density /= static_cast<double>(shell_count);
  shell_density_sq /= static_cast<double>(shell_count);
  radial_velocity_shell /= static_cast<double>(shell_count);
  const double shell_density_variance = shell_density_sq - shell_density * shell_density;

  ahv::requireOrThrow(central_pressure > outer_pressure, "Sedov: central pressure excess vanished");
  ahv::requireOrThrow(radial_velocity_shell > 1.0e-5, "Sedov: shell did not move outward");
  ahv::requireOrThrow(
      shell_density_variance / std::max(shell_density * shell_density, 1.0e-14) < 0.35,
      "Sedov: refined blast shell is too asymmetric for CI guard");
}

}  // namespace

int main() {
  testAmrSedovRefinement();
  return 0;
}
