#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "amr_hydro_validation_helpers.hpp"

namespace {

namespace ahv = cosmosim::tests::amr_hydro_validation;

[[nodiscard]] std::array<cosmosim::amr::PatchDescriptor, 2> shockTubePatches() {
  return {
      cosmosim::amr::PatchDescriptor{
          .patch_id = 101,
          .parent_patch_id = 0,
          .level = 0,
          .morton_key = 101,
          .origin_comov = {0.0, 0.0, 0.0},
          .extent_comov = {0.5, 1.0, 1.0},
          .cell_dims = {4, 1, 1}},
      cosmosim::amr::PatchDescriptor{
          .patch_id = 102,
          .parent_patch_id = 0,
          .level = 0,
          .morton_key = 102,
          .origin_comov = {0.5, 0.0, 0.0},
          .extent_comov = {0.5, 1.0, 1.0},
          .cell_dims = {4, 1, 1}}};
}

[[nodiscard]] cosmosim::hydro::HydroPrimitiveState sodPrimitive(const ahv::CellCenter& center) {
  if (center.x_comoving < 0.5) {
    return cosmosim::hydro::HydroPrimitiveState{
        .rho_comoving = 1.0,
        .pressure_comoving = 1.0};
  }
  return cosmosim::hydro::HydroPrimitiveState{
      .rho_comoving = 0.125,
      .pressure_comoving = 0.1};
}

void testAmrShockTubeRefinement() {
  const auto patches = shockTubePatches();
  cosmosim::core::SimulationState state = ahv::makePatchState(
      {patches[0], patches[1]},
      [](const cosmosim::amr::PatchDescriptor&, std::size_t, const ahv::CellCenter& center) {
        return sodPrimitive(center);
      });
  const auto before = ahv::totalState(state, cosmosim::amr::buildProductionAmrPatchDescriptors(state));

  const auto refine = cosmosim::amr::refineProductionPatchInSimulationState(state, patches[0], 1000, 20000);
  ahv::requireOrThrow(refine.refined_patch_count == 1U, "shock tube: discontinuity-adjacent patch was not refined");
  ahv::requireOrThrow(refine.created_gas_cell_count == 32U, "shock tube: unexpected child gas-cell count");
  ahv::requireFinitePositiveState(state, "shock tube after refine");

  const auto diagnostics = ahv::advanceProductionHydroSteps(state, 8U, 2.5e-5);
  ahv::requireOrThrow(diagnostics.ghost_fill.coarse_to_fine_ghosts_filled > 0U, "shock tube: no coarse-fine ghosts filled");
  ahv::requireOrThrow(diagnostics.flux_register_entry_count > 0U, "shock tube: no flux registers generated");
  ahv::requireOrThrow(diagnostics.reflux.complete_register_count > 0U, "shock tube: no complete reflux register");
  ahv::requireOrThrow(diagnostics.reflux.corrected_cells > 0U, "shock tube: reflux did not correct a coarse owner");
  ahv::requireFinitePositiveState(state, "shock tube after evolution");

  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  double left_density = 0.0;
  double right_density = 0.0;
  double left_pressure = 0.0;
  double right_pressure = 0.0;
  double max_velocity_x = -1.0e30;
  std::size_t left_count = 0;
  std::size_t right_count = 0;
  for (const auto& patch : descriptors) {
    for (const std::uint32_t row : state.gas_cell_identity.rowsForPatch(patch.patch_id)) {
      if (state.cells.center_x_comoving[row] < 0.5) {
        left_density += state.gas_cells.density_code[row];
        left_pressure += state.gas_cells.pressure_code[row];
        ++left_count;
      } else {
        right_density += state.gas_cells.density_code[row];
        right_pressure += state.gas_cells.pressure_code[row];
        ++right_count;
      }
      max_velocity_x = std::max(max_velocity_x, state.gas_cells.velocity_x_peculiar[row]);
    }
  }
  left_density /= static_cast<double>(left_count);
  right_density /= static_cast<double>(right_count);
  left_pressure /= static_cast<double>(left_count);
  right_pressure /= static_cast<double>(right_count);

  ahv::requireOrThrow(left_density > right_density, "shock tube: density jump trend was erased");
  ahv::requireOrThrow(left_pressure > right_pressure, "shock tube: pressure jump trend was erased");
  ahv::requireOrThrow(max_velocity_x > 1.0e-5, "shock tube: shock/contact did not move in +x");
  const auto after = ahv::totalState(state, descriptors);
  ahv::requireOrThrow(ahv::relativeDifference(after.mass, before.mass) < 5.0e-3, "shock tube: mass drift");
  ahv::requireOrThrow(
      ahv::relativeDifference(after.total_energy, before.total_energy) < 5.0e-3,
      "shock tube: total-energy drift");
}

}  // namespace

int main() {
  testAmrShockTubeRefinement();
  return 0;
}
