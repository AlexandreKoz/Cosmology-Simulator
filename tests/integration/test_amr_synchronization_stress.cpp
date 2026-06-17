#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "amr_hydro_validation_helpers.hpp"

namespace {

namespace ahv = cosmosim::tests::amr_hydro_validation;

[[nodiscard]] std::array<cosmosim::amr::PatchDescriptor, 2> stressPatches() {
  return {
      cosmosim::amr::PatchDescriptor{
          .patch_id = 501,
          .parent_patch_id = 0,
          .level = 0,
          .morton_key = 501,
          .origin_comov = {0.0, 0.0, 0.0},
          .extent_comov = {0.5, 1.0, 1.0},
          .cell_dims = {2, 2, 1}},
      cosmosim::amr::PatchDescriptor{
          .patch_id = 502,
          .parent_patch_id = 0,
          .level = 0,
          .morton_key = 502,
          .origin_comov = {0.5, 0.0, 0.0},
          .extent_comov = {0.5, 1.0, 1.0},
          .cell_dims = {2, 2, 1}}};
}

[[nodiscard]] cosmosim::hydro::HydroPrimitiveState stressPrimitive(const ahv::CellCenter& center) {
  const double wave = std::sin(2.0 * 3.14159265358979323846 * center.y_comoving);
  return cosmosim::hydro::HydroPrimitiveState{
      .rho_comoving = center.x_comoving < 0.5 ? 1.0 + 0.1 * wave : 0.75 - 0.05 * wave,
      .vel_x_peculiar = center.x_comoving < 0.5 ? 0.06 : -0.04,
      .vel_y_peculiar = 0.015 * wave,
      .vel_z_peculiar = 0.0,
      .pressure_comoving = center.x_comoving < 0.5 ? 1.0 : 0.7};
}

void testRefineDerefineSynchronizationStress() {
  const auto patches = stressPatches();
  cosmosim::core::SimulationState state = ahv::makePatchState(
      {patches[0], patches[1]},
      [](const cosmosim::amr::PatchDescriptor&, std::size_t, const ahv::CellCenter& center) {
        return stressPrimitive(center);
      });
  const auto initial = ahv::totalState(state, cosmosim::amr::buildProductionAmrPatchDescriptors(state));

  auto diagnostics = ahv::advanceProductionHydroSteps(state, 3U, 1.0e-5);
  ahv::requireOrThrow(diagnostics.advanced_patch_count == 2U, "sync stress: initial coarse patches were not advanced");
  ahv::requireFinitePositiveState(state, "sync stress before refine");

  const auto refine = cosmosim::amr::refineProductionPatchInSimulationState(state, patches[0], 700, 40000);
  ahv::requireOrThrow(refine.refined_patch_count == 1U, "sync stress: refine did not run");
  ahv::requireOrThrow(refine.created_gas_cell_count == 32U, "sync stress: unexpected refined cell count");
  ahv::requireOrThrow(
      ahv::relativeDifference(refine.conserved_total_energy_after, refine.conserved_total_energy_before) < 1.0e-12,
      "sync stress: conservative prolongation changed energy");
  ahv::requireFinitePositiveState(state, "sync stress after refine");

  diagnostics = ahv::advanceProductionHydroSteps(state, 6U, 1.0e-5);
  ahv::requireOrThrow(diagnostics.ghost_fill.coarse_to_fine_ghosts_filled > 0U, "sync stress: no coarse-fine ghosts");
  ahv::requireOrThrow(diagnostics.flux_register_entry_count > 0U, "sync stress: no flux registers");
  ahv::requireOrThrow(diagnostics.reflux.complete_register_count > 0U, "sync stress: no complete reflux register");
  ahv::requireOrThrow(diagnostics.reflux.corrected_cells > 0U, "sync stress: reflux did not correct state");
  ahv::requireFinitePositiveState(state, "sync stress after refined evolution");

  const auto derefine = cosmosim::amr::derefineProductionPatchInSimulationState(state, patches[0], 50000);
  ahv::requireOrThrow(derefine.derefined_patch_count == 1U, "sync stress: derefine did not run");
  ahv::requireOrThrow(derefine.created_gas_cell_count == 4U, "sync stress: unexpected restricted parent cell count");
  ahv::requireOrThrow(derefine.retired_gas_cell_count == 32U, "sync stress: unexpected retired child count");
  ahv::requireOrThrow(
      ahv::relativeDifference(derefine.conserved_mass_after, derefine.conserved_mass_before) < 1.0e-12,
      "sync stress: conservative restriction changed mass");
  ahv::requireOrThrow(
      ahv::relativeDifference(derefine.conserved_total_energy_after, derefine.conserved_total_energy_before) < 1.0e-12,
      "sync stress: conservative restriction changed energy");
  ahv::requireFinitePositiveState(state, "sync stress after derefine");

  diagnostics = ahv::advanceProductionHydroSteps(state, 3U, 1.0e-5);
  ahv::requireOrThrow(diagnostics.advanced_patch_count == 2U, "sync stress: final coarse patches were not advanced");
  ahv::requireFinitePositiveState(state, "sync stress final state");

  const auto final = ahv::totalState(state, cosmosim::amr::buildProductionAmrPatchDescriptors(state));
  ahv::requireConservedClose(final, initial, 1.0e-2, "sync stress");
}

}  // namespace

int main() {
  testRefineDerefineSynchronizationStress();
  return 0;
}
