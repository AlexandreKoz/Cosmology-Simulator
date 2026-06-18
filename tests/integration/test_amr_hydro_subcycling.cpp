#include <cassert>
#include <cmath>
#include <cstdint>
#include <optional>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace {

constexpr double k_gamma = 1.4;

void setPatch(
    cosmosim::core::SimulationState& state,
    std::size_t patch_index,
    const cosmosim::amr::PatchDescriptor& descriptor,
    std::uint32_t first_cell,
    std::uint32_t cell_count) {
  cosmosim::amr::writePatchDescriptorToStateRow(state, patch_index, descriptor);
  state.patches.first_cell[patch_index] = first_cell;
  state.patches.cell_count[patch_index] = cell_count;
  state.patches.owning_rank[patch_index] = 0;
}

void setCell(
    cosmosim::core::SimulationState& state,
    std::uint32_t row,
    double x,
    double rho,
    double velocity_x,
    double pressure,
    std::uint32_t patch_index,
    std::uint64_t patch_id,
    std::uint64_t gas_cell_id,
    std::vector<cosmosim::core::GasCellIdentityRecord>& records) {
  state.cells.center_x_comoving[row] = x;
  state.cells.center_y_comoving[row] = 0.5;
  state.cells.center_z_comoving[row] = 0.5;
  state.cells.patch_index[row] = patch_index;
  state.cells.time_bin[row] = 0;
  state.cells.mass_code[row] = rho;
  state.gas_cells.gas_cell_id[row] = gas_cell_id;
  state.gas_cells.parent_particle_id[row] = 0;
  state.gas_cells.density_code[row] = rho;
  state.gas_cells.pressure_code[row] = pressure;
  state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
  state.gas_cells.velocity_x_peculiar[row] = velocity_x;
  state.gas_cells.velocity_y_peculiar[row] = 0.0;
  state.gas_cells.velocity_z_peculiar[row] = 0.0;
  state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
  state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * pressure / rho);
  records.push_back(cosmosim::core::GasCellIdentityRecord{
      .gas_cell_id = gas_cell_id,
      .parent_particle_id = std::nullopt,
      .owning_patch_id = patch_id,
      .local_cell_row = row});
}

[[nodiscard]] cosmosim::core::SimulationState makeCoarseFineState() {
  cosmosim::core::SimulationState state;
  state.resizeCells(4);
  state.resizePatches(2);
  setPatch(state, 0, cosmosim::amr::PatchDescriptor{
      .patch_id = 101,
      .parent_patch_id = 0,
      .level = 0,
      .morton_key = 101,
      .origin_comov = {0.0, 0.0, 0.0},
      .extent_comov = {1.0, 1.0, 1.0},
      .cell_dims = {2, 1, 1}}, 0, 2);
  setPatch(state, 1, cosmosim::amr::PatchDescriptor{
      .patch_id = 201,
      .parent_patch_id = 101,
      .level = 1,
      .morton_key = 201,
      .origin_comov = {1.0, 0.0, 0.0},
      .extent_comov = {0.5, 1.0, 1.0},
      .cell_dims = {2, 1, 1}}, 2, 2);
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  setCell(state, 0, 0.25, 1.0, 0.20, 1.0, 0, 101, 9001, records);
  setCell(state, 1, 0.75, 0.9, 0.35, 0.9, 0, 101, 9002, records);
  setCell(state, 2, 1.125, 0.8, -0.25, 0.8, 1, 201, 9101, records);
  setCell(state, 3, 1.375, 0.7, -0.10, 0.7, 1, 201, 9102, records);
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

void testTwoLevelSubcyclingAdvancesFineTwice() {
  auto state = makeCoarseFineState();
  assert(cosmosim::amr::hasProductionAmrHydroCoverage(state));
  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::HlleRiemannSolver riemann;
  const cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-5, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::amr::ProductionAmrHydroOptions options;
  options.adiabatic_index = k_gamma;
  options.density_floor = 1.0e-12;
  options.pressure_floor = 1.0e-12;
  options.refinement_ratio = 2;
  const auto diagnostics = cosmosim::amr::advanceProductionAmrHydroSubcycled(
      state,
      {},
      update,
      source_context,
      solver,
      riemann,
      {},
      options);
  assert(diagnostics.levels_advanced == 2U);
  assert(diagnostics.substeps_by_level.size() >= 2U);
  assert(diagnostics.substeps_by_level[0] == 1U);
  assert(diagnostics.substeps_by_level[1] == 2U);
  assert(diagnostics.advanced_patch_count == 3U);
  assert(diagnostics.flux_register_entry_count > 0U);
  assert(diagnostics.pending_register_applied_count > 0U);
  assert(state.pending_flux_registers.empty());
  assert(state.gas_cell_identity.coversDenseLocalRows(state.cells.size()));
  assert(state.gasCellIdentityMapMatchesSidecarLanes());
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    assert(std::isfinite(state.gas_cells.density_code[row]));
    assert(std::isfinite(state.gas_cells.pressure_code[row]));
    assert(std::isfinite(state.gas_cells.internal_energy_code[row]));
    assert(state.gas_cells.density_code[row] > 0.0);
    assert(state.gas_cells.pressure_code[row] > 0.0);
    assert(state.gas_cells.internal_energy_code[row] > 0.0);
  }
}

void testIncompleteSubstepCoverageDoesNotRefluxEarly() {
  auto state = makeCoarseFineState();
  std::vector<cosmosim::amr::FluxRegisterEntry> entries;
  cosmosim::amr::FluxRegisterEntry coarse;
  coarse.register_key = 123;
  coarse.coarse_patch_id = 101;
  coarse.coarse_gas_cell_id = 9002;
  coarse.coarse_cell_index = 1;
  coarse.level = 0;
  coarse.axis = cosmosim::hydro::HydroFaceAxis::kX;
  coarse.orientation = cosmosim::hydro::HydroFaceSide::kUpper;
  coarse.coarse_face_flux_code.mass_code = 0.1;
  coarse.face_area_comov = 0.5;
  coarse.coarse_area_comov = 0.5;
  coarse.dt_code = 1.0e-5;
  coarse.coarse_face_count = 1;
  entries.push_back(coarse);
  cosmosim::amr::ProductionAmrHydroOptions options;
  options.adiabatic_index = k_gamma;
  options.persist_incomplete_flux_registers = true;
  options.expected_fine_substeps = 2;
  options.reflux_coarse_dt_code = 1.0e-5;
  (void)cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, options);
  entries.clear();
  cosmosim::amr::FluxRegisterEntry fine = coarse;
  fine.coarse_face_count = 0;
  fine.coarse_area_comov = 0.0;
  fine.fine_face_count = 1;
  fine.fine_area_comov = 0.5;
  fine.fine_face_flux_code.mass_code = 0.4;
  fine.dt_code = 5.0e-6;
  entries.push_back(fine);
  options.fine_substep_index = 0;
  (void)cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, options);
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  const auto diagnostics = cosmosim::amr::applyCompletePendingFluxRegistersToSimulationState(state, descriptors, k_gamma);
  assert(diagnostics.complete_register_count == 0U);
  assert(diagnostics.skipped_incomplete_register_count == 1U);
  assert(state.pending_flux_registers.size() == 1U);
}

}  // namespace

int main() {
  testTwoLevelSubcyclingAdvancesFineTwice();
  testIncompleteSubstepCoverageDoesNotRefluxEarly();
  return 0;
}
