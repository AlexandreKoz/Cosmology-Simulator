#include <cassert>
#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"

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
    double pressure,
    std::uint32_t patch_index,
    std::uint64_t patch_id,
    std::uint64_t gas_cell_id,
    std::vector<cosmosim::core::GasCellIdentityRecord>& records) {
  state.cells.center_x_comoving[row] = x;
  state.cells.center_y_comoving[row] = 0.5;
  state.cells.center_z_comoving[row] = 0.5;
  state.cells.patch_index[row] = patch_index;
  state.cells.mass_code[row] = rho;
  state.cells.time_bin[row] = 0;
  state.gas_cells.gas_cell_id[row] = gas_cell_id;
  state.gas_cells.parent_particle_id[row] = 0;
  state.gas_cells.density_code[row] = rho;
  state.gas_cells.pressure_code[row] = pressure;
  state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
  state.gas_cells.velocity_x_peculiar[row] = 0.0;
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

[[nodiscard]] cosmosim::core::SimulationState makeState() {
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
  setCell(state, 0, 0.25, 1.0, 1.0, 0, 101, 9001, records);
  setCell(state, 1, 0.75, 0.9, 0.9, 0, 101, 9002, records);
  setCell(state, 2, 1.125, 0.8, 0.8, 1, 201, 9101, records);
  setCell(state, 3, 1.375, 0.7, 0.7, 1, 201, 9102, records);
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

[[nodiscard]] cosmosim::amr::FluxRegisterEntry makeEntry(
    bool coarse,
    bool fine,
    double fine_area = 0.5,
    std::uint64_t coarse_gas_cell_id = 9002) {
  cosmosim::amr::FluxRegisterEntry entry;
  entry.register_key = 77;
  entry.coarse_patch_id = 101;
  entry.coarse_gas_cell_id = coarse_gas_cell_id;
  entry.coarse_cell_index = 1;
  entry.level = 0;
  entry.axis = cosmosim::hydro::HydroFaceAxis::kX;
  entry.orientation = cosmosim::hydro::HydroFaceSide::kUpper;
  entry.face_area_comov = 0.5;
  entry.dt_code = coarse ? 1.0e-4 : 5.0e-5;
  if (coarse) {
    entry.coarse_area_comov = 0.5;
    entry.coarse_face_count = 1;
    entry.coarse_face_flux_code.mass_code = 0.1;
    entry.coarse_face_flux_code.momentum_x_code = 0.2;
    entry.coarse_face_flux_code.total_energy_code = 0.3;
  }
  if (fine) {
    entry.fine_area_comov = fine_area;
    entry.fine_face_count = 1;
    entry.fine_face_flux_code.mass_code = 0.4;
    entry.fine_face_flux_code.momentum_x_code = 0.5;
    entry.fine_face_flux_code.total_energy_code = 0.6;
  }
  return entry;
}

[[nodiscard]] cosmosim::amr::ProductionAmrHydroOptions pendingOptions(std::uint32_t substep) {
  cosmosim::amr::ProductionAmrHydroOptions options;
  options.adiabatic_index = k_gamma;
  options.persist_incomplete_flux_registers = true;
  options.reflux_coarse_dt_code = 1.0e-4;
  options.reflux_interval_start_code = 0.0;
  options.reflux_interval_end_code = 1.0e-4;
  options.expected_fine_substeps = 2;
  options.fine_substep_index = substep;
  return options;
}

void testPendingAppliesOnlyAfterFinalFineSubstep() {
  auto state = makeState();
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  const double density_before = state.gas_cells.density_code[1];

  std::vector<cosmosim::amr::FluxRegisterEntry> entries{makeEntry(true, false)};
  assert(cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, pendingOptions(0)) == 1U);
  auto diagnostics = cosmosim::amr::applyCompletePendingFluxRegistersToSimulationState(state, descriptors, k_gamma);
  assert(diagnostics.complete_register_count == 0U);
  assert(state.pending_flux_registers.size() == 1U);
  assert(state.gas_cells.density_code[1] == density_before);

  entries = {makeEntry(false, true)};
  assert(cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, pendingOptions(0)) == 0U);
  diagnostics = cosmosim::amr::applyCompletePendingFluxRegistersToSimulationState(state, descriptors, k_gamma);
  assert(diagnostics.complete_register_count == 0U);
  assert(state.pending_flux_registers.size() == 1U);

  assert(cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, pendingOptions(1)) == 0U);
  diagnostics = cosmosim::amr::applyCompletePendingFluxRegistersToSimulationState(state, descriptors, k_gamma);
  assert(diagnostics.complete_register_count == 1U);
  assert(diagnostics.corrected_cells == 1U);
  assert(state.pending_flux_registers.empty());
  assert(state.gas_cells.density_code[1] != density_before);
}

void testAreaMismatchAndMissingTargetStayRejected() {
  auto state = makeState();
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  std::vector<cosmosim::amr::FluxRegisterEntry> entries{makeEntry(true, false)};
  (void)cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, pendingOptions(0));
  entries = {makeEntry(false, true, 0.25)};
  (void)cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, pendingOptions(0));
  (void)cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, pendingOptions(1));
  auto diagnostics = cosmosim::amr::applyCompletePendingFluxRegistersToSimulationState(state, descriptors, k_gamma);
  assert(diagnostics.complete_register_count == 0U);
  assert(diagnostics.skipped_area_mismatch_count == 1U);
  assert(state.pending_flux_registers.size() == 1U);

  state = makeState();
  entries = {makeEntry(true, false, 0.5, 0)};
  (void)cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, pendingOptions(0));
  entries = {makeEntry(false, true, 0.5, 0)};
  (void)cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, pendingOptions(0));
  (void)cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, pendingOptions(1));
  diagnostics = cosmosim::amr::applyCompletePendingFluxRegistersToSimulationState(state, descriptors, k_gamma);
  assert(diagnostics.complete_register_count == 0U);
  assert(diagnostics.skipped_missing_target_count == 1U);
}

void testStaleGenerationRejected() {
  auto state = makeState();
  auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  std::vector<cosmosim::amr::FluxRegisterEntry> entries{makeEntry(true, false)};
  (void)cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, pendingOptions(0));
  entries = {makeEntry(false, true)};
  (void)cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, pendingOptions(0));
  (void)cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, pendingOptions(1));
  std::vector<cosmosim::core::GasCellIdentityRecord> rebuilt(
      state.gas_cell_identity.records().begin(), state.gas_cell_identity.records().end());
  state.gas_cell_identity.assign(std::move(rebuilt));
  const auto diagnostics = cosmosim::amr::applyCompletePendingFluxRegistersToSimulationState(state, descriptors, k_gamma);
  assert(diagnostics.complete_register_count == 0U);
  assert(diagnostics.skipped_missing_target_count == 1U);
  assert(state.pending_flux_registers.size() == 1U);
}

}  // namespace

int main() {
  testPendingAppliesOnlyAfterFinalFineSubstep();
  testAreaMismatchAndMissingTargetStayRejected();
  testStaleGenerationRejected();
  return 0;
}
