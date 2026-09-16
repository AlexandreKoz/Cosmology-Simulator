#include <cassert>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"

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
      .patch_id = 101, .parent_patch_id = 0, .level = 0, .morton_key = 101,
      .origin_comov = {0.0, 0.0, 0.0}, .extent_comov = {1.0, 1.0, 1.0},
      .cell_dims = {2, 1, 1}}, 0, 2);
  setPatch(state, 1, cosmosim::amr::PatchDescriptor{
      .patch_id = 201, .parent_patch_id = 101, .level = 1, .morton_key = 201,
      .origin_comov = {1.0, 0.0, 0.0}, .extent_comov = {0.5, 1.0, 1.0},
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
    std::uint64_t key, std::uint64_t gas_id, std::size_t cell_index) {
  cosmosim::amr::FluxRegisterEntry entry;
  entry.register_key = key;
  entry.coarse_patch_id = 101;
  entry.coarse_gas_cell_id = gas_id;
  entry.coarse_cell_index = cell_index;
  entry.level = 0;
  entry.axis = cosmosim::hydro::HydroFaceAxis::kX;
  entry.orientation = cosmosim::hydro::HydroFaceSide::kUpper;
  entry.face_area_comov = 0.5;
  entry.coarse_area_comov = 0.5;
  entry.fine_area_comov = 0.5;
  entry.dt_code = 1.0e-4;
  entry.coarse_face_count = 1;
  entry.fine_face_count = 1;
  entry.coarse_face_flux_code.mass_code = 0.1;
  entry.fine_face_flux_code.mass_code = 0.4;
  return entry;
}

[[nodiscard]] cosmosim::amr::ProductionAmrHydroOptions pendingOptions(std::uint32_t substep) {
  cosmosim::amr::ProductionAmrHydroOptions options;
  options.adiabatic_index = k_gamma;
  options.persist_incomplete_flux_registers = true;
  options.reflux_coarse_dt_code = 1.0e-4;
  options.reflux_interval_start_code = 0.0;
  options.reflux_interval_end_code = 1.0e-4;
  options.expected_fine_substeps = 1;
  options.fine_substep_index = substep;
  return options;
}

void testSparseAllActiveEmptyRows() {
  auto state = makeState();
  cosmosim::amr::AmrActiveLevelRowWorkspace workspace;
  cosmosim::core::MemoryGovernor governor;
  // Empty requested = all active.
  cosmosim::amr::prepareAmrActiveLevelRows(state, {}, 0, 1, &governor, workspace);
  assert(workspace.coarse_rows.size() == 2U);
  assert(workspace.fine_rows.size() == 2U);
  // Sparse subset: only row 1 (coarse) requested.
  std::vector<std::uint32_t> subset = {1U};
  cosmosim::amr::AmrActiveLevelRowWorkspace sparse;
  cosmosim::amr::prepareAmrActiveLevelRows(
      state, std::span<const std::uint32_t>(subset.data(), subset.size()), 0, 1,
      &governor, sparse);
  assert(sparse.coarse_rows.size() == 1U);
  assert(sparse.fine_rows.empty());
  // Empty level: level 7 has no rows.
  std::vector<std::uint32_t> empty_out;
  cosmosim::amr::fillActiveRowsForLevelInto(state, 7, {}, empty_out);
  assert(empty_out.empty());
}

void testCoexistenceAdmissionAndStability() {
  auto state = makeState();
  // Exact needs for the fixture: 2 coarse + 2 fine rows coexist.
  const std::uint64_t exact =
      cosmosim::amr::activeLevelRowStagingBytes(2U, 2U);
  assert(exact == 4U * sizeof(std::uint32_t));
  cosmosim::core::MemoryGovernor tight({.hard_limit_bytes = exact - 1U});
  tight.setBaselineOwnedBytes(0U);
  cosmosim::amr::AmrActiveLevelRowWorkspace workspace;
  bool rejected = false;
  try {
    cosmosim::amr::prepareAmrActiveLevelRows(state, {}, 0, 1, &tight, workspace);
  } catch (const cosmosim::core::MemoryAdmissionError&) {
    rejected = true;
  }
  assert(rejected);
  // Retry with sufficient headroom succeeds and stabilizes.
  cosmosim::core::MemoryGovernor ample({.hard_limit_bytes = exact + 1024U});
  ample.setBaselineOwnedBytes(0U);
  cosmosim::amr::prepareAmrActiveLevelRows(state, {}, 0, 1, &ample, workspace);
  const auto cap_coarse = workspace.coarse_rows.capacity();
  const auto cap_fine = workspace.fine_rows.capacity();
  assert(cap_coarse > 0U && cap_fine > 0U);
  cosmosim::amr::prepareAmrActiveLevelRows(state, {}, 0, 1, &ample, workspace);
  assert(workspace.coarse_rows.capacity() == cap_coarse);
  assert(workspace.fine_rows.capacity() == cap_fine);
}

void testInvalidLaterRecordCausesZeroMutation() {
  auto state = makeState();
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  const double density0_before = state.gas_cells.density_code[0];
  const double density1_before = state.gas_cells.density_code[1];
  // Two complete pending records: first valid (targets row 1 / gas 9002),
  // second stale (wrong cell index triggers the stale-mapping throw).
  std::vector<cosmosim::amr::FluxRegisterEntry> entries{
      makeEntry(77U, 9002U, 1U), makeEntry(78U, 9002U, 0U)};
  auto options = pendingOptions(0);
  cosmosim::core::MemoryGovernor governor;
  governor.setBaselineOwnedBytes(cosmosim::core::memoryReportBaselineOwnedBytes(
      cosmosim::core::collectSimulationMemoryReport(state)));
  options.regrid_memory_governor = &governor;
  (void)cosmosim::amr::mergeFluxRegistersIntoPendingStore(state, entries, options);
  assert(state.pending_flux_registers.size() == 2U);
  bool threw = false;
  try {
    (void)cosmosim::amr::applyCompletePendingFluxRegistersToSimulationState(
        state, descriptors, k_gamma);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
  // Canonical conserved state unchanged and pending state unchanged (retry possible).
  assert(state.gas_cells.density_code[0] == density0_before);
  assert(state.gas_cells.density_code[1] == density1_before);
  assert(state.pending_flux_registers.size() == 2U);
}

}  // namespace

int main() {
  testSparseAllActiveEmptyRows();
  testCoexistenceAdmissionAndStability();
  testInvalidLaterRecordCausesZeroMutation();
  return 0;
}
