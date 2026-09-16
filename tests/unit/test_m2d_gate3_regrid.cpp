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
constexpr double k_tol = 1.0e-12;

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

[[nodiscard]] cosmosim::core::SimulationState makeRefineState() {
  cosmosim::core::SimulationState state;
  state.resizeCells(8);
  state.resizePatches(1);
  setPatch(state, 0, cosmosim::amr::PatchDescriptor{
      .patch_id = 301, .parent_patch_id = 0, .level = 0, .morton_key = 301,
      .origin_comov = {0.0, 0.0, 0.0}, .extent_comov = {1.0, 1.0, 1.0},
      .cell_dims = {2, 2, 2}}, 0, 8);
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(8);
  for (std::uint32_t row = 0; row < 8U; ++row) {
    const double rho = 1.0 + 0.05 * static_cast<double>(row);
    const double pressure = 1.0 + 0.02 * static_cast<double>(row);
    state.cells.center_x_comoving[row] = 0.25 + 0.5 * static_cast<double>(row % 2U);
    state.cells.center_y_comoving[row] = 0.25 + 0.5 * static_cast<double>((row / 2U) % 2U);
    state.cells.center_z_comoving[row] = 0.25 + 0.5 * static_cast<double>(row / 4U);
    state.cells.patch_index[row] = 0;
    state.cells.time_bin[row] = 0;
    state.cells.mass_code[row] = rho * 0.125;
    state.gas_cells.gas_cell_id[row] = 9201 + row;
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
        .gas_cell_id = 9201 + row,
        .parent_particle_id = std::nullopt,
        .owning_patch_id = 301,
        .local_cell_row = row});
  }
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

[[nodiscard]] cosmosim::amr::ProductionAmrHydroOptions regridOptions(
    const cosmosim::core::SimulationState& state, cosmosim::core::MemoryGovernor& governor) {
  const auto baseline = cosmosim::core::memoryReportBaselineOwnedBytes(
      cosmosim::core::collectSimulationMemoryReport(state));
  governor.setBaselineOwnedBytes(baseline);
  cosmosim::amr::ProductionAmrHydroOptions options;
  options.adiabatic_index = k_gamma;
  options.regrid_memory_governor = &governor;
  return options;
}

void testRefineHeadroomRejectionAndRetry() {
  auto state = makeRefineState();
  const auto parent = cosmosim::amr::buildProductionAmrPatchDescriptors(state).front();
  const std::size_t cells_before = state.cells.size();
  const std::size_t patches_before = state.patches.size();
  cosmosim::core::MemoryGovernor tight({.hard_limit_bytes = 1U});
  tight.setBaselineOwnedBytes(UINT64_MAX - 1U);
  // Force admission failure by using a nearly-exhausted governor.
  cosmosim::amr::ProductionAmrHydroOptions tight_options;
  tight_options.adiabatic_index = k_gamma;
  tight_options.regrid_memory_governor = &tight;
  bool rejected = false;
  try {
    (void)cosmosim::amr::refineProductionPatchInSimulationState(
        state, parent, 400, 20000, tight_options);
  } catch (const std::exception&) {
    rejected = true;
  }
  assert(rejected);
  assert(state.cells.size() == cells_before);
  assert(state.patches.size() == patches_before);
  // Retry with sufficient headroom preserves conservation.
  cosmosim::core::MemoryGovernor ample;
  auto options = regridOptions(state, ample);
  const auto diag = cosmosim::amr::refineProductionPatchInSimulationState(
      state, parent, 400, 20000, options);
  assert(diag.refined_patch_count == 1U);
  assert(std::abs(diag.conserved_mass_before - diag.conserved_mass_after) < k_tol);
}

void testDerefineHeadroomRejectionAndRetry() {
  auto state = makeRefineState();
  cosmosim::core::MemoryGovernor ample;
  auto options = regridOptions(state, ample);
  const auto parent = cosmosim::amr::buildProductionAmrPatchDescriptors(state).front();
  const auto refine = cosmosim::amr::refineProductionPatchInSimulationState(
      state, parent, 400, 20000, options);
  assert(refine.refined_patch_count == 1U);
  const std::size_t cells_after_refine = state.cells.size();
  cosmosim::core::MemoryGovernor tight({.hard_limit_bytes = 1U});
  tight.setBaselineOwnedBytes(UINT64_MAX - 1U);
  cosmosim::amr::ProductionAmrHydroOptions tight_options;
  tight_options.adiabatic_index = k_gamma;
  tight_options.regrid_memory_governor = &tight;
  bool rejected = false;
  try {
    (void)cosmosim::amr::derefineProductionPatchInSimulationState(
        state, parent, 30000, tight_options);
  } catch (const std::exception&) {
    rejected = true;
  }
  assert(rejected);
  assert(state.cells.size() == cells_after_refine);
  auto options2 = regridOptions(state, ample);
  const auto derefine = cosmosim::amr::derefineProductionPatchInSimulationState(
      state, parent, 30000, options2);
  assert(derefine.derefined_patch_count == 1U);
  assert(std::abs(derefine.conserved_mass_before - derefine.conserved_mass_after) < k_tol);
}

}  // namespace

int main() {
  testRefineHeadroomRejectionAndRetry();
  testDerefineHeadroomRejectionAndRetry();
  return 0;
}
