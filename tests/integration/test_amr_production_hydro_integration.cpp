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
constexpr double k_tol = 1.0e-9;

struct Totals {
  double mass = 0.0;
  double momentum_x = 0.0;
  double momentum_y = 0.0;
  double momentum_z = 0.0;
  double total_energy = 0.0;
};

[[nodiscard]] cosmosim::hydro::HydroPrimitiveState primitiveForRow(
    const cosmosim::core::SimulationState& state,
    std::uint32_t row) {
  return cosmosim::hydro::HydroPrimitiveState{
      .rho_comoving = state.gas_cells.density_code[row],
      .vel_x_peculiar = state.gas_cells.velocity_x_peculiar[row],
      .vel_y_peculiar = state.gas_cells.velocity_y_peculiar[row],
      .vel_z_peculiar = state.gas_cells.velocity_z_peculiar[row],
      .pressure_comoving = state.gas_cells.pressure_code[row]};
}

[[nodiscard]] Totals totalState(
    const cosmosim::core::SimulationState& state,
    const std::vector<cosmosim::amr::PatchDescriptor>& descriptors) {
  Totals totals;
  for (const cosmosim::amr::PatchDescriptor& patch : descriptors) {
    const double volume = patch.extent_comov[0] * patch.extent_comov[1] * patch.extent_comov[2] /
        static_cast<double>(static_cast<std::size_t>(patch.cell_dims[0]) * patch.cell_dims[1] * patch.cell_dims[2]);
    const std::vector<std::uint32_t> rows = state.gas_cell_identity.rowsForPatch(patch.patch_id);
    for (const std::uint32_t row : rows) {
      const auto conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(
          primitiveForRow(state, row), k_gamma);
      totals.mass += conserved.mass_density_comoving * volume;
      totals.momentum_x += conserved.momentum_density_x_comoving * volume;
      totals.momentum_y += conserved.momentum_density_y_comoving * volume;
      totals.momentum_z += conserved.momentum_density_z_comoving * volume;
      totals.total_energy += conserved.total_energy_density_comoving * volume;
    }
  }
  return totals;
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
  state.patches.patch_id[0] = 101;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 2;
  state.patches.owning_rank[0] = 0;
  state.patches.patch_id[1] = 201;
  state.patches.level[1] = 1;
  state.patches.first_cell[1] = 2;
  state.patches.cell_count[1] = 2;
  state.patches.owning_rank[1] = 0;

  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(4);
  setCell(state, 0, 0.25, 1.0, 0.20, 1.0, 0, 101, 9001, records);
  setCell(state, 1, 0.75, 0.9, 0.35, 0.9, 0, 101, 9002, records);
  setCell(state, 2, 1.125, 0.8, -0.25, 0.8, 1, 201, 9101, records);
  setCell(state, 3, 1.375, 0.7, -0.10, 0.7, 1, 201, 9102, records);
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

void testProductionAmrGeometryScatterAndRefluxPath() {
  cosmosim::core::SimulationState state = makeCoarseFineState();
  assert(cosmosim::amr::hasProductionAmrHydroCoverage(state));

  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  assert(descriptors.size() == 2U);
  const Totals before = totalState(state, descriptors);

  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::HlleRiemannSolver riemann;
  const cosmosim::hydro::HydroUpdateContext update{
      .dt_code = 1.0e-4,
      .scale_factor = 1.0,
      .hubble_rate_code = 0.0};
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  const std::vector<std::uint32_t> active_rows{0, 1, 2, 3};
  const cosmosim::amr::ProductionAmrHydroOptions options{
      .physical_boundary_kind = cosmosim::hydro::HydroBoundaryKind::kOpen,
      .adiabatic_index = k_gamma,
      .density_floor = 1.0e-12,
      .pressure_floor = 1.0e-12};

  const auto diagnostics = cosmosim::amr::advanceProductionAmrHydro(
      state,
      active_rows,
      update,
      source_context,
      solver,
      riemann,
      {},
      options);

  assert(diagnostics.patch_count == 2U);
  assert(diagnostics.advanced_patch_count == 2U);
  assert(diagnostics.ghost_fill.same_level_ghosts_filled + diagnostics.ghost_fill.coarse_to_fine_ghosts_filled + diagnostics.ghost_fill.fine_to_coarse_ghosts_filled > 0U);
  assert(diagnostics.flux_register_entry_count > 0U);
  assert(diagnostics.reflux.corrected_cells > 0U);
  assert(diagnostics.reflux.corrected_mass_code >= 0.0);
  assert(diagnostics.reflux.corrected_momentum_x_code >= 0.0);
  assert(diagnostics.reflux.corrected_energy_code >= 0.0);

  assert(state.gas_cell_identity.coversDenseLocalRows(state.cells.size()));
  assert(state.gasCellIdentityMapMatchesSidecarLanes());
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    assert(state.gas_cells.gas_cell_id[row] != 0U);
    assert(state.gas_cells.density_code[row] > 0.0);
    assert(state.gas_cells.pressure_code[row] > 0.0);
    assert(std::isfinite(state.gas_cells.velocity_x_peculiar[row]));
    assert(std::isfinite(state.gas_cells.internal_energy_code[row]));
  }

  const Totals after = totalState(state, cosmosim::amr::buildProductionAmrPatchDescriptors(state));
  assert(std::isfinite(after.mass));
  assert(std::isfinite(after.total_energy));
  assert(std::abs(after.mass - before.mass) < 1.0e-2);
}

[[nodiscard]] cosmosim::core::SimulationState makeRefineState() {
  cosmosim::core::SimulationState state;
  state.resizeCells(8);
  state.resizePatches(1);
  state.patches.patch_id[0] = 301;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 8;
  state.patches.owning_rank[0] = 0;
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(8);
  std::uint32_t row = 0;
  for (std::uint32_t k = 0; k < 2; ++k) {
    for (std::uint32_t j = 0; j < 2; ++j) {
      for (std::uint32_t i = 0; i < 2; ++i) {
        const double rho = 1.0 + 0.05 * static_cast<double>(row);
        const double pressure = 1.0 + 0.02 * static_cast<double>(row);
        state.cells.center_x_comoving[row] = 0.25 + 0.5 * static_cast<double>(i);
        state.cells.center_y_comoving[row] = 0.25 + 0.5 * static_cast<double>(j);
        state.cells.center_z_comoving[row] = 0.25 + 0.5 * static_cast<double>(k);
        state.cells.patch_index[row] = 0;
        state.cells.time_bin[row] = 0;
        state.cells.mass_code[row] = rho * 0.125;
        state.gas_cells.gas_cell_id[row] = 9201 + row;
        state.gas_cells.parent_particle_id[row] = 0;
        state.gas_cells.density_code[row] = rho;
        state.gas_cells.pressure_code[row] = pressure;
        state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
        state.gas_cells.velocity_x_peculiar[row] = 0.1 * static_cast<double>(i);
        state.gas_cells.velocity_y_peculiar[row] = -0.05 * static_cast<double>(j);
        state.gas_cells.velocity_z_peculiar[row] = 0.02 * static_cast<double>(k);
        state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
        state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * pressure / rho);
        records.push_back(cosmosim::core::GasCellIdentityRecord{
            .gas_cell_id = 9201 + row,
            .parent_particle_id = std::nullopt,
            .owning_patch_id = 301,
            .local_cell_row = row});
        ++row;
      }
    }
  }
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

void assertClose(double lhs, double rhs) {
  assert(std::abs(lhs - rhs) < k_tol);
}

void testProductionRefineAndDerefineConserveSimulationState() {
  cosmosim::core::SimulationState state = makeRefineState();
  const auto parent = cosmosim::amr::buildProductionAmrPatchDescriptors(state).front();
  const auto refine = cosmosim::amr::refineProductionPatchInSimulationState(state, parent, 400, 20000);
  assert(refine.refined_patch_count == 1U);
  assert(refine.created_gas_cell_count == 64U);
  assert(refine.retired_gas_cell_count == 8U);
  assert(state.patches.size() == 8U);
  assert(state.cells.size() == 64U);
  assert(state.gas_cell_identity.coversDenseLocalRows(state.cells.size()));
  assert(state.gasCellIdentityMapMatchesSidecarLanes());
  assertClose(refine.conserved_mass_before, refine.conserved_mass_after);
  assertClose(refine.conserved_momentum_x_before, refine.conserved_momentum_x_after);
  assertClose(refine.conserved_momentum_y_before, refine.conserved_momentum_y_after);
  assertClose(refine.conserved_momentum_z_before, refine.conserved_momentum_z_after);
  assertClose(refine.conserved_total_energy_before, refine.conserved_total_energy_after);

  const auto derefine = cosmosim::amr::derefineProductionPatchInSimulationState(state, parent, 30000);
  assert(derefine.derefined_patch_count == 1U);
  assert(derefine.created_gas_cell_count == 8U);
  assert(derefine.retired_gas_cell_count == 64U);
  assert(state.patches.size() == 1U);
  assert(state.cells.size() == 8U);
  assert(state.gas_cell_identity.coversDenseLocalRows(state.cells.size()));
  assert(state.gasCellIdentityMapMatchesSidecarLanes());
  assertClose(derefine.conserved_mass_before, derefine.conserved_mass_after);
  assertClose(derefine.conserved_momentum_x_before, derefine.conserved_momentum_x_after);
  assertClose(derefine.conserved_momentum_y_before, derefine.conserved_momentum_y_after);
  assertClose(derefine.conserved_momentum_z_before, derefine.conserved_momentum_z_after);
  assertClose(derefine.conserved_total_energy_before, derefine.conserved_total_energy_after);
}

}  // namespace

int main() {
  testProductionAmrGeometryScatterAndRefluxPath();
  testProductionRefineAndDerefineConserveSimulationState();
  return 0;
}
