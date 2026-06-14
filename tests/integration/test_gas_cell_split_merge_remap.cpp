#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <optional>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace {

constexpr double k_tol = 1.0e-12;

struct Totals {
  double mass = 0.0;
  double momentum_x = 0.0;
  double momentum_y = 0.0;
  double momentum_z = 0.0;
  double energy = 0.0;
};

Totals totals(const cosmosim::hydro::HydroConservedStateSoa& conserved) {
  Totals total;
  for (std::size_t row = 0; row < conserved.size(); ++row) {
    const auto cell = conserved.loadCell(row);
    total.mass += cell.mass_density_comoving;
    total.momentum_x += cell.momentum_density_x_comoving;
    total.momentum_y += cell.momentum_density_y_comoving;
    total.momentum_z += cell.momentum_density_z_comoving;
    total.energy += cell.total_energy_density_comoving;
  }
  return total;
}

void requireSameTotals(const Totals& lhs, const Totals& rhs) {
  assert(std::abs(lhs.mass - rhs.mass) < k_tol);
  assert(std::abs(lhs.momentum_x - rhs.momentum_x) < k_tol);
  assert(std::abs(lhs.momentum_y - rhs.momentum_y) < k_tol);
  assert(std::abs(lhs.momentum_z - rhs.momentum_z) < k_tol);
  assert(std::abs(lhs.energy - rhs.energy) < k_tol);
}

cosmosim::hydro::HydroConservedState conserved(
    double mass,
    double momentum_x,
    double momentum_y,
    double momentum_z,
    double energy) {
  return cosmosim::hydro::HydroConservedState{
      .mass_density_comoving = mass,
      .momentum_density_x_comoving = momentum_x,
      .momentum_density_y_comoving = momentum_y,
      .momentum_density_z_comoving = momentum_z,
      .total_energy_density_comoving = energy,
  };
}

void mirrorIdentityToSidecar(cosmosim::core::SimulationState& state) {
  for (const auto& record : state.gas_cell_identity.records()) {
    const std::uint32_t row = record.local_cell_row;
    state.gas_cells.gas_cell_id[row] = record.gas_cell_id;
    state.gas_cells.parent_particle_id[row] =
        record.parent_particle_id.has_value() ? *record.parent_particle_id : 0U;
  }
}

void test_split_keeps_parent_lineage_nonunique_and_conserves() {
  cosmosim::core::SimulationState old_state;
  old_state.resizeCells(1);
  old_state.gas_cell_identity.assign({
      {.gas_cell_id = 1001, .parent_particle_id = 42, .owning_patch_id = 0, .local_cell_row = 0},
  });
  mirrorIdentityToSidecar(old_state);
  old_state.gas_cells.velocity_x_peculiar[0] = 2.0;
  old_state.gas_cells.velocity_y_peculiar[0] = -1.0;
  old_state.gas_cells.velocity_z_peculiar[0] = 0.5;
  assert(old_state.validateOwnershipInvariants());

  cosmosim::hydro::HydroConservedStateSoa before(1);
  before.storeCell(0, conserved(12.0, 24.0, -12.0, 6.0, 90.0));

  cosmosim::core::SimulationState split_state;
  split_state.resizeCells(2);
  split_state.gas_cell_identity.assign({
      {.gas_cell_id = 1001, .parent_particle_id = 42, .owning_patch_id = 0, .local_cell_row = 0},
      {.gas_cell_id = 1002, .parent_particle_id = 42, .owning_patch_id = 0, .local_cell_row = 1},
  });
  mirrorIdentityToSidecar(split_state);
  assert(split_state.gas_cell_identity.rowsForParentParticleId(42) == std::vector<std::uint32_t>({0, 1}));
  assert(split_state.gasCellIdentityMapMatchesSidecarLanes());
  assert(split_state.validateOwnershipInvariants());

  const auto new_to_old = cosmosim::core::buildGasCellNewToOldRowMap(
      old_state.gas_cell_identity, split_state.gas_cell_identity);
  assert(new_to_old[0] == 0U);
  assert(new_to_old[1] == cosmosim::core::kInvalidGasCellRow);

  cosmosim::hydro::HydroConservedStateSoa after(2);
  after.storeCell(0, conserved(5.0, 10.0, -5.0, 2.5, 37.5));
  after.storeCell(1, conserved(7.0, 14.0, -7.0, 3.5, 52.5));
  requireSameTotals(totals(before), totals(after));
}

void test_merge_parentless_result_conserves() {
  cosmosim::core::SimulationState merged_state;
  merged_state.resizeCells(1);
  merged_state.gas_cell_identity.assign({
      {.gas_cell_id = 3001, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 0},
  });
  mirrorIdentityToSidecar(merged_state);
  assert(!merged_state.parentParticleIdForGasCellRow(0).has_value());
  assert(merged_state.gas_cells.parent_particle_id[0] == 0U);
  assert(merged_state.validateOwnershipInvariants());

  cosmosim::hydro::HydroConservedStateSoa before(2);
  before.storeCell(0, conserved(3.0, 1.0, 2.0, 3.0, 11.0));
  before.storeCell(1, conserved(4.0, 5.0, 7.0, 9.0, 23.0));

  cosmosim::hydro::HydroConservedStateSoa after(1);
  after.storeCell(0, conserved(7.0, 6.0, 9.0, 12.0, 34.0));
  requireSameTotals(totals(before), totals(after));
}

void test_row_reorder_scatter_by_gas_cell_id() {
  cosmosim::core::SimulationState state;
  state.resizeCells(3);
  state.gas_cell_identity.assign({
      {.gas_cell_id = 7101, .parent_particle_id = 90, .owning_patch_id = 0, .local_cell_row = 0},
      {.gas_cell_id = 7102, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 1},
      {.gas_cell_id = 7103, .parent_particle_id = 90, .owning_patch_id = 0, .local_cell_row = 2},
  });
  mirrorIdentityToSidecar(state);
  for (std::size_t row = 0; row < state.cells.size(); ++row) {
    state.gas_cells.density_code[row] = 100.0 + static_cast<double>(row);
    state.gas_cells.pressure_code[row] = 200.0 + static_cast<double>(row);
  }

  cosmosim::core::TransientStepWorkspace workspace;
  const std::array<std::uint32_t, 1> active{1};
  auto view = cosmosim::core::buildHydroCellKernelView(state, active, workspace);
  assert(view.gas_cell_id[0] == 7102U);

  state.gas_cell_identity.assign({
      {.gas_cell_id = 7103, .parent_particle_id = 90, .owning_patch_id = 0, .local_cell_row = 0},
      {.gas_cell_id = 7101, .parent_particle_id = 90, .owning_patch_id = 0, .local_cell_row = 1},
      {.gas_cell_id = 7102, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 2},
  });
  mirrorIdentityToSidecar(state);
  state.gas_cells.density_code = {102.0, 100.0, 101.0};
  state.gas_cells.pressure_code = {202.0, 200.0, 201.0};

  view.source_gas_cell_identity_generation = state.gasCellIdentityGeneration();
  view.density_code[0] = 444.0;
  view.pressure_code[0] = 555.0;
  cosmosim::core::scatterHydroCellKernelView(view, state);

  assert(state.gas_cells.density_code[0] == 102.0);
  assert(state.gas_cells.density_code[1] == 100.0);
  assert(state.gas_cells.density_code[2] == 444.0);
  assert(state.gas_cells.pressure_code[2] == 555.0);
  assert(state.validateOwnershipInvariants());
}

}  // namespace

int main() {
  test_split_keeps_parent_lineage_nonunique_and_conserves();
  test_merge_parentless_result_conserves();
  test_row_reorder_scatter_by_gas_cell_id();
  return 0;
}
