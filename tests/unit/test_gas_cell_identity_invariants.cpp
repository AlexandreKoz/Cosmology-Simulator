#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"

namespace {

cosmosim::core::SimulationState makeGasContractState() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(4);
  state.resizeCells(3);
  state.species.count_by_species = {1, 3, 0, 0, 0};

  state.particle_sidecar.particle_id = {100, 200, 101, 102};
  state.particle_sidecar.species_tag = {
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas)};
  state.particle_sidecar.sfc_key = {30, 0, 20, 10};

  for (std::size_t cell = 0; cell < state.cells.size(); ++cell) {
    state.cells.center_x_comoving[cell] = 0.1 * static_cast<double>(cell + 1);
    state.cells.mass_code[cell] = 10.0 + static_cast<double>(cell);
    state.gas_cells.density_code[cell] = 100.0 + static_cast<double>(cell);
    state.gas_cells.pressure_code[cell] = 200.0 + static_cast<double>(cell);
    state.gas_cells.internal_energy_code[cell] = 300.0 + static_cast<double>(cell);
    state.gas_cells.sound_speed_code[cell] = 0.5 + static_cast<double>(cell);
  }

  state.rebuildSpeciesIndex();
  state.refreshGasCellIdentityFromParticleOrder();
  return state;
}

std::unordered_map<std::uint64_t, double> densityByGasParticleId(const cosmosim::core::SimulationState& state) {
  std::unordered_map<std::uint64_t, double> by_id;
  cosmosim::core::requireParticleBoundGasCellContract(state, "densityByGasParticleId");
  for (std::size_t cell = 0; cell < state.cells.size(); ++cell) {
    by_id.emplace(
        cosmosim::core::parentParticleIdForGasCellRow(state, static_cast<std::uint32_t>(cell)).value(),
        state.gas_cells.density_code[cell]);
  }
  return by_id;
}

void test_gas_cell_identity_invariants() {
  cosmosim::core::SimulationState state = makeGasContractState();
  cosmosim::core::requireParticleBoundGasCellContract(state, "test_gas_cell_identity_invariants");
  assert(state.gasCellIdentityMapMatchesParticleBoundState());
  state.requireGasCellIdentityMapCoversDenseRows("test_gas_cell_identity_invariants");
  assert(state.gas_cell_identity.coversDenseLocalRows(state.cells.size()));
  assert(state.gas_cell_identity.gasCellIdForLocalRow(0).has_value());
  assert(*state.gas_cell_identity.gasCellIdForLocalRow(0) == 100);
  assert(state.gas_cell_identity.rowForGasCellId(102).has_value());
  assert(*state.gas_cell_identity.rowForGasCellId(102) == 2);
  assert(state.parentParticleIdForGasCellRow(0).value() == 100);
  assert(state.parentParticleIdForGasCellRow(1).value() == 101);
  assert(state.gasCellRowForParticleId(102) == 2);
  assert(state.gasParticleIndexForCellRow(1) == 2);

  const auto baseline_density = densityByGasParticleId(state);
  assert(baseline_density.at(100) == 100.0);
  assert(baseline_density.at(101) == 101.0);
  assert(baseline_density.at(102) == 102.0);

  cosmosim::core::TransientStepWorkspace workspace;
  const std::array<std::uint32_t, 2> active_cells{0, 2};
  auto hydro_view = cosmosim::core::buildHydroCellKernelView(state, active_cells, workspace);
  hydro_view.density_code[0] = 555.0;
  hydro_view.pressure_code[1] = 777.0;
  cosmosim::core::scatterHydroCellKernelView(hydro_view, state);

  assert(state.gas_cells.density_code[0] == 555.0);
  assert(state.gas_cells.pressure_code[2] == 777.0);
  // Reconstruction scratch lanes are persistent sidecar fields and are not part of HydroCellKernelView.
  assert(state.gas_cells.sound_speed_code[0] == 0.5);
  assert(state.gas_cells.internal_energy_code[2] == 302.0);
}

void test_simulation_state_identity_map_materialization_and_drift_rejection() {
  cosmosim::core::SimulationState state = makeGasContractState();
  const std::uint64_t initial_generation = state.gasCellIdentityGeneration();
  assert(initial_generation > 0);

  state.refreshGasCellIdentityMapFromParticleBoundState();
  assert(state.gasCellIdentityGeneration() == initial_generation + 1);
  assert(state.gasCellIdentityMapMatchesParticleBoundState());
  state.requireGasCellIdentityMapFresh(state.gasCellIdentityGeneration(), "fresh generation");

  bool stale_generation_threw = false;
  try {
    state.requireGasCellIdentityMapFresh(initial_generation, "stale generation");
  } catch (const std::runtime_error&) {
    stale_generation_threw = true;
  }
  assert(stale_generation_threw);

  state.gas_cells.gas_cell_id[1] = 999999;
  assert(!state.gasCellIdentityMapMatchesParticleBoundState());
  bool contract_threw = false;
  try {
    cosmosim::core::requireParticleBoundGasCellContract(state, "drift rejection");
  } catch (const std::runtime_error&) {
    contract_threw = true;
  }
  assert(contract_threw);

  state.gas_cells.gas_cell_id[1] = 101;
  assert(state.gasCellIdentityMapMatchesParticleBoundState());

  state.gas_cells.parent_particle_id[0] = 0;
  bool zero_parent_threw = false;
  try {
    state.refreshGasCellIdentityMapFromParticleBoundState();
  } catch (const std::runtime_error&) {
    zero_parent_threw = true;
  }
  assert(zero_parent_threw);
}

void test_sidecar_identity_refresh_synchronizes_patch_ownership() {
  cosmosim::core::SimulationState state;
  state.resizeCells(2);
  state.resizePatches(1);
  state.patches.patch_id[0] = 42;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 2;
  state.patches.owning_rank[0] = 0;
  state.cells.patch_index = {0, 0};
  state.gas_cells.gas_cell_id = {7001, 7002};
  state.gas_cells.parent_particle_id = {0, 9001};

  state.gas_cell_identity.assign({
      {.gas_cell_id = 7001, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 0},
      {.gas_cell_id = 7002, .parent_particle_id = 9001, .owning_patch_id = 0, .local_cell_row = 1},
  });
  assert(!state.gasCellIdentityMapMatchesSidecarLanes());

  state.refreshGasCellIdentityMapFromSidecarLanes();
  assert(state.gasCellIdentityMapMatchesSidecarLanes());
  const auto* first = state.gas_cell_identity.findByLocalRow(0);
  const auto* second = state.gas_cell_identity.findByLocalRow(1);
  assert(first != nullptr && second != nullptr);
  assert(first->owning_patch_id == 42);
  assert(second->owning_patch_id == 42);
  assert(!first->parent_particle_id.has_value());
  assert(second->parent_particle_id.has_value());
  assert(*second->parent_particle_id == 9001);
}

void test_hydro_view_rejects_stale_gas_identity_generation() {
  cosmosim::core::SimulationState state = makeGasContractState();
  cosmosim::core::TransientStepWorkspace workspace;
  const std::array<std::uint32_t, 1> active_cell{1};
  auto view = cosmosim::core::buildHydroCellKernelView(state, active_cell, workspace);

  state.refreshGasCellIdentityMapFromParticleBoundState();
  assert(view.source_cell_index_generation == state.cellIndexGeneration());
  assert(view.source_gas_cell_identity_generation != state.gasCellIdentityGeneration());

  bool stale_identity_threw = false;
  try {
    cosmosim::core::scatterHydroCellKernelView(view, state);
  } catch (const std::runtime_error&) {
    stale_identity_threw = true;
  }
  assert(stale_identity_threw);
}

void test_hydro_view_scatters_by_stable_gas_cell_id_without_parent() {
  cosmosim::core::SimulationState state;
  state.resizeCells(3);
  for (std::size_t row = 0; row < state.cells.size(); ++row) {
    state.cells.center_x_comoving[row] = 10.0 + static_cast<double>(row);
    state.cells.mass_code[row] = 1.0 + static_cast<double>(row);
    state.gas_cells.density_code[row] = 100.0 + static_cast<double>(row);
    state.gas_cells.pressure_code[row] = 200.0 + static_cast<double>(row);
  }
  state.gas_cell_identity.assign({
      {.gas_cell_id = 7003, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 0},
      {.gas_cell_id = 7001, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 1},
      {.gas_cell_id = 7002, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 2},
  });

  cosmosim::core::TransientStepWorkspace workspace;
  const std::array<std::uint32_t, 1> active_cell{1};
  auto view = cosmosim::core::buildHydroCellKernelView(state, active_cell, workspace);
  assert(view.gas_cell_id[0] == 7001);
  assert(view.local_cell_row[0] == 1);

  view.cell_index[0] = 0;
  view.local_cell_row[0] = 0;
  view.density_code[0] = 444.0;
  view.pressure_code[0] = 555.0;
  cosmosim::core::scatterHydroCellKernelView(view, state);

  assert(state.gas_cells.density_code[0] == 100.0);
  assert(state.gas_cells.pressure_code[0] == 200.0);
  assert(state.gas_cells.density_code[1] == 444.0);
  assert(state.gas_cells.pressure_code[1] == 555.0);
}

void test_gas_cell_reorder_resize_invariants() {
  cosmosim::core::SimulationState state = makeGasContractState();
  const auto baseline_density = densityByGasParticleId(state);

  // Allowed reorder: gas relative order unchanged.
  const auto by_species = cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kBySpecies);
  cosmosim::core::reorderParticles(state, by_species);
  assert(densityByGasParticleId(state) == baseline_density);

  // Forbidden reorder: gas relative order changes without gas-cell rebuild.
  const auto by_sfc = cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kBySfcKey);
  bool reorder_threw = false;
  try {
    cosmosim::core::reorderParticles(state, by_sfc);
  } catch (const std::invalid_argument&) {
    reorder_threw = true;
  }
  assert(reorder_threw);

  cosmosim::core::TransientStepWorkspace workspace;
  const std::array<std::uint32_t, 1> active_cell{0};
  auto view = cosmosim::core::buildHydroCellKernelView(state, active_cell, workspace);
  state.resizeCells(2);

  bool stale_scatter_threw = false;
  try {
    cosmosim::core::scatterHydroCellKernelView(view, state);
  } catch (const std::runtime_error&) {
    stale_scatter_threw = true;
  }
  assert(stale_scatter_threw);

  bool contract_threw = false;
  try {
    cosmosim::core::requireParticleBoundGasCellContract(state, "test resized gas state");
  } catch (const std::runtime_error&) {
    contract_threw = true;
  }
  assert(contract_threw);

  bool active_view_threw = false;
  try {
    (void)cosmosim::core::buildHydroCellKernelView(state, active_cell, workspace);
  } catch (const std::runtime_error&) {
    active_view_threw = true;
  }
  assert(active_view_threw);
}


void test_decoupled_gas_cell_identity_map_api_shape() {
  cosmosim::core::GasCellIdentityMap map;
  assert(map.generation() == 0);

  map.assign({
      {.gas_cell_id = 9001, .parent_particle_id = 100, .owning_patch_id = 7, .local_cell_row = 0},
      {.gas_cell_id = 9002, .parent_particle_id = 100, .owning_patch_id = 7, .local_cell_row = 1},
      {.gas_cell_id = 9100, .parent_particle_id = std::nullopt, .owning_patch_id = 8, .local_cell_row = 2},
  });

  assert(map.size() == 3);
  assert(map.generation() == 1);
  assert(map.isConsistent());
  assert(map.coversDenseLocalRows(3));
  map.requireCoversDenseLocalRows(3, "test_decoupled_gas_cell_identity_map_api_shape");

  const auto* first_child = map.findByGasCellId(9001);
  assert(first_child != nullptr);
  assert(first_child->parent_particle_id.has_value());
  assert(*first_child->parent_particle_id == 100);
  assert(first_child->owning_patch_id == 7);
  assert(first_child->local_cell_row == 0);
  assert(map.rowForGasCellId(9002).has_value());
  assert(*map.rowForGasCellId(9002) == 1);

  const auto* second_child = map.findByGasCellId(9002);
  assert(second_child != nullptr);
  assert(second_child->parent_particle_id.has_value());
  assert(*second_child->parent_particle_id == 100);
  assert(second_child->local_cell_row == 1);

  const auto* parentless = map.findByLocalRow(2);
  assert(parentless != nullptr);
  assert(parentless->gas_cell_id == 9100);
  assert(!parentless->parent_particle_id.has_value());
  assert(parentless->owning_patch_id == 8);
  assert(map.gasCellIdForLocalRow(2).has_value());
  assert(*map.gasCellIdForLocalRow(2) == 9100);

  const auto parent_rows = map.rowsForParentParticleId(100);
  assert((parent_rows == std::vector<std::uint32_t>{0, 1}));
  const auto patch_rows = map.rowsForPatch(7);
  assert((patch_rows == std::vector<std::uint32_t>{0, 1}));

  bool duplicate_cell_id_threw = false;
  try {
    map.assign({
        {.gas_cell_id = 9001, .parent_particle_id = 100, .owning_patch_id = 7, .local_cell_row = 0},
        {.gas_cell_id = 9001, .parent_particle_id = 101, .owning_patch_id = 7, .local_cell_row = 1},
    });
  } catch (const std::invalid_argument&) {
    duplicate_cell_id_threw = true;
  }
  assert(duplicate_cell_id_threw);
  assert(map.generation() == 1);

  bool duplicate_local_row_threw = false;
  try {
    map.assign({
        {.gas_cell_id = 9201, .parent_particle_id = 100, .owning_patch_id = 7, .local_cell_row = 0},
        {.gas_cell_id = 9202, .parent_particle_id = 101, .owning_patch_id = 8, .local_cell_row = 0},
    });
  } catch (const std::invalid_argument&) {
    duplicate_local_row_threw = true;
  }
  assert(duplicate_local_row_threw);

  bool zero_cell_id_threw = false;
  try {
    map.assign({
        {.gas_cell_id = 0, .parent_particle_id = 100, .owning_patch_id = 7, .local_cell_row = 0},
    });
  } catch (const std::invalid_argument&) {
    zero_cell_id_threw = true;
  }
  assert(zero_cell_id_threw);

  cosmosim::core::GasCellIdentityMap sparse_rows;
  sparse_rows.assign({
      {.gas_cell_id = 9301, .parent_particle_id = 100, .owning_patch_id = 7, .local_cell_row = 0},
      {.gas_cell_id = 9302, .parent_particle_id = 101, .owning_patch_id = 7, .local_cell_row = 2},
  });
  assert(!sparse_rows.coversDenseLocalRows(2));
  bool dense_rows_threw = false;
  try {
    sparse_rows.requireCoversDenseLocalRows(2, "sparse_rows");
  } catch (const std::runtime_error&) {
    dense_rows_threw = true;
  }
  assert(dense_rows_threw);

  map.clear();
  assert(map.empty());
  assert(map.generation() == 2);
}

void test_gas_cell_identity_map_row_remap_by_stable_id() {
  cosmosim::core::GasCellIdentityMap old_map;
  old_map.assign({
      {.gas_cell_id = 7001, .parent_particle_id = 11, .owning_patch_id = 4, .local_cell_row = 0},
      {.gas_cell_id = 7002, .parent_particle_id = 11, .owning_patch_id = 4, .local_cell_row = 1},
      {.gas_cell_id = 7003, .parent_particle_id = 12, .owning_patch_id = 5, .local_cell_row = 2},
  });

  cosmosim::core::GasCellIdentityMap new_map;
  new_map.assign({
      {.gas_cell_id = 7003, .parent_particle_id = 12, .owning_patch_id = 5, .local_cell_row = 0},
      {.gas_cell_id = 7001, .parent_particle_id = 11, .owning_patch_id = 4, .local_cell_row = 1},
      {.gas_cell_id = 7999, .parent_particle_id = std::nullopt, .owning_patch_id = 6, .local_cell_row = 2},
  });

  const auto new_to_old = cosmosim::core::buildGasCellNewToOldRowMap(old_map, new_map);
  assert(new_to_old.size() == 3);
  assert(new_to_old[0] == 2);
  assert(new_to_old[1] == 0);
  assert(new_to_old[2] == cosmosim::core::kInvalidGasCellRow);
}

}  // namespace

int main() {
  test_gas_cell_identity_invariants();
  test_simulation_state_identity_map_materialization_and_drift_rejection();
  test_sidecar_identity_refresh_synchronizes_patch_ownership();
  test_hydro_view_rejects_stale_gas_identity_generation();
  test_hydro_view_scatters_by_stable_gas_cell_id_without_parent();
  test_gas_cell_reorder_resize_invariants();
  test_decoupled_gas_cell_identity_map_api_shape();
  test_gas_cell_identity_map_row_remap_by_stable_id();
  return 0;
}
