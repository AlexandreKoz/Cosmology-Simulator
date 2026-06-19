#include <cassert>
#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

#include "cosmosim/amr/amr_ghost_fill.hpp"
#include "cosmosim/amr/amr_hydro_orchestrator.hpp"

namespace {

constexpr double k_gamma = 1.4;
constexpr double k_tol = 1.0e-12;

void setPatch(
    cosmosim::core::SimulationState& state,
    std::size_t patch_row,
    const cosmosim::amr::PatchDescriptor& patch,
    std::uint32_t first_cell) {
  cosmosim::amr::writePatchDescriptorToStateRow(state, patch_row, patch);
  state.patches.first_cell[patch_row] = first_cell;
  state.patches.cell_count[patch_row] = static_cast<std::uint32_t>(
      static_cast<std::size_t>(patch.cell_dims[0]) * patch.cell_dims[1] * patch.cell_dims[2]);
  state.patches.owning_rank[patch_row] = 0;
}

void setCell(
    cosmosim::core::SimulationState& state,
    std::uint32_t row,
    std::uint32_t patch_row,
    std::uint64_t patch_id,
    std::uint64_t gas_cell_id,
    double x,
    double rho,
    double pressure,
    std::vector<cosmosim::core::GasCellIdentityRecord>& identities) {
  state.cells.patch_index[row] = patch_row;
  state.cells.center_x_comoving[row] = x;
  state.cells.center_y_comoving[row] = 0.5;
  state.cells.center_z_comoving[row] = 0.5;
  state.cells.mass_code[row] = rho;
  state.gas_cells.gas_cell_id[row] = gas_cell_id;
  state.gas_cells.parent_particle_id[row] = 0;
  state.gas_cells.density_code[row] = rho;
  state.gas_cells.pressure_code[row] = pressure;
  state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
  state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
  state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * pressure / rho);
  state.gas_cells.velocity_x_peculiar[row] = 0.0;
  state.gas_cells.velocity_y_peculiar[row] = 0.0;
  state.gas_cells.velocity_z_peculiar[row] = 0.0;
  identities.push_back(cosmosim::core::GasCellIdentityRecord{
      .gas_cell_id = gas_cell_id,
      .parent_particle_id = std::nullopt,
      .owning_patch_id = patch_id,
      .local_cell_row = row});
}

[[nodiscard]] cosmosim::core::SimulationState makeState() {
  cosmosim::core::SimulationState state;
  state.resizeCells(4);
  state.resizePatches(2);
  const cosmosim::amr::PatchDescriptor coarse{
      .patch_id = 101,
      .parent_patch_id = 0,
      .level = 0,
      .morton_key = 101,
      .origin_comov = {0.0, 0.0, 0.0},
      .extent_comov = {1.0, 1.0, 1.0},
      .cell_dims = {2, 1, 1}};
  const cosmosim::amr::PatchDescriptor fine{
      .patch_id = 201,
      .parent_patch_id = 101,
      .level = 1,
      .morton_key = 201,
      .origin_comov = {1.0, 0.0, 0.0},
      .extent_comov = {0.5, 1.0, 1.0},
      .cell_dims = {2, 1, 1}};
  setPatch(state, 0, coarse, 0);
  setPatch(state, 1, fine, 2);
  std::vector<cosmosim::core::GasCellIdentityRecord> identities;
  setCell(state, 0, 0, 101, 9001, 0.25, 1.0, 1.0, identities);
  setCell(state, 1, 0, 101, 9002, 0.75, 1.0, 1.0, identities);
  setCell(state, 2, 1, 201, 9101, 1.125, 0.8, 0.8, identities);
  setCell(state, 3, 1, 201, 9102, 1.375, 0.7, 0.7, identities);
  state.gas_cell_identity.assign(std::move(identities));
  return state;
}

[[nodiscard]] cosmosim::amr::AmrHydroGeometryOptions coarseOptions() {
  cosmosim::amr::AmrHydroGeometryOptions options;
  options.boundary_classes[1] = cosmosim::amr::AmrHydroBoundaryClass::kCoarseFine;
  return options;
}

[[nodiscard]] cosmosim::amr::AmrHydroGeometryOptions fineOptions() {
  cosmosim::amr::AmrHydroGeometryOptions options;
  options.boundary_classes[0] = cosmosim::amr::AmrHydroBoundaryClass::kCoarseFine;
  return options;
}

[[nodiscard]] cosmosim::amr::AmrHydroGhostDescriptor& lowerXGhost(
    cosmosim::amr::AmrHydroPatchGeometry& geometry) {
  for (auto& descriptor : geometry.ghosts) {
    const auto& ghost = geometry.geometry.ghost_cells[descriptor.ghost_slot];
    if (ghost.axis == cosmosim::hydro::HydroFaceAxis::kX &&
        ghost.side == cosmosim::hydro::HydroFaceSide::kLower &&
        ghost.owner_real_cell == 0U) {
      return descriptor;
    }
  }
  assert(false);
  return geometry.ghosts.front();
}

[[nodiscard]] cosmosim::amr::AmrHydroGhostDescriptor& upperXGhost(
    cosmosim::amr::AmrHydroPatchGeometry& geometry) {
  for (auto& descriptor : geometry.ghosts) {
    const auto& ghost = geometry.geometry.ghost_cells[descriptor.ghost_slot];
    if (ghost.axis == cosmosim::hydro::HydroFaceAxis::kX &&
        ghost.side == cosmosim::hydro::HydroFaceSide::kUpper &&
        ghost.owner_real_cell == 1U) {
      return descriptor;
    }
  }
  assert(false);
  return geometry.ghosts.front();
}

struct Fixture {
  cosmosim::core::SimulationState state = makeState();
  std::vector<cosmosim::amr::PatchDescriptor> descriptors =
      cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  cosmosim::amr::AmrHydroPatchGeometry coarse =
      cosmosim::amr::buildAmrHydroPatchGeometry(state, descriptors[0], coarseOptions());
  cosmosim::amr::AmrHydroPatchGeometry fine =
      cosmosim::amr::buildAmrHydroPatchGeometry(state, descriptors[1], fineOptions());
  cosmosim::hydro::HydroConservedStateSoa coarse_conserved =
      cosmosim::amr::loadAmrHydroConservedState(state, coarse, k_gamma);
  cosmosim::hydro::HydroConservedStateSoa fine_conserved =
      cosmosim::amr::loadAmrHydroConservedState(state, fine, k_gamma);

  void captureInterval() {
    cosmosim::amr::captureAmrTemporalBoundaryHistoryStart(state, descriptors, 0.0, k_gamma);
    state.gas_cells.density_code[1] = 3.0;
    state.gas_cells.pressure_code[1] = 1.0;
    state.gas_cells.internal_energy_code[1] = 1.0 / ((k_gamma - 1.0) * 3.0);
    state.gas_cells.temperature_code[1] = state.gas_cells.internal_energy_code[1];
    state.gas_cells.sound_speed_code[1] = std::sqrt(k_gamma / 3.0);
    cosmosim::amr::captureAmrTemporalBoundaryHistoryEnd(state, descriptors, 1.0, k_gamma);
  }

  [[nodiscard]] std::vector<cosmosim::amr::AmrHydroGhostFillPatch> views(double fill_time) {
    return {
        {.geometry = &coarse,
         .conserved = &coarse_conserved,
         .target_state_time_code = 0.0,
         .ghost_fill_time_code = fill_time,
         .source_current_state_time_code = 1.0,
         .temporal_boundary_history = &state.amr_temporal_boundary_history,
         .enable_temporal_coarse_to_fine = true},
        {.geometry = &fine,
         .conserved = &fine_conserved,
         .target_state_time_code = fill_time,
         .ghost_fill_time_code = fill_time,
         .source_current_state_time_code = fill_time,
         .temporal_boundary_history = &state.amr_temporal_boundary_history,
         .enable_temporal_coarse_to_fine = true}};
  }
};

void testTemporalEndpointsAndMidpoint() {
  Fixture fixture;
  fixture.captureInterval();

  auto views = fixture.views(0.0);
  auto diagnostics = cosmosim::amr::fillAmrHydroGhostCells(views, k_gamma);
  auto& ghost = lowerXGhost(fixture.fine);
  assert(ghost.fill_status == cosmosim::amr::AmrHydroGhostFillStatus::kFilledCoarseToFineTemporal);
  const auto start = fixture.fine_conserved.loadCell(ghost.ghost_cell);
  assert(std::abs(start.mass_density_comoving - 1.0) < k_tol);
  assert(diagnostics.temporal_endpoint_ghosts_filled > 0U);

  views = fixture.views(0.5);
  diagnostics = cosmosim::amr::fillAmrHydroGhostCells(views, k_gamma);
  const auto midpoint = fixture.fine_conserved.loadCell(ghost.ghost_cell);
  assert(std::abs(midpoint.mass_density_comoving - 2.0) < k_tol);
  assert(std::abs(midpoint.total_energy_density_comoving - 2.5) < k_tol);
  const auto primitive = cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(midpoint, k_gamma);
  assert(std::isfinite(primitive.pressure_comoving));
  assert(primitive.rho_comoving > 0.0 && primitive.pressure_comoving > 0.0);
  assert(diagnostics.temporal_coarse_to_fine_ghosts_filled > 0U);

  views = fixture.views(1.0);
  diagnostics = cosmosim::amr::fillAmrHydroGhostCells(views, k_gamma);
  const auto end = fixture.fine_conserved.loadCell(ghost.ghost_cell);
  assert(std::abs(end.mass_density_comoving - 3.0) < k_tol);
  assert(diagnostics.temporal_endpoint_ghosts_filled > 0U);
}

void testRejectsOutOfRangeAndStaleHistory() {
  Fixture fixture;
  fixture.captureInterval();
  auto views = fixture.views(1.1);
  bool threw = false;
  try {
    (void)cosmosim::amr::fillAmrHydroGhostCells(views, k_gamma);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
  assert(lowerXGhost(fixture.fine).fill_status ==
      cosmosim::amr::AmrHydroGhostFillStatus::kRejectedTemporalTimeOutOfRange);

  Fixture geometry_fixture;
  geometry_fixture.captureInterval();
  geometry_fixture.coarse.patch.extent_comov[0] = 0.9;
  views = geometry_fixture.views(0.5);
  threw = false;
  try {
    (void)cosmosim::amr::fillAmrHydroGhostCells(views, k_gamma);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
  assert(lowerXGhost(geometry_fixture.fine).fill_status ==
      cosmosim::amr::AmrHydroGhostFillStatus::kRejectedTemporalHistoryInvalid);

  Fixture identity_fixture;
  identity_fixture.captureInterval();
  ++identity_fixture.coarse.source_gas_cell_identity_generation;
  views = identity_fixture.views(0.5);
  threw = false;
  try {
    (void)cosmosim::amr::fillAmrHydroGhostCells(views, k_gamma);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
  assert(lowerXGhost(identity_fixture.fine).fill_status ==
      cosmosim::amr::AmrHydroGhostFillStatus::kRejectedTemporalHistoryInvalid);

}

void testRejectsUnsynchronizedFineToCoarse() {
  Fixture fixture;
  fixture.captureInterval();
  auto views = fixture.views(0.0);
  views[0].target_state_time_code = 0.0;
  views[0].ghost_fill_time_code = 0.0;
  views[0].source_current_state_time_code = 0.0;
  views[1].source_current_state_time_code = 0.5;
  bool threw = false;
  try {
    (void)cosmosim::amr::fillAmrHydroGhostCells(views, k_gamma);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
  assert(upperXGhost(fixture.coarse).fill_status ==
      cosmosim::amr::AmrHydroGhostFillStatus::kRejectedTemporalFineToCoarseMismatch);
}

void testActiveHistoryBlocksTopologyMutation() {
  Fixture fixture;
  fixture.captureInterval();
  bool threw = false;
  try {
    (void)cosmosim::amr::refineProductionPatchInSimulationState(
        fixture.state, fixture.descriptors[0], 1000, 20000);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}


}  // namespace

int main() {
  testTemporalEndpointsAndMidpoint();
  testRejectsOutOfRangeAndStaleHistory();
  testRejectsUnsynchronizedFineToCoarse();
  testActiveHistoryBlocksTopologyMutation();
  return 0;
}
