#include <cassert>
#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

#include "cosmosim/amr/amr_hydro_geometry.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace {

constexpr double k_tol = 1.0e-12;

cosmosim::core::SimulationState makePatchState() {
  cosmosim::core::SimulationState state;
  state.resizeCells(8);
  state.resizePatches(1);
  state.patches.patch_id[0] = 31;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 8;
  state.patches.owning_rank[0] = 0;

  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(8);
  for (std::uint32_t row = 0; row < 8; ++row) {
    state.cells.patch_index[row] = 0;
    state.cells.mass_code[row] = 1.0;
    state.gas_cells.gas_cell_id[row] = 9001U + row;
    state.gas_cells.parent_particle_id[row] = 0;
    state.gas_cells.density_code[row] = 1.0 + 0.01 * static_cast<double>(row);
    state.gas_cells.pressure_code[row] = 1.0;
    state.gas_cells.velocity_x_peculiar[row] = 0.01 * static_cast<double>(row);
    state.gas_cells.velocity_y_peculiar[row] = 0.0;
    state.gas_cells.velocity_z_peculiar[row] = 0.0;
    records.push_back(cosmosim::core::GasCellIdentityRecord{
        .gas_cell_id = 9001U + row,
        .parent_particle_id = std::nullopt,
        .owning_patch_id = 31,
        .local_cell_row = row});
  }
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

cosmosim::amr::PatchDescriptor makeDescriptor() {
  cosmosim::amr::PatchDescriptor patch;
  patch.patch_id = 31;
  patch.level = 0;
  patch.origin_comov = {1.0, 2.0, 3.0};
  patch.extent_comov = {1.0, 2.0, 4.0};
  patch.cell_dims = {2, 2, 2};
  return patch;
}

void testPatchGeometryCountsAndIdentityCoverage() {
  const cosmosim::core::SimulationState state = makePatchState();
  const auto patch = makeDescriptor();

  cosmosim::amr::AmrHydroGeometryOptions options;
  options.boundary_classes = {
      cosmosim::amr::AmrHydroBoundaryClass::kPhysical,
      cosmosim::amr::AmrHydroBoundaryClass::kSameLevel,
      cosmosim::amr::AmrHydroBoundaryClass::kCoarseFine,
      cosmosim::amr::AmrHydroBoundaryClass::kPhysical,
      cosmosim::amr::AmrHydroBoundaryClass::kSameLevel,
      cosmosim::amr::AmrHydroBoundaryClass::kCoarseFine,
  };

  const cosmosim::amr::AmrHydroPatchGeometry view =
      cosmosim::amr::buildAmrHydroPatchGeometry(state, patch, options);

  assert(view.geometry.cellCount() == 8U);
  assert(view.real_cells.size() == 8U);
  assert(view.gasCellIds().size() == 8U);
  assert(view.gas_cell_ids.front() == 9001U);
  assert(view.gas_cell_ids.back() == 9008U);
  assert(view.local_cell_rows.front() == 0U);
  assert(view.local_cell_rows.back() == 7U);
  assert(std::abs(view.geometry.cell_width_x_comoving - 0.5) < k_tol);
  assert(std::abs(view.geometry.cell_width_y_comoving - 1.0) < k_tol);
  assert(std::abs(view.geometry.cell_width_z_comoving - 2.0) < k_tol);
  assert(std::abs(view.geometry.cell_volume_comoving - 1.0) < k_tol);

  assert(view.internalFaceIndices().size() == 12U);
  assert(view.ghosts.size() == 24U);
  assert(view.geometry.ghost_cells.size() == 24U);
  assert(view.geometry.faces.size() == 36U);
  assert(view.faces.size() == view.geometry.faces.size());

  std::size_t x_faces = 0;
  std::size_t y_faces = 0;
  std::size_t z_faces = 0;
  std::size_t physical_ghosts = 0;
  std::size_t same_level_ghosts = 0;
  std::size_t coarse_fine_ghosts = 0;
  for (const cosmosim::amr::AmrHydroFaceDescriptor& face : view.faces) {
    if (face.axis == cosmosim::hydro::HydroFaceAxis::kX) {
      ++x_faces;
      assert(std::abs(face.area_comoving - 2.0) < k_tol);
      assert(std::abs(std::abs(face.normal_x) - 1.0) < k_tol);
      assert(std::abs(face.normal_y) < k_tol);
      assert(std::abs(face.normal_z) < k_tol);
    } else if (face.axis == cosmosim::hydro::HydroFaceAxis::kY) {
      ++y_faces;
      assert(std::abs(face.area_comoving - 1.0) < k_tol);
      assert(std::abs(face.normal_x) < k_tol);
      assert(std::abs(std::abs(face.normal_y) - 1.0) < k_tol);
      assert(std::abs(face.normal_z) < k_tol);
    } else {
      ++z_faces;
      assert(std::abs(face.area_comoving - 0.5) < k_tol);
      assert(std::abs(face.normal_x) < k_tol);
      assert(std::abs(face.normal_y) < k_tol);
      assert(std::abs(std::abs(face.normal_z) - 1.0) < k_tol);
    }
    assert(face.patch_id == 31U);
    assert(face.face_id != 0U);
    assert(face.owner_gas_cell_id >= 9001U && face.owner_gas_cell_id <= 9008U);
  }
  assert(x_faces == 12U);
  assert(y_faces == 12U);
  assert(z_faces == 12U);

  for (const cosmosim::amr::AmrHydroGhostDescriptor& ghost : view.ghosts) {
    assert(ghost.patch_id == 31U);
    assert(ghost.owner_gas_cell_id >= 9001U && ghost.owner_gas_cell_id <= 9008U);
    if (ghost.boundary_class == cosmosim::amr::AmrHydroBoundaryClass::kPhysical) {
      ++physical_ghosts;
      assert(ghost.fill_status == cosmosim::amr::AmrHydroGhostFillStatus::kUnfilledPhysicalBoundary);
    } else if (ghost.boundary_class == cosmosim::amr::AmrHydroBoundaryClass::kSameLevel) {
      ++same_level_ghosts;
      assert(ghost.fill_status == cosmosim::amr::AmrHydroGhostFillStatus::kUnfilledSameLevel);
    } else {
      ++coarse_fine_ghosts;
      assert(ghost.fill_status == cosmosim::amr::AmrHydroGhostFillStatus::kUnfilledCoarseFine);
    }
  }
  assert(physical_ghosts == 8U);
  assert(same_level_ghosts == 8U);
  assert(coarse_fine_ghosts == 8U);
}

void testRejectsStaleIdentityGeneration() {
  cosmosim::core::SimulationState state = makePatchState();
  const auto view = cosmosim::amr::buildAmrHydroPatchGeometry(state, makeDescriptor());

  state.gas_cell_identity.assign({
      {.gas_cell_id = 9101, .parent_particle_id = std::nullopt, .owning_patch_id = 31, .local_cell_row = 0},
      {.gas_cell_id = 9102, .parent_particle_id = std::nullopt, .owning_patch_id = 31, .local_cell_row = 1},
      {.gas_cell_id = 9103, .parent_particle_id = std::nullopt, .owning_patch_id = 31, .local_cell_row = 2},
      {.gas_cell_id = 9104, .parent_particle_id = std::nullopt, .owning_patch_id = 31, .local_cell_row = 3},
      {.gas_cell_id = 9105, .parent_particle_id = std::nullopt, .owning_patch_id = 31, .local_cell_row = 4},
      {.gas_cell_id = 9106, .parent_particle_id = std::nullopt, .owning_patch_id = 31, .local_cell_row = 5},
      {.gas_cell_id = 9107, .parent_particle_id = std::nullopt, .owning_patch_id = 31, .local_cell_row = 6},
      {.gas_cell_id = 9108, .parent_particle_id = std::nullopt, .owning_patch_id = 31, .local_cell_row = 7},
  });

  bool threw = false;
  try {
    (void)cosmosim::amr::loadAmrHydroConservedState(state, view, 1.4);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}

void testHydroSolverAdvancesPatchThroughAdapter() {
  const cosmosim::core::SimulationState state = makePatchState();
  const auto view = cosmosim::amr::buildAmrHydroPatchGeometry(state, makeDescriptor());
  cosmosim::hydro::HydroConservedStateSoa conserved =
      cosmosim::amr::loadAmrHydroConservedState(state, view, 1.4);
  const cosmosim::hydro::HydroConservedState before = conserved.loadCell(0);

  const std::vector<std::size_t> active_cells{0U, 1U, 2U, 3U, 4U, 5U, 6U, 7U};
  const std::vector<std::size_t> active_faces = view.internalFaceIndices();
  cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-4, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::HydroCoreSolver solver(1.4);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code / view.geometry.cell_width_x_comoving,
      .dt_over_cell_width_code = {
          update.dt_code / view.geometry.cell_width_x_comoving,
          update.dt_code / view.geometry.cell_width_y_comoving,
          update.dt_code / view.geometry.cell_width_z_comoving},
      .rho_floor = 1.0e-12,
      .pressure_floor = 1.0e-12,
      .enable_muscl_hancock_predictor = true,
      .adiabatic_index = 1.4});
  cosmosim::hydro::HllcRiemannSolver riemann;
  cosmosim::hydro::HydroScratchBuffers scratch;
  cosmosim::hydro::HydroPrimitiveCacheSoa cache(conserved.size());
  cosmosim::hydro::HydroProfileEvent profile;

  solver.advancePatchActiveSetWithScratch(
      conserved,
      view.geometry,
      cosmosim::hydro::HydroActiveSetView{.active_cells = active_cells, .active_faces = active_faces},
      update,
      reconstruction,
      riemann,
      {},
      source_context,
      scratch,
      &cache,
      &profile);

  const cosmosim::hydro::HydroConservedState after = conserved.loadCell(0);
  assert(profile.face_count == active_faces.size());
  assert(std::abs(after.mass_density_comoving - before.mass_density_comoving) > 0.0);
}

}  // namespace

int main() {
  testPatchGeometryCountsAndIdentityCoverage();
  testRejectsStaleIdentityGeneration();
  testHydroSolverAdvancesPatchThroughAdapter();
  return 0;
}
