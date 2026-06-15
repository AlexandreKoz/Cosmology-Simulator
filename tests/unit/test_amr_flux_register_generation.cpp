#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/amr/amr_framework.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace {

constexpr double k_gamma = 1.4;
constexpr double k_tol = 1.0e-12;
constexpr std::uint64_t k_register_key = 0x4d2U;
constexpr std::uint64_t k_coarse_patch_id = 77U;
constexpr std::size_t k_coarse_cell_index = 3U;

[[nodiscard]] cosmosim::hydro::HydroPrimitiveState primitive(
    double rho,
    double velocity_x,
    double pressure) {
  return cosmosim::hydro::HydroPrimitiveState{
      .rho_comoving = rho,
      .vel_x_peculiar = velocity_x,
      .vel_y_peculiar = 0.0,
      .vel_z_peculiar = 0.0,
      .pressure_comoving = pressure};
}

[[nodiscard]] cosmosim::hydro::HydroPatchGeometry makeOneFaceGeometry(
    double area_comoving,
    double normal_x,
    cosmosim::hydro::HydroFaceSide side,
    cosmosim::hydro::HydroFluxRegisterFaceRole role,
    double coarse_orientation_sign) {
  cosmosim::hydro::HydroPatchGeometry geometry;
  geometry.nx = 1;
  geometry.ny = 1;
  geometry.nz = 1;
  geometry.cell_width_x_comoving = 1.0;
  geometry.cell_width_y_comoving = 1.0;
  geometry.cell_width_z_comoving = 1.0;
  geometry.cell_volume_comoving = 1.0;
  geometry.ghost_cells.push_back(cosmosim::hydro::HydroGhostCell{
      .owner_real_cell = 0,
      .source_real_cell = cosmosim::hydro::k_invalid_cell_index,
      .ghost_cell = 1,
      .ghost_slot = 0,
      .boundary_kind = cosmosim::hydro::HydroBoundaryKind::kImportedMpi,
      .axis = cosmosim::hydro::HydroFaceAxis::kX,
      .side = side,
      .mutation_rights = cosmosim::hydro::HydroGhostMutationRights::kReadOnlyImported});
  geometry.faces.push_back(cosmosim::hydro::HydroFace{
      .owner_cell = 0,
      .neighbor_cell = 1,
      .owner_minus_cell = cosmosim::hydro::k_invalid_cell_index,
      .neighbor_plus_cell = cosmosim::hydro::k_invalid_cell_index,
      .ghost_cell_slot = 0,
      .area_comoving = area_comoving,
      .normal_x = normal_x,
      .normal_y = 0.0,
      .normal_z = 0.0,
      .axis = cosmosim::hydro::HydroFaceAxis::kX});
  geometry.flux_register_faces.push_back(cosmosim::hydro::HydroFluxRegisterFace{
      .role = role,
      .register_key = k_register_key,
      .coarse_patch_id = k_coarse_patch_id,
      .coarse_cell_index = k_coarse_cell_index,
      .level = 0,
      .axis = cosmosim::hydro::HydroFaceAxis::kX,
      .orientation = cosmosim::hydro::HydroFaceSide::kUpper,
      .coarse_orientation_sign = coarse_orientation_sign});
  return geometry;
}

[[nodiscard]] cosmosim::hydro::HydroConservedStateSoa makeConserved(
    const cosmosim::hydro::HydroPrimitiveState& owner,
    const cosmosim::hydro::HydroPrimitiveState& ghost) {
  cosmosim::hydro::HydroConservedStateSoa conserved(2);
  conserved.storeCell(0, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(owner, k_gamma));
  conserved.storeCell(1, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(ghost, k_gamma));
  return conserved;
}

void runSweep(
    cosmosim::hydro::HydroConservedStateSoa& conserved,
    const cosmosim::hydro::HydroPatchGeometry& geometry,
    cosmosim::amr::FluxRegisterAccumulator& accumulator) {
  const std::vector<std::size_t> active_cells{0U};
  const std::vector<std::size_t> active_faces{0U};
  const cosmosim::hydro::HydroUpdateContext update{
      .dt_code = 0.125,
      .scale_factor = 1.0,
      .hubble_rate_code = 0.0};
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::PiecewiseConstantReconstruction reconstruction;
  cosmosim::hydro::HlleRiemannSolver riemann;
  cosmosim::hydro::HydroScratchBuffers scratch;
  cosmosim::hydro::HydroPrimitiveCacheSoa cache(conserved.size());
  cosmosim::hydro::HydroProfileEvent profile;

  solver.advancePatchActiveSetWithScratch(
      conserved,
      geometry,
      cosmosim::hydro::HydroActiveSetView{.active_cells = active_cells, .active_faces = active_faces},
      update,
      reconstruction,
      riemann,
      {},
      source_context,
      scratch,
      &cache,
      &profile,
      &accumulator);
  assert(profile.face_count == 1U);
}

void assertNear(double actual, double expected) {
  assert(std::abs(actual - expected) < k_tol);
}

void assertConservedNear(
    const cosmosim::amr::ConservedState& actual,
    const cosmosim::hydro::HydroConservedState& expected) {
  assertNear(actual.mass_code, expected.mass_density_comoving);
  assertNear(actual.momentum_x_code, expected.momentum_density_x_comoving);
  assertNear(actual.momentum_y_code, expected.momentum_density_y_comoving);
  assertNear(actual.momentum_z_code, expected.momentum_density_z_comoving);
  assertNear(actual.total_energy_code, expected.total_energy_density_comoving);
}

void testActualSweepFluxesGenerateRegisterEntry() {
  const auto coarse_owner = primitive(1.0, 0.7, 1.1);
  const auto coarse_ghost = primitive(0.8, -0.2, 0.9);
  const auto fine_owner_a = primitive(1.2, -0.3, 1.3);
  const auto fine_owner_b = primitive(0.9, -0.1, 0.7);
  const auto fine_coarse_ghost = primitive(1.0, 0.7, 1.1);

  const auto coarse_geometry = makeOneFaceGeometry(
      1.0,
      1.0,
      cosmosim::hydro::HydroFaceSide::kUpper,
      cosmosim::hydro::HydroFluxRegisterFaceRole::kCoarse,
      1.0);
  const auto fine_geometry = makeOneFaceGeometry(
      0.5,
      -1.0,
      cosmosim::hydro::HydroFaceSide::kLower,
      cosmosim::hydro::HydroFluxRegisterFaceRole::kFine,
      -1.0);

  cosmosim::hydro::HlleRiemannSolver expected_riemann;
  const auto expected_coarse_flux = expected_riemann.computeFlux(
      coarse_owner,
      coarse_ghost,
      coarse_geometry.faces.front(),
      k_gamma);
  const auto expected_fine_flux_a =
      -1.0 * expected_riemann.computeFlux(fine_owner_a, fine_coarse_ghost, fine_geometry.faces.front(), k_gamma);
  const auto expected_fine_flux_b =
      -1.0 * expected_riemann.computeFlux(fine_owner_b, fine_coarse_ghost, fine_geometry.faces.front(), k_gamma);
  const auto expected_fine_average = 0.5 * (expected_fine_flux_a + expected_fine_flux_b);

  cosmosim::amr::FluxRegisterAccumulator accumulator;
  auto coarse_conserved = makeConserved(coarse_owner, coarse_ghost);
  auto fine_conserved_a = makeConserved(fine_owner_a, fine_coarse_ghost);
  auto fine_conserved_b = makeConserved(fine_owner_b, fine_coarse_ghost);
  runSweep(coarse_conserved, coarse_geometry, accumulator);
  runSweep(fine_conserved_a, fine_geometry, accumulator);
  runSweep(fine_conserved_b, fine_geometry, accumulator);

  const auto entries = accumulator.entries();
  assert(entries.size() == 1U);
  const auto& entry = entries.front();
  assert(entry.register_key == k_register_key);
  assert(entry.coarse_patch_id == k_coarse_patch_id);
  assert(entry.coarse_cell_index == k_coarse_cell_index);
  assert(entry.level == 0U);
  assert(entry.axis == cosmosim::hydro::HydroFaceAxis::kX);
  assert(entry.orientation == cosmosim::hydro::HydroFaceSide::kUpper);
  assert(entry.coarse_face_count == 1U);
  assert(entry.fine_face_count == 2U);
  assertNear(entry.face_area_comov, 1.0);
  assertNear(entry.dt_code, 0.125);
  assertConservedNear(entry.coarse_face_flux_code, expected_coarse_flux);
  assertConservedNear(entry.fine_face_flux_code, expected_fine_average);
}

void testSameLevelFacesDoNotGenerateRegisterEntries() {
  auto geometry = makeOneFaceGeometry(
      1.0,
      1.0,
      cosmosim::hydro::HydroFaceSide::kUpper,
      cosmosim::hydro::HydroFluxRegisterFaceRole::kCoarse,
      1.0);
  geometry.flux_register_faces.clear();
  auto conserved = makeConserved(primitive(1.0, 0.1, 1.0), primitive(1.0, -0.1, 1.0));
  cosmosim::amr::FluxRegisterAccumulator accumulator;
  runSweep(conserved, geometry, accumulator);
  assert(accumulator.entries().empty());
}

}  // namespace

int main() {
  testActualSweepFluxesGenerateRegisterEntry();
  testSameLevelFacesDoNotGenerateRegisterEntries();
  return 0;
}
