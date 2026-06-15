#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/amr/amr_framework.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace {

constexpr double k_tolerance = 1.0e-12;
constexpr double k_gamma = 1.4;
constexpr std::uint64_t k_register_key = 401U;

void seedRefinementHotspot(cosmosim::amr::AmrPatch& patch) {
  auto metrics = patch.metricsView();
  for (std::size_t i = 0; i < metrics.size(); ++i) {
    metrics[i].density_code = (i == 0) ? 16.0 : 0.01;
    metrics[i].sound_speed_code = 1.0;
    metrics[i].gradient_indicator = (i == 0) ? 2.0 : 0.0;
    metrics[i].particle_count = (i == 0) ? 16 : 0;
  }
}

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

[[nodiscard]] cosmosim::hydro::HydroPatchGeometry makeFluxRegisterGeometry(
    double face_area_comov,
    double normal_x,
    cosmosim::hydro::HydroFaceSide side,
    cosmosim::hydro::HydroFluxRegisterFaceRole role,
    double coarse_orientation_sign,
    std::uint64_t coarse_patch_id,
    std::size_t coarse_cell_index) {
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
      .area_comoving = face_area_comov,
      .normal_x = normal_x,
      .normal_y = 0.0,
      .normal_z = 0.0,
      .axis = cosmosim::hydro::HydroFaceAxis::kX});
  geometry.flux_register_faces.push_back(cosmosim::hydro::HydroFluxRegisterFace{
      .role = role,
      .register_key = k_register_key,
      .coarse_patch_id = coarse_patch_id,
      .coarse_cell_index = coarse_cell_index,
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

void recordSweep(
    cosmosim::hydro::HydroConservedStateSoa& conserved,
    const cosmosim::hydro::HydroPatchGeometry& geometry,
    cosmosim::amr::FluxRegisterAccumulator& accumulator) {
  const std::vector<std::size_t> active_cells{0U};
  const std::vector<std::size_t> active_faces{0U};
  const cosmosim::hydro::HydroUpdateContext update{.dt_code = 0.1, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::PiecewiseConstantReconstruction reconstruction;
  cosmosim::hydro::HlleRiemannSolver riemann;
  cosmosim::hydro::HydroScratchBuffers scratch;
  cosmosim::hydro::HydroPrimitiveCacheSoa cache(conserved.size());

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
      nullptr,
      &accumulator);
}

[[nodiscard]] std::vector<cosmosim::amr::FluxRegisterEntry> generateFluxEntriesFromHydroSweeps(
    std::uint64_t coarse_patch_id,
    std::size_t coarse_cell_index) {
  const auto coarse_geometry = makeFluxRegisterGeometry(
      0.25,
      1.0,
      cosmosim::hydro::HydroFaceSide::kUpper,
      cosmosim::hydro::HydroFluxRegisterFaceRole::kCoarse,
      1.0,
      coarse_patch_id,
      coarse_cell_index);
  const auto fine_geometry = makeFluxRegisterGeometry(
      0.125,
      -1.0,
      cosmosim::hydro::HydroFaceSide::kLower,
      cosmosim::hydro::HydroFluxRegisterFaceRole::kFine,
      -1.0,
      coarse_patch_id,
      coarse_cell_index);

  cosmosim::amr::FluxRegisterAccumulator accumulator;
  auto coarse_conserved = makeConserved(primitive(1.0, 0.4, 1.0), primitive(0.7, -0.1, 0.8));
  auto fine_conserved_a = makeConserved(primitive(1.1, -0.2, 1.2), primitive(1.0, 0.4, 1.0));
  auto fine_conserved_b = makeConserved(primitive(0.9, -0.05, 0.9), primitive(1.0, 0.4, 1.0));
  recordSweep(coarse_conserved, coarse_geometry, accumulator);
  recordSweep(fine_conserved_a, fine_geometry, accumulator);
  recordSweep(fine_conserved_b, fine_geometry, accumulator);
  return accumulator.entries();
}

}  // namespace

int main() {
  cosmosim::amr::PatchHierarchy hierarchy;
  cosmosim::amr::PatchDescriptor root;
  root.cell_dims = {4, 4, 4};
  const std::uint64_t root_id = hierarchy.createRootPatch(root);

  auto* root_patch = hierarchy.findPatch(root_id);
  assert(root_patch != nullptr);

  auto root_conserved = root_patch->conservedView();
  for (auto& cell : root_conserved) {
    cell.mass_code = 1.0;
    cell.total_energy_code = 2.0;
  }

  seedRefinementHotspot(*root_patch);

  cosmosim::amr::RefinementCriteria criteria;
  criteria.mass_threshold_code = 0.5;
  criteria.gradient_threshold = 1.0;
  criteria.particle_threshold = 8;
  criteria.jeans_resolution_cells = 0.5;

  cosmosim::amr::RefineDerefineManager manager(criteria);
  const auto regrid_diag = manager.regrid(hierarchy);
  assert(regrid_diag.refined_patch_count == 1);
  assert(hierarchy.levelView(1).size() == 8);

  const auto flux_entries = generateFluxEntriesFromHydroSweeps(root_id, 0U);
  assert(flux_entries.size() == 1U);
  assert(flux_entries.front().coarse_face_count == 1U);
  assert(flux_entries.front().fine_face_count == 2U);

  const auto reflux_diag = cosmosim::amr::RefluxSynchronizer::apply(hierarchy, flux_entries);
  assert(reflux_diag.corrected_cells == 1);

  root_patch = hierarchy.findPatch(root_id);
  assert(root_patch != nullptr);
  root_conserved = root_patch->conservedView();

  const auto& entry = flux_entries.front();
  const double expected_mass = 1.0 -
      ((entry.fine_face_flux_code.mass_code - entry.coarse_face_flux_code.mass_code) *
       entry.face_area_comov * entry.dt_code / root_patch->cellVolumeComov());
  const double expected_energy = 2.0 -
      ((entry.fine_face_flux_code.total_energy_code - entry.coarse_face_flux_code.total_energy_code) *
       entry.face_area_comov * entry.dt_code / root_patch->cellVolumeComov());

  assert(std::abs(root_conserved[0].mass_code - expected_mass) < k_tolerance);
  assert(std::abs(root_conserved[0].total_energy_code - expected_energy) < k_tolerance);

  return 0;
}
