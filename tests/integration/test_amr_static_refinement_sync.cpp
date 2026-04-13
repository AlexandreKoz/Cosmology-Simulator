#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/amr/amr_framework.hpp"

namespace {

constexpr double k_tolerance = 1.0e-12;

void seedRefinementHotspot(cosmosim::amr::AmrPatch& patch) {
  auto metrics = patch.metricsView();
  for (std::size_t i = 0; i < metrics.size(); ++i) {
    metrics[i].density_code = (i == 0) ? 16.0 : 0.01;
    metrics[i].sound_speed_code = 1.0;
    metrics[i].gradient_indicator = (i == 0) ? 2.0 : 0.0;
    metrics[i].particle_count = (i == 0) ? 16 : 0;
  }
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

  // Synchronize one coarse-fine interface through explicit flux register entries.
  std::vector<cosmosim::amr::FluxRegisterEntry> flux_entries;
  flux_entries.push_back(cosmosim::amr::FluxRegisterEntry{
      .coarse_patch_id = root_id,
      .coarse_cell_index = 0,
      .coarse_face_flux_code = {.mass_code = 2.0, .total_energy_code = 4.0},
      .fine_face_flux_code = {.mass_code = 2.5, .total_energy_code = 4.8},
      .face_area_comov = 0.25,
      .dt_code = 0.1,
  });

  const auto reflux_diag = cosmosim::amr::RefluxSynchronizer::apply(hierarchy, flux_entries);
  assert(reflux_diag.corrected_cells == 1);

  root_patch = hierarchy.findPatch(root_id);
  assert(root_patch != nullptr);
  root_conserved = root_patch->conservedView();

  const double expected_mass = 1.0 - ((2.5 - 2.0) * 0.25 * 0.1 / root_patch->cellVolumeComov());
  const double expected_energy = 2.0 - ((4.8 - 4.0) * 0.25 * 0.1 / root_patch->cellVolumeComov());

  assert(std::abs(root_conserved[0].mass_code - expected_mass) < k_tolerance);
  assert(std::abs(root_conserved[0].total_energy_code - expected_energy) < k_tolerance);

  return 0;
}
