#include <cassert>
#include <cmath>

#include "cosmosim/amr/amr_framework.hpp"

namespace {

constexpr double k_tolerance = 1.0e-12;

void testRefineAndDerefineLifecycle() {
  cosmosim::amr::PatchHierarchy hierarchy;
  cosmosim::amr::PatchDescriptor root;
  root.cell_dims = {2, 2, 2};
  const std::uint64_t root_id = hierarchy.createRootPatch(root);

  auto* root_patch = hierarchy.findPatch(root_id);
  assert(root_patch != nullptr);

  for (auto& metric : root_patch->metricsView()) {
    metric.density_code = 8.0;
    metric.sound_speed_code = 0.5;
    metric.gradient_indicator = 0.0;
    metric.particle_count = 1;
  }

  cosmosim::amr::RefinementCriteria criteria;
  criteria.mass_threshold_code = 0.2;
  criteria.gradient_threshold = 1.0;
  criteria.particle_threshold = 64;
  criteria.jeans_resolution_cells = 0.5;

  cosmosim::amr::RefineDerefineManager manager(criteria);
  const auto refine_diag = manager.regrid(hierarchy);
  assert(refine_diag.refined_patch_count == 1);
  assert(hierarchy.levelCount() == 2);
  assert(hierarchy.levelView(1).size() == 8);

  for (auto& child_patch : hierarchy.levelView(1)) {
    for (auto& metric : child_patch.metricsView()) {
      metric.density_code = 0.01;
      metric.sound_speed_code = 10.0;
      metric.gradient_indicator = 0.0;
      metric.particle_count = 0;
    }
  }

  const auto derefine_diag = manager.regrid(hierarchy);
  assert(derefine_diag.derefined_patch_count == 1);
  assert(hierarchy.levelView(1).empty());
}

void testConservativeProlongRestrict() {
  const cosmosim::amr::ConservedState coarse{
      .mass_code = 8.0,
      .momentum_x_code = 4.0,
      .momentum_y_code = -2.0,
      .momentum_z_code = 1.5,
      .total_energy_code = 16.0,
  };

  const auto fine = cosmosim::amr::ConservativeTransfer::prolongateFromCoarse(coarse, 8);
  assert(fine.size() == 8);

  cosmosim::amr::ConservedState total_fine;
  for (const auto& cell : fine) {
    total_fine += cell;
  }

  assert(std::abs(total_fine.mass_code - coarse.mass_code) < k_tolerance);
  assert(std::abs(total_fine.total_energy_code - coarse.total_energy_code) < k_tolerance);

  const auto restricted = cosmosim::amr::ConservativeTransfer::restrictToCoarse(fine);
  assert(std::abs(restricted.mass_code - coarse.mass_code) < k_tolerance);
  assert(std::abs(restricted.momentum_x_code - coarse.momentum_x_code) < k_tolerance);
  assert(std::abs(restricted.total_energy_code - coarse.total_energy_code) < k_tolerance);
}

void testRefluxCorrection() {
  cosmosim::amr::PatchHierarchy hierarchy;
  cosmosim::amr::PatchDescriptor root;
  root.cell_dims = {2, 2, 2};
  const std::uint64_t root_id = hierarchy.createRootPatch(root);

  auto* root_patch = hierarchy.findPatch(root_id);
  assert(root_patch != nullptr);

  auto conserved = root_patch->conservedView();
  conserved[0].mass_code = 10.0;
  conserved[0].total_energy_code = 20.0;

  const double coarse_volume = root_patch->cellVolumeComov();
  const double dt_code = 0.25;
  const double area = 0.5;

  cosmosim::amr::FluxRegisterEntry entry;
  entry.coarse_patch_id = root_id;
  entry.coarse_cell_index = 0;
  entry.coarse_face_flux_code.mass_code = 2.0;
  entry.fine_face_flux_code.mass_code = 3.0;
  entry.coarse_face_flux_code.total_energy_code = 4.0;
  entry.fine_face_flux_code.total_energy_code = 6.0;
  entry.face_area_comov = area;
  entry.dt_code = dt_code;

  const auto diag = cosmosim::amr::RefluxSynchronizer::apply(hierarchy, {&entry, 1});
  assert(diag.corrected_cells == 1);

  const double expected_mass = 10.0 - ((3.0 - 2.0) * area * dt_code / coarse_volume);
  const double expected_energy = 20.0 - ((6.0 - 4.0) * area * dt_code / coarse_volume);

  conserved = root_patch->conservedView();
  assert(std::abs(conserved[0].mass_code - expected_mass) < k_tolerance);
  assert(std::abs(conserved[0].total_energy_code - expected_energy) < k_tolerance);
}

}  // namespace

int main() {
  testRefineAndDerefineLifecycle();
  testConservativeProlongRestrict();
  testRefluxCorrection();
  return 0;
}
