#include <cassert>
#include <cmath>
#include <cstdint>

#include "cosmosim/amr/amr_framework.hpp"

namespace {

constexpr double k_tolerance = 1.0e-12;

void assertConservedClose(
    const cosmosim::amr::ConservedState& actual,
    const cosmosim::amr::ConservedState& expected) {
  assert(std::abs(actual.mass_code - expected.mass_code) < k_tolerance);
  assert(std::abs(actual.momentum_x_code - expected.momentum_x_code) < k_tolerance);
  assert(std::abs(actual.momentum_y_code - expected.momentum_y_code) < k_tolerance);
  assert(std::abs(actual.momentum_z_code - expected.momentum_z_code) < k_tolerance);
  assert(std::abs(actual.total_energy_code - expected.total_energy_code) < k_tolerance);
}

void seedRefinementState(cosmosim::amr::AmrPatch& patch) {
  auto conserved = patch.conservedView();
  auto metrics = patch.metricsView();
  for (std::size_t i = 0; i < conserved.size(); ++i) {
    const double row = static_cast<double>(i + 1);
    conserved[i].mass_code = 0.5 + row;
    conserved[i].momentum_x_code = 0.125 * row;
    conserved[i].momentum_y_code = -0.25 * row;
    conserved[i].momentum_z_code = 0.375 * row;
    conserved[i].total_energy_code = 4.0 + 2.0 * row;

    metrics[i].density_code = (i == 0) ? 32.0 : 0.5 + row;
    metrics[i].pressure_code = 0.25 + row;
    metrics[i].sound_speed_code = 1.0 + 0.1 * row;
    metrics[i].gradient_indicator = (i == 0) ? 3.0 : 0.0;
    metrics[i].particle_count = static_cast<std::uint32_t>((i == 0) ? 17U : (i % 5U));
  }
}

}  // namespace

int main() {
  cosmosim::amr::PatchHierarchy hierarchy;
  cosmosim::amr::PatchDescriptor root;
  root.cell_dims = {4, 3, 2};
  const std::uint64_t root_id = hierarchy.createRootPatch(root);

  auto* root_patch = hierarchy.findPatch(root_id);
  assert(root_patch != nullptr);
  seedRefinementState(*root_patch);

  const auto parent_total_before_refine = root_patch->totalConserved();

  cosmosim::amr::RefinementCriteria criteria;
  criteria.mass_threshold_code = 0.1;
  criteria.gradient_threshold = 1.0;
  criteria.particle_threshold = 8;
  criteria.jeans_resolution_cells = 0.5;

  cosmosim::amr::RefineDerefineManager manager(criteria);
  const auto regrid_diag = manager.regrid(hierarchy);
  assert(regrid_diag.refined_patch_count == 1);
  assert(hierarchy.levelView(1).size() == 8);

  cosmosim::amr::ConservedState child_total_after_refine;
  for (const auto& child_patch : hierarchy.levelView(1)) {
    assert(child_patch.descriptor().parent_patch_id == root_id);
    for (const auto& cell : child_patch.conservedView()) {
      child_total_after_refine += cell;
      assert(cell.mass_code > 0.0);
      assert(cell.total_energy_code > 0.0);
    }
  }

  assertConservedClose(child_total_after_refine, parent_total_before_refine);

  return 0;
}
