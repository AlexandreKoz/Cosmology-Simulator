#include <cassert>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_set>
#include <vector>

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

  for (auto& cell : root_patch->conservedView()) {
    cell.mass_code = 1.0;
    cell.momentum_x_code = 0.5;
    cell.momentum_y_code = -0.25;
    cell.momentum_z_code = 0.125;
    cell.total_energy_code = 2.0;
  }
  const auto parent_total_before_refine = root_patch->totalConserved();

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

  cosmosim::amr::ConservedState child_total_after_refine;
  for (const auto& child_patch : hierarchy.levelView(1)) {
    child_total_after_refine += child_patch.totalConserved();
    for (const auto& cell : child_patch.conservedView()) {
      assert(cell.mass_code > 0.0);
      assert(cell.total_energy_code > 0.0);
    }
  }
  assertConservedClose(child_total_after_refine, parent_total_before_refine);

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

  assertConservedClose(total_fine, coarse);

  const auto restricted = cosmosim::amr::ConservativeTransfer::restrictToCoarse(fine);
  assertConservedClose(restricted, coarse);
}

void testRefinePatchInitializesChildrenConservatively() {
  cosmosim::amr::PatchHierarchy hierarchy;
  cosmosim::amr::PatchDescriptor root;
  root.cell_dims = {3, 2, 2};
  const std::uint64_t root_id = hierarchy.createRootPatch(root);

  auto* root_patch = hierarchy.findPatch(root_id);
  assert(root_patch != nullptr);

  auto root_conserved = root_patch->conservedView();
  auto root_metrics = root_patch->metricsView();
  for (std::size_t i = 0; i < root_conserved.size(); ++i) {
    const double row = static_cast<double>(i + 1);
    root_conserved[i].mass_code = row;
    root_conserved[i].momentum_x_code = 0.25 * row;
    root_conserved[i].momentum_y_code = -0.5 * row;
    root_conserved[i].momentum_z_code = 0.75 * row;
    root_conserved[i].total_energy_code = 2.0 * row;

    root_metrics[i].density_code = 3.0 + row;
    root_metrics[i].pressure_code = 1.5 + row;
    root_metrics[i].sound_speed_code = 0.25 + row;
    root_metrics[i].gradient_indicator = 0.1 * row;
    root_metrics[i].particle_count = static_cast<std::uint32_t>(9U + i);
  }

  const auto parent_total_before_refine = root_patch->totalConserved();
  const std::uint32_t parent_particle_count_before_refine = [&root_metrics]() {
    std::uint32_t total = 0;
    for (const auto& metrics : root_metrics) {
      total += metrics.particle_count;
    }
    return total;
  }();

  const auto child_ids = hierarchy.refinePatch(root_id);
  assert(hierarchy.levelView(1).size() == child_ids.size());

  cosmosim::amr::ConservedState child_total_after_refine;
  std::uint32_t child_particle_count_after_refine = 0;
  for (const auto child_id : child_ids) {
    const auto* child_patch = hierarchy.findPatch(child_id);
    assert(child_patch != nullptr);
    assert(child_patch->descriptor().parent_patch_id == root_id);
    for (const auto& cell : child_patch->conservedView()) {
      child_total_after_refine += cell;
      assert(cell.mass_code > 0.0);
      assert(cell.total_energy_code > 0.0);
    }
    for (const auto& metrics : child_patch->metricsView()) {
      child_particle_count_after_refine += metrics.particle_count;
      assert(metrics.density_code > 0.0);
      assert(metrics.pressure_code > 0.0);
      assert(metrics.sound_speed_code > 0.0);
    }
  }

  assertConservedClose(child_total_after_refine, parent_total_before_refine);
  assert(child_particle_count_after_refine == parent_particle_count_before_refine);
}

void testDerefineRestrictsModifiedChildrenAndArchivesIds() {
  cosmosim::amr::PatchHierarchy hierarchy;
  cosmosim::amr::PatchDescriptor root;
  root.cell_dims = {2, 2, 2};
  const std::uint64_t root_id = hierarchy.createRootPatch(root);

  auto* root_patch = hierarchy.findPatch(root_id);
  assert(root_patch != nullptr);

  const std::vector<std::uint64_t> parent_gas_cell_ids(
      root_patch->gasCellIdView().begin(),
      root_patch->gasCellIdView().end());
  assert(parent_gas_cell_ids.size() == root_patch->cellCount());

  const auto child_ids = hierarchy.refinePatch(root_id);
  assert(hierarchy.levelView(1).size() == child_ids.size());

  std::unordered_set<std::uint64_t> child_gas_cell_ids;
  cosmosim::amr::ConservedState expected_restricted_total;
  std::uint32_t expected_particle_count = 0;
  for (std::size_t child_patch_index = 0; child_patch_index < child_ids.size(); ++child_patch_index) {
    auto* child_patch = hierarchy.findPatch(child_ids[child_patch_index]);
    assert(child_patch != nullptr);

    auto child_conserved = child_patch->conservedView();
    auto child_metrics = child_patch->metricsView();
    auto child_ids_view = child_patch->gasCellIdView();
    for (std::size_t cell_index = 0; cell_index < child_conserved.size(); ++cell_index) {
      const double row =
          static_cast<double>((child_patch_index + 1U) * 100U + static_cast<unsigned>(cell_index + 1U));
      child_conserved[cell_index].mass_code = 1.0 + row;
      child_conserved[cell_index].momentum_x_code = 0.5 * row;
      child_conserved[cell_index].momentum_y_code = -0.25 * row;
      child_conserved[cell_index].momentum_z_code = 0.125 * row;
      child_conserved[cell_index].total_energy_code = 3.0 + 0.75 * row;
      expected_restricted_total += child_conserved[cell_index];

      child_metrics[cell_index].particle_count = static_cast<std::uint32_t>((cell_index % 3U) + 1U);
      expected_particle_count += child_metrics[cell_index].particle_count;
      assert(child_gas_cell_ids.insert(child_ids_view[cell_index]).second);
      assert(std::find(
                 parent_gas_cell_ids.begin(),
                 parent_gas_cell_ids.end(),
                 child_ids_view[cell_index]) == parent_gas_cell_ids.end());
    }
  }

  assert(hierarchy.derefinePatch(root_id));
  root_patch = hierarchy.findPatch(root_id);
  assert(root_patch != nullptr);
  assert(root_patch->isLeaf());
  assert(hierarchy.levelView(1).empty());

  assertConservedClose(root_patch->totalConserved(), expected_restricted_total);
  assert(std::vector<std::uint64_t>(
             root_patch->gasCellIdView().begin(),
             root_patch->gasCellIdView().end()) == parent_gas_cell_ids);

  std::uint32_t parent_particle_count = 0;
  for (const auto& metrics : root_patch->metricsView()) {
    parent_particle_count += metrics.particle_count;
    assert(metrics.density_code > 0.0);
  }
  assert(parent_particle_count == expected_particle_count);

  const auto retired_ids = hierarchy.retiredGasCellIds();
  assert(retired_ids.size() == child_gas_cell_ids.size());
  for (const auto retired_id : retired_ids) {
    assert(child_gas_cell_ids.contains(retired_id));
  }
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
  testRefinePatchInitializesChildrenConservatively();
  testDerefineRestrictsModifiedChildrenAndArchivesIds();
  testRefluxCorrection();
  return 0;
}
