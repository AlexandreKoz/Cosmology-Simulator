#include <cassert>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_set>
#include <vector>

#include "cosmosim/amr/amr_framework.hpp"

namespace {

constexpr double k_tolerance = 1.0e-12;

[[nodiscard]] std::size_t cellIndex(
    std::array<std::uint16_t, 3> dims,
    std::uint16_t x,
    std::uint16_t y,
    std::uint16_t z) {
  return (static_cast<std::size_t>(z) * dims[1] + y) * dims[0] + x;
}

[[nodiscard]] std::size_t parentCellIndexForChildCell(
    const cosmosim::amr::PatchDescriptor& child,
    std::uint16_t x,
    std::uint16_t y,
    std::uint16_t z) {
  const std::uint8_t octant = static_cast<std::uint8_t>(child.morton_key & 7U);
  const std::uint16_t fine_x =
      static_cast<std::uint16_t>(x + (((octant & 1U) != 0U) ? child.cell_dims[0] : 0U));
  const std::uint16_t fine_y =
      static_cast<std::uint16_t>(y + (((octant & 2U) != 0U) ? child.cell_dims[1] : 0U));
  const std::uint16_t fine_z =
      static_cast<std::uint16_t>(z + (((octant & 4U) != 0U) ? child.cell_dims[2] : 0U));

  return cellIndex(
      child.cell_dims,
      static_cast<std::uint16_t>(fine_x / 2U),
      static_cast<std::uint16_t>(fine_y / 2U),
      static_cast<std::uint16_t>(fine_z / 2U));
}

void assertConservedClose(
    const cosmosim::amr::ConservedState& actual,
    const cosmosim::amr::ConservedState& expected) {
  assert(std::abs(actual.mass_code - expected.mass_code) < k_tolerance);
  assert(std::abs(actual.momentum_x_code - expected.momentum_x_code) < k_tolerance);
  assert(std::abs(actual.momentum_y_code - expected.momentum_y_code) < k_tolerance);
  assert(std::abs(actual.momentum_z_code - expected.momentum_z_code) < k_tolerance);
  assert(std::abs(actual.total_energy_code - expected.total_energy_code) < k_tolerance);
}

}  // namespace

int main() {
  cosmosim::amr::PatchHierarchy hierarchy;
  cosmosim::amr::PatchDescriptor root;
  root.cell_dims = {3, 2, 2};
  const std::uint64_t root_id = hierarchy.createRootPatch(root);

  auto* root_patch = hierarchy.findPatch(root_id);
  assert(root_patch != nullptr);
  const std::vector<std::uint64_t> parent_gas_cell_ids_before(
      root_patch->gasCellIdView().begin(),
      root_patch->gasCellIdView().end());

  const auto child_patch_ids = hierarchy.refinePatch(root_id);
  assert(child_patch_ids.size() == 8);

  std::vector<cosmosim::amr::ConservedState> expected_parent(root_patch->cellCount());
  std::vector<std::uint32_t> expected_particle_count(root_patch->cellCount(), 0U);
  std::unordered_set<std::uint64_t> child_gas_cell_ids;

  for (std::size_t child_ordinal = 0; child_ordinal < child_patch_ids.size(); ++child_ordinal) {
    auto* child_patch = hierarchy.findPatch(child_patch_ids[child_ordinal]);
    assert(child_patch != nullptr);

    const auto& descriptor = child_patch->descriptor();
    auto child_conserved = child_patch->conservedView();
    auto child_metrics = child_patch->metricsView();
    const auto child_gas_ids = child_patch->gasCellIdView();

    for (std::uint16_t z = 0; z < descriptor.cell_dims[2]; ++z) {
      for (std::uint16_t y = 0; y < descriptor.cell_dims[1]; ++y) {
        for (std::uint16_t x = 0; x < descriptor.cell_dims[0]; ++x) {
          const std::size_t child_index = cellIndex(descriptor.cell_dims, x, y, z);
          const std::size_t parent_index = parentCellIndexForChildCell(descriptor, x, y, z);
          const double stamp = static_cast<double>(
              1000U * static_cast<unsigned>(child_ordinal + 1U) +
              static_cast<unsigned>(child_index + 1U));

          child_conserved[child_index] = {
              .mass_code = 0.25 + stamp,
              .momentum_x_code = 0.5 + 0.125 * stamp,
              .momentum_y_code = -0.75 - 0.0625 * stamp,
              .momentum_z_code = 1.25 + 0.03125 * stamp,
              .total_energy_code = 2.0 + 0.5 * stamp,
          };
          expected_parent[parent_index] += child_conserved[child_index];

          child_metrics[child_index].density_code = 100.0 + stamp;
          child_metrics[child_index].pressure_code = 10.0 + 0.01 * stamp;
          child_metrics[child_index].sound_speed_code = 1.0 + 0.001 * stamp;
          child_metrics[child_index].gradient_indicator = 0.001 * stamp;
          child_metrics[child_index].particle_count =
              static_cast<std::uint32_t>((child_index % 5U) + 1U);
          expected_particle_count[parent_index] += child_metrics[child_index].particle_count;

          assert(child_gas_cell_ids.insert(child_gas_ids[child_index]).second);
        }
      }
    }
  }

  assert(hierarchy.derefinePatch(root_id));
  root_patch = hierarchy.findPatch(root_id);
  assert(root_patch != nullptr);
  assert(root_patch->isLeaf());
  assert(hierarchy.levelView(1).empty());

  const auto parent_conserved_after = root_patch->conservedView();
  const auto parent_metrics_after = root_patch->metricsView();
  for (std::size_t i = 0; i < parent_conserved_after.size(); ++i) {
    assertConservedClose(parent_conserved_after[i], expected_parent[i]);
    assert(parent_metrics_after[i].particle_count == expected_particle_count[i]);
  }

  const std::vector<std::uint64_t> parent_gas_cell_ids_after(
      root_patch->gasCellIdView().begin(),
      root_patch->gasCellIdView().end());
  assert(parent_gas_cell_ids_after == parent_gas_cell_ids_before);

  const auto retired_gas_cell_ids = hierarchy.retiredGasCellIds();
  assert(retired_gas_cell_ids.size() == child_gas_cell_ids.size());
  for (const auto retired_id : retired_gas_cell_ids) {
    assert(child_gas_cell_ids.contains(retired_id));
    assert(std::find(
               parent_gas_cell_ids_after.begin(),
               parent_gas_cell_ids_after.end(),
               retired_id) == parent_gas_cell_ids_after.end());
  }

  return 0;
}
