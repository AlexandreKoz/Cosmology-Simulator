#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/amr/amr_framework.hpp"

int main() {
  constexpr std::size_t k_patch_cells_per_dim = 8;
  constexpr std::size_t k_leaf_patches = 64;
  constexpr std::size_t k_regrid_steps = 50;

  cosmosim::amr::PatchHierarchy hierarchy;
  cosmosim::amr::PatchDescriptor root;
  root.cell_dims = {
      static_cast<std::uint16_t>(k_patch_cells_per_dim),
      static_cast<std::uint16_t>(k_patch_cells_per_dim),
      static_cast<std::uint16_t>(k_patch_cells_per_dim)};

  const std::uint64_t root_id = hierarchy.createRootPatch(root);
  auto* root_patch = hierarchy.findPatch(root_id);
  if (root_patch == nullptr) {
    return 1;
  }

  auto root_metrics = root_patch->metricsView();
  for (std::size_t i = 0; i < root_metrics.size(); ++i) {
    root_metrics[i].density_code = (i % 17 == 0) ? 12.0 : 0.05;
    root_metrics[i].sound_speed_code = 1.0;
    root_metrics[i].gradient_indicator = (i % 13 == 0) ? 2.0 : 0.0;
    root_metrics[i].particle_count = (i % 11 == 0) ? 20 : 0;
  }

  cosmosim::amr::RefinementCriteria criteria;
  criteria.mass_threshold_code = 0.6;
  criteria.gradient_threshold = 1.0;
  criteria.particle_threshold = 8;
  criteria.jeans_resolution_cells = 0.5;

  cosmosim::amr::RefineDerefineManager manager(criteria);

  const auto setup_begin = std::chrono::steady_clock::now();
  const auto initial_diag = manager.regrid(hierarchy);
  const auto setup_end = std::chrono::steady_clock::now();

  std::vector<cosmosim::amr::FluxRegisterEntry> flux_entries;
  flux_entries.reserve(k_leaf_patches);
  for (std::size_t i = 0; i < k_leaf_patches; ++i) {
    flux_entries.push_back(cosmosim::amr::FluxRegisterEntry{
        .coarse_patch_id = root_id,
        .coarse_cell_index = i % root_patch->cellCount(),
        .coarse_face_flux_code = {.mass_code = 1.0, .total_energy_code = 1.5},
        .fine_face_flux_code = {.mass_code = 1.2, .total_energy_code = 1.8},
        .face_area_comov = 0.05,
        .dt_code = 0.01,
    });
  }

  std::uint64_t total_refined = initial_diag.refined_patch_count;
  std::uint64_t total_derefined = initial_diag.derefined_patch_count;
  std::uint64_t total_reflux_cells = 0;

  const auto run_begin = std::chrono::steady_clock::now();
  for (std::size_t step = 0; step < k_regrid_steps; ++step) {
    const auto diag = manager.regrid(hierarchy);
    total_refined += diag.refined_patch_count;
    total_derefined += diag.derefined_patch_count;

    const auto reflux_diag = cosmosim::amr::RefluxSynchronizer::apply(hierarchy, flux_entries);
    total_reflux_cells += reflux_diag.corrected_cells;
  }
  const auto run_end = std::chrono::steady_clock::now();

  const auto setup_us = std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_begin).count();
  const auto run_us = std::chrono::duration_cast<std::chrono::microseconds>(run_end - run_begin).count();

  std::cout << "build_type="
#ifdef NDEBUG
            << "release";
#else
            << "debug";
#endif
  std::cout << '\n';
  std::cout << "hardware_hint=single_node_cpu" << '\n';
  std::cout << "threads=1" << '\n';
  std::cout << "enabled_features=amr" << '\n';
  std::cout << "setup_us=" << setup_us << '\n';
  std::cout << "run_us=" << run_us << '\n';
  std::cout << "regrid_steps=" << k_regrid_steps << '\n';
  std::cout << "initial_refined=" << initial_diag.refined_patch_count << '\n';
  std::cout << "total_refined=" << total_refined << '\n';
  std::cout << "total_derefined=" << total_derefined << '\n';
  std::cout << "total_reflux_cells=" << total_reflux_cells << '\n';
  std::cout << "reflux_cells_per_us="
            << (run_us > 0 ? static_cast<double>(total_reflux_cells) / static_cast<double>(run_us) : 0.0)
            << '\n';

  return 0;
}
