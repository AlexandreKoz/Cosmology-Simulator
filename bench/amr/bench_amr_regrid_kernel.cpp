#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "bench/reporting/bench_report.hpp"
#include "cosmosim/amr/amr_framework.hpp"

namespace {

void seedPatchMetrics(cosmosim::amr::AmrPatch& patch) {
  auto metrics = patch.metricsView();
  for (std::size_t i = 0; i < metrics.size(); ++i) {
    metrics[i].density_code = (i % 11U == 0) ? 8.0 : 0.1;
    metrics[i].gradient_indicator = (i % 7U == 0) ? 2.5 : 0.1;
    metrics[i].sound_speed_code = 1.0;
    metrics[i].particle_count = (i % 13U == 0) ? 12U : 0U;
  }
}

}  // namespace

int main() {
  constexpr std::size_t k_cells_per_dim = 8;
  constexpr std::size_t k_steps = 64;

  const auto execution = cosmosim::bench::defaultExecutionConfig(2, k_steps);

  cosmosim::amr::PatchHierarchy hierarchy;
  cosmosim::amr::PatchDescriptor root;
  root.cell_dims = {
      static_cast<std::uint16_t>(k_cells_per_dim),
      static_cast<std::uint16_t>(k_cells_per_dim),
      static_cast<std::uint16_t>(k_cells_per_dim)};

  const std::uint64_t root_id = hierarchy.createRootPatch(root);
  auto* root_patch = hierarchy.findPatch(root_id);
  if (root_patch == nullptr) {
    return 1;
  }
  seedPatchMetrics(*root_patch);

  cosmosim::amr::RefinementCriteria criteria;
  criteria.mass_threshold_code = 0.55;
  criteria.gradient_threshold = 1.0;
  criteria.particle_threshold = 8;
  criteria.jeans_resolution_cells = 0.5;

  cosmosim::amr::RefineDerefineManager manager(criteria);

  std::vector<cosmosim::amr::FluxRegisterEntry> reflux_entries;
  reflux_entries.reserve(64);
  for (std::size_t i = 0; i < 64; ++i) {
    reflux_entries.push_back(cosmosim::amr::FluxRegisterEntry{
        .coarse_patch_id = root_id,
        .coarse_cell_index = i % root_patch->cellCount(),
        .coarse_face_flux_code = {.mass_code = 1.0, .total_energy_code = 0.8},
        .fine_face_flux_code = {.mass_code = 1.1, .total_energy_code = 0.9},
        .face_area_comov = 0.1,
        .dt_code = 0.02,
    });
  }

  for (std::size_t iter = 0; iter < execution.warmup_iterations; ++iter) {
    const auto diag = manager.regrid(hierarchy);
    (void)diag;
    const auto reflux_diag = cosmosim::amr::RefluxSynchronizer::apply(hierarchy, reflux_entries);
    (void)reflux_diag;
  }

  std::size_t total_refined = 0;
  std::size_t total_derefined = 0;
  std::size_t total_reflux_cells = 0;

  const auto begin = cosmosim::bench::BenchmarkClock::now();
  for (std::size_t iter = 0; iter < execution.measurement_iterations; ++iter) {
    const auto diag = manager.regrid(hierarchy);
    total_refined += diag.refined_patch_count;
    total_derefined += diag.derefined_patch_count;

    const auto reflux_diag = cosmosim::amr::RefluxSynchronizer::apply(hierarchy, reflux_entries);
    total_reflux_cells += reflux_diag.corrected_cells;
  }
  const auto end = cosmosim::bench::BenchmarkClock::now();

  const double measurement_ms = cosmosim::bench::BenchmarkClock::millisecondsBetween(begin, end);
  const double reflux_cells_per_second =
      measurement_ms > 0.0 ? static_cast<double>(total_reflux_cells) / (measurement_ms * 1.0e-3) : 0.0;

  const std::uint64_t bytes_moved_proxy = static_cast<std::uint64_t>(total_reflux_cells) * sizeof(cosmosim::amr::ConservedState);

  cosmosim::bench::BenchmarkReporter reporter("bench_amr_regrid_reflux_kernel");
  cosmosim::bench::addExecutionFields(reporter, execution);
  reporter.addField("root_cells_per_dim", k_cells_per_dim);
  reporter.addField("measurement_ms", measurement_ms);
  reporter.addField("refined_patch_count", total_refined);
  reporter.addField("derefined_patch_count", total_derefined);
  reporter.addField("reflux_corrected_cells", total_reflux_cells);
  reporter.addField("reflux_cells_per_second", reflux_cells_per_second);
  reporter.addField("patch_count_final", hierarchy.patchCount());
  cosmosim::bench::addBandwidthFields(reporter, bytes_moved_proxy, measurement_ms, "bytes_moved_proxy", "effective_bandwidth_proxy_gb_s");
  reporter.write();

  return 0;
}
