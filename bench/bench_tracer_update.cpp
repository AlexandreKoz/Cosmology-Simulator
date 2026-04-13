#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/physics/tracer_support.hpp"

int main() {
  constexpr std::size_t k_cells = 1u << 18;
  constexpr std::size_t k_tracers = 1u << 18;
  constexpr std::size_t k_iterations = 10;

  cosmosim::core::SimulationState baseline_state;
  baseline_state.resizeCells(k_cells);
  for (std::size_t i = 0; i < k_cells; ++i) {
    baseline_state.cells.mass_code[i] = 1.0 + static_cast<double>(i % 23);
  }
  std::vector<std::uint32_t> active_cells(k_cells);
  for (std::size_t i = 0; i < k_cells; ++i) {
    active_cells[i] = static_cast<std::uint32_t>(i);
  }

  const auto baseline_start = std::chrono::steady_clock::now();
  double baseline_sink = 0.0;
  for (std::size_t it = 0; it < k_iterations; ++it) {
    for (std::size_t i = 0; i < k_cells; ++i) {
      baseline_sink += baseline_state.cells.mass_code[i] * 1.0e-12;
    }
  }
  const auto baseline_end = std::chrono::steady_clock::now();

  cosmosim::core::SimulationState tracer_state;
  tracer_state.resizeCells(k_cells);
  tracer_state.resizeParticles(k_tracers);
  tracer_state.tracers.resize(k_tracers);

  for (std::size_t i = 0; i < k_cells; ++i) {
    tracer_state.cells.mass_code[i] = 1.0 + static_cast<double>(i % 23);
  }
  for (std::size_t i = 0; i < k_tracers; ++i) {
    tracer_state.particle_sidecar.species_tag[i] =
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kTracer);
    tracer_state.tracers.particle_index[i] = static_cast<std::uint32_t>(i);
    tracer_state.tracers.host_cell_index[i] = static_cast<std::uint32_t>(i % k_cells);
    tracer_state.tracers.mass_fraction_of_host[i] = 0.01;
    tracer_state.tracers.last_host_mass_code[i] = tracer_state.cells.mass_code[i % k_cells];
    tracer_state.tracers.cumulative_exchanged_mass_code[i] = 0.0;
  }

  cosmosim::physics::TracerConfig tracer_config;
  tracer_config.enabled = true;
  tracer_config.track_mass = true;
  cosmosim::physics::TracerModel tracer_model(tracer_config);

  const auto setup_start = std::chrono::steady_clock::now();
  const auto warm = tracer_model.updateMassFromHostCells(tracer_state, active_cells);
  const auto setup_end = std::chrono::steady_clock::now();

  std::uint64_t updated_total = warm.updated_tracers;
  double delta_total = warm.cumulative_absolute_mass_delta_code;

  const auto steady_start = std::chrono::steady_clock::now();
  for (std::size_t it = 0; it < k_iterations; ++it) {
    const auto counters = tracer_model.updateMassFromHostCells(tracer_state, active_cells);
    updated_total += counters.updated_tracers;
    delta_total += counters.cumulative_absolute_mass_delta_code;
  }
  const auto steady_end = std::chrono::steady_clock::now();

  const double baseline_ms =
      std::chrono::duration<double, std::milli>(baseline_end - baseline_start).count();
  const double setup_ms =
      std::chrono::duration<double, std::milli>(setup_end - setup_start).count();
  const double steady_ms =
      std::chrono::duration<double, std::milli>(steady_end - steady_start).count();
  const double steady_s = std::max(steady_ms * 1.0e-3, 1.0e-12);
  const double updates = static_cast<double>(k_tracers * k_iterations);

  std::cout << "bench_tracer_update"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " hardware=cpu"
            << " threads=1"
            << " features=tracer_update"
            << " setup_ms=" << setup_ms
            << " steady_ms=" << steady_ms
            << " baseline_ms=" << baseline_ms
            << " tracers=" << k_tracers
            << " cells=" << k_cells
            << " iterations=" << k_iterations
            << " tracer_updates_per_s=" << (updates / steady_s)
            << " estimated_effective_bandwidth_gb_s="
            << ((updates * 6.0 * sizeof(double)) / steady_s * 1.0e-9)
            << " baseline_sink=" << baseline_sink
            << " updated_tracers_total=" << updated_total
            << " cumulative_absolute_mass_delta_code=" << delta_total
            << '\n';

  return 0;
}
