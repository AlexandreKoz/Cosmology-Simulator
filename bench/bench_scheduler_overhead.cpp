#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/core/time_integration.hpp"

int main() {
  constexpr std::size_t k_particle_count = 200000;
  constexpr std::size_t k_cell_count = 80000;
  constexpr int k_step_count = 20000;

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 1.0e-4;

  std::vector<std::uint32_t> active_particles;
  std::vector<std::uint32_t> active_cells;
  active_particles.reserve(k_particle_count / 8);
  active_cells.reserve(k_cell_count / 8);

  for (std::uint32_t i = 0; i < k_particle_count; i += 8) {
    active_particles.push_back(i);
  }
  for (std::uint32_t i = 0; i < k_cell_count; i += 8) {
    active_cells.push_back(i);
  }

  cosmosim::core::ActiveSetDescriptor active_set{
      .particle_indices = active_particles,
      .cell_indices = active_cells,
      .particles_are_subset = true,
      .cells_are_subset = true,
  };

  cosmosim::core::StageScheduler scheduler;

  std::uint64_t checksum = 0;
  const auto start = std::chrono::steady_clock::now();
  for (int step = 0; step < k_step_count; ++step) {
    integrator_state.step_index = static_cast<std::uint64_t>(step);
    const auto stages = scheduler.schedule(integrator_state, active_set);
    checksum += stages.size();
    checksum += static_cast<std::uint64_t>(active_set.hasParticleSubset(k_particle_count));
    checksum += static_cast<std::uint64_t>(active_set.hasCellSubset(k_cell_count));
  }
  const auto end = std::chrono::steady_clock::now();

  const auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

  std::cout << "steps=" << k_step_count << '\n';
  std::cout << "active_particles=" << active_particles.size() << '\n';
  std::cout << "active_cells=" << active_cells.size() << '\n';
  std::cout << "elapsed_us=" << elapsed_us << '\n';
  std::cout << "us_per_step=" << (static_cast<double>(elapsed_us) / static_cast<double>(k_step_count))
            << '\n';
  std::cout << "checksum=" << checksum << '\n';

  return 0;
}
