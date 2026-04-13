#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/star_formation.hpp"

int main() {
  constexpr std::size_t k_cells = 1u << 19;
  constexpr std::size_t k_iterations = 8;

  cosmosim::core::SimulationState state;
  state.resizeCells(k_cells);
  for (std::size_t i = 0; i < k_cells; ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.center_y_comoving[i] = 0.0;
    state.cells.center_z_comoving[i] = 0.0;
    state.cells.mass_code[i] = 2.5;
    state.gas_cells.density_code[i] = 25.0;
    state.gas_cells.temperature_code[i] = 5.0e3;
  }

  std::vector<std::uint32_t> active(k_cells);
  std::vector<double> div_v(k_cells, -0.2);
  std::vector<double> metallicity(k_cells, 0.01);
  for (std::size_t i = 0; i < active.size(); ++i) {
    active[i] = static_cast<std::uint32_t>(i);
  }

  cosmosim::physics::StarFormationConfig config;
  config.epsilon_ff = 0.01;
  config.stochastic_spawning = true;
  config.min_star_particle_mass_code = 0.05;
  cosmosim::physics::StarFormationModel model(config);

  const auto setup_start = std::chrono::steady_clock::now();
  const auto warm = model.apply(state, active, div_v, metallicity, 1.0e7, 0.95, 0, 0);
  const auto setup_end = std::chrono::steady_clock::now();

  std::uint64_t spawned_particles = warm.counters.spawned_particles;
  double spawned_mass = warm.counters.spawned_mass_code;

  const auto steady_start = std::chrono::steady_clock::now();
  for (std::size_t iter = 0; iter < k_iterations; ++iter) {
    const auto report = model.apply(state, active, div_v, metallicity, 1.0e7, 0.95, iter + 1, 0);
    spawned_particles += report.counters.spawned_particles;
    spawned_mass += report.counters.spawned_mass_code;
  }
  const auto steady_end = std::chrono::steady_clock::now();

  const double setup_ms = std::chrono::duration<double, std::milli>(setup_end - setup_start).count();
  const double steady_ms = std::chrono::duration<double, std::milli>(steady_end - steady_start).count();
  const double steady_s = std::max(steady_ms * 1.0e-3, 1.0e-12);
  const double cell_updates = static_cast<double>(k_cells * k_iterations);

  std::cout << "bench_star_formation_spawn"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " hardware=cpu"
            << " threads=1"
            << " features=star_formation+stochastic_spawning"
            << " setup_ms=" << setup_ms
            << " steady_ms=" << steady_ms
            << " cells=" << k_cells
            << " iterations=" << k_iterations
            << " cell_updates_per_s=" << (cell_updates / steady_s)
            << " effective_read_write_bandwidth_gb_s="
            << ((cell_updates * 6.0 * sizeof(double)) / steady_s * 1.0e-9)
            << " spawned_particles=" << spawned_particles
            << " spawned_mass_code=" << spawned_mass
            << '\n';

  return 0;
}
