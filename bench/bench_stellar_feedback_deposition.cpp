#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"

int main() {
  constexpr std::size_t k_star_count = 1u << 10;
  constexpr std::size_t k_cell_count = 1u << 11;
  constexpr std::size_t k_iterations = 3;

  cosmosim::core::SimulationState state;
  state.resizeCells(k_cell_count);
  for (std::size_t cell_index = 0; cell_index < k_cell_count; ++cell_index) {
    state.cells.center_x_comoving[cell_index] = static_cast<double>(cell_index & 127u);
    state.cells.center_y_comoving[cell_index] = static_cast<double>((cell_index >> 7u) & 127u);
    state.cells.center_z_comoving[cell_index] = static_cast<double>((cell_index >> 14u) & 7u);
    state.cells.mass_code[cell_index] = 1.0;
    state.gas_cells.density_code[cell_index] = 1.0;
    state.gas_cells.internal_energy_code[cell_index] = 1.0;
  }

  state.resizeParticles(k_star_count);
  state.star_particles.resize(k_star_count);
  for (std::size_t star_index = 0; star_index < k_star_count; ++star_index) {
    state.particles.position_x_comoving[star_index] = static_cast<double>(star_index & 127u) + 0.5;
    state.particles.position_y_comoving[star_index] = static_cast<double>((star_index >> 7u) & 127u) + 0.5;
    state.particles.position_z_comoving[star_index] = 0.5;
    state.star_particles.particle_index[star_index] = static_cast<std::uint32_t>(star_index);
    state.star_particles.birth_mass_code[star_index] = 1.0;
  }

  std::vector<std::uint32_t> active_stars(k_star_count);
  std::vector<double> returned_mass_delta(k_star_count, 0.02);
  std::vector<double> returned_metals_delta(k_star_count, 0.002);
  for (std::size_t i = 0; i < k_star_count; ++i) {
    active_stars[i] = static_cast<std::uint32_t>(i);
  }

  cosmosim::physics::StellarFeedbackConfig config;
  config.mode = cosmosim::physics::StellarFeedbackMode::kThermalKineticMomentum;
  config.variant = cosmosim::physics::StellarFeedbackVariant::kNone;
  config.neighbor_count = 16;
  cosmosim::physics::StellarFeedbackModel model(config);
  cosmosim::physics::StellarFeedbackModuleState module_state;

  const auto setup_begin = std::chrono::steady_clock::now();
  auto warm = model.apply(state, module_state, active_stars, returned_mass_delta, returned_metals_delta, 1.0e-3);
  const auto setup_end = std::chrono::steady_clock::now();

  double deposited_mass_total = warm.counters.deposited_mass_code;
  const auto steady_begin = std::chrono::steady_clock::now();
  for (std::size_t i = 0; i < k_iterations; ++i) {
    const auto report = model.apply(state, module_state, active_stars, returned_mass_delta, returned_metals_delta, 1.0e-3);
    deposited_mass_total += report.counters.deposited_mass_code;
  }
  const auto steady_end = std::chrono::steady_clock::now();

  const double setup_ms = std::chrono::duration<double, std::milli>(setup_end - setup_begin).count();
  const double steady_ms = std::chrono::duration<double, std::milli>(steady_end - steady_begin).count();
  const double steady_s = std::max(steady_ms * 1.0e-3, 1.0e-12);
  const double star_updates = static_cast<double>(k_star_count * k_iterations);

  std::cout << "bench_stellar_feedback_deposition"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " hardware=cpu"
            << " threads=1"
            << " features=feedback_deposition"
            << " setup_ms=" << setup_ms
            << " steady_ms=" << steady_ms
            << " stars=" << k_star_count
            << " cells=" << k_cell_count
            << " iterations=" << k_iterations
            << " star_updates_per_s=" << (star_updates / steady_s)
            << " approx_effective_bandwidth_gb_s="
            << ((star_updates * 16.0 * sizeof(double)) / steady_s * 1.0e-9)
            << " deposited_mass_code=" << deposited_mass_total << '\n';

  return 0;
}
