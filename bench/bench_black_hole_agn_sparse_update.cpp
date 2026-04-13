#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/black_hole_agn.hpp"

int main() {
  constexpr std::size_t k_cells = 1u << 20;
  constexpr std::size_t k_bh = 2048;
  constexpr std::size_t k_iterations = 32;

  cosmosim::core::SimulationState state;
  state.resizeCells(k_cells);
  for (std::size_t i = 0; i < k_cells; ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.center_y_comoving[i] = 0.0;
    state.cells.center_z_comoving[i] = 0.0;
    state.cells.time_bin[i] = 0;
    state.gas_cells.density_code[i] = 5.0;
    state.gas_cells.sound_speed_code[i] = 10.0;
    state.gas_cells.internal_energy_code[i] = 0.0;
  }

  cosmosim::physics::BlackHoleAgnConfig config;
  config.enabled = true;
  config.seed_halo_mass_threshold_code = 1.0;
  config.seed_mass_code = 5.0;
  cosmosim::physics::BlackHoleAgnModel model(config);

  std::vector<cosmosim::physics::BlackHoleSeedCandidate> seeds;
  seeds.reserve(k_bh);
  for (std::size_t i = 0; i < k_bh; ++i) {
    seeds.push_back(cosmosim::physics::BlackHoleSeedCandidate{static_cast<std::uint32_t>((i * 503) % k_cells), 10.0, 0});
  }

  const auto setup_start = std::chrono::steady_clock::now();
  auto report = model.apply(state, seeds, 1.0, 0);
  const auto setup_end = std::chrono::steady_clock::now();

  double accreted_mass = report.counters.integrated_accreted_mass_code;
  double feedback_energy = report.counters.integrated_feedback_energy_code;

  const auto steady_start = std::chrono::steady_clock::now();
  for (std::size_t iter = 0; iter < k_iterations; ++iter) {
    report = model.apply(state, std::vector<cosmosim::physics::BlackHoleSeedCandidate>{}, 0.2, iter + 1);
    accreted_mass += report.counters.integrated_accreted_mass_code;
    feedback_energy += report.counters.integrated_feedback_energy_code;
  }
  const auto steady_end = std::chrono::steady_clock::now();

  const double setup_ms = std::chrono::duration<double, std::milli>(setup_end - setup_start).count();
  const double steady_ms = std::chrono::duration<double, std::milli>(steady_end - steady_start).count();
  const double steady_s = std::max(steady_ms * 1.0e-3, 1.0e-12);
  const double bh_updates = static_cast<double>(k_bh * k_iterations);

  std::cout << "bench_black_hole_agn_sparse_update"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " hardware=cpu"
            << " threads=1"
            << " features=black_hole_agn+eddington_cap"
            << " setup_ms=" << setup_ms
            << " steady_ms=" << steady_ms
            << " bh_count=" << k_bh
            << " iterations=" << k_iterations
            << " bh_updates_per_s=" << (bh_updates / steady_s)
            << " effective_read_write_bandwidth_gb_s="
            << ((bh_updates * 12.0 * sizeof(double)) / steady_s * 1.0e-9)
            << " integrated_accreted_mass_code=" << accreted_mass
            << " integrated_feedback_energy_code=" << feedback_energy
            << '\n';

  return 0;
}
