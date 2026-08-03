#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"

int main() {
  constexpr std::size_t kStarCount = 1U << 17;
  constexpr std::size_t kIterations = 8;
  cosmosim::core::SimulationState state;
  state.resizeParticles(kStarCount);
  state.star_particles.resize(kStarCount);
  std::vector<std::uint32_t> active(kStarCount);
  for (std::size_t i = 0; i < kStarCount; ++i) {
    state.particles.mass_code[i] = 1.0;
    state.star_particles.particle_index[i] = static_cast<std::uint32_t>(i);
    state.star_particles.birth_mass_code[i] = 1.0;
    state.star_particles.metallicity_mass_fraction[i] =
        0.03 * static_cast<double>(i % 1024U) / 1023.0;
    active[i] = static_cast<std::uint32_t>(i);
  }
  const auto table = cosmosim::physics::StellarEvolutionTable::loadFromTextFile(
      std::string(COSMOSIM_SOURCE_DIR) +
      "/resources/stellar_evolution/test_synthetic_v2.txt");
  cosmosim::physics::StellarEvolutionBookkeeper bookkeeper({}, table);
  const auto start = std::chrono::steady_clock::now();
  double returned = 0.0;
  for (std::size_t iteration = 0; iteration < kIterations; ++iteration) {
    const auto report = bookkeeper.evaluateElapsedYears(state, active, 1.0e7);
    returned += report.counters.returned_mass_code;
    bookkeeper.commitBudgets(state, report);
  }
  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  const double updates = static_cast<double>(kStarCount * kIterations);
  std::cout << "bench_metal_enrichment"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " stars=" << kStarCount
            << " iterations=" << kIterations
            << " seconds=" << seconds
            << " star_intervals_per_s=" << updates / std::max(seconds, 1.0e-12)
            << " new_persistent_bytes_per_star=" << 8U * sizeof(double)
            << " returned_mass_code=" << returned
            << '\n';
}
