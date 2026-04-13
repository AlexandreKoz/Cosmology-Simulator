#include <chrono>
#include <cstddef>
#include <iostream>

#include "cosmosim/core/simulation_state.hpp"

int main() {
  constexpr std::size_t particle_count = 3'000'000;
  constexpr std::size_t cell_count = 1'200'000;
  constexpr std::size_t patch_count = 40'000;

  cosmosim::core::SimulationState state;

  const auto start = std::chrono::steady_clock::now();
  state.resizeParticles(particle_count);
  state.resizeCells(cell_count);
  state.resizePatches(patch_count);
  state.species.count_by_species = {particle_count, 0, 0, 0, 0};
  const auto stop = std::chrono::steady_clock::now();

  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

  std::cout << "bench_simulation_state_resize"
            << " particles=" << state.particles.size()
            << " cells=" << state.cells.size()
            << " patches=" << state.patches.size()
            << " elapsed_ms=" << elapsed_ms << '\n';

  return state.validateOwnershipInvariants() ? 0 : 2;
}
