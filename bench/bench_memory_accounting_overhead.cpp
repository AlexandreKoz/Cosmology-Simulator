#include <chrono>
#include <cstdint>
#include <iostream>

#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/simulation_state.hpp"

int main() {
  constexpr std::size_t k_particle_count = 200000;
  constexpr std::size_t k_cell_count = 100000;
  constexpr int k_iterations = 200;

  cosmosim::core::SimulationState state;
  state.resizeParticles(k_particle_count);
  state.resizeCells(k_cell_count);
  for (std::size_t i = 0; i < k_particle_count; ++i) {
    state.particle_sidecar.particle_id[i] = 1000000 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[i] = 0;
    state.particles.mass_code[i] = 1.0;
  }
  for (std::size_t i = 0; i < k_cell_count; ++i) {
    state.cells.mass_code[i] = 1.0;
    state.gas_cells.density_code[i] = 1.0;
  }

  cosmosim::core::TransientStepWorkspace workspace;
  workspace.gravity_particle_index.reserve(k_particle_count / 8);
  workspace.particle_position_x_comoving.reserve(k_particle_count / 8);
  workspace.hydro_cell_index.reserve(k_cell_count / 8);
  workspace.hydro_recon_gradient_x.reserve(k_cell_count / 8);
  static_cast<void>(workspace.scratch.allocateBytes(1U << 20U, alignof(double)));

  volatile std::uint64_t checksum = 0;
  const auto begin = std::chrono::steady_clock::now();
  for (int iter = 0; iter < k_iterations; ++iter) {
    const auto report = cosmosim::core::collectSimulationMemoryReport(state, &workspace);
    checksum += report.totals.persistent_total_bytes;
    checksum += report.totals.transient_total_bytes;
    checksum += report.entries.size();
  }
  const auto end = std::chrono::steady_clock::now();
  const double elapsed_ms = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - begin).count();

  std::cout << "benchmark=memory_accounting_overhead"
            << " particles=" << k_particle_count
            << " cells=" << k_cell_count
            << " iterations=" << k_iterations
            << " elapsed_ms=" << elapsed_ms
            << " report_us=" << (elapsed_ms * 1000.0 / static_cast<double>(k_iterations))
            << " checksum=" << checksum
            << '\n';
  return 0;
}
