#include <chrono>
#include <cstdint>
#include <iostream>

#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "core_benchmark_process_memory.hpp"

namespace {

std::uint64_t knownOwnedBytes(const cosmosim::core::MemoryReport& report) {
  std::uint64_t total = 0;
  for (const auto& entry : report.entries) {
    if (!entry.estimate_only && entry.lifetime != cosmosim::core::MemoryLifetime::kUnknown) {
      total += entry.owned_capacity_bytes;
    }
  }
  return total;
}

std::uint64_t estimatedBytes(const cosmosim::core::MemoryReport& report) {
  std::uint64_t total = 0;
  for (const auto& entry : report.entries) {
    if (entry.estimate_only) {
      total += entry.owned_capacity_bytes;
    }
  }
  return total;
}

std::size_t unknownExternalEntryCount(const cosmosim::core::MemoryReport& report) {
  std::size_t count = 0;
  for (const auto& entry : report.entries) {
    if (entry.lifetime == cosmosim::core::MemoryLifetime::kUnknown) {
      ++count;
    }
  }
  return count;
}

}  // namespace

int main() {
  constexpr std::size_t k_particle_count = 200000;
  constexpr std::size_t k_cell_count = 100000;
  constexpr int k_iterations = 200;

  const auto baseline_memory = cosmosim::bench::sampleProcessMemory();

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

  const auto allocated_memory = cosmosim::bench::sampleProcessMemory();
  const auto calibration_report = cosmosim::core::collectSimulationMemoryReport(state, &workspace);
  const std::uint64_t known_owned_bytes = knownOwnedBytes(calibration_report);
  const std::uint64_t estimated_bytes = estimatedBytes(calibration_report);
  const std::size_t unknown_entries = unknownExternalEntryCount(calibration_report);

  volatile std::uint64_t checksum = 0;
  const auto begin = std::chrono::steady_clock::now();
  for (int iter = 0; iter < k_iterations; ++iter) {
    const auto report = cosmosim::core::collectSimulationMemoryReport(state, &workspace);
    checksum += report.totals.persistent_total_bytes;
    checksum += report.totals.transient_total_bytes;
    checksum += report.entries.size();
  }
  const auto end = std::chrono::steady_clock::now();
  const auto final_memory = cosmosim::bench::sampleProcessMemory();
  const double elapsed_ms = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - begin).count();

  std::cout << "benchmark=memory_accounting_overhead"
            << " particles=" << k_particle_count
            << " cells=" << k_cell_count
            << " iterations=" << k_iterations
            << " elapsed_ms=" << elapsed_ms
            << " report_us=" << (elapsed_ms * 1000.0 / static_cast<double>(k_iterations))
            << " known_owned_bytes=" << known_owned_bytes
            << " estimated_bytes=" << estimated_bytes
            << " unknown_external_entries=" << unknown_entries;
  if (baseline_memory.current_rss_bytes && allocated_memory.current_rss_bytes) {
    const std::uint64_t baseline = *baseline_memory.current_rss_bytes;
    const std::uint64_t allocated = *allocated_memory.current_rss_bytes;
    const std::uint64_t rss_delta = allocated >= baseline ? allocated - baseline : 0U;
    std::cout << " baseline_rss_bytes=" << baseline
              << " allocated_rss_bytes=" << allocated
              << " rss_delta_bytes=" << rss_delta
              << " known_owned_to_rss_delta_ratio="
              << (rss_delta > 0U ? static_cast<double>(known_owned_bytes) / static_cast<double>(rss_delta) : 0.0);
  } else {
    std::cout << " rss_current=unavailable";
  }
  if (final_memory.peak_rss_bytes) {
    std::cout << " peak_rss_bytes=" << *final_memory.peak_rss_bytes;
  } else {
    std::cout << " peak_rss=unavailable";
  }
  std::cout << " checksum=" << checksum << '\n';
  return 0;
}
