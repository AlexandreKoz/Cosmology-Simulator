#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"

int main() {
  constexpr std::size_t k_particle_count = 500000;
  constexpr std::size_t k_active_stride = 4;

  cosmosim::core::SimulationState state;
  state.resizeParticles(k_particle_count);

  for (std::size_t i = 0; i < k_particle_count; ++i) {
    state.particle_sidecar.particle_id[i] = 100000 + i;
    state.particle_sidecar.sfc_key[i] = static_cast<std::uint64_t>(i * 17U);
    state.particle_sidecar.species_tag[i] = 0;
    state.particles.position_x_comoving[i] = static_cast<double>(i) * 0.01;
    state.particles.position_y_comoving[i] = static_cast<double>(i) * 0.02;
    state.particles.position_z_comoving[i] = static_cast<double>(i) * 0.03;
    state.particles.velocity_x_peculiar[i] = 0.1;
    state.particles.velocity_y_peculiar[i] = 0.2;
    state.particles.velocity_z_peculiar[i] = 0.3;
    state.particles.mass_code[i] = 1.0;
    state.particles.time_bin[i] = static_cast<std::uint8_t>(i % 4U);
  }

  std::vector<std::uint32_t> active_indices;
  active_indices.reserve(k_particle_count / k_active_stride);
  for (std::uint32_t i = 0; i < k_particle_count; i += k_active_stride) {
    active_indices.push_back(i);
  }

  volatile double full_checksum = 0.0;
  const auto full_begin = std::chrono::steady_clock::now();
  for (int repeat = 0; repeat < 30; ++repeat) {
    for (const auto index : active_indices) {
      full_checksum += state.particles.position_x_comoving[index] + state.particles.mass_code[index];
    }
  }
  const auto full_end = std::chrono::steady_clock::now();

  cosmosim::core::TransientStepWorkspace workspace;
  auto compact_view = cosmosim::core::buildGravityParticleKernelView(state, active_indices, workspace);

  volatile double compact_checksum = 0.0;
  const auto compact_begin = std::chrono::steady_clock::now();
  for (int repeat = 0; repeat < 30; ++repeat) {
    for (std::size_t i = 0; i < compact_view.size(); ++i) {
      compact_checksum += compact_view.position_x_comoving[i] + compact_view.mass_code[i];
    }
  }
  const auto compact_end = std::chrono::steady_clock::now();

  const auto full_us =
      std::chrono::duration_cast<std::chrono::microseconds>(full_end - full_begin).count();
  const auto compact_us =
      std::chrono::duration_cast<std::chrono::microseconds>(compact_end - compact_begin).count();

  const std::size_t full_bytes_per_particle = sizeof(double) * 8U + sizeof(std::uint8_t);
  const std::size_t compact_bytes_per_particle = sizeof(double) * 7U + sizeof(std::uint32_t);

  std::cout << "active_count=" << active_indices.size() << '\n';
  std::cout << "full_sweep_us=" << full_us << '\n';
  std::cout << "compact_sweep_us=" << compact_us << '\n';
  std::cout << "full_bytes_touched=" << (active_indices.size() * full_bytes_per_particle) << '\n';
  std::cout << "compact_bytes_touched=" << (active_indices.size() * compact_bytes_per_particle) << '\n';
  std::cout << "checksums=" << full_checksum << "," << compact_checksum << '\n';

  return 0;
}
