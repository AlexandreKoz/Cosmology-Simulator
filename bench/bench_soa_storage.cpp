#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
#include <utility>

#include "cosmosim/core/soa_storage.hpp"

int main() {
  constexpr std::size_t particle_count = 400000;
  cosmosim::core::ParticleSoaStorage storage;
  storage.resize(particle_count);

  auto pos_x = storage.span<cosmosim::core::ParticleSoaField::kPosX>();
  auto pos_y = storage.span<cosmosim::core::ParticleSoaField::kPosY>();
  auto pos_z = storage.span<cosmosim::core::ParticleSoaField::kPosZ>();
  auto rho = storage.span<cosmosim::core::ParticleSoaField::kRho>();
  auto u_int = storage.span<cosmosim::core::ParticleSoaField::kUInt>();

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = static_cast<double>(i) * 0.01;
    pos_y[i] = static_cast<double>(i) * 0.02;
    pos_z[i] = static_cast<double>(i) * 0.03;
    rho[i] = 1.0 + static_cast<double>(i % 17U) * 0.1;
    u_int[i] = 0.5 + static_cast<double>(i % 13U) * 0.03;
  }

  const auto sweep_start = std::chrono::steady_clock::now();
  double sweep_checksum = 0.0;
  for (std::size_t i = 0; i < particle_count; ++i) {
    sweep_checksum += pos_x[i] + pos_y[i] + pos_z[i];
  }
  const auto sweep_stop = std::chrono::steady_clock::now();

  std::vector<std::uint32_t> active_indices;
  active_indices.reserve(particle_count / 4U);
  for (std::uint32_t i = 0; i < particle_count; i += 4U) {
    active_indices.push_back(i);
  }

  std::vector<double> rho_active(active_indices.size());
  std::vector<double> u_int_active(active_indices.size());
  const auto gather_start = std::chrono::steady_clock::now();
  cosmosim::core::gatherSpan<double>(
      std::as_const(storage).span<cosmosim::core::ParticleSoaField::kRho>(), active_indices, rho_active);
  cosmosim::core::gatherSpan<double>(
      std::as_const(storage).span<cosmosim::core::ParticleSoaField::kUInt>(), active_indices, u_int_active);
  const auto gather_stop = std::chrono::steady_clock::now();

  double gather_checksum = 0.0;
  for (std::size_t i = 0; i < rho_active.size(); ++i) {
    gather_checksum += rho_active[i] * u_int_active[i];
  }

  std::vector<std::uint8_t> keep_flags(particle_count, 0U);
  for (std::size_t i = 0; i < particle_count; ++i) {
    keep_flags[i] = static_cast<std::uint8_t>((i % 3U) != 0U);
  }

  const auto compact_start = std::chrono::steady_clock::now();
  const std::size_t compact_size = storage.stableCompact(keep_flags);
  const auto compact_stop = std::chrono::steady_clock::now();

  const auto sweep_us =
      std::chrono::duration_cast<std::chrono::microseconds>(sweep_stop - sweep_start).count();
  const auto gather_us =
      std::chrono::duration_cast<std::chrono::microseconds>(gather_stop - gather_start).count();
  const auto compact_us =
      std::chrono::duration_cast<std::chrono::microseconds>(compact_stop - compact_start).count();

  std::cout << "bench_soa_storage sweep_checksum=" << sweep_checksum << " sweep_us=" << sweep_us
            << " gather_checksum=" << gather_checksum << " gather_us=" << gather_us
            << " compact_size=" << compact_size << " compact_us=" << compact_us << '\n';
  return 0;
}
