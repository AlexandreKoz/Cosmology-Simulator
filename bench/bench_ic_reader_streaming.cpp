#include <chrono>
#include <cstddef>
#include <iostream>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/ic_reader.hpp"

int main() {
  cosmosim::core::SimulationConfig config;
  config.output.run_name = "bench_ic_reader";

  constexpr std::size_t iterations = 40;
  constexpr std::size_t dm_count = 1u << 18;
  constexpr std::size_t gas_count = 1u << 16;

  const auto setup_start = std::chrono::steady_clock::now();
  const auto warmup = cosmosim::io::buildGeneratedIsolatedIc(config, dm_count, gas_count);
  const auto setup_end = std::chrono::steady_clock::now();

  std::size_t checksum = warmup.state.particles.size();
  const auto steady_start = std::chrono::steady_clock::now();
  for (std::size_t i = 0; i < iterations; ++i) {
    const auto result = cosmosim::io::buildGeneratedIsolatedIc(config, dm_count, gas_count, 10 + i * 1000);
    checksum += result.state.particle_sidecar.particle_id[0];
  }
  const auto steady_end = std::chrono::steady_clock::now();

  const auto setup_us =
      std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start).count();
  const auto steady_us =
      std::chrono::duration_cast<std::chrono::microseconds>(steady_end - steady_start).count();

  const double particles_loaded = static_cast<double>(iterations) * static_cast<double>(dm_count + gas_count);
  const double particles_per_second = particles_loaded / (static_cast<double>(steady_us) * 1.0e-6);

  std::cout << "bench_ic_reader_streaming build=unknown"
            << " hdf5_enabled=" << COSMOSIM_ENABLE_HDF5 << " threads=1"
            << " setup_us=" << setup_us << " steady_us=" << steady_us
            << " particles_per_second=" << particles_per_second << " checksum=" << checksum << '\n';

  return 0;
}
