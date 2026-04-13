#include <chrono>
#include <cstddef>
#include <iostream>

#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/units.hpp"

int main() {
  cosmosim::core::CosmologyBackgroundConfig cfg;
  cfg.hubble_param = 0.674;
  cfg.omega_matter = 0.315;
  cfg.omega_lambda = 0.685;
  cosmosim::core::LambdaCdmBackground background(cfg);
  const auto units = cosmosim::core::makeUnitSystem("mpc", "msun", "km_s");

  constexpr std::size_t iterations = 500000;
  double checksum = 0.0;

  const auto start = std::chrono::steady_clock::now();
  for (std::size_t i = 0; i < iterations; ++i) {
    const double a = 0.05 + (0.95 * static_cast<double>(i % 1000) / 1000.0);
    checksum += background.hubbleSi(a);
    checksum += background.criticalDensitySi(a);
    checksum += units.densityCodeToSi(1.0 + a);
    checksum += cosmosim::core::comovingToPhysicalLength(10.0, a);
  }
  const auto stop = std::chrono::steady_clock::now();

  const auto elapsed_us =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

  std::cout << "bench_units_cosmology_smoke iterations=" << iterations
            << " checksum=" << checksum << " elapsed_us=" << elapsed_us << '\n';
  return 0;
}
