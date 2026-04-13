#include <chrono>
#include <cstddef>
#include <iostream>
#include <string>

#include "cosmosim/core/config.hpp"

int main() {
  const std::string config_text = R"(
[mode]
mode = zoom_in
ic_file = ics_zoom.hdf5
zoom_high_res_region = true
zoom_region_file = region_zoom.hdf5

[cosmology]
omega_matter = 0.31
omega_lambda = 0.69
box_size = 50 mpc

[numerics]
time_begin_code = 0.0
time_end_code = 1.0
max_global_steps = 1024
gravity_softening = 1.0 kpc

[output]
run_name = bench_parser
output_stem = snapshot
restart_stem = restart
)";

  constexpr std::size_t iterations = 5000;
  std::size_t checksum = 0;
  const auto start = std::chrono::steady_clock::now();
  for (std::size_t i = 0; i < iterations; ++i) {
    const auto frozen = cosmosim::core::loadFrozenConfigFromString(config_text, "bench");
    checksum += frozen.provenance.config_hash_hex.size();
  }
  const auto stop = std::chrono::steady_clock::now();
  const auto elapsed_us =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

  std::cout << "bench_config_parser iterations=" << iterations << " checksum=" << checksum
            << " elapsed_us=" << elapsed_us << '\n';
  return 0;
}
