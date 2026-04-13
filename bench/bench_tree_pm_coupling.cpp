#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"

int main() {
  constexpr std::size_t particle_count =
#if COSMOSIM_ENABLE_FFTW
      180000;
#else
      12000;
#endif
  constexpr std::size_t active_count =
#if COSMOSIM_ENABLE_FFTW
      65536;
#else
      4096;
#endif
  constexpr double box_size_comoving = 1.0;

  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 1.0);

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod(static_cast<double>((53U * i + 11U) % 104729U) * 9.31e-5, box_size_comoving);
    pos_y[i] = std::fmod(static_cast<double>((67U * i + 13U) % 130363U) * 7.73e-5, box_size_comoving);
    pos_z[i] = std::fmod(static_cast<double>((79U * i + 17U) % 156007U) * 6.29e-5, box_size_comoving);
    mass[i] = 0.7 + static_cast<double>(i % 13U) * 0.015;
  }

  std::vector<std::uint32_t> active(active_count, 0U);
  for (std::size_t i = 0; i < active_count; ++i) {
    active[i] = static_cast<std::uint32_t>((19U * i) % particle_count);
  }

  std::vector<double> accel_x(active_count, 0.0);
  std::vector<double> accel_y(active_count, 0.0);
  std::vector<double> accel_z(active_count, 0.0);

  cosmosim::gravity::TreePmForceAccumulatorView accumulator{active, accel_x, accel_y, accel_z};

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = box_size_comoving;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.enable_window_deconvolution = true;
  options.tree_options.opening_theta = 0.7;
  options.tree_options.max_leaf_size = 16;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 3.5e-4;
  options.split_policy.split_scale_comoving = 0.03;

  const cosmosim::gravity::PmGridShape pm_shape =
#if COSMOSIM_ENABLE_FFTW
      {128, 128, 128};
#else
      {16, 16, 16};
#endif
  cosmosim::gravity::TreePmCoordinator coordinator(pm_shape);

  cosmosim::gravity::TreePmProfileEvent profile;
  cosmosim::gravity::TreePmDiagnostics diagnostics;
  coordinator.solveActiveSet(pos_x, pos_y, pos_z, mass, accumulator, options, &profile, &diagnostics);

  const double active_mpart_s = static_cast<double>(active_count) / (profile.coupling_overhead_ms * 1.0e3) * 1.0e-6;
  const double bytes_proxy = static_cast<double>(profile.pm_profile.bytes_moved) /
      std::max(profile.pm_profile.assign_ms + profile.pm_profile.interpolate_ms, 1.0e-9) / 1.0e6;

  std::cout << "bench_tree_pm_coupling"
            << " build_type=manual"
            << " threads=" << std::max(1U, std::thread::hardware_concurrency())
            << " features=pm+tree+gaussian_split"
            << " particle_count=" << particle_count
            << " active_count=" << active_count
            << " assign_ms=" << profile.pm_profile.assign_ms
            << " fft_forward_ms=" << profile.pm_profile.fft_forward_ms
            << " poisson_ms=" << profile.pm_profile.poisson_ms
            << " gradient_ms=" << profile.pm_profile.gradient_ms
            << " fft_inverse_ms=" << profile.pm_profile.fft_inverse_ms
            << " interpolate_ms=" << profile.pm_profile.interpolate_ms
            << " tree_build_ms=" << profile.tree_profile.build_ms
            << " tree_multipole_ms=" << profile.tree_profile.multipole_ms
            << " tree_traversal_ms=" << profile.tree_profile.traversal_ms
            << " coupling_overhead_ms=" << profile.coupling_overhead_ms
            << " split_scale_comoving=" << diagnostics.split_scale_comoving
            << " split_composition_error=" << diagnostics.composition_error_at_split
            << " bytes_proxy_mb_per_s=" << bytes_proxy
            << " active_mpart_s=" << active_mpart_s
            << '\n';

  return 0;
}
