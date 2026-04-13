#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/pm_solver.hpp"

int main() {
  if (!cosmosim::gravity::PmSolver::cudaBackendAvailable()) {
    std::cout << "bench_pm_gpu_transfer"
              << " backend=unavailable"
              << " build_type=" << COSMOSIM_BUILD_TYPE
              << " threads=1"
              << " note=no_cuda_device"
              << '\n';
    return 0;
  }

  const cosmosim::gravity::PmGridShape shape{64, 64, 64};
  cosmosim::gravity::PmGridStorage grid(shape);
  cosmosim::gravity::PmSolver solver(shape);

  const std::size_t particle_count = 200000;
  std::vector<double> pos_x(particle_count);
  std::vector<double> pos_y(particle_count);
  std::vector<double> pos_z(particle_count);
  std::vector<double> mass(particle_count, 1.0);
  std::vector<double> accel_x(particle_count, 0.0);
  std::vector<double> accel_y(particle_count, 0.0);
  std::vector<double> accel_z(particle_count, 0.0);

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod(static_cast<double>((37U * i) % 100000U) * 0.00037, 1.0);
    pos_y[i] = std::fmod(static_cast<double>((53U * i) % 100000U) * 0.00029, 1.0);
    pos_z[i] = std::fmod(static_cast<double>((97U * i) % 100000U) * 0.00019, 1.0);
  }

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.execution_policy = cosmosim::core::ExecutionPolicy::kCuda;
  options.data_residency = cosmosim::gravity::PmDataResidencyPolicy::kPreferDevice;

  cosmosim::gravity::PmProfileEvent profile;
  const auto t0 = std::chrono::steady_clock::now();
  solver.solveForParticles(grid, pos_x, pos_y, pos_z, mass, accel_x, accel_y, accel_z, options, &profile);
  const auto t1 = std::chrono::steady_clock::now();

  const double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  const double transfer_total_ms = profile.transfer_h2d_ms + profile.transfer_d2h_ms;
  const double kernel_total_ms = profile.device_kernel_ms;
  const double transfer_bytes = static_cast<double>(profile.bytes_moved);
  const double transfer_bw_gb_s = transfer_bytes / (std::max(transfer_total_ms, 1.0e-9) * 1.0e-3) * 1.0e-9;

  std::cout << "bench_pm_gpu_transfer"
            << " backend=" << cosmosim::gravity::PmSolver::fftBackendName() << "+cuda"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " threads=1"
            << " particle_count=" << particle_count
            << " setup_plus_steady_ms=" << total_ms
            << " transfer_h2d_ms=" << profile.transfer_h2d_ms
            << " transfer_d2h_ms=" << profile.transfer_d2h_ms
            << " transfer_total_ms=" << transfer_total_ms
            << " device_kernel_ms=" << kernel_total_ms
            << " fft_forward_ms=" << profile.fft_forward_ms
            << " fft_inverse_ms=" << profile.fft_inverse_ms
            << " bytes_moved=" << profile.bytes_moved
            << " transfer_effective_bandwidth_gb_s=" << transfer_bw_gb_s
            << '\n';

  return 0;
}
