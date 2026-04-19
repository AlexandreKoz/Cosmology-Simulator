#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/pm_solver.hpp"

int main() {
  const bool fft_backend_fast = cosmosim::gravity::PmSolver::fftBackendAvailable();
  const cosmosim::gravity::PmGridShape shape = fft_backend_fast
      ? cosmosim::gravity::PmGridShape{64, 64, 64}
      : cosmosim::gravity::PmGridShape{12, 12, 12};
  cosmosim::gravity::PmGridStorage grid(shape);
  cosmosim::gravity::PmSolver solver(shape);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 32.0;
  options.scale_factor = 1.0;
  options.gravitational_constant_code = 1.0;

  const std::size_t particle_count = fft_backend_fast ? 256000 : 4096;
  std::vector<double> pos_x(particle_count);
  std::vector<double> pos_y(particle_count);
  std::vector<double> pos_z(particle_count);
  std::vector<double> mass(particle_count, 1.0);
  std::vector<double> accel_x(particle_count, 0.0);
  std::vector<double> accel_y(particle_count, 0.0);
  std::vector<double> accel_z(particle_count, 0.0);

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod(static_cast<double>((37U * i) % 100000U) * 0.00037, options.box_size_mpc_comoving);
    pos_y[i] = std::fmod(static_cast<double>((53U * i) % 100000U) * 0.00029, options.box_size_mpc_comoving);
    pos_z[i] = std::fmod(static_cast<double>((97U * i) % 100000U) * 0.00019, options.box_size_mpc_comoving);
  }

  cosmosim::gravity::PmProfileEvent profile;
  constexpr std::size_t warmup_solves = 1;
  constexpr std::size_t measured_solves = 8;

  const auto setup_start = std::chrono::steady_clock::now();
  for (std::size_t i = 0; i < warmup_solves; ++i) {
    solver.assignDensity(grid, pos_x, pos_y, pos_z, mass, options, &profile);
    solver.solvePoissonPeriodic(grid, options, &profile);
    solver.interpolateForces(grid, pos_x, pos_y, pos_z, accel_x, accel_y, accel_z, options, &profile);
  }
  const auto setup_stop = std::chrono::steady_clock::now();

  cosmosim::gravity::PmProfileEvent measured_profile;
  const auto steady_start = std::chrono::steady_clock::now();
  for (std::size_t i = 0; i < measured_solves; ++i) {
    solver.assignDensity(grid, pos_x, pos_y, pos_z, mass, options, &measured_profile);
    solver.solvePoissonPeriodic(grid, options, &measured_profile);
    solver.interpolateForces(grid, pos_x, pos_y, pos_z, accel_x, accel_y, accel_z, options, &measured_profile);
  }
  const auto steady_stop = std::chrono::steady_clock::now();

  const auto setup_ms = std::chrono::duration<double, std::milli>(setup_stop - setup_start).count();
  const auto steady_ms = std::chrono::duration<double, std::milli>(steady_stop - steady_start).count();

  const double assign_throughput_mpart_s =
      static_cast<double>(particle_count * measured_solves) / (measured_profile.assign_ms * 1.0e3);
  const double interp_throughput_mpart_s =
      static_cast<double>(particle_count * measured_solves) / (measured_profile.interpolate_ms * 1.0e3);
  const double bytes_per_s = static_cast<double>(measured_profile.bytes_moved) / (steady_ms * 1.0e-3);

  std::cout << "bench_pm_solver"
            << " backend=" << cosmosim::gravity::PmSolver::fftBackendName()
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " threads=1"
            << " setup_ms=" << setup_ms
            << " steady_ms=" << steady_ms
            << " warmup_solves=" << warmup_solves
            << " measured_solves=" << measured_solves
            << " assign_ms=" << measured_profile.assign_ms
            << " fft_forward_ms=" << measured_profile.fft_forward_ms
            << " poisson_ms=" << measured_profile.poisson_ms
            << " gradient_ms=" << measured_profile.gradient_ms
            << " fft_inverse_ms=" << measured_profile.fft_inverse_ms
            << " interpolate_ms=" << measured_profile.interpolate_ms
            << " transfer_h2d_ms=" << measured_profile.transfer_h2d_ms
            << " transfer_d2h_ms=" << measured_profile.transfer_d2h_ms
            << " device_kernel_ms=" << measured_profile.device_kernel_ms
            << " plan_cache_size=" << solver.cachedPlanCount()
            << " plan_build_count=" << solver.planBuildCount()
            << " assign_mpart_s=" << assign_throughput_mpart_s * 1.0e-6
            << " interp_mpart_s=" << interp_throughput_mpart_s * 1.0e-6
            << " bytes_moved=" << measured_profile.bytes_moved
            << " effective_bandwidth_gb_s=" << bytes_per_s * 1.0e-9
            << '\n';

  return 0;
}
