#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/gravity/pm_solver.hpp"

namespace {

constexpr double k_tolerance = 1.0e-8;
constexpr double k_pi = 3.141592653589793238462643383279502884;

void testCicMassConservation() {
  const cosmosim::gravity::PmGridShape shape{16, 8, 4};
  cosmosim::gravity::PmGridStorage grid(shape);
  cosmosim::gravity::PmSolver solver(shape);

  const std::vector<double> pos_x{0.1, 2.1, 3.9, 7.2};
  const std::vector<double> pos_y{0.2, 1.2, 2.2, 3.2};
  const std::vector<double> pos_z{0.3, 0.4, 0.5, 0.6};
  const std::vector<double> mass{2.0, 3.0, 4.0, 5.0};

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 8.0;
  options.scale_factor = 1.0;

  solver.assignDensity(grid, pos_x, pos_y, pos_z, mass, options, nullptr);

  const double total_density = std::accumulate(grid.density().begin(), grid.density().end(), 0.0);
  const double total_mass = std::accumulate(mass.begin(), mass.end(), 0.0);
  const double cell_volume = std::pow(options.box_size_mpc_comoving, 3.0) / static_cast<double>(shape.cellCount());
  assert(std::abs(total_density * cell_volume - total_mass) < k_tolerance);
}

void testPoissonAnalyticMode() {
  const cosmosim::gravity::PmGridShape shape{16, 8, 8};
  cosmosim::gravity::PmGridStorage grid(shape);
  cosmosim::gravity::PmSolver solver(shape);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 0.5;
  options.gravitational_constant_code = 1.2;

  const double amplitude = 0.25;
  const double kx = 2.0 * k_pi / options.box_size_mpc_comoving;
  for (std::size_t ix = 0; ix < shape.nx; ++ix) {
    const double x = (static_cast<double>(ix) + 0.5) / static_cast<double>(shape.nx) * options.box_size_mpc_comoving;
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        grid.density()[grid.linearIndex(ix, iy, iz)] = amplitude * std::sin(kx * x);
      }
    }
  }

  solver.solvePoissonPeriodic(grid, options, nullptr);

  const double expected_amp = 4.0 * k_pi * options.gravitational_constant_code *
      options.scale_factor * options.scale_factor * amplitude / kx;

  double corr = 0.0;
  double norm_expected = 0.0;
  double norm_got = 0.0;
  for (std::size_t ix = 0; ix < shape.nx; ++ix) {
    const double x = (static_cast<double>(ix) + 0.5) / static_cast<double>(shape.nx) * options.box_size_mpc_comoving;
    const double expected = expected_amp * std::cos(kx * x);
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        const double got = grid.force_x()[grid.linearIndex(ix, iy, iz)];
        corr += expected * got;
        norm_expected += expected * expected;
        norm_got += got * got;
      }
    }
  }

  const double cosine_similarity = corr / std::sqrt(std::max(norm_expected * norm_got, 1.0e-20));
#if COSMOSIM_ENABLE_FFTW
  assert(cosine_similarity > 0.98);
#else
  assert(std::isfinite(cosine_similarity));
  assert(norm_got > 0.0);
#endif
}

void testTreePmBuildGate() {
  assert(cosmosim::gravity::treePmSupportedByBuild());
  cosmosim::gravity::requireTreePmSupportOrThrow(cosmosim::core::GravitySolver::kTreePm);
  const std::string backend = cosmosim::gravity::PmSolver::fftBackendName();
  assert(!backend.empty());
}

void testExecutionPolicyValidation() {
  const cosmosim::gravity::PmGridShape shape{8, 8, 8};
  cosmosim::gravity::PmGridStorage grid(shape);
  cosmosim::gravity::PmSolver solver(shape);

  const std::vector<double> pos_x{0.1};
  const std::vector<double> pos_y{0.2};
  const std::vector<double> pos_z{0.3};
  const std::vector<double> mass{1.0};
  std::vector<double> accel_x(1, 0.0);
  std::vector<double> accel_y(1, 0.0);
  std::vector<double> accel_z(1, 0.0);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.execution_policy = cosmosim::core::ExecutionPolicy::kCuda;
  options.data_residency = cosmosim::gravity::PmDataResidencyPolicy::kHostOnly;

  bool threw = false;
  try {
    solver.solveForParticles(grid, pos_x, pos_y, pos_z, mass, accel_x, accel_y, accel_z, options, nullptr);
  } catch (const std::invalid_argument&) {
    threw = true;
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}

void testDeviceCpuAgreementWhenCudaAvailable() {
  if (!cosmosim::gravity::PmSolver::cudaBackendAvailable()) {
    return;
  }

  const cosmosim::gravity::PmGridShape shape{16, 8, 8};
  cosmosim::gravity::PmGridStorage grid_cpu(shape);
  cosmosim::gravity::PmGridStorage grid_gpu(shape);
  cosmosim::gravity::PmSolver solver(shape);

  const std::size_t particle_count = 64;
  std::vector<double> pos_x(particle_count);
  std::vector<double> pos_y(particle_count);
  std::vector<double> pos_z(particle_count);
  std::vector<double> mass(particle_count, 1.0);
  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod(0.091 * static_cast<double>(i), 1.0);
    pos_y[i] = std::fmod(0.173 * static_cast<double>(i), 1.0);
    pos_z[i] = std::fmod(0.287 * static_cast<double>(i), 1.0);
  }

  std::vector<double> accel_x_cpu(particle_count, 0.0);
  std::vector<double> accel_y_cpu(particle_count, 0.0);
  std::vector<double> accel_z_cpu(particle_count, 0.0);
  std::vector<double> accel_x_gpu(particle_count, 0.0);
  std::vector<double> accel_y_gpu(particle_count, 0.0);
  std::vector<double> accel_z_gpu(particle_count, 0.0);

  cosmosim::gravity::PmSolveOptions cpu_options;
  cpu_options.box_size_mpc_comoving = 1.0;
  cpu_options.execution_policy = cosmosim::core::ExecutionPolicy::kHostSerial;

  cosmosim::gravity::PmSolveOptions gpu_options = cpu_options;
  gpu_options.execution_policy = cosmosim::core::ExecutionPolicy::kCuda;
  gpu_options.data_residency = cosmosim::gravity::PmDataResidencyPolicy::kPreferDevice;

  solver.solveForParticles(
      grid_cpu, pos_x, pos_y, pos_z, mass, accel_x_cpu, accel_y_cpu, accel_z_cpu, cpu_options, nullptr);
  solver.solveForParticles(
      grid_gpu, pos_x, pos_y, pos_z, mass, accel_x_gpu, accel_y_gpu, accel_z_gpu, gpu_options, nullptr);

  constexpr double agreement_tol = 1.0e-6;
  for (std::size_t i = 0; i < particle_count; ++i) {
    assert(std::abs(accel_x_cpu[i] - accel_x_gpu[i]) < agreement_tol);
    assert(std::abs(accel_y_cpu[i] - accel_y_gpu[i]) < agreement_tol);
    assert(std::abs(accel_z_cpu[i] - accel_z_gpu[i]) < agreement_tol);
  }
}

}  // namespace

int main() {
  testCicMassConservation();
  testPoissonAnalyticMode();
  testTreePmBuildGate();
  testExecutionPolicyValidation();
  testDeviceCpuAgreementWhenCudaAvailable();
  return 0;
}
