#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/gravity/pm_solver.hpp"

namespace {

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void testPmCudaSmoke() {
  if (!cosmosim::gravity::PmSolver::cudaBackendAvailable()) {
    return;
  }

  const cosmosim::gravity::PmGridShape shape{24, 12, 12};
  cosmosim::gravity::PmGridStorage grid(shape);
  cosmosim::gravity::PmSolver solver(shape);

  const std::size_t particle_count = 192;
  std::vector<double> pos_x(particle_count);
  std::vector<double> pos_y(particle_count);
  std::vector<double> pos_z(particle_count);
  std::vector<double> mass(particle_count, 1.0);
  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod(0.137 * static_cast<double>(i), 1.0);
    pos_y[i] = std::fmod(0.211 * static_cast<double>(i), 1.0);
    pos_z[i] = std::fmod(0.331 * static_cast<double>(i), 1.0);
    mass[i] = 1.0 + 0.2 * std::sin(6.283185307179586 * pos_x[i]);
  }

  std::vector<double> accel_x(particle_count, 0.0);
  std::vector<double> accel_y(particle_count, 0.0);
  std::vector<double> accel_z(particle_count, 0.0);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.execution_policy = cosmosim::core::ExecutionPolicy::kCuda;
  options.data_residency = cosmosim::gravity::PmDataResidencyPolicy::kPreferDevice;

  cosmosim::gravity::PmProfileEvent profile;
  solver.solveForParticles(grid, pos_x, pos_y, pos_z, mass, accel_x, accel_y, accel_z, options, &profile);

  double signal_norm = 0.0;
  for (std::size_t i = 0; i < particle_count; ++i) {
    signal_norm += accel_x[i] * accel_x[i] + accel_y[i] * accel_y[i] + accel_z[i] * accel_z[i];
  }

  requireOrThrow(std::isfinite(signal_norm), "CUDA PM smoke produced non-finite accelerations");
  requireOrThrow(signal_norm > 0.0, "CUDA PM smoke produced zero acceleration signal");
  requireOrThrow(profile.transfer_h2d_ms >= 0.0, "CUDA PM smoke transfer_h2d_ms should be non-negative");
  requireOrThrow(profile.transfer_d2h_ms >= 0.0, "CUDA PM smoke transfer_d2h_ms should be non-negative");
  requireOrThrow(profile.device_kernel_ms >= 0.0, "CUDA PM smoke device_kernel_ms should be non-negative");
}

}  // namespace

int main() {
  testPmCudaSmoke();
  return 0;
}
