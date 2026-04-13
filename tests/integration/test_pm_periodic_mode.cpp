#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/pm_solver.hpp"

namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void testPeriodicSinusoidalForceResponse() {
  const cosmosim::gravity::PmGridShape shape{32, 8, 8};
  cosmosim::gravity::PmGridStorage grid(shape);
  cosmosim::gravity::PmSolver solver(shape);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 0.8;
  options.gravitational_constant_code = 1.0;

  const std::size_t particle_count = 32;
  std::vector<double> pos_x(particle_count);
  std::vector<double> pos_y(particle_count, 0.25);
  std::vector<double> pos_z(particle_count, 0.75);
  std::vector<double> mass(particle_count, 1.0);
  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = (static_cast<double>(i) + 0.5) / static_cast<double>(particle_count) * options.box_size_mpc_comoving;
    mass[i] = 1.0 + 0.1 * std::sin(2.0 * k_pi * pos_x[i] / options.box_size_mpc_comoving);
  }

  std::vector<double> accel_x(particle_count, 0.0);
  std::vector<double> accel_y(particle_count, 0.0);
  std::vector<double> accel_z(particle_count, 0.0);
  cosmosim::gravity::PmProfileEvent profile;

  solver.solveForParticles(grid, pos_x, pos_y, pos_z, mass, accel_x, accel_y, accel_z, options, &profile);

  const double kx = 2.0 * k_pi / options.box_size_mpc_comoving;
  const double expected_amp = 4.0 * k_pi * options.gravitational_constant_code *
      options.scale_factor * options.scale_factor * 0.1 / kx;

  double corr = 0.0;
  double norm_expected = 0.0;
  double norm_got = 0.0;
  double max_transverse = 0.0;
  for (std::size_t i = 0; i < particle_count; ++i) {
    const double expected = expected_amp * std::cos(kx * pos_x[i]);
    corr += expected * accel_x[i];
    norm_expected += expected * expected;
    norm_got += accel_x[i] * accel_x[i];
    max_transverse = std::max(max_transverse, std::max(std::abs(accel_y[i]), std::abs(accel_z[i])));
  }

  const double cosine_similarity = corr / std::sqrt(std::max(norm_expected * norm_got, 1.0e-20));
#if COSMOSIM_ENABLE_FFTW
  const double min_cosine_similarity = 0.9;
#else
  const double min_cosine_similarity = 0.05;
#endif

  std::ostringstream diag;
  diag << "PM periodic sinusoidal response validation failed: backend='" << cosmosim::gravity::PmSolver::fftBackendName()
       << "', cosine_similarity=" << cosine_similarity
       << " (required >= " << min_cosine_similarity << ")"
       << ", transverse_max=" << max_transverse
       << ", signal_norm=" << std::sqrt(norm_got)
       << ", build_flag.COSMOSIM_ENABLE_FFTW=" << (COSMOSIM_ENABLE_FFTW ? "ON" : "OFF");
  requireOrThrow(std::isfinite(cosine_similarity), diag.str());
  requireOrThrow(norm_got > 0.0, diag.str());
  requireOrThrow(std::abs(cosine_similarity) > min_cosine_similarity, diag.str());
  requireOrThrow(max_transverse < 5.0e-2, diag.str());

  requireOrThrow(profile.assign_ms >= 0.0, "PM profile.assign_ms must be non-negative");
  requireOrThrow(profile.fft_forward_ms >= 0.0, "PM profile.fft_forward_ms must be non-negative");
  requireOrThrow(profile.fft_inverse_ms >= 0.0, "PM profile.fft_inverse_ms must be non-negative");
  requireOrThrow(profile.bytes_moved > 0, "PM profile.bytes_moved must be positive");
}

}  // namespace

int main() {
  testPeriodicSinusoidalForceResponse();
  return 0;
}
