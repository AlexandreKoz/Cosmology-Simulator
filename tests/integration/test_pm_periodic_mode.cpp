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

void testPeriodicSinusoidalForceResponseForScheme(cosmosim::gravity::PmAssignmentScheme scheme) {
  const cosmosim::gravity::PmGridShape shape{32, 8, 8};
  cosmosim::gravity::PmGridStorage grid(shape);
  cosmosim::gravity::PmSolver solver(shape);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 0.8;
  options.gravitational_constant_code = 1.0;
  options.assignment_scheme = scheme;

  const std::size_t particle_count = shape.nx * shape.ny * shape.nz;
  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 1.0);
  std::vector<double> accel_x(particle_count, 0.0);
  std::vector<double> accel_y(particle_count, 0.0);
  std::vector<double> accel_z(particle_count, 0.0);
  std::vector<double> phi_particles(particle_count, 0.0);

  const double amplitude = 0.1;
  const double kx = 2.0 * k_pi / options.box_size_mpc_comoving;
  std::size_t slot = 0;
  for (std::size_t ix = 0; ix < shape.nx; ++ix) {
    const double x = (static_cast<double>(ix) + 0.5) / static_cast<double>(shape.nx) * options.box_size_mpc_comoving;
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      const double y = (static_cast<double>(iy) + 0.5) / static_cast<double>(shape.ny) * options.box_size_mpc_comoving;
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        const double z = (static_cast<double>(iz) + 0.5) / static_cast<double>(shape.nz) * options.box_size_mpc_comoving;
        pos_x[slot] = x;
        pos_y[slot] = y;
        pos_z[slot] = z;
        mass[slot] = 1.0 + amplitude * std::sin(kx * x);
        ++slot;
      }
    }
  }

  cosmosim::gravity::PmProfileEvent profile;
  solver.solveForParticles(grid, pos_x, pos_y, pos_z, mass, accel_x, accel_y, accel_z, options, &profile);
  solver.interpolatePotential(grid, pos_x, pos_y, pos_z, phi_particles, options, nullptr);

  const double expected_force_amp = 4.0 * k_pi * options.gravitational_constant_code *
      options.scale_factor * options.scale_factor * amplitude / kx;
  const double expected_phi_amp = -4.0 * k_pi * options.gravitational_constant_code *
      options.scale_factor * options.scale_factor * amplitude / (kx * kx);

  double corr_force = 0.0;
  double norm_force_expected = 0.0;
  double norm_force_got = 0.0;
  double corr_phi = 0.0;
  double norm_phi_expected = 0.0;
  double norm_phi_got = 0.0;
  double consistency_rms = 0.0;
  double consistency_ref_rms = 0.0;
  double max_transverse = 0.0;
  const double dx = options.box_size_mpc_comoving / static_cast<double>(shape.nx);
  for (std::size_t ix = 0; ix < shape.nx; ++ix) {
    const double x = (static_cast<double>(ix) + 0.5) / static_cast<double>(shape.nx) * options.box_size_mpc_comoving;
    const double expected_force = expected_force_amp * std::cos(kx * x);
    const double expected_phi = expected_phi_amp * std::sin(kx * x);
    const std::size_t ix_prev = (ix + shape.nx - 1U) % shape.nx;
    const std::size_t ix_next = (ix + 1U) % shape.nx;
    const std::size_t sample = (ix * shape.ny + 0U) * shape.nz + 0U;
    const std::size_t grid_index = grid.linearIndex(ix, 0, 0);
    const std::size_t prev = (ix_prev * shape.ny + 0U) * shape.nz + 0U;
    const std::size_t next = (ix_next * shape.ny + 0U) * shape.nz + 0U;
    const double accel_from_phi = -(phi_particles[next] - phi_particles[prev]) / (2.0 * dx);

    corr_force += expected_force * accel_x[sample];
    norm_force_expected += expected_force * expected_force;
    norm_force_got += accel_x[sample] * accel_x[sample];
    corr_phi += expected_phi * phi_particles[sample];
    norm_phi_expected += expected_phi * expected_phi;
    norm_phi_got += phi_particles[sample] * phi_particles[sample];
    const double consistency_force = grid.force_x()[grid_index];
    const double consistency_error = accel_from_phi - consistency_force;
    consistency_rms += consistency_error * consistency_error;
    consistency_ref_rms += consistency_force * consistency_force;
    max_transverse = std::max(max_transverse, std::max(std::abs(accel_y[sample]), std::abs(accel_z[sample])));
  }

  const double force_cosine_similarity = corr_force / std::sqrt(std::max(norm_force_expected * norm_force_got, 1.0e-20));
  const double phi_cosine_similarity = corr_phi / std::sqrt(std::max(norm_phi_expected * norm_phi_got, 1.0e-20));
  const double force_consistency_rel = std::sqrt(consistency_rms / std::max(consistency_ref_rms, 1.0e-20));
#if COSMOSIM_ENABLE_FFTW
  const double min_cosine_similarity = 0.9;
  const double max_consistency_rel = 0.12;
#else
  const double min_cosine_similarity = 0.05;
  const double max_consistency_rel = 5.0;
#endif

  std::ostringstream diag;
  diag << "PM periodic sinusoidal response validation failed: backend='" << cosmosim::gravity::PmSolver::fftBackendName()
       << "', force_cosine_similarity=" << force_cosine_similarity
       << " (required >= " << min_cosine_similarity << ")"
       << ", potential_cosine_similarity=" << phi_cosine_similarity
       << ", force_from_potential_rel_l2=" << force_consistency_rel
       << " (required <= " << max_consistency_rel << ")"
       << ", transverse_max=" << max_transverse
       << ", force_signal_norm=" << std::sqrt(norm_force_got)
       << ", potential_signal_norm=" << std::sqrt(norm_phi_got)
       << ", build_flag.COSMOSIM_ENABLE_FFTW=" << (COSMOSIM_ENABLE_FFTW ? "ON" : "OFF");
  requireOrThrow(std::isfinite(force_cosine_similarity), diag.str());
  requireOrThrow(std::isfinite(phi_cosine_similarity), diag.str());
  requireOrThrow(std::isfinite(force_consistency_rel), diag.str());
  requireOrThrow(norm_force_got > 0.0, diag.str());
  requireOrThrow(norm_phi_got > 0.0, diag.str());
  requireOrThrow(std::abs(force_cosine_similarity) > min_cosine_similarity, diag.str());
  requireOrThrow(std::abs(phi_cosine_similarity) > min_cosine_similarity, diag.str());
  requireOrThrow(force_consistency_rel <= max_consistency_rel, diag.str());
  requireOrThrow(max_transverse < 5.0e-2, diag.str());

  requireOrThrow(profile.assign_ms >= 0.0, "PM profile.assign_ms must be non-negative");
  requireOrThrow(profile.fft_forward_ms >= 0.0, "PM profile.fft_forward_ms must be non-negative");
  requireOrThrow(profile.fft_inverse_ms >= 0.0, "PM profile.fft_inverse_ms must be non-negative");
  requireOrThrow(profile.bytes_moved > 0, "PM profile.bytes_moved must be positive");
}

}  // namespace

int main() {
  testPeriodicSinusoidalForceResponseForScheme(cosmosim::gravity::PmAssignmentScheme::kCic);
  testPeriodicSinusoidalForceResponseForScheme(cosmosim::gravity::PmAssignmentScheme::kTsc);
  return 0;
}
