#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "bench/reporting/bench_report.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/gravity/tree_gravity.hpp"

namespace {

struct ParticleBatch {
  std::vector<double> pos_x_comov;
  std::vector<double> pos_y_comov;
  std::vector<double> pos_z_comov;
  std::vector<double> mass_code;
};

[[nodiscard]] ParticleBatch makeParticles(std::size_t particle_count) {
  ParticleBatch particles;
  particles.pos_x_comov.resize(particle_count);
  particles.pos_y_comov.resize(particle_count);
  particles.pos_z_comov.resize(particle_count);
  particles.mass_code.resize(particle_count);

  for (std::size_t i = 0; i < particle_count; ++i) {
    particles.pos_x_comov[i] = std::fmod(static_cast<double>((37U * i + 13U) % 104729U) * 7.13e-5, 1.0);
    particles.pos_y_comov[i] = std::fmod(static_cast<double>((73U * i + 17U) % 130363U) * 5.11e-5, 1.0);
    particles.pos_z_comov[i] = std::fmod(static_cast<double>((97U * i + 19U) % 156007U) * 4.73e-5, 1.0);
    particles.mass_code[i] = 0.9 + static_cast<double>(i % 7U) * 0.05;
  }

  return particles;
}

}  // namespace

int main() {
  constexpr std::size_t k_particle_count = 65536;
  constexpr std::size_t k_active_count = 16384;
  constexpr cosmosim::gravity::PmGridShape k_grid_shape{64, 64, 64};

  const auto execution = cosmosim::bench::defaultExecutionConfig(1, 4);
  const auto particles = makeParticles(k_particle_count);

  std::vector<std::uint32_t> active_index(k_active_count);
  for (std::size_t i = 0; i < k_active_count; ++i) {
    active_index[i] = static_cast<std::uint32_t>((19U * i) % k_particle_count);
  }

  std::vector<double> tree_accel_x(k_active_count, 0.0);
  std::vector<double> tree_accel_y(k_active_count, 0.0);
  std::vector<double> tree_accel_z(k_active_count, 0.0);
  std::vector<double> pm_accel_x(k_particle_count, 0.0);
  std::vector<double> pm_accel_y(k_particle_count, 0.0);
  std::vector<double> pm_accel_z(k_particle_count, 0.0);

  cosmosim::gravity::TreeGravityOptions tree_options;
  tree_options.opening_theta = 0.7;
  tree_options.max_leaf_size = 24;
  tree_options.softening.epsilon_comoving = 3.0e-4;

  cosmosim::gravity::TreeGravitySolver tree_solver;
  cosmosim::gravity::TreeGravityProfile tree_profile{};

  for (std::size_t iter = 0; iter < execution.warmup_iterations; ++iter) {
    tree_solver.build(
        particles.pos_x_comov, particles.pos_y_comov, particles.pos_z_comov, particles.mass_code, tree_options, nullptr);
    tree_solver.evaluateActiveSet(
        particles.pos_x_comov,
        particles.pos_y_comov,
        particles.pos_z_comov,
        particles.mass_code,
        active_index,
        tree_accel_x,
        tree_accel_y,
        tree_accel_z,
        tree_options,
        nullptr);
  }

  const auto tree_begin = cosmosim::bench::BenchmarkClock::now();
  for (std::size_t iter = 0; iter < execution.measurement_iterations; ++iter) {
    tree_solver.build(
        particles.pos_x_comov, particles.pos_y_comov, particles.pos_z_comov, particles.mass_code, tree_options, &tree_profile);
    tree_solver.evaluateActiveSet(
        particles.pos_x_comov,
        particles.pos_y_comov,
        particles.pos_z_comov,
        particles.mass_code,
        active_index,
        tree_accel_x,
        tree_accel_y,
        tree_accel_z,
        tree_options,
        &tree_profile);
  }
  const auto tree_end = cosmosim::bench::BenchmarkClock::now();

  const double tree_ms = cosmosim::bench::BenchmarkClock::millisecondsBetween(tree_begin, tree_end);
  const double tree_targets_per_second =
      tree_ms > 0.0 ? static_cast<double>(k_active_count * execution.measurement_iterations) / (tree_ms * 1.0e-3) : 0.0;

  cosmosim::bench::BenchmarkReporter tree_report("bench_gravity_tree_kernel");
  cosmosim::bench::addExecutionFields(tree_report, execution);
  tree_report.addField("particle_count", k_particle_count);
  tree_report.addField("active_count", k_active_count);
  tree_report.addField("measurement_ms", tree_ms);
  tree_report.addField("build_ms", tree_profile.build_ms);
  tree_report.addField("multipole_ms", tree_profile.multipole_ms);
  tree_report.addField("traversal_ms", tree_profile.traversal_ms);
  tree_report.addField("visited_nodes", tree_profile.visited_nodes);
  tree_report.addField("accepted_nodes", tree_profile.accepted_nodes);
  tree_report.addField("particle_particle_interactions", tree_profile.particle_particle_interactions);
  tree_report.addField("targets_per_second", tree_targets_per_second);
  tree_report.write();

  cosmosim::bench::BenchmarkReporter pm_report("bench_gravity_pm_kernel");
  cosmosim::bench::addExecutionFields(pm_report, execution);
  pm_report.addField("particle_count", k_particle_count);
  pm_report.addField("grid_nx", k_grid_shape.nx);
  pm_report.addField("grid_ny", k_grid_shape.ny);
  pm_report.addField("grid_nz", k_grid_shape.nz);

  if (!cosmosim::gravity::PmSolver::fftBackendAvailable()) {
    pm_report.addField("skipped", "fft_backend_unavailable");
    pm_report.addField("fft_backend", cosmosim::gravity::PmSolver::fftBackendName());
    pm_report.write();
    return 0;
  }

  cosmosim::gravity::PmSolver pm_solver(k_grid_shape);
  cosmosim::gravity::PmGridStorage pm_grid(k_grid_shape);

  cosmosim::gravity::PmSolveOptions pm_options;
  pm_options.box_size_mpc_comoving = 10.0;
  pm_options.scale_factor = 1.0;
  pm_options.gravitational_constant_code = 1.0;
  pm_options.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kCic;

  cosmosim::gravity::PmProfileEvent pm_profile{};

  for (std::size_t iter = 0; iter < execution.warmup_iterations; ++iter) {
    pm_solver.solveForParticles(
        pm_grid,
        particles.pos_x_comov,
        particles.pos_y_comov,
        particles.pos_z_comov,
        particles.mass_code,
        pm_accel_x,
        pm_accel_y,
        pm_accel_z,
        pm_options,
        nullptr);
  }

  const auto pm_begin = cosmosim::bench::BenchmarkClock::now();
  for (std::size_t iter = 0; iter < execution.measurement_iterations; ++iter) {
    pm_solver.solveForParticles(
        pm_grid,
        particles.pos_x_comov,
        particles.pos_y_comov,
        particles.pos_z_comov,
        particles.mass_code,
        pm_accel_x,
        pm_accel_y,
        pm_accel_z,
        pm_options,
        &pm_profile);
  }
  const auto pm_end = cosmosim::bench::BenchmarkClock::now();

  const double pm_ms = cosmosim::bench::BenchmarkClock::millisecondsBetween(pm_begin, pm_end);
  const double particle_updates_per_second =
      pm_ms > 0.0 ? static_cast<double>(k_particle_count * execution.measurement_iterations) / (pm_ms * 1.0e-3) : 0.0;

  pm_report.addField("measurement_ms", pm_ms);
  pm_report.addField("assign_ms", pm_profile.assign_ms);
  pm_report.addField("fft_forward_ms", pm_profile.fft_forward_ms);
  pm_report.addField("poisson_ms", pm_profile.poisson_ms);
  pm_report.addField("gradient_ms", pm_profile.gradient_ms);
  pm_report.addField("fft_inverse_ms", pm_profile.fft_inverse_ms);
  pm_report.addField("interpolate_ms", pm_profile.interpolate_ms);
  pm_report.addField("particle_updates_per_second", particle_updates_per_second);
  cosmosim::bench::addBandwidthFields(pm_report, pm_profile.bytes_moved, pm_ms);
  pm_report.write();

  return 0;
}
