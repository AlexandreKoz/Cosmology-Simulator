#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "validation_tolerance.hpp"

namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void directSumAcceleration(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    std::span<const std::uint32_t> active,
    const cosmosim::gravity::TreeGravityOptions& options,
    std::span<double> ax,
    std::span<double> ay,
    std::span<double> az) {
  for (std::size_t i = 0; i < active.size(); ++i) {
    const std::uint32_t target = active[i];
    const double tx = pos_x[target];
    const double ty = pos_y[target];
    const double tz = pos_z[target];
    double acc_x = 0.0;
    double acc_y = 0.0;
    double acc_z = 0.0;
    for (std::size_t j = 0; j < mass.size(); ++j) {
      if (j == target) {
        continue;
      }
      const double dx = pos_x[j] - tx;
      const double dy = pos_y[j] - ty;
      const double dz = pos_z[j] - tz;
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double factor = options.gravitational_constant_code * mass[j] *
          cosmosim::gravity::softenedInvR3(r2, options.softening);
      acc_x += factor * dx;
      acc_y += factor * dy;
      acc_z += factor * dz;
    }
    ax[i] = acc_x;
    ay[i] = acc_y;
    az[i] = acc_z;
  }
}

[[nodiscard]] double meanTreeRelativeError(const cosmosim::gravity::TreeGravityOptions& options) {
  constexpr std::size_t particle_count = 128;
  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 0.0);

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod(static_cast<double>((29U * i + 7U) % 997U) * 0.013, 1.0);
    pos_y[i] = std::fmod(static_cast<double>((47U * i + 11U) % 991U) * 0.017, 1.0);
    pos_z[i] = std::fmod(static_cast<double>((71U * i + 3U) % 983U) * 0.019, 1.0);
    mass[i] = 0.5 + static_cast<double>((13U * i) % 9U) * 0.05;
  }

  std::vector<std::uint32_t> active(32);
  for (std::size_t i = 0; i < active.size(); ++i) {
    active[i] = static_cast<std::uint32_t>(i * 2U);
  }

  std::vector<double> tree_ax(active.size(), 0.0);
  std::vector<double> tree_ay(active.size(), 0.0);
  std::vector<double> tree_az(active.size(), 0.0);

  cosmosim::gravity::TreeGravitySolver solver;
  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);
  solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, tree_ax, tree_ay, tree_az, options, nullptr);

  std::vector<double> direct_ax(active.size(), 0.0);
  std::vector<double> direct_ay(active.size(), 0.0);
  std::vector<double> direct_az(active.size(), 0.0);
  directSumAcceleration(pos_x, pos_y, pos_z, mass, active, options, direct_ax, direct_ay, direct_az);

  double mean_relative_error = 0.0;
  for (std::size_t i = 0; i < active.size(); ++i) {
    const double direct_norm = std::sqrt(
        direct_ax[i] * direct_ax[i] + direct_ay[i] * direct_ay[i] + direct_az[i] * direct_az[i]);
    const double diff_norm = std::sqrt(
        (tree_ax[i] - direct_ax[i]) * (tree_ax[i] - direct_ax[i]) +
        (tree_ay[i] - direct_ay[i]) * (tree_ay[i] - direct_ay[i]) +
        (tree_az[i] - direct_az[i]) * (tree_az[i] - direct_az[i]));
    mean_relative_error += diff_norm / std::max(direct_norm, 1.0e-12);
  }

  return mean_relative_error / static_cast<double>(active.size());
}

cosmosim::hydro::HydroPatchGeometry makePeriodic1dGeometry(std::size_t cell_count) {
  cosmosim::hydro::HydroPatchGeometry geometry;
  geometry.cell_volume_comoving = 1.0 / static_cast<double>(cell_count);
  geometry.faces.reserve(cell_count);
  for (std::size_t i = 0; i < cell_count; ++i) {
    geometry.faces.push_back(cosmosim::hydro::HydroFace{
        .owner_cell = i,
        .neighbor_cell = (i + 1U) % cell_count,
        .area_comoving = 1.0,
        .normal_x = 1.0,
        .normal_y = 0.0,
        .normal_z = 0.0});
  }
  return geometry;
}

[[nodiscard]] double runHydroSineWaveL1(std::size_t cell_count) {
  constexpr double k_gamma = 1.4;
  constexpr double k_final_time_code = 0.05;
  constexpr double k_velocity = 0.3;

  cosmosim::hydro::HydroConservedStateSoa conserved(cell_count);

  for (std::size_t i = 0; i < cell_count; ++i) {
    const double x = (static_cast<double>(i) + 0.5) / static_cast<double>(cell_count);
    cosmosim::hydro::HydroPrimitiveState primitive{};
    primitive.rho_comoving = 1.0 + 0.2 * std::sin(2.0 * k_pi * x);
    primitive.vel_x_peculiar = k_velocity;
    primitive.pressure_comoving = 1.0;
    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
  }

  const auto geometry = makePeriodic1dGeometry(cell_count);

  const double dx = 1.0 / static_cast<double>(cell_count);
  cosmosim::hydro::HydroUpdateContext update{};
  update.dt_code = 0.2 * dx;
  update.scale_factor = 1.0;
  update.hubble_rate_code = 0.0;

  cosmosim::hydro::HydroSourceContext source_context{};
  source_context.update = update;

  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kVanLeer,
      .dt_over_dx_code = update.dt_code / dx,
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = true});
  cosmosim::hydro::HllcRiemannSolver riemann;

  const std::size_t step_count = static_cast<std::size_t>(std::ceil(k_final_time_code / update.dt_code));
  update.dt_code = k_final_time_code / static_cast<double>(step_count);
  source_context.update = update;

  for (std::size_t step = 0; step < step_count; ++step) {
    solver.advancePatch(conserved, geometry, update, reconstruction, riemann, {}, source_context, nullptr);
  }

  const double shift = std::fmod(k_velocity * k_final_time_code, 1.0);
  double l1 = 0.0;
  for (std::size_t i = 0; i < cell_count; ++i) {
    const double x = (static_cast<double>(i) + 0.5) / static_cast<double>(cell_count);
    const double exact = 1.0 + 0.2 * std::sin(2.0 * k_pi * (x - shift));
    l1 += std::abs(conserved.massDensityComoving()[i] - exact) * geometry.cell_volume_comoving;
  }

  return l1;
}

void testGravityOpeningConvergence(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  cosmosim::gravity::TreeGravityOptions options;
  options.gravitational_constant_code = 1.0;
  options.max_leaf_size = 4;
  options.softening.epsilon_comoving = 1.0e-3;

  options.opening_theta = 0.8;
  const double err_08 = meanTreeRelativeError(options);

  options.opening_theta = 0.6;
  const double err_06 = meanTreeRelativeError(options);

  options.opening_theta = 0.4;
  const double err_04 = meanTreeRelativeError(options);

  requireOrThrow(err_08 <= tolerances.require("gravity_tree_opening.mean_relative_error_theta_0p8"), "tree theta=0.8 tolerance");
  requireOrThrow(err_06 <= tolerances.require("gravity_tree_opening.mean_relative_error_theta_0p6"), "tree theta=0.6 tolerance");
  requireOrThrow(err_04 <= tolerances.require("gravity_tree_opening.mean_relative_error_theta_0p4"), "tree theta=0.4 tolerance");
  requireOrThrow(err_08 >= err_06 && err_06 >= err_04, "tree opening convergence must be monotonic");
}

void testHydroResolutionConvergence(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  const double l1_n64 = runHydroSineWaveL1(64);
  const double l1_n128 = runHydroSineWaveL1(128);

  requireOrThrow(l1_n64 <= tolerances.require("hydro_sine_wave.l1_density_error_n64"), "hydro sine n64 tolerance");
  requireOrThrow(l1_n128 <= tolerances.require("hydro_sine_wave.l1_density_error_n128"), "hydro sine n128 tolerance");

  const double observed_order = std::log(l1_n64 / std::max(l1_n128, 1.0e-16)) / std::log(2.0);
  requireOrThrow(
      observed_order >= tolerances.require("hydro_sine_wave.observed_order_min"),
      "hydro sine observed order below minimum");
}

}  // namespace

int main() {
  const auto tolerances = cosmosim::validation::ValidationToleranceTable::loadFromFile(
      std::string(COSMOSIM_SOURCE_DIR) + "/validation/reference/validation_tolerances_v1.txt");

  testGravityOpeningConvergence(tolerances);
  testHydroResolutionConvergence(tolerances);
  return 0;
}
