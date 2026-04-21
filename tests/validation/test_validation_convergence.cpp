#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "validation_tolerance.hpp"

namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

struct TwoBodyEnergyDriftResult {
  double dt_code = 0.0;
  int bin_index = 0;
  double relative_energy_drift = 0.0;
};

struct HaloProfileRow {
  double radius_code = 0.0;
  double tree_radial_force_code = 0.0;
  double reference_radial_force_code = 0.0;
  double tree_potential_code_estimate = 0.0;
  double reference_potential_code = 0.0;
  double radial_force_relative_error = 0.0;
};

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

double softenedDirectPotentialAt(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    std::uint32_t active_index,
    const cosmosim::gravity::TreeGravityOptions& options) {
  const double tx = pos_x[active_index];
  const double ty = pos_y[active_index];
  const double tz = pos_z[active_index];
  double phi = 0.0;
  const double eps2 = options.softening.epsilon_comoving * options.softening.epsilon_comoving;
  for (std::size_t j = 0; j < mass.size(); ++j) {
    if (j == active_index) {
      continue;
    }
    const double dx = pos_x[j] - tx;
    const double dy = pos_y[j] - ty;
    const double dz = pos_z[j] - tz;
    const double r2 = dx * dx + dy * dy + dz * dz;
    phi += -options.gravitational_constant_code * mass[j] / std::sqrt(r2 + eps2);
  }
  return phi;
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

TwoBodyEnergyDriftResult runTwoBodyOrbitEnergyDrift(double dt, int bin_index) {
  constexpr std::size_t steps = 4000;
  constexpr double g = 1.0;
  constexpr double m0 = 1.0;
  constexpr double m1 = 1.0;

  cosmosim::gravity::TreeGravityOptions options;
  options.opening_theta = 0.2;
  options.max_leaf_size = 1;
  options.gravitational_constant_code = g;
  options.softening.epsilon_comoving = 1.0e-3;

  std::vector<double> px = {-0.5, 0.5};
  std::vector<double> py = {0.0, 0.0};
  std::vector<double> pz = {0.0, 0.0};
  const double v = std::sqrt(g * (m0 + m1) / 1.0) * 0.5;
  std::vector<double> vx = {0.0, 0.0};
  std::vector<double> vy = {v, -v};
  std::vector<double> vz = {0.0, 0.0};
  std::vector<double> mass = {m0, m1};
  std::vector<std::uint32_t> active = {0U, 1U};

  auto totalEnergy = [&]() {
    const double dx = px[1] - px[0];
    const double dy = py[1] - py[0];
    const double dz = pz[1] - pz[0];
    const double r2 = dx * dx + dy * dy + dz * dz;
    const double kinetic = 0.5 * m0 * (vx[0] * vx[0] + vy[0] * vy[0] + vz[0] * vz[0]) +
        0.5 * m1 * (vx[1] * vx[1] + vy[1] * vy[1] + vz[1] * vz[1]);
    const double potential = -g * m0 * m1 / std::sqrt(r2 + options.softening.epsilon_comoving * options.softening.epsilon_comoving);
    return kinetic + potential;
  };

  std::vector<double> ax(2, 0.0), ay(2, 0.0), az(2, 0.0);
  cosmosim::gravity::TreeGravitySolver solver;
  solver.build(px, py, pz, mass, options, nullptr);
  solver.evaluateActiveSet(px, py, pz, mass, active, ax, ay, az, options, nullptr);

  const double e0 = totalEnergy();
  for (std::size_t step = 0; step < steps; ++step) {
    for (std::size_t i = 0; i < 2; ++i) {
      vx[i] += 0.5 * dt * ax[i];
      vy[i] += 0.5 * dt * ay[i];
      vz[i] += 0.5 * dt * az[i];
      px[i] += dt * vx[i];
      py[i] += dt * vy[i];
      pz[i] += dt * vz[i];
    }
    solver.build(px, py, pz, mass, options, nullptr);
    solver.evaluateActiveSet(px, py, pz, mass, active, ax, ay, az, options, nullptr);
    for (std::size_t i = 0; i < 2; ++i) {
      vx[i] += 0.5 * dt * ax[i];
      vy[i] += 0.5 * dt * ay[i];
      vz[i] += 0.5 * dt * az[i];
    }
  }

  const double e1 = totalEnergy();
  const double rel_drift = std::abs(e1 - e0) / std::max(std::abs(e0), 1.0e-20);
  return TwoBodyEnergyDriftResult{.dt_code = dt, .bin_index = bin_index, .relative_energy_drift = rel_drift};
}

std::vector<TwoBodyEnergyDriftResult> testTwoBodyOrbitEnergyDrift(
    const cosmosim::validation::ValidationToleranceTable& tolerances) {
  const auto baseline = runTwoBodyOrbitEnergyDrift(1.0e-3, -1);
  requireOrThrow(
      baseline.relative_energy_drift <= tolerances.require("gravity_two_body_orbit.max_relative_energy_drift"),
      "gravity_two_body_orbit failed: relative energy drift above tolerance");

  const cosmosim::core::TimeStepLimits limits{
      .min_dt_time_code = 1.0e-3,
      .max_dt_time_code = 8.0e-3,
      .max_bin = 3,
  };

  std::vector<TwoBodyEnergyDriftResult> rows;
  rows.push_back(baseline);
  for (int bin = 0; bin <= 2; ++bin) {
    const double dt = cosmosim::core::binIndexToDt(bin, limits);
    rows.push_back(runTwoBodyOrbitEnergyDrift(dt, bin));
  }
  return rows;
}

std::vector<HaloProfileRow> testStaticHaloRadialProfile(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  constexpr std::size_t shell_count = 6;
  constexpr std::size_t particle_count = 96;
  constexpr double g = 1.0;

  std::vector<double> px(particle_count, 0.0);
  std::vector<double> py(particle_count, 0.0);
  std::vector<double> pz(particle_count, 0.0);
  std::vector<double> mass(particle_count, 1.0 / static_cast<double>(particle_count));

  for (std::size_t i = 0; i < particle_count; ++i) {
    const double u = (static_cast<double>(i) + 0.5) / static_cast<double>(particle_count);
    const double v = std::fmod(0.61803398875 * static_cast<double>(i), 1.0);
    const double w = std::fmod(0.41421356237 * static_cast<double>(i), 1.0);
    const double r = 0.05 + 0.45 * std::cbrt(u);
    const double cos_t = 1.0 - 2.0 * v;
    const double sin_t = std::sqrt(std::max(0.0, 1.0 - cos_t * cos_t));
    const double phi = 2.0 * k_pi * w;
    px[i] = r * sin_t * std::cos(phi);
    py[i] = r * sin_t * std::sin(phi);
    pz[i] = r * cos_t;
  }

  std::vector<double> test_x(shell_count, 0.0);
  std::vector<double> test_y(shell_count, 0.0);
  std::vector<double> test_z(shell_count, 0.0);
  std::vector<double> test_m(shell_count, 1.0e-6);
  std::vector<std::uint32_t> active(shell_count, 0U);
  for (std::size_t i = 0; i < shell_count; ++i) {
    const double r = 0.08 + 0.06 * static_cast<double>(i);
    test_x[i] = r;
    active[i] = static_cast<std::uint32_t>(i);
  }

  std::vector<double> all_x = test_x;
  std::vector<double> all_y = test_y;
  std::vector<double> all_z = test_z;
  std::vector<double> all_m = test_m;
  all_x.insert(all_x.end(), px.begin(), px.end());
  all_y.insert(all_y.end(), py.begin(), py.end());
  all_z.insert(all_z.end(), pz.begin(), pz.end());
  all_m.insert(all_m.end(), mass.begin(), mass.end());

  cosmosim::gravity::TreeGravityOptions options;
  options.opening_theta = 0.35;
  options.max_leaf_size = 4;
  options.gravitational_constant_code = g;
  options.softening.epsilon_comoving = 5.0e-3;

  cosmosim::gravity::TreeGravitySolver solver;
  solver.build(all_x, all_y, all_z, all_m, options, nullptr);
  std::vector<double> ax(shell_count, 0.0), ay(shell_count, 0.0), az(shell_count, 0.0);
  solver.evaluateActiveSet(all_x, all_y, all_z, all_m, active, ax, ay, az, options, nullptr);

  std::vector<double> ref_ax(shell_count, 0.0), ref_ay(shell_count, 0.0), ref_az(shell_count, 0.0);
  directSumAcceleration(all_x, all_y, all_z, all_m, active, options, ref_ax, ref_ay, ref_az);

  std::vector<HaloProfileRow> rows;
  rows.reserve(shell_count);
  double max_radial_rel = 0.0;
  for (std::size_t i = 0; i < shell_count; ++i) {
    const double ar_tree = ax[i];
    const double ar_ref = ref_ax[i];
    const double rel = std::abs(ar_tree - ar_ref) / std::max(std::abs(ar_ref), 1.0e-12);
    max_radial_rel = std::max(max_radial_rel, rel);
    rows.push_back(HaloProfileRow{
        .radius_code = test_x[i],
        .tree_radial_force_code = ar_tree,
        .reference_radial_force_code = ar_ref,
        .tree_potential_code_estimate = 0.0,
        .reference_potential_code = softenedDirectPotentialAt(all_x, all_y, all_z, all_m, active[i], options),
        .radial_force_relative_error = rel,
    });
  }

  for (std::size_t i = rows.size(); i > 1; --i) {
    const std::size_t inner = i - 2;
    const std::size_t outer = i - 1;
    const double dr = rows[outer].radius_code - rows[inner].radius_code;
    rows[inner].tree_potential_code_estimate =
        rows[outer].tree_potential_code_estimate + 0.5 * (rows[outer].tree_radial_force_code + rows[inner].tree_radial_force_code) * dr;
  }

  requireOrThrow(
      max_radial_rel <= tolerances.require("gravity_static_halo.max_radial_relative_force_error"),
      "gravity_static_halo failed: radial profile mismatch vs softened direct-sum reference");
  return rows;
}

void writePhase3Artifacts(
    std::span<const TwoBodyEnergyDriftResult> drift_rows,
    std::span<const HaloProfileRow> halo_rows,
    const cosmosim::validation::ValidationToleranceTable& tolerances) {
  const std::filesystem::path root =
      std::filesystem::path(COSMOSIM_SOURCE_DIR) / "validation" / "artifacts" / "research_grade" / "phase3";
  const std::filesystem::path force_dir = root / "force_accuracy";
  const std::filesystem::path time_dir = root / "time_integration";
  std::filesystem::create_directories(force_dir);
  std::filesystem::create_directories(time_dir);

  {
    std::ofstream out(force_dir / "halo_force_potential_profile.csv");
    out << std::setprecision(12);
    out << "observable,reference_target,tolerance_envelope,radius_code,tree_force_code,reference_force_code,tree_potential_estimate_code,reference_potential_code,force_relative_error\n";
    for (const auto& row : halo_rows) {
      out << "halo_radial_profile,softened_direct_sum,"
          << tolerances.require("gravity_static_halo.max_radial_relative_force_error")
          << ',' << row.radius_code
          << ',' << row.tree_radial_force_code
          << ',' << row.reference_radial_force_code
          << ',' << row.tree_potential_code_estimate
          << ',' << row.reference_potential_code
          << ',' << row.radial_force_relative_error
          << '\n';
    }
  }

  {
    std::ofstream out(time_dir / "hierarchical_time_integration_accuracy.csv");
    out << std::setprecision(12);
    out << "observable,reference_target,tolerance_envelope,time_bin_index,dt_code,relative_energy_drift\n";
    for (const auto& row : drift_rows) {
      out << "two_body_energy_drift,leapfrog_softened_orbit,"
          << tolerances.require("gravity_two_body_orbit.max_relative_energy_drift")
          << ',' << row.bin_index
          << ',' << row.dt_code
          << ',' << row.relative_energy_drift
          << '\n';
    }
  }
}

}  // namespace

int main() {
  const auto tolerances = cosmosim::validation::ValidationToleranceTable::loadFromFile(
      std::string(COSMOSIM_SOURCE_DIR) + "/validation/reference/validation_tolerances_v1.txt");

  testGravityOpeningConvergence(tolerances);
  testHydroResolutionConvergence(tolerances);
  const auto drift_rows = testTwoBodyOrbitEnergyDrift(tolerances);
  const auto halo_rows = testStaticHaloRadialProfile(tolerances);
  writePhase3Artifacts(drift_rows, halo_rows, tolerances);
  return 0;
}
