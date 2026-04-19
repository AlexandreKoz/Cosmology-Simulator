#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/physics/cooling_heating.hpp"
#include "validation_tolerance.hpp"

namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
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

void fillSodLikeInitialState(cosmosim::hydro::HydroConservedStateSoa& conserved, double gamma) {
  for (std::size_t i = 0; i < conserved.size(); ++i) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    if (i < conserved.size() / 2U) {
      primitive.rho_comoving = 1.0;
      primitive.pressure_comoving = 1.0;
    } else {
      primitive.rho_comoving = 0.125;
      primitive.pressure_comoving = 0.1;
    }
    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }
}

void testPmSingleMode(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  const cosmosim::gravity::PmGridShape shape{32, 8, 8};
  cosmosim::gravity::PmGridStorage grid(shape);
  cosmosim::gravity::PmSolver solver(shape);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 1.0;
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
  std::vector<double> phi(particle_count, 0.0);

  solver.assignDensity(grid, pos_x, pos_y, pos_z, mass, options, nullptr);
  solver.solvePoissonPeriodic(grid, options, nullptr);
  solver.interpolateForces(grid, pos_x, pos_y, pos_z, accel_x, accel_y, accel_z, options, nullptr);
  solver.interpolatePotential(grid, pos_x, pos_y, pos_z, phi, options, nullptr);

  const double kx = 2.0 * k_pi / options.box_size_mpc_comoving;
  const double expected_amp = 4.0 * k_pi * options.gravitational_constant_code *
      options.scale_factor * options.scale_factor * 0.1 / kx;
  const double expected_phi_amp = -4.0 * k_pi * options.gravitational_constant_code *
      options.scale_factor * options.scale_factor * 0.1 / (kx * kx);

  double corr = 0.0;
  double norm_expected = 0.0;
  double norm_got = 0.0;
  double corr_phi = 0.0;
  double norm_phi_expected = 0.0;
  double norm_phi_got = 0.0;
  for (std::size_t i = 0; i < particle_count; ++i) {
    const double expected = expected_amp * std::cos(kx * pos_x[i]);
    const double expected_phi = expected_phi_amp * std::sin(kx * pos_x[i]);
    corr += expected * accel_x[i];
    norm_expected += expected * expected;
    norm_got += accel_x[i] * accel_x[i];
    corr_phi += expected_phi * phi[i];
    norm_phi_expected += expected_phi * expected_phi;
    norm_phi_got += phi[i] * phi[i];
  }

  const double cosine_similarity = corr / std::sqrt(std::max(norm_expected * norm_got, 1.0e-20));
  const double cosine_similarity_phi = corr_phi / std::sqrt(std::max(norm_phi_expected * norm_phi_got, 1.0e-20));
  requireOrThrow(
      std::abs(cosine_similarity) >= tolerances.require("gravity_pm_single_mode.min_cosine_similarity"),
      "gravity_pm_single_mode failed: cosine similarity below tolerance");
  requireOrThrow(
      std::abs(cosine_similarity_phi) >= tolerances.require("gravity_pm_single_mode.min_cosine_similarity"),
      "gravity_pm_single_mode failed: potential cosine similarity below tolerance");
}

void testPmUniformDensityCancellation(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  const cosmosim::gravity::PmGridShape shape{16, 16, 16};
  cosmosim::gravity::PmGridStorage grid(shape);
  cosmosim::gravity::PmSolver solver(shape);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 1.0;
  options.gravitational_constant_code = 1.0;

  auto computeRms = [&](double jitter_amplitude) {
    const std::size_t nside = 8;
    std::vector<double> pos_x;
    std::vector<double> pos_y;
    std::vector<double> pos_z;
    std::vector<double> mass;
    pos_x.reserve(nside * nside * nside);
    pos_y.reserve(nside * nside * nside);
    pos_z.reserve(nside * nside * nside);
    mass.reserve(nside * nside * nside);
    const double cell = 1.0 / static_cast<double>(nside);
    for (std::size_t ix = 0; ix < nside; ++ix) {
      for (std::size_t iy = 0; iy < nside; ++iy) {
        for (std::size_t iz = 0; iz < nside; ++iz) {
          const double base_x = (static_cast<double>(ix) + 0.5) * cell;
          const double base_y = (static_cast<double>(iy) + 0.5) * cell;
          const double base_z = (static_cast<double>(iz) + 0.5) * cell;
          const double phase = static_cast<double>(37U * ix + 57U * iy + 73U * iz + 11U);
          const double jitter_x = jitter_amplitude * cell * std::sin(2.0 * k_pi * (0.61803398875 * phase));
          const double jitter_y = jitter_amplitude * cell * std::sin(2.0 * k_pi * (0.41421356237 * phase));
          const double jitter_z = jitter_amplitude * cell * std::sin(2.0 * k_pi * (0.73205080757 * phase));
          pos_x.push_back(std::fmod(base_x + jitter_x + 1.0, 1.0));
          pos_y.push_back(std::fmod(base_y + jitter_y + 1.0, 1.0));
          pos_z.push_back(std::fmod(base_z + jitter_z + 1.0, 1.0));
          mass.push_back(1.0);
        }
      }
    }

    std::vector<double> ax(pos_x.size(), 0.0);
    std::vector<double> ay(pos_x.size(), 0.0);
    std::vector<double> az(pos_x.size(), 0.0);
    solver.assignDensity(grid, pos_x, pos_y, pos_z, mass, options, nullptr);
    solver.solvePoissonPeriodic(grid, options, nullptr);
    solver.interpolateForces(grid, pos_x, pos_y, pos_z, ax, ay, az, options, nullptr);

    double rms = 0.0;
    for (std::size_t i = 0; i < ax.size(); ++i) {
      rms += ax[i] * ax[i] + ay[i] * ay[i] + az[i] * az[i];
    }
    return std::sqrt(rms / std::max<std::size_t>(ax.size(), 1));
  };

  const double lattice_rms = computeRms(0.0);
  const double glass_like_rms = computeRms(0.18);
  requireOrThrow(
      lattice_rms <= tolerances.require("gravity_pm_uniform_density.max_rms_accel"),
      "gravity_pm_uniform_density failed: PM acceleration non-zero for uniform lattice");
  requireOrThrow(
      glass_like_rms <= tolerances.require("gravity_pm_glass_like.max_rms_accel"),
      "gravity_pm_glass_like failed: PM acceleration too large for deterministic jittered-lattice proxy");
}

struct ForceField {
  std::vector<double> ax;
  std::vector<double> ay;
  std::vector<double> az;
};

[[nodiscard]] double minimumImageDelta(double delta, double box_size_comoving) {
  return delta - box_size_comoving * std::nearbyint(delta / box_size_comoving);
}

ForceField solveTreePm(
    std::size_t pm_grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    double asmth_cells,
    double rcut_cells,
    double theta,
    std::size_t max_leaf) {
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.enable_window_deconvolution = true;
  options.pm_options.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc;
  options.tree_options.opening_theta = theta;
  options.tree_options.max_leaf_size = max_leaf;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 0.01;
  const double mesh_spacing = 1.0 / static_cast<double>(pm_grid);
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(asmth_cells, rcut_cells, mesh_spacing);

  std::vector<std::uint32_t> active(pos_x.size(), 0U);
  for (std::size_t i = 0; i < active.size(); ++i) {
    active[i] = static_cast<std::uint32_t>(i);
  }

  ForceField f{
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0)};
  cosmosim::gravity::TreePmForceAccumulatorView acc{active, f.ax, f.ay, f.az};
  cosmosim::gravity::TreePmCoordinator coordinator({pm_grid, pm_grid, pm_grid});
  coordinator.solveActiveSet(pos_x, pos_y, pos_z, mass, acc, options, nullptr, nullptr);
  return f;
}

ForceField solveFinePmLongRangeReference(
    std::size_t pm_grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    double split_scale_comoving) {
  cosmosim::gravity::PmGridStorage grid({pm_grid, pm_grid, pm_grid});
  cosmosim::gravity::PmSolver solver({pm_grid, pm_grid, pm_grid});
  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 1.0;
  options.gravitational_constant_code = 1.0;
  options.enable_window_deconvolution = true;
  options.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc;
  options.tree_pm_split_scale_comoving = split_scale_comoving;

  ForceField field{
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0)};
  solver.assignDensity(grid, pos_x, pos_y, pos_z, mass, options, nullptr);
  solver.solvePoissonPeriodic(grid, options, nullptr);
  solver.interpolateForces(grid, pos_x, pos_y, pos_z, field.ax, field.ay, field.az, options, nullptr);
  return field;
}

ForceField computeDirectShortRangeResidual(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    const cosmosim::gravity::TreePmOptions& options) {
  ForceField field{
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0)};
  const double cutoff2 = options.split_policy.cutoff_radius_comoving * options.split_policy.cutoff_radius_comoving;
  for (std::size_t target = 0; target < pos_x.size(); ++target) {
    for (std::size_t source = 0; source < pos_x.size(); ++source) {
      if (target == source) {
        continue;
      }
      const double dx = minimumImageDelta(pos_x[source] - pos_x[target], options.pm_options.box_size_mpc_comoving);
      const double dy = minimumImageDelta(pos_y[source] - pos_y[target], options.pm_options.box_size_mpc_comoving);
      const double dz = minimumImageDelta(pos_z[source] - pos_z[target], options.pm_options.box_size_mpc_comoving);
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 > cutoff2) {
        continue;
      }
      const double r = std::sqrt(std::max(r2, 1.0e-30));
      const double split_factor = cosmosim::gravity::treePmGaussianShortRangeForceFactor(
          r, options.split_policy.split_scale_comoving);
      const double softened_factor = cosmosim::gravity::softenedInvR3(r2, options.tree_options.softening) *
          split_factor * options.tree_options.gravitational_constant_code;
      field.ax[target] += softened_factor * mass[source] * dx;
      field.ay[target] += softened_factor * mass[source] * dy;
      field.az[target] += softened_factor * mass[source] * dz;
    }
  }
  return field;
}

ForceField addFields(const ForceField& lhs, const ForceField& rhs) {
  ForceField out{
      std::vector<double>(lhs.ax.size(), 0.0),
      std::vector<double>(lhs.ay.size(), 0.0),
      std::vector<double>(lhs.az.size(), 0.0)};
  for (std::size_t i = 0; i < lhs.ax.size(); ++i) {
    out.ax[i] = lhs.ax[i] + rhs.ax[i];
    out.ay[i] = lhs.ay[i] + rhs.ay[i];
    out.az[i] = lhs.az[i] + rhs.az[i];
  }
  return out;
}

double relativeL2(const ForceField& lhs, const ForceField& rhs) {
  double ref2 = 0.0;
  double err2 = 0.0;
  for (std::size_t i = 0; i < lhs.ax.size(); ++i) {
    const double dx = lhs.ax[i] - rhs.ax[i];
    const double dy = lhs.ay[i] - rhs.ay[i];
    const double dz = lhs.az[i] - rhs.az[i];
    err2 += dx * dx + dy * dy + dz * dz;
    ref2 += rhs.ax[i] * rhs.ax[i] + rhs.ay[i] * rhs.ay[i] + rhs.az[i] * rhs.az[i];
  }
  return std::sqrt(err2 / std::max(ref2, 1.0e-24));
}

void testTreePmPeriodicReferenceAndConsistency(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  constexpr std::size_t particle_count = 24;
  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 1.0);

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod((13.0 * static_cast<double>(i) + 1.0) * 0.047, 1.0);
    pos_y[i] = std::fmod((19.0 * static_cast<double>(i) + 2.0) * 0.031, 1.0);
    pos_z[i] = std::fmod((23.0 * static_cast<double>(i) + 3.0) * 0.029, 1.0);
    mass[i] = 0.9 + 0.05 * static_cast<double>(i % 5U);
  }

  const ForceField tree_only = solveTreePm(32, pos_x, pos_y, pos_z, mass, 100.0, 40.0, 0.55, 8);
  const ForceField pm_only = solveTreePm(32, pos_x, pos_y, pos_z, mass, 0.2, 4.5, 0.55, 8);
  const ForceField split = solveTreePm(32, pos_x, pos_y, pos_z, mass, 1.25, 4.5, 0.55, 8);

  cosmosim::gravity::TreePmOptions ref_options;
  ref_options.pm_options.box_size_mpc_comoving = 1.0;
  ref_options.pm_options.scale_factor = 1.0;
  ref_options.pm_options.gravitational_constant_code = 1.0;
  ref_options.tree_options.gravitational_constant_code = 1.0;
  ref_options.tree_options.softening.epsilon_comoving = 0.01;
  ref_options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 4.5, 1.0 / 32.0);

  const ForceField periodic_spectral_reference = solveFinePmLongRangeReference(
      96, pos_x, pos_y, pos_z, mass, ref_options.split_policy.split_scale_comoving);
  const ForceField direct_short_range_reference = computeDirectShortRangeResidual(pos_x, pos_y, pos_z, mass, ref_options);
  const ForceField periodic_proxy_reference = addFields(periodic_spectral_reference, direct_short_range_reference);

  const double tree_only_rel = relativeL2(tree_only, periodic_proxy_reference);
  const double pm_only_rel = relativeL2(pm_only, periodic_proxy_reference);
  const double split_rel = relativeL2(split, periodic_proxy_reference);

  std::ostringstream msg;
  msg << "tree-pm periodic consistency failure (reference = fine spectral PM + exact pairwise short-range residual; not Ewald exact): tree_only_rel="
      << tree_only_rel << ", split_rel=" << split_rel << ", pm_only_rel=" << pm_only_rel;

  requireOrThrow(tree_only_rel <= tolerances.require("gravity_tree_pm_periodic_proxy.tree_only_rel_l2_max"), msg.str());
  requireOrThrow(pm_only_rel <= tolerances.require("gravity_tree_pm_periodic_proxy.pm_only_rel_l2_max"), msg.str());
  requireOrThrow(split_rel <= tolerances.require("gravity_tree_pm_periodic_proxy.split_rel_l2_max"), msg.str());
  requireOrThrow(split_rel <= pm_only_rel + 1.0e-9, msg.str());
}

void testHydroSodMassConservation(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  constexpr std::size_t k_cell_count = 64;
  constexpr std::size_t k_step_count = 24;
  constexpr double k_gamma = 1.4;

  cosmosim::hydro::HydroConservedStateSoa conserved(k_cell_count);
  fillSodLikeInitialState(conserved, k_gamma);

  const auto geometry = makePeriodic1dGeometry(k_cell_count);

  cosmosim::hydro::HydroUpdateContext update{};
  update.dt_code = 6.0e-4;
  update.scale_factor = 1.0;

  cosmosim::hydro::HydroSourceContext source_context{};
  source_context.update = update;

  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code * static_cast<double>(k_cell_count),
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = true});
  cosmosim::hydro::HllcRiemannSolver riemann;

  double initial_mass = 0.0;
  for (double rho : conserved.massDensityComoving()) {
    initial_mass += rho * geometry.cell_volume_comoving;
  }

  for (std::size_t step = 0; step < k_step_count; ++step) {
    solver.advancePatch(conserved, geometry, update, reconstruction, riemann, {}, source_context, nullptr);
  }

  double final_mass = 0.0;
  for (double rho : conserved.massDensityComoving()) {
    final_mass += rho * geometry.cell_volume_comoving;
  }

  const double mass_error = std::abs(final_mass - initial_mass);
  requireOrThrow(
      mass_error <= tolerances.require("hydro_sod_like.mass_abs_error"),
      "hydro_sod_like failed: mass conservation error above tolerance");
}

void testCoolingEnergyMonotonicity(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  constexpr double k_gamma = 5.0 / 3.0;

  cosmosim::hydro::HydroPrimitiveState primitive{};
  primitive.rho_comoving = 1.0;
  primitive.pressure_comoving = 1.0;
  auto conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma);

  cosmosim::physics::CoolingModelConfig model{};
  model.temperature_floor_k = 100.0;
  model.max_subcycles = 32;
  model.uv_background_model = cosmosim::core::UvBackgroundModel::kNone;

  cosmosim::physics::CoolingRateProvider provider(model);
  cosmosim::physics::CoolingSourceIntegrator integrator(1.0e-8);
  cosmosim::physics::CoolingHeatingSource cooling(provider, integrator);

  std::vector<double> nh(1, 5.0);
  std::vector<double> metal(1, 0.02);
  std::vector<double> temp(1, 2.0e6);
  cosmosim::hydro::HydroSourceContext source_context{};
  source_context.update.dt_code = 0.005;
  source_context.update.scale_factor = 1.0;
  source_context.hydrogen_number_density_cgs = nh;
  source_context.metallicity_mass_fraction = metal;
  source_context.temperature_k = temp;

  const double initial_total_energy = conserved.total_energy_density_comoving;
  for (int step = 0; step < 16; ++step) {
    const auto source = cooling.sourceForCell(0, conserved, primitive, source_context);
    conserved.total_energy_density_comoving += source.total_energy_density_comoving * source_context.update.dt_code;
  }

  const double allowed_increase = tolerances.require("cooling_box.total_energy_nonincrease_eps");
  requireOrThrow(
      conserved.total_energy_density_comoving <= initial_total_energy + allowed_increase,
      "cooling_box failed: thermal energy increased beyond allowed epsilon");
}

}  // namespace

int main() {
  const auto tolerances = cosmosim::validation::ValidationToleranceTable::loadFromFile(
      std::string(COSMOSIM_SOURCE_DIR) + "/validation/reference/validation_tolerances_v1.txt");

  testPmSingleMode(tolerances);
  testPmUniformDensityCancellation(tolerances);
  testTreePmPeriodicReferenceAndConsistency(tolerances);
  testHydroSodMassConservation(tolerances);
  testCoolingEnergyMonotonicity(tolerances);
  return 0;
}
