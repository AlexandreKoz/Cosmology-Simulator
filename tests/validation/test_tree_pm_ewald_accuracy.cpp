#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "periodic_ewald_reference.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

constexpr double k_required_relative_l2 = 1.0e-2;
constexpr double k_required_p99_normalized = 5.0e-2;
constexpr double k_normalized_error_floor_fraction = 1.0e-3;
constexpr double k_required_net_force_fraction = 1.0e-12;
constexpr double k_required_periodic_translation_drift = 1.0e-11;
constexpr double k_asmth_cells = 1.25;
#if COSMOSIM_ENABLE_FFTW
constexpr double k_rcut_cells = 6.25;
#else
// The naive-DFT meshes are intentionally small. Their non-certifying profile
// uses a bounded cutoff that remains strictly inside half the shortest box
// axis; FFTW coverage above retains the normalized production value.
constexpr double k_rcut_cells = 3.9;
#endif

// The FFTW-backed TSC path is the accuracy-certified profile. CIC is still
// evaluated against the identical targets, but remains diagnostic because the
// documented split and mesh choices exceed the targets on the uniform fixture.
// The reduced naive-DFT meshes are diagnostics too;
// they are deliberately small enough for non-FFTW developer builds.

struct ForceField {
  std::vector<double> ax;
  std::vector<double> ay;
  std::vector<double> az;
};

struct ParticleFixture {
  std::string label;
  std::string topology_label = "interior";
  std::string separation_bin = "aggregate";
  std::string mesh_direction_label = "mixed";
  cosmosim::test_support::PeriodicEwaldBox box{};
  cosmosim::gravity::PmGridShape pm_shape{};
  std::vector<double> pos_x;
  std::vector<double> pos_y;
  std::vector<double> pos_z;
  std::vector<double> mass;
  double classified_separation = std::numeric_limits<double>::quiet_NaN();
  double classified_split_scale = std::numeric_limits<double>::quiet_NaN();
  double classified_cutoff = std::numeric_limits<double>::quiet_NaN();
};

struct IndependentTargetFixture {
  ParticleFixture sources{};
  std::vector<double> target_pos_x;
  std::vector<double> target_pos_y;
  std::vector<double> target_pos_z;
};

struct SolverProfile {
  std::string label;
  cosmosim::gravity::PmAssignmentScheme assignment_scheme =
      cosmosim::gravity::PmAssignmentScheme::kTsc;
  bool enable_window_deconvolution = true;
  cosmosim::gravity::TreeMultipoleOrder multipole_order =
      cosmosim::gravity::TreeMultipoleOrder::kQuadrupole;
  cosmosim::gravity::TreeOpeningCriterion opening_criterion =
      cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance;
  double opening_theta = 0.7;
  double relative_force_tolerance = 0.005;
  double relative_force_acceleration_floor_code = 1.0e-30;
  std::size_t max_leaf_size = 16;
  bool supply_converged_reference_force_scale = false;
  bool accuracy_certified = false;
};

struct ErrorMetrics {
  double absolute_rms = 0.0;
  double relative_l2 = 0.0;
  double normalized_floor = 0.0;
  double median_normalized = 0.0;
  double p90_normalized = 0.0;
  double p95_normalized = 0.0;
  double p99_normalized = 0.0;
  double max_normalized = 0.0;
};

struct MatchedMacResult {
  SolverProfile solver_profile{};
  ErrorMetrics error{};
  cosmosim::gravity::TreePmProfileEvent profile{};
};

struct AccuracyResult {
  std::string fixture_label;
  std::string topology_label;
  std::string separation_bin;
  std::string mesh_direction_label;
  SolverProfile solver_profile{};
  ErrorMetrics error{};
  double reference_force_l2 = 0.0;
  double net_force_fraction = 0.0;
  double reference_net_force_fraction = 0.0;
  double translation_relative_l2 = std::numeric_limits<double>::quiet_NaN();
  double classified_separation = std::numeric_limits<double>::quiet_NaN();
  double classified_split_scale = std::numeric_limits<double>::quiet_NaN();
  double classified_cutoff = std::numeric_limits<double>::quiet_NaN();
  cosmosim::gravity::PmGridShape pm_shape{};
  cosmosim::gravity::TreePmDiagnostics diagnostics{};
  cosmosim::gravity::TreePmProfileEvent profile{};
};

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

[[nodiscard]] double squaredNorm(double x, double y, double z) {
  return x * x + y * y + z * z;
}

[[nodiscard]] double vectorNorm(double x, double y, double z) {
  return std::sqrt(squaredNorm(x, y, z));
}

[[nodiscard]] double forceL2Norm(const ForceField& force) {
  double sum_squared = 0.0;
  for (std::size_t i = 0; i < force.ax.size(); ++i) {
    sum_squared += squaredNorm(force.ax[i], force.ay[i], force.az[i]);
  }
  return std::sqrt(sum_squared);
}

[[nodiscard]] double wrapPeriodic(double coordinate, double length) {
  double wrapped = std::fmod(coordinate, length);
  if (wrapped < 0.0) {
    wrapped += length;
  }
  return wrapped == length ? 0.0 : wrapped;
}

[[nodiscard]] cosmosim::gravity::PmGridShape cubicPmShape() {
#if COSMOSIM_ENABLE_FFTW
  return {32U, 32U, 32U};
#else
  // The fallback is the same spectral operator evaluated by a naive DFT. A
  // 12^3 mesh keeps this validation runnable without treating that backend as
  // a production-performance path.
  return {12U, 12U, 12U};
#endif
}

[[nodiscard]] cosmosim::gravity::PmGridShape rectangularPmShape() {
#if COSMOSIM_ENABLE_FFTW
  return {32U, 40U, 48U};
#else
  return {8U, 10U, 12U};
#endif
}

void validateFixture(const ParticleFixture& fixture) {
  const std::size_t count = fixture.mass.size();
  requireOrThrow(count > 0, fixture.label + ": fixture must contain sources");
  requireOrThrow(
      fixture.pos_x.size() == count && fixture.pos_y.size() == count && fixture.pos_z.size() == count,
      fixture.label + ": position/mass extents differ");
  requireOrThrow(fixture.pm_shape.isValid(), fixture.label + ": PM shape is invalid");
  requireOrThrow(
      fixture.box.length_x > 0.0 && fixture.box.length_y > 0.0 && fixture.box.length_z > 0.0,
      fixture.label + ": box lengths must be positive");
  for (std::size_t i = 0; i < count; ++i) {
    requireOrThrow(
        std::isfinite(fixture.pos_x[i]) && std::isfinite(fixture.pos_y[i]) &&
            std::isfinite(fixture.pos_z[i]) && std::isfinite(fixture.mass[i]) && fixture.mass[i] > 0.0,
        fixture.label + ": fixture contains non-finite state or non-positive mass");
  }
  if (std::isfinite(fixture.classified_separation)) {
    requireOrThrow(
        fixture.classified_separation > 0.0 &&
            std::isfinite(fixture.classified_split_scale) && fixture.classified_split_scale > 0.0 &&
            std::isfinite(fixture.classified_cutoff) && fixture.classified_cutoff > 0.0,
        fixture.label + ": classified separation metadata is invalid");
  }
}

[[nodiscard]] ParticleFixture makeQuasiUniformFixture() {
  ParticleFixture fixture;
  fixture.label = "quasi_uniform_cubic";
  fixture.box = {.length_x = 1.0, .length_y = 1.0, .length_z = 1.0};
  fixture.pm_shape = cubicPmShape();
  constexpr std::size_t particle_count = 24;
  fixture.pos_x.resize(particle_count);
  fixture.pos_y.resize(particle_count);
  fixture.pos_z.resize(particle_count);
  fixture.mass.resize(particle_count);
  for (std::size_t i = 0; i < particle_count; ++i) {
    const double x_fraction = (static_cast<double>((17U * i + 5U) % 97U) + 0.5) / 97.0;
    const double y_fraction = (static_cast<double>((37U * i + 11U) % 101U) + 0.5) / 101.0;
    const double z_fraction = (static_cast<double>((53U * i + 19U) % 103U) + 0.5) / 103.0;
    fixture.pos_x[i] = x_fraction;
    fixture.pos_y[i] = y_fraction;
    fixture.pos_z[i] = z_fraction;
    fixture.mass[i] = 0.8 + 0.04 * static_cast<double>((7U * i + 3U) % 7U);
  }
  return fixture;
}

[[nodiscard]] ParticleFixture makeDmoZeldovichLatticeFixture() {
  ParticleFixture fixture;
  fixture.label = "dmo_zeldovich_4cubed_initial";
  fixture.topology_label = "cancellation_dominated_lattice";
  fixture.box = {.length_x = 1.0, .length_y = 1.0, .length_z = 1.0};
  fixture.pm_shape = {16U, 16U, 16U};

  constexpr std::size_t lattice_n = 4U;
  constexpr std::size_t particle_count = lattice_n * lattice_n * lattice_n;
  constexpr double linear_density_amplitude = 5.0e-3;
  constexpr double hubble0_code = 67.4;
  constexpr double omega_matter = 0.315;
  constexpr double pi = 3.141592653589793238462643383279502884;
  constexpr double wave_number = 2.0 * pi;
  constexpr double displacement_amplitude = linear_density_amplitude / wave_number;

  // The workflow derives M_box=3 H0^2 Omega_m V/(8 pi G).  Store G*M_i
  // here and solve with G=1 in both TreePM and Ewald, factoring out the
  // physical code-unit G without changing the production force scale.
  constexpr double effective_particle_mass =
      3.0 * hubble0_code * hubble0_code * omega_matter /
      (8.0 * pi * static_cast<double>(particle_count));

  fixture.pos_x.reserve(particle_count);
  fixture.pos_y.reserve(particle_count);
  fixture.pos_z.reserve(particle_count);
  fixture.mass.assign(particle_count, effective_particle_mass);
  constexpr double spacing = 1.0 / static_cast<double>(lattice_n);
  for (std::size_t ix = 0; ix < lattice_n; ++ix) {
    for (std::size_t iy = 0; iy < lattice_n; ++iy) {
      for (std::size_t iz = 0; iz < lattice_n; ++iz) {
        const double qx = (static_cast<double>(ix) + 0.5) * spacing;
        fixture.pos_x.push_back(wrapPeriodic(
            qx + displacement_amplitude * std::sin(wave_number * qx), 1.0));
        fixture.pos_y.push_back((static_cast<double>(iy) + 0.5) * spacing);
        fixture.pos_z.push_back((static_cast<double>(iz) + 0.5) * spacing);
      }
    }
  }
  return fixture;
}

[[nodiscard]] ParticleFixture makeClusteredFixture() {
  ParticleFixture fixture;
  fixture.label = "two_cluster_cubic";
  fixture.box = {.length_x = 1.0, .length_y = 1.0, .length_z = 1.0};
  fixture.pm_shape = cubicPmShape();
  constexpr std::size_t particle_count = 20;
  fixture.pos_x.resize(particle_count);
  fixture.pos_y.resize(particle_count);
  fixture.pos_z.resize(particle_count);
  fixture.mass.resize(particle_count);
  for (std::size_t i = 0; i < particle_count; ++i) {
    const bool lower_cluster = i < particle_count / 2U;
    const double phase = static_cast<double>(i % (particle_count / 2U));
    const double center_x = lower_cluster ? 0.27 : 0.73;
    const double center_y = lower_cluster ? 0.31 : 0.69;
    const double center_z = lower_cluster ? 0.24 : 0.76;
    fixture.pos_x[i] = center_x + 0.043 * std::sin(0.73 * phase + 0.17);
    fixture.pos_y[i] = center_y + 0.037 * std::cos(0.61 * phase + 0.29);
    fixture.pos_z[i] = center_z + 0.041 * std::sin(0.47 * phase + 0.83);
    fixture.mass[i] = 0.65 + 0.05 * static_cast<double>((5U * i + 1U) % 6U);
  }
  return fixture;
}

[[nodiscard]] ParticleFixture makeMacComparisonFixture() {
  ParticleFixture fixture;
  fixture.label = "dense_irregular_mac_comparison";
  fixture.box = {.length_x = 1.0, .length_y = 1.0, .length_z = 1.0};
  fixture.pm_shape = cubicPmShape();
  constexpr std::size_t particle_count = 48;
  fixture.pos_x.resize(particle_count);
  fixture.pos_y.resize(particle_count);
  fixture.pos_z.resize(particle_count);
  fixture.mass.resize(particle_count);
  for (std::size_t i = 0; i < particle_count; ++i) {
    const double x_fraction = (static_cast<double>((17U * i + 5U) % 53U) + 0.5) / 53.0 - 0.5;
    const double y_fraction = (static_cast<double>((29U * i + 7U) % 59U) + 0.5) / 59.0 - 0.5;
    const double z_fraction = (static_cast<double>((43U * i + 11U) % 61U) + 0.5) / 61.0 - 0.5;
    fixture.pos_x[i] = 0.5 + 0.09 * x_fraction;
    fixture.pos_y[i] = 0.5 + 0.09 * y_fraction;
    fixture.pos_z[i] = 0.5 + 0.09 * z_fraction;
    fixture.mass[i] = 0.7 + 0.025 * static_cast<double>((13U * i + 3U) % 11U);
  }
  return fixture;
}

[[nodiscard]] IndependentTargetFixture makeMatchedMacFixture() {
  IndependentTargetFixture fixture;
  fixture.sources.label = "compact_sources_independent_targets";
  fixture.sources.box = {.length_x = 1.0, .length_y = 1.0, .length_z = 1.0};
  fixture.sources.pm_shape = cubicPmShape();
  constexpr std::size_t source_count = 64U;
  fixture.sources.pos_x.resize(source_count);
  fixture.sources.pos_y.resize(source_count);
  fixture.sources.pos_z.resize(source_count);
  fixture.sources.mass.resize(source_count);
  for (std::size_t i = 0; i < source_count; ++i) {
    const double dx = (static_cast<double>((17U * i + 3U) % 67U) / 66.0 - 0.5) * 0.012;
    const double dy = (static_cast<double>((29U * i + 7U) % 71U) / 70.0 - 0.5) * 0.010;
    const double dz = (static_cast<double>((43U * i + 11U) % 73U) / 72.0 - 0.5) * 0.011;
    fixture.sources.pos_x[i] = 0.34 + dx;
    fixture.sources.pos_y[i] = 0.50 + dy;
    fixture.sources.pos_z[i] = 0.50 + dz;
    fixture.sources.mass[i] = 0.75 + 0.02 * static_cast<double>((11U * i + 5U) % 13U);
  }

  constexpr std::size_t target_count = 16U;
  fixture.target_pos_x.resize(target_count);
  fixture.target_pos_y.resize(target_count);
  fixture.target_pos_z.resize(target_count);
  for (std::size_t i = 0; i < target_count; ++i) {
    fixture.target_pos_x[i] = 0.505 + 0.00125 * static_cast<double>(i);
    fixture.target_pos_y[i] = 0.497 + 0.0004 * static_cast<double>((5U * i + 1U) % 13U);
    fixture.target_pos_z[i] = 0.498 + 0.00035 * static_cast<double>((7U * i + 2U) % 11U);
  }
  return fixture;
}

[[nodiscard]] ParticleFixture makeSeamFixture() {
  ParticleFixture fixture;
  fixture.label = "xyz_seam_rectangular";
  fixture.topology_label = "periodic_seam_xyz";
  fixture.box = {.length_x = 1.0, .length_y = 1.25, .length_z = 1.5};
  fixture.pm_shape = rectangularPmShape();
  constexpr std::size_t particle_count = 18;
  fixture.pos_x.resize(particle_count);
  fixture.pos_y.resize(particle_count);
  fixture.pos_z.resize(particle_count);
  fixture.mass.resize(particle_count);
  for (std::size_t i = 0; i < particle_count; ++i) {
    const double phase = static_cast<double>(i);
    const double local_x = 0.046 * std::sin(0.71 * phase + 0.11);
    const double local_y = 0.052 * std::cos(0.59 * phase + 0.37);
    const double local_z = 0.057 * std::sin(0.43 * phase + 0.91);
    fixture.pos_x[i] = wrapPeriodic(local_x, fixture.box.length_x);
    fixture.pos_y[i] = wrapPeriodic(local_y, fixture.box.length_y);
    fixture.pos_z[i] = wrapPeriodic(local_z, fixture.box.length_z);
    fixture.mass[i] = 0.7 + 0.03 * static_cast<double>((11U * i + 2U) % 9U);
  }
  return fixture;
}

[[nodiscard]] ParticleFixture makeTwoBodySeparationFixture(
    std::string label,
    std::string separation_bin,
    std::array<double, 3> direction,
    std::string mesh_direction_label,
    double separation) {
  ParticleFixture fixture;
  fixture.label = std::move(label);
  fixture.separation_bin = std::move(separation_bin);
  fixture.mesh_direction_label = std::move(mesh_direction_label);
  fixture.box = {.length_x = 1.0, .length_y = 1.0, .length_z = 1.0};
  fixture.pm_shape = cubicPmShape();
  const double mesh_spacing = 1.0 / static_cast<double>(fixture.pm_shape.nx);
  fixture.classified_separation = separation;
  fixture.classified_split_scale = k_asmth_cells * mesh_spacing;
  fixture.classified_cutoff = k_rcut_cells * mesh_spacing;
  const double direction_norm = vectorNorm(direction[0], direction[1], direction[2]);
  requireOrThrow(direction_norm > 0.0, fixture.label + ": separation direction must be nonzero");
  for (double& component : direction) {
    component /= direction_norm;
  }
  fixture.pos_x = {0.5 - 0.5 * separation * direction[0], 0.5 + 0.5 * separation * direction[0]};
  fixture.pos_y = {0.5 - 0.5 * separation * direction[1], 0.5 + 0.5 * separation * direction[1]};
  fixture.pos_z = {0.5 - 0.5 * separation * direction[2], 0.5 + 0.5 * separation * direction[2]};
  fixture.mass = {1.0, 1.0};
  return fixture;
}

[[nodiscard]] ParticleFixture makeTwoBodySplitFixture() {
  const cosmosim::gravity::PmGridShape pm_shape = cubicPmShape();
  const double mesh_spacing = 1.0 / static_cast<double>(pm_shape.nx);
  const double split_scale = k_asmth_cells * mesh_spacing;
  return makeTwoBodySeparationFixture(
      "two_body_at_split_transition",
      "gaussian_transition_q_1",
      {0.7427813527082074, 0.5570860145311556, 0.3713906763541037},
      "oblique",
      2.0 * split_scale);
}

[[nodiscard]] std::vector<ParticleFixture> makeSeparationClassificationFixtures() {
  const cosmosim::gravity::PmGridShape pm_shape = cubicPmShape();
  const double mesh_spacing = 1.0 / static_cast<double>(pm_shape.nx);
  const double split_scale = k_asmth_cells * mesh_spacing;
  const double cutoff = k_rcut_cells * mesh_spacing;
  constexpr std::array<double, 3> mesh_axis_x{1.0, 0.0, 0.0};
  constexpr double inverse_sqrt_three = 0.577350269189625764509148780501957456;
  constexpr std::array<double, 3> body_diagonal{
      inverse_sqrt_three,
      inverse_sqrt_three,
      inverse_sqrt_three,
  };

  // q = r/(2*r_s) is the dimensionless argument of the Gaussian split
  // kernel. The compact set brackets q=1 and the hard short-range cutoff in
  // both a mesh-aligned and an oblique direction without multiplying every
  // other diagnostic dimension into this already-expensive FFT/Ewald test.
  std::vector<ParticleFixture> fixtures;
  fixtures.push_back(makeTwoBodySeparationFixture(
      "pair_split_q_0p8_axis_x", "below_gaussian_transition", mesh_axis_x, "axis_x", 1.6 * split_scale));
  fixtures.push_back(makeTwoBodySeparationFixture(
      "pair_split_q_1_axis_x", "at_gaussian_transition", mesh_axis_x, "axis_x", 2.0 * split_scale));
  fixtures.push_back(makeTwoBodySeparationFixture(
      "pair_split_q_1_oblique", "at_gaussian_transition", body_diagonal, "body_diagonal", 2.0 * split_scale));
  fixtures.push_back(makeTwoBodySeparationFixture(
      "pair_split_q_1p2_oblique", "above_gaussian_transition", body_diagonal, "body_diagonal", 2.4 * split_scale));
  fixtures.push_back(makeTwoBodySeparationFixture(
      "pair_cutoff_0p9_axis_x", "below_short_range_cutoff", mesh_axis_x, "axis_x", 0.9 * cutoff));
  fixtures.push_back(makeTwoBodySeparationFixture(
      "pair_cutoff_0p999_axis_x", "immediately_below_short_range_cutoff", mesh_axis_x, "axis_x", 0.999 * cutoff));
  fixtures.push_back(makeTwoBodySeparationFixture(
      "pair_cutoff_1_axis_x", "at_short_range_cutoff", mesh_axis_x, "axis_x", cutoff));
  fixtures.push_back(makeTwoBodySeparationFixture(
      "pair_cutoff_1_oblique", "at_short_range_cutoff", body_diagonal, "body_diagonal", cutoff));
  fixtures.push_back(makeTwoBodySeparationFixture(
      "pair_cutoff_1p001_axis_x", "immediately_above_short_range_cutoff", mesh_axis_x, "axis_x", 1.001 * cutoff));
  fixtures.push_back(makeTwoBodySeparationFixture(
      "pair_cutoff_1p1_oblique", "above_short_range_cutoff", body_diagonal, "body_diagonal", 1.1 * cutoff));
  return fixtures;
}

[[nodiscard]] cosmosim::gravity::TreePmOptions makeProductionOptions(
    const ParticleFixture& fixture,
    const SolverProfile& solver_profile) {
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_x_mpc_comoving = fixture.box.length_x;
  options.pm_options.box_size_y_mpc_comoving = fixture.box.length_y;
  options.pm_options.box_size_z_mpc_comoving = fixture.box.length_z;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.assignment_scheme = solver_profile.assignment_scheme;
  options.pm_options.enable_window_deconvolution = solver_profile.enable_window_deconvolution;
  options.pm_options.boundary_condition = cosmosim::gravity::PmBoundaryCondition::kPeriodic;

  // Exercise the documented production defaults on the FFTW certification
  // path. The small naive-DFT fallback differs only in its explicitly
  // non-certifying bounded cutoff.
  options.tree_options.opening_theta = solver_profile.opening_theta;
  options.tree_options.opening_criterion = solver_profile.opening_criterion;
  options.tree_options.multipole_order = solver_profile.multipole_order;
  options.tree_options.relative_force_tolerance = solver_profile.relative_force_tolerance;
  options.tree_options.relative_force_acceleration_floor_code =
      solver_profile.relative_force_acceleration_floor_code;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.max_leaf_size = solver_profile.max_leaf_size;
  options.tree_options.softening.epsilon_comoving = 0.0;

  const double spacing_x = fixture.box.length_x / static_cast<double>(fixture.pm_shape.nx);
  const double spacing_y = fixture.box.length_y / static_cast<double>(fixture.pm_shape.ny);
  const double spacing_z = fixture.box.length_z / static_cast<double>(fixture.pm_shape.nz);
  const double representative_spacing = std::cbrt(spacing_x * spacing_y * spacing_z);
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      k_asmth_cells, k_rcut_cells, representative_spacing);
  return options;
}

[[nodiscard]] ForceField solveProduction(
    const ParticleFixture& fixture,
    const SolverProfile& solver_profile,
    const ForceField* converged_reference_force_scale,
    cosmosim::gravity::TreePmProfileEvent* profile,
    cosmosim::gravity::TreePmDiagnostics* diagnostics) {
  std::vector<std::uint32_t> active(fixture.mass.size());
  std::iota(active.begin(), active.end(), 0U);
  std::vector<double> previous_acceleration_magnitude;
  if (solver_profile.supply_converged_reference_force_scale) {
    requireOrThrow(
        converged_reference_force_scale != nullptr,
        solver_profile.label + ": relative MAC profile requires a converged reference force scale");
    previous_acceleration_magnitude.resize(active.size());
    for (std::size_t i = 0; i < active.size(); ++i) {
      previous_acceleration_magnitude[i] = vectorNorm(
          converged_reference_force_scale->ax[i],
          converged_reference_force_scale->ay[i],
          converged_reference_force_scale->az[i]);
    }
  }
  ForceField force{
      std::vector<double>(active.size(), 0.0),
      std::vector<double>(active.size(), 0.0),
      std::vector<double>(active.size(), 0.0),
  };
  const cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      .active_particle_index = active,
      .accel_x_comoving = force.ax,
      .accel_y_comoving = force.ay,
      .accel_z_comoving = force.az,
      .previous_acceleration_magnitude_code = previous_acceleration_magnitude,
  };
  cosmosim::gravity::TreePmCoordinator coordinator(fixture.pm_shape);
  const cosmosim::gravity::TreePmOptions options = makeProductionOptions(fixture, solver_profile);
  coordinator.solveActiveSet(
      fixture.pos_x,
      fixture.pos_y,
      fixture.pos_z,
      fixture.mass,
      accumulator,
      options,
      profile,
      diagnostics);
  return force;
}

[[nodiscard]] ForceField solveProductionAtIndependentTargets(
    const IndependentTargetFixture& fixture,
    const SolverProfile& solver_profile,
    const ForceField& converged_reference_force_scale,
    cosmosim::gravity::TreePmProfileEvent* profile) {
  const std::size_t target_count = fixture.target_pos_x.size();
  requireOrThrow(
      fixture.target_pos_y.size() == target_count && fixture.target_pos_z.size() == target_count,
      fixture.sources.label + ": independent target coordinate extents differ");
  requireOrThrow(
      converged_reference_force_scale.ax.size() == target_count &&
          converged_reference_force_scale.ay.size() == target_count &&
          converged_reference_force_scale.az.size() == target_count,
      fixture.sources.label + ": independent target reference extent differs");

  std::vector<std::uint32_t> active(
      target_count, std::numeric_limits<std::uint32_t>::max());
  std::vector<double> previous_acceleration_magnitude;
  if (solver_profile.supply_converged_reference_force_scale) {
    previous_acceleration_magnitude.resize(target_count);
    for (std::size_t i = 0; i < target_count; ++i) {
      previous_acceleration_magnitude[i] = vectorNorm(
          converged_reference_force_scale.ax[i],
          converged_reference_force_scale.ay[i],
          converged_reference_force_scale.az[i]);
    }
  }
  ForceField force{
      std::vector<double>(target_count, 0.0),
      std::vector<double>(target_count, 0.0),
      std::vector<double>(target_count, 0.0),
  };
  const cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      .active_particle_index = active,
      .accel_x_comoving = force.ax,
      .accel_y_comoving = force.ay,
      .accel_z_comoving = force.az,
      .previous_acceleration_magnitude_code = previous_acceleration_magnitude,
      .target_pos_x_comoving = fixture.target_pos_x,
      .target_pos_y_comoving = fixture.target_pos_y,
      .target_pos_z_comoving = fixture.target_pos_z,
  };
  cosmosim::gravity::TreePmCoordinator coordinator(fixture.sources.pm_shape);
  const cosmosim::gravity::TreePmOptions options =
      makeProductionOptions(fixture.sources, solver_profile);
  coordinator.solveActiveSet(
      fixture.sources.pos_x,
      fixture.sources.pos_y,
      fixture.sources.pos_z,
      fixture.sources.mass,
      accumulator,
      options,
      profile,
      nullptr);
  return force;
}

[[nodiscard]] ForceField solveEwaldReferenceAtIndependentTargets(
    const IndependentTargetFixture& fixture) {
  std::vector<cosmosim::test_support::PeriodicEwaldSource> sources;
  sources.reserve(fixture.sources.mass.size());
  for (std::size_t i = 0; i < fixture.sources.mass.size(); ++i) {
    sources.push_back({
        .position = {
            .x = fixture.sources.pos_x[i],
            .y = fixture.sources.pos_y[i],
            .z = fixture.sources.pos_z[i],
        },
        .mass = fixture.sources.mass[i],
    });
  }
  std::vector<cosmosim::test_support::PeriodicEwaldTarget> targets;
  targets.reserve(fixture.target_pos_x.size());
  for (std::size_t i = 0; i < fixture.target_pos_x.size(); ++i) {
    targets.push_back({
        .position = {
            .x = fixture.target_pos_x[i],
            .y = fixture.target_pos_y[i],
            .z = fixture.target_pos_z[i],
        },
    });
  }
  const cosmosim::test_support::PeriodicEwaldOptions options{
      .gravitational_constant = 1.0,
      .alpha_inverse_length = 2.5,
      .real_image_limits = {3, 3, 3},
      .reciprocal_mode_limits = {10, 10, 10},
  };
  const auto acceleration = cosmosim::test_support::periodicEwaldAccelerations(
      sources, targets, fixture.sources.box, options);
  ForceField reference{
      std::vector<double>(acceleration.size()),
      std::vector<double>(acceleration.size()),
      std::vector<double>(acceleration.size()),
  };
  for (std::size_t i = 0; i < acceleration.size(); ++i) {
    reference.ax[i] = acceleration[i].x;
    reference.ay[i] = acceleration[i].y;
    reference.az[i] = acceleration[i].z;
  }
  return reference;
}

[[nodiscard]] ForceField solveEwaldReference(const ParticleFixture& fixture) {
  std::vector<cosmosim::test_support::PeriodicEwaldSource> sources;
  sources.reserve(fixture.mass.size());
  for (std::size_t i = 0; i < fixture.mass.size(); ++i) {
    sources.push_back({
        .position = {.x = fixture.pos_x[i], .y = fixture.pos_y[i], .z = fixture.pos_z[i]},
        .mass = fixture.mass[i],
    });
  }

  const double minimum_length = std::min({
      fixture.box.length_x,
      fixture.box.length_y,
      fixture.box.length_z,
  });
  const auto mode_limit = [minimum_length](double length) {
    return static_cast<int>(std::ceil(10.0 * length / minimum_length));
  };
  const cosmosim::test_support::PeriodicEwaldOptions options{
      .gravitational_constant = 1.0,
      .alpha_inverse_length = 2.5 / minimum_length,
      .real_image_limits = {3, 3, 3},
      .reciprocal_mode_limits = {
          mode_limit(fixture.box.length_x),
          mode_limit(fixture.box.length_y),
          mode_limit(fixture.box.length_z),
      },
  };
  const auto acceleration = cosmosim::test_support::periodicEwaldAccelerationsAtSources(
      sources, fixture.box, options);
  ForceField reference{
      std::vector<double>(acceleration.size()),
      std::vector<double>(acceleration.size()),
      std::vector<double>(acceleration.size()),
  };
  for (std::size_t i = 0; i < acceleration.size(); ++i) {
    reference.ax[i] = acceleration[i].x;
    reference.ay[i] = acceleration[i].y;
    reference.az[i] = acceleration[i].z;
  }
  return reference;
}

[[nodiscard]] double percentileNearestRank(std::vector<double> values, double probability) {
  requireOrThrow(!values.empty(), "percentile requires at least one value");
  requireOrThrow(probability >= 0.0 && probability <= 1.0, "percentile probability is outside [0,1]");
  std::sort(values.begin(), values.end());
  if (probability == 0.0) {
    return values.front();
  }
  const std::size_t rank = static_cast<std::size_t>(
      std::ceil(probability * static_cast<double>(values.size())));
  return values[std::min(rank - 1U, values.size() - 1U)];
}

[[nodiscard]] ErrorMetrics computeErrorMetrics(const ForceField& value, const ForceField& reference) {
  requireOrThrow(value.ax.size() == reference.ax.size(), "force/reference size mismatch");
  const std::size_t count = value.ax.size();
  double error_sum_squared = 0.0;
  double reference_sum_squared = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    error_sum_squared += squaredNorm(
        value.ax[i] - reference.ax[i],
        value.ay[i] - reference.ay[i],
        value.az[i] - reference.az[i]);
    reference_sum_squared += squaredNorm(reference.ax[i], reference.ay[i], reference.az[i]);
  }
  const double reference_rms = std::sqrt(reference_sum_squared / static_cast<double>(count));
  ErrorMetrics metrics;
  metrics.absolute_rms = std::sqrt(error_sum_squared / static_cast<double>(count));
  metrics.relative_l2 = std::sqrt(error_sum_squared / std::max(reference_sum_squared, 1.0e-300));
  metrics.normalized_floor = k_normalized_error_floor_fraction * reference_rms;

  std::vector<double> normalized_errors;
  normalized_errors.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    const double error = vectorNorm(
        value.ax[i] - reference.ax[i],
        value.ay[i] - reference.ay[i],
        value.az[i] - reference.az[i]);
    const double reference_norm = vectorNorm(reference.ax[i], reference.ay[i], reference.az[i]);
    normalized_errors.push_back(error / std::max(reference_norm, metrics.normalized_floor));
  }
  metrics.median_normalized = percentileNearestRank(normalized_errors, 0.5);
  metrics.p90_normalized = percentileNearestRank(normalized_errors, 0.9);
  metrics.p95_normalized = percentileNearestRank(normalized_errors, 0.95);
  metrics.p99_normalized = percentileNearestRank(normalized_errors, 0.99);
  metrics.max_normalized = *std::max_element(normalized_errors.begin(), normalized_errors.end());
  return metrics;
}

[[nodiscard]] double massWeightedNetForceFraction(
    const ParticleFixture& fixture,
    const ForceField& force) {
  double net_x = 0.0;
  double net_y = 0.0;
  double net_z = 0.0;
  double force_scale = 0.0;
  for (std::size_t i = 0; i < fixture.mass.size(); ++i) {
    net_x += fixture.mass[i] * force.ax[i];
    net_y += fixture.mass[i] * force.ay[i];
    net_z += fixture.mass[i] * force.az[i];
    force_scale += fixture.mass[i] * vectorNorm(force.ax[i], force.ay[i], force.az[i]);
  }
  return vectorNorm(net_x, net_y, net_z) / std::max(force_scale, 1.0e-300);
}

[[nodiscard]] double relativeL2Difference(const ForceField& value, const ForceField& reference) {
  double difference_squared = 0.0;
  double reference_squared = 0.0;
  for (std::size_t i = 0; i < value.ax.size(); ++i) {
    difference_squared += squaredNorm(
        value.ax[i] - reference.ax[i],
        value.ay[i] - reference.ay[i],
        value.az[i] - reference.az[i]);
    reference_squared += squaredNorm(reference.ax[i], reference.ay[i], reference.az[i]);
  }
  return std::sqrt(difference_squared / std::max(reference_squared, 1.0e-300));
}

[[nodiscard]] ParticleFixture translateByIndependentBoxImages(const ParticleFixture& fixture) {
  // Independent integer-image shifts are stronger than shifting the complete
  // fixture by one common box vector: every source representation changes while
  // the physical periodic system remains identical.
  ParticleFixture translated = fixture;
  for (std::size_t i = 0; i < fixture.mass.size(); ++i) {
    const int shift_x = static_cast<int>((3U * i + 1U) % 7U) - 3;
    const int shift_y = static_cast<int>((5U * i + 2U) % 9U) - 4;
    const int shift_z = static_cast<int>((7U * i + 3U) % 11U) - 5;
    translated.pos_x[i] += static_cast<double>(shift_x) * fixture.box.length_x;
    translated.pos_y[i] += static_cast<double>(shift_y) * fixture.box.length_y;
    translated.pos_z[i] += static_cast<double>(shift_z) * fixture.box.length_z;
  }
  return translated;
}

[[nodiscard]] const char* assignmentName(cosmosim::gravity::PmAssignmentScheme scheme) {
  return scheme == cosmosim::gravity::PmAssignmentScheme::kCic ? "cic" : "tsc";
}

[[nodiscard]] const char* multipoleName(cosmosim::gravity::TreeMultipoleOrder order) {
  return order == cosmosim::gravity::TreeMultipoleOrder::kMonopole ? "monopole" : "quadrupole";
}

[[nodiscard]] const char* openingCriterionName(cosmosim::gravity::TreeOpeningCriterion criterion) {
  switch (criterion) {
    case cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric:
      return "geometric";
    case cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance:
      return "com_distance";
    case cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError:
      return "relative_force_error";
  }
  return "unknown";
}

[[nodiscard]] AccuracyResult runAccuracyCase(
    const ParticleFixture& fixture,
    const ForceField& ewald_reference,
    const SolverProfile& solver_profile,
    bool measure_translation) {
  cosmosim::gravity::TreePmDiagnostics diagnostics;
  cosmosim::gravity::TreePmProfileEvent profile;
  const ForceField production = solveProduction(
      fixture, solver_profile, &ewald_reference, &profile, &diagnostics);
  AccuracyResult result;
  result.fixture_label = fixture.label;
  result.topology_label = fixture.topology_label;
  result.separation_bin = fixture.separation_bin;
  result.mesh_direction_label = fixture.mesh_direction_label;
  result.solver_profile = solver_profile;
  result.error = computeErrorMetrics(production, ewald_reference);
  result.reference_force_l2 = forceL2Norm(ewald_reference);
  result.net_force_fraction = massWeightedNetForceFraction(fixture, production);
  result.reference_net_force_fraction = massWeightedNetForceFraction(fixture, ewald_reference);
  result.classified_separation = fixture.classified_separation;
  result.classified_split_scale = fixture.classified_split_scale;
  result.classified_cutoff = fixture.classified_cutoff;
  result.pm_shape = fixture.pm_shape;
  result.diagnostics = diagnostics;
  result.profile = profile;
  if (measure_translation) {
    const ParticleFixture translated = translateByIndependentBoxImages(fixture);
    cosmosim::gravity::TreePmDiagnostics translated_diagnostics;
    cosmosim::gravity::TreePmProfileEvent translated_profile;
    const ForceField translated_force = solveProduction(
        translated,
        solver_profile,
        &ewald_reference,
        &translated_profile,
        &translated_diagnostics);
    result.translation_relative_l2 = relativeL2Difference(translated_force, production);
  }
  return result;
}

[[nodiscard]] double requiredRelativeL2(const AccuracyResult& result) {
  static_cast<void>(result);
  return k_required_relative_l2;
}

[[nodiscard]] double requiredP99Normalized(const AccuracyResult& result) {
  static_cast<void>(result);
  return k_required_p99_normalized;
}

[[nodiscard]] bool meetsAccuracyTarget(const AccuracyResult& result) {
  return result.error.relative_l2 <= requiredRelativeL2(result) &&
         result.error.p99_normalized <= requiredP99Normalized(result);
}

[[nodiscard]] bool isCertifiedProfile(const AccuracyResult& result) {
#if COSMOSIM_ENABLE_FFTW
  return result.solver_profile.accuracy_certified;
#else
  static_cast<void>(result);
  return false;
#endif
}

[[nodiscard]] const char* profileName(const AccuracyResult& result) {
  if (isCertifiedProfile(result)) {
    return "certified";
  }
#if COSMOSIM_ENABLE_FFTW
  return "diagnostic";
#else
  return "fallback_diagnostic";
#endif
}

void reportResult(const AccuracyResult& result) {
  const auto& solver_profile = result.solver_profile;
  std::cout << std::scientific << std::setprecision(9)
            << "treepm_ewald fixture=" << result.fixture_label
            << " topology=" << result.topology_label
            << " separation_bin=" << result.separation_bin
            << " mesh_direction=" << result.mesh_direction_label
            << " solver_profile=" << solver_profile.label
            << " assignment=" << assignmentName(solver_profile.assignment_scheme)
            << " deconvolution=" << (solver_profile.enable_window_deconvolution ? "on" : "off")
            << " multipole=" << multipoleName(solver_profile.multipole_order)
            << " opening_criterion=" << openingCriterionName(solver_profile.opening_criterion)
            << " theta=" << solver_profile.opening_theta
            << " relative_alpha=" << solver_profile.relative_force_tolerance
            << " relative_floor=" << solver_profile.relative_force_acceleration_floor_code
            << " previous_force_scale="
            << (solver_profile.supply_converged_reference_force_scale ? "converged_ewald" : "none")
            << " max_leaf_size=" << solver_profile.max_leaf_size
            << " profile=" << profileName(result)
            << " target_result=" << (meetsAccuracyTarget(result) ? "pass" : "fail")
            << " abs_rms=" << result.error.absolute_rms
            << " rel_l2=" << result.error.relative_l2
            << " reference_force_l2=" << result.reference_force_l2
            << " median=" << result.error.median_normalized
            << " p90=" << result.error.p90_normalized
            << " p95=" << result.error.p95_normalized
            << " p99=" << result.error.p99_normalized
            << " max=" << result.error.max_normalized
            << " denominator_floor=" << result.error.normalized_floor
            << " net_force_fraction=" << result.net_force_fraction
            << " reference_net_force_fraction=" << result.reference_net_force_fraction;
  if (std::isfinite(result.classified_separation)) {
    std::cout << " separation=" << result.classified_separation
              << " separation_over_split="
              << result.classified_separation / result.classified_split_scale
              << " gaussian_q="
              << result.classified_separation / (2.0 * result.classified_split_scale)
              << " separation_over_cutoff="
              << result.classified_separation / result.classified_cutoff;
  } else {
    std::cout << " separation=aggregate separation_over_split=aggregate"
              << " gaussian_q=aggregate separation_over_cutoff=aggregate";
  }
  if (std::isfinite(result.translation_relative_l2)) {
    std::cout << " translation_drift_rel_l2=" << result.translation_relative_l2;
  } else {
    std::cout << " translation_drift_rel_l2=not_measured";
  }
  std::cout << " pm_grid=" << result.pm_shape.nx << 'x' << result.pm_shape.ny << 'x' << result.pm_shape.nz
            << " mesh_spacing=" << result.diagnostics.mesh_spacing_comoving
            << " split_scale=" << result.diagnostics.split_scale_comoving
            << " cutoff=" << result.diagnostics.cutoff_radius_comoving
            << " pm_force_l2=" << result.diagnostics.force_l2_pm_global
            << " tree_force_l2=" << result.diagnostics.force_l2_tree_short_range
            << " total_force_l2=" << result.diagnostics.force_l2_total
            << " tree_visited=" << result.profile.tree_profile.visited_nodes
            << " tree_accepted=" << result.profile.tree_profile.accepted_nodes
            << " tree_particle_pairs=" << result.profile.tree_profile.particle_particle_interactions
            << " tree_build_ms=" << result.profile.tree_profile.build_ms
            << " tree_multipole_ms=" << result.profile.tree_profile.multipole_ms
            << " tree_traversal_ms=" << result.profile.tree_profile.traversal_ms
            << std::endl;
}

void requirePhysicalInvariants(const AccuracyResult& result) {
  std::ostringstream message;
  message << std::setprecision(17)
          << result.fixture_label << '/' << result.solver_profile.label
          << " TreePM-vs-Ewald validation failed: rel_l2=" << result.error.relative_l2
          << " required<=" << requiredRelativeL2(result)
          << ", p99=" << result.error.p99_normalized
          << " required<=" << requiredP99Normalized(result)
          << ", max=" << result.error.max_normalized
          << ", denominator_floor=" << result.error.normalized_floor
          << ", epsilon=0, G=1, a=1, asmth=" << k_asmth_cells
          << ", rcut=" << k_rcut_cells
          << ", theta=" << result.solver_profile.opening_theta
          << ", multipole=" << multipoleName(result.solver_profile.multipole_order)
          << ", MAC=" << openingCriterionName(result.solver_profile.opening_criterion);
  requireOrThrow(
      std::isfinite(result.error.absolute_rms) && std::isfinite(result.error.relative_l2) &&
          std::isfinite(result.error.median_normalized) &&
          std::isfinite(result.error.p90_normalized) && std::isfinite(result.error.p95_normalized) &&
          std::isfinite(result.error.p99_normalized) && std::isfinite(result.error.max_normalized),
      message.str());
  requireOrThrow(
      std::isfinite(result.net_force_fraction),
      message.str());
  requireOrThrow(
      std::isfinite(result.reference_net_force_fraction) &&
          result.reference_net_force_fraction <= k_required_net_force_fraction,
      message.str());
  if (std::isfinite(result.translation_relative_l2)) {
    requireOrThrow(
        result.translation_relative_l2 <= k_required_periodic_translation_drift,
        message.str());
  }
}

void requireCertifiedAccuracy(const AccuracyResult& result) {
  requirePhysicalInvariants(result);
  if (!isCertifiedProfile(result)) {
    return;
  }
  requireOrThrow(
      result.net_force_fraction <= k_required_net_force_fraction,
      result.fixture_label + '/' + result.solver_profile.label +
          " certified profile violated the mass-weighted net-force tolerance");
  std::ostringstream message;
  message << std::setprecision(17)
          << result.fixture_label << '/' << result.solver_profile.label
          << " certified TreePM-vs-Ewald accuracy failed: rel_l2=" << result.error.relative_l2
          << " required<=" << requiredRelativeL2(result)
          << ", p99=" << result.error.p99_normalized
          << " required<=" << requiredP99Normalized(result)
          << ", epsilon=0, G=1, a=1, asmth=" << k_asmth_cells
          << ", rcut=" << k_rcut_cells
          << ", theta=" << result.solver_profile.opening_theta
          << ", multipole=" << multipoleName(result.solver_profile.multipole_order)
          << ", MAC=" << openingCriterionName(result.solver_profile.opening_criterion);
  requireOrThrow(meetsAccuracyTarget(result), message.str());
}

[[nodiscard]] SolverProfile certifiedDefaultProfile() {
  const double configured_default_theta =
      cosmosim::core::makeUnvalidatedSimulationConfigForTests()
          .numerics.treepm_tree_opening_theta;
  return SolverProfile{
      .label = "default_tsc_deconv_quadrupole_com",
      .assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc,
      .enable_window_deconvolution = true,
      .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole,
      .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance,
      .opening_theta = configured_default_theta,
      .relative_force_tolerance = 0.005,
      .relative_force_acceleration_floor_code = 1.0e-30,
      .max_leaf_size = 16,
      .supply_converged_reference_force_scale = false,
      .accuracy_certified = true,
  };
}

void runAndRecord(
    const ParticleFixture& fixture,
    const ForceField& ewald_reference,
    const SolverProfile& solver_profile,
    bool measure_translation,
    std::vector<AccuracyResult>* results) {
  requireOrThrow(results != nullptr, "accuracy result destination is null");
  results->push_back(runAccuracyCase(
      fixture, ewald_reference, solver_profile, measure_translation));
  reportResult(results->back());
}

void testTreePmAgainstIndependentEwald() {
  std::vector<ParticleFixture> fixtures;
  fixtures.push_back(makeQuasiUniformFixture());
  fixtures.push_back(makeDmoZeldovichLatticeFixture());
  fixtures.push_back(makeClusteredFixture());
  fixtures.push_back(makeSeamFixture());
  fixtures.push_back(makeTwoBodySplitFixture());
  for (const auto& fixture : fixtures) {
    validateFixture(fixture);
  }

  std::vector<AccuracyResult> results;
  const SolverProfile certified_profile = certifiedDefaultProfile();

  // Accuracy certification: these fixtures span quasi-uniform, the exact
  // cancellation-dominated 4^3 DMO initial lattice, clustered,
  // periodic-seam/rectangular, and Gaussian split-transition configurations.
  // The DMO mesh is 16^3 so its r_cut=6.25/16 remains strictly inside the
  // half-box minimum-image domain used by the production short-range tree.
  // Only this exact FFTW/TSC/deconvolved/quadrupole/COM-distance profile is a
  // release gate. All following matrix variants are deliberately diagnostic:
  // they print the same metrics and pass/fail classification without moving
  // the calibrated production threshold to accommodate a weaker profile.
  for (const auto& fixture : fixtures) {
    const ForceField ewald_reference = solveEwaldReference(fixture);
    const bool measure_translation =
        fixture.label == "xyz_seam_rectangular" ||
        fixture.label == "dmo_zeldovich_4cubed_initial";
    runAndRecord(fixture, ewald_reference, certified_profile, measure_translation, &results);
    if (fixture.label == "dmo_zeldovich_4cubed_initial") {
      requireOrThrow(
          results.back().diagnostics.cutoff_radius_comoving <
              0.5 * fixture.box.length_x,
          "DMO lattice fixture violated the periodic minimum-image cutoff invariant");
    }
  }

  // PM assignment/window matrix. The certified TSC+deconvolution-on corner is
  // already emitted above for the same quasi-uniform fixture, leaving three
  // non-duplicated diagnostic corners here.
  const ParticleFixture assignment_fixture = makeQuasiUniformFixture();
  const ForceField assignment_reference = solveEwaldReference(assignment_fixture);
  for (const auto [assignment_scheme, deconvolution] : {
           std::pair{cosmosim::gravity::PmAssignmentScheme::kCic, false},
           std::pair{cosmosim::gravity::PmAssignmentScheme::kCic, true},
           std::pair{cosmosim::gravity::PmAssignmentScheme::kTsc, false},
       }) {
    SolverProfile profile = certified_profile;
    profile.label = std::string("assignment_window_") + assignmentName(assignment_scheme) +
                    (deconvolution ? "_deconv_on" : "_deconv_off");
    profile.assignment_scheme = assignment_scheme;
    profile.enable_window_deconvolution = deconvolution;
    profile.accuracy_certified = false;
    runAndRecord(assignment_fixture, assignment_reference, profile, false, &results);
  }

  // Multipole comparison on source targets. TreePM's distributed-equivalence
  // envelope deliberately drives monopoles to leaves, so this is not used to
  // pretend that a theta value controls monopole work.
  const ParticleFixture mac_fixture = makeMacComparisonFixture();
  validateFixture(mac_fixture);
  const ForceField mac_reference = solveEwaldReference(mac_fixture);
  for (const auto multipole_order : {
           cosmosim::gravity::TreeMultipoleOrder::kMonopole,
           cosmosim::gravity::TreeMultipoleOrder::kQuadrupole,
       }) {
    SolverProfile profile = certified_profile;
    profile.label = std::string("multipole_matrix_") + multipoleName(multipole_order);
    profile.multipole_order = multipole_order;
    profile.max_leaf_size = 4;
    profile.accuracy_certified = false;
    runAndRecord(mac_fixture, mac_reference, profile, false, &results);
  }

  // Calibrated MAC comparison. Independent targets sit outside a compact
  // source node but inside r_cut, avoiding target-containing-node recursion.
  // For the relative criterion with a converged force scale,
  // (width/r)^2 <= alpha. Therefore theta=sqrt(alpha) gives all three MACs
  // the same leading geometric acceptance scale. The value remains below the
  // decomposition-stability width/r envelope, so the configured MAC -- not
  // that safety envelope -- controls internal-node acceptance.
  const IndependentTargetFixture matched_fixture = makeMatchedMacFixture();
  validateFixture(matched_fixture.sources);
  const ForceField matched_reference =
      solveEwaldReferenceAtIndependentTargets(matched_fixture);
  constexpr double k_matched_relative_alpha = 0.005;
  constexpr double k_matched_opening_theta = 0.07071067811865475;
  std::vector<MatchedMacResult> matched_results;
  for (const auto opening_criterion : {
           cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric,
           cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance,
           cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError,
       }) {
    SolverProfile profile = certified_profile;
    profile.label = std::string("matched_mac_") + openingCriterionName(opening_criterion);
    profile.opening_criterion = opening_criterion;
    profile.opening_theta = k_matched_opening_theta;
    profile.relative_force_tolerance = k_matched_relative_alpha;
    profile.max_leaf_size = 4;
    profile.supply_converged_reference_force_scale =
        opening_criterion == cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError;
    profile.accuracy_certified = false;

    MatchedMacResult result;
    result.solver_profile = profile;
    const ForceField production = solveProductionAtIndependentTargets(
        matched_fixture, profile, matched_reference, &result.profile);
    result.error = computeErrorMetrics(production, matched_reference);
    std::cout << std::scientific << std::setprecision(9)
              << "treepm_ewald fixture=" << matched_fixture.sources.label
              << " solver_profile=" << profile.label
              << " multipole=" << multipoleName(profile.multipole_order)
              << " opening_criterion=" << openingCriterionName(profile.opening_criterion)
              << " theta=" << profile.opening_theta
              << " relative_alpha=" << profile.relative_force_tolerance
              << " calibration=theta_squared_equals_alpha"
              << " rel_l2=" << result.error.relative_l2
              << " p99=" << result.error.p99_normalized
              << " tree_visited=" << result.profile.tree_profile.visited_nodes
              << " tree_accepted=" << result.profile.tree_profile.accepted_nodes
              << " tree_opened=" << result.profile.tree_profile.opened_nodes
              << " tree_particle_pairs=" << result.profile.tree_profile.particle_particle_interactions
              << " tree_traversal_ms=" << result.profile.tree_profile.traversal_ms
              << std::endl;

    const std::uint64_t direct_pair_count = static_cast<std::uint64_t>(
        matched_fixture.sources.mass.size() * matched_fixture.target_pos_x.size());
    std::ostringstream context;
    context << profile.label << " calibrated MAC validation failed: rel_l2="
            << result.error.relative_l2 << ", p99=" << result.error.p99_normalized
            << ", visited=" << result.profile.tree_profile.visited_nodes
            << ", accepted=" << result.profile.tree_profile.accepted_nodes
            << ", opened=" << result.profile.tree_profile.opened_nodes
            << ", particle_pairs=" << result.profile.tree_profile.particle_particle_interactions
            << ", direct_pairs=" << direct_pair_count;
    requireOrThrow(meetsAccuracyTarget(AccuracyResult{.error = result.error}), context.str());
    requireOrThrow(result.profile.tree_profile.visited_nodes > 0U, context.str());
    requireOrThrow(result.profile.tree_profile.accepted_nodes > 0U, context.str());
    requireOrThrow(result.profile.tree_profile.opened_nodes > 0U, context.str());
    requireOrThrow(
        result.profile.tree_profile.particle_particle_interactions < direct_pair_count,
        context.str());
    matched_results.push_back(std::move(result));
  }

  const auto error_bounds = std::minmax_element(
      matched_results.begin(), matched_results.end(),
      [](const MatchedMacResult& lhs, const MatchedMacResult& rhs) {
        return lhs.error.relative_l2 < rhs.error.relative_l2;
      });
  const auto p99_bounds = std::minmax_element(
      matched_results.begin(), matched_results.end(),
      [](const MatchedMacResult& lhs, const MatchedMacResult& rhs) {
        return lhs.error.p99_normalized < rhs.error.p99_normalized;
      });
  const auto visited_bounds = std::minmax_element(
      matched_results.begin(), matched_results.end(),
      [](const MatchedMacResult& lhs, const MatchedMacResult& rhs) {
        return lhs.profile.tree_profile.visited_nodes < rhs.profile.tree_profile.visited_nodes;
      });
  requireOrThrow(
      error_bounds.second->error.relative_l2 <=
          1.25 * std::max(error_bounds.first->error.relative_l2, 1.0e-12),
      "calibrated geometric/COM/relative MAC relative-L2 errors are not approximately matched");
  requireOrThrow(
      p99_bounds.second->error.p99_normalized <=
          1.25 * std::max(p99_bounds.first->error.p99_normalized, 1.0e-12),
      "calibrated geometric/COM/relative MAC p99 errors are not approximately matched");
  requireOrThrow(
      visited_bounds.second->profile.tree_profile.visited_nodes * 2U <=
          visited_bounds.first->profile.tree_profile.visited_nodes * 3U,
      "calibrated geometric/COM/relative MAC traversal work differs by more than 1.5x");

  // Separation/direction classification uses the certified numerical settings
  // but remains diagnostic because two-body lattice anisotropy near a hard
  // cutoff is not the aggregate calibration fixture. The cases bracket q=1
  // for the Gaussian split and r/r_cut=1 for both axis-aligned and diagonal
  // directions.
  for (const ParticleFixture& separation_fixture : makeSeparationClassificationFixtures()) {
    validateFixture(separation_fixture);
    const ForceField separation_reference = solveEwaldReference(separation_fixture);
    SolverProfile profile = certified_profile;
    profile.label = "separation_direction_classification";
    profile.accuracy_certified = false;
    runAndRecord(separation_fixture, separation_reference, profile, false, &results);
  }

  for (const auto& result : results) {
    requireCertifiedAccuracy(result);
  }
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
#endif
  testTreePmAgainstIndependentEwald();
#if COSMOSIM_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
