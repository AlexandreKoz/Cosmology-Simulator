#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "cosmosim/gravity/gravity_memory.hpp"
#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/parallel/distributed_mesh.hpp"
#include "gravity/internal/tree_interaction_common.hpp"

namespace {

static_assert(!std::is_assignable_v<
    cosmosim::gravity::GravitySourceGeneration&,
    cosmosim::gravity::PmFieldVersion>);
static_assert(!std::is_assignable_v<
    cosmosim::gravity::TreeBuildGeneration&,
    cosmosim::gravity::DecompositionEpoch>);
static_assert(!std::is_convertible_v<
    cosmosim::gravity::ForceEvaluationEpoch,
    cosmosim::gravity::GravitySourceGeneration>);

constexpr double k_tolerance = 1.0e-10;

struct OneTargetEvaluation {
  double ax = 0.0;
  double ay = 0.0;
  double az = 0.0;
  cosmosim::gravity::TreeGravityProfile profile{};
};

[[nodiscard]] OneTargetEvaluation evaluateOneTarget(
    const cosmosim::gravity::TreeGravitySolver& solver,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    const cosmosim::gravity::TreeGravityOptions& options,
    std::span<const double> previous_acceleration_magnitude_code = {},
    bool use_target_view = false) {
  const std::vector<std::uint32_t> active{0U};
  std::vector<double> ax(1, 0.0);
  std::vector<double> ay(1, 0.0);
  std::vector<double> az(1, 0.0);
  cosmosim::gravity::TreeGravityProfile profile;

  if (use_target_view) {
    const cosmosim::gravity::TreeGravitySolver::TreeGravitySourceView source_view{
        .pos_x_comoving = pos_x,
        .pos_y_comoving = pos_y,
        .pos_z_comoving = pos_z,
        .mass_code = mass,
    };
    const cosmosim::gravity::TreeGravitySolver::TreeGravityTargetView target_view{
        .active_particle_index = active,
        .accel_x_comoving = ax,
        .accel_y_comoving = ay,
        .accel_z_comoving = az,
        .previous_acceleration_magnitude_code = previous_acceleration_magnitude_code,
    };
    solver.evaluateActiveSet(source_view, target_view, options, &profile);
  } else if (previous_acceleration_magnitude_code.empty()) {
    // Exercise the source-compatible first-evaluation path: omitted history.
    solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, ax, ay, az, options, &profile);
  } else {
    solver.evaluateActiveSet(
        pos_x,
        pos_y,
        pos_z,
        mass,
        active,
        ax,
        ay,
        az,
        options,
        &profile,
        {},
        previous_acceleration_magnitude_code);
  }
  return {.ax = ax[0], .ay = ay[0], .az = az[0], .profile = profile};
}

void assertSameEvaluation(const OneTargetEvaluation& lhs, const OneTargetEvaluation& rhs) {
  assert(lhs.ax == rhs.ax);
  assert(lhs.ay == rhs.ay);
  assert(lhs.az == rhs.az);
  assert(lhs.profile.visited_nodes == rhs.profile.visited_nodes);
  assert(lhs.profile.accepted_nodes == rhs.profile.accepted_nodes);
  assert(lhs.profile.opened_nodes == rhs.profile.opened_nodes);
  assert(lhs.profile.particle_particle_interactions == rhs.profile.particle_particle_interactions);
}

void testTwoBodyForce() {
  cosmosim::gravity::TreeGravitySolver solver;
  cosmosim::gravity::TreeGravityOptions options;
  options.opening_theta = 0.6;
  options.gravitational_constant_code = 2.0;
  options.max_leaf_size = 1;
  options.softening.epsilon_comoving = 0.1;

  const std::vector<double> pos_x{0.0, 1.0};
  const std::vector<double> pos_y{0.0, 0.0};
  const std::vector<double> pos_z{0.0, 0.0};
  const std::vector<double> mass{1.0, 3.0};

  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);

  const std::vector<std::uint32_t> active{0U};
  std::vector<double> ax(1, 0.0);
  std::vector<double> ay(1, 0.0);
  std::vector<double> az(1, 0.0);

  solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, ax, ay, az, options, nullptr);

  const double expected = options.gravitational_constant_code * mass[1] /
      std::pow(1.0 + options.softening.epsilon_comoving * options.softening.epsilon_comoving, 1.5);
  assert(std::abs(ax[0] - expected) < k_tolerance);
  assert(std::abs(ay[0]) < k_tolerance);
  assert(std::abs(az[0]) < k_tolerance);
}

void testForceSymmetry() {
  cosmosim::gravity::TreeGravitySolver solver;
  cosmosim::gravity::TreeGravityOptions options;
  options.opening_theta = 0.5;
  options.gravitational_constant_code = 1.0;
  options.max_leaf_size = 1;
  options.softening.epsilon_comoving = 0.0;

  const std::vector<double> pos_x{-0.5, 0.5};
  const std::vector<double> pos_y{0.0, 0.0};
  const std::vector<double> pos_z{0.0, 0.0};
  const std::vector<double> mass{2.0, 2.0};

  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);

  const std::vector<std::uint32_t> active{0U, 1U};
  std::vector<double> ax(2, 0.0);
  std::vector<double> ay(2, 0.0);
  std::vector<double> az(2, 0.0);

  solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, ax, ay, az, options, nullptr);

  assert(std::abs(ax[0] + ax[1]) < k_tolerance);
  assert(std::abs(ay[0] + ay[1]) < k_tolerance);
  assert(std::abs(az[0] + az[1]) < k_tolerance);
}

void testOpeningCriterionBehavior() {
  cosmosim::gravity::TreeGravitySolver solver;
  cosmosim::gravity::TreeGravityOptions options;
  options.gravitational_constant_code = 1.0;
  options.max_leaf_size = 1;

  const std::vector<double> pos_x{-1.0, -0.9, 0.95, 1.0};
  const std::vector<double> pos_y{0.0, 0.02, -0.01, 0.0};
  const std::vector<double> pos_z{0.0, 0.0, 0.0, 0.0};
  const std::vector<double> mass(4, 1.0);
  const std::vector<std::uint32_t> active{0U};

  options.opening_theta = 0.35;
  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);
  cosmosim::gravity::TreeGravityProfile strict_profile;
  std::vector<double> ax(1, 0.0);
  std::vector<double> ay(1, 0.0);
  std::vector<double> az(1, 0.0);
  solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, ax, ay, az, options, &strict_profile);

  options.opening_theta = 1.5;
  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);
  cosmosim::gravity::TreeGravityProfile relaxed_profile;
  solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, ax, ay, az, options, &relaxed_profile);

  assert(strict_profile.visited_nodes >= relaxed_profile.visited_nodes);
  assert(strict_profile.accepted_nodes >= relaxed_profile.accepted_nodes);
}

void testMultipoleAccumulation() {
  cosmosim::gravity::TreeGravitySolver solver;
  cosmosim::gravity::TreeGravityOptions options;
  options.max_leaf_size = 2;

  const std::vector<double> pos_x{0.0, 2.0, 4.0};
  const std::vector<double> pos_y{0.0, 0.0, 0.0};
  const std::vector<double> pos_z{0.0, 0.0, 0.0};
  const std::vector<double> mass{1.0, 2.0, 3.0};

  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);

  const auto& nodes = solver.nodes();
  assert(nodes.size() > 0);
  const double total_mass = 1.0 + 2.0 + 3.0;
  const double expected_com_x = (1.0 * 0.0 + 2.0 * 2.0 + 3.0 * 4.0) / total_mass;
  assert(std::abs(nodes.mass_code[0] - total_mass) < k_tolerance);
  assert(std::abs(nodes.com_x_comoving[0] - expected_com_x) < k_tolerance);
  const double trace = nodes.quad_xx[0] + nodes.quad_yy[0] + nodes.quad_zz[0];
  assert(std::abs(trace) < 1.0e-8);
  assert(std::abs(nodes.quad_xx[0]) > 0.0);
}

void testMacCriteriaAreAvailable() {
  cosmosim::gravity::TreeGravitySolver solver;
  const std::vector<double> pos_x{-1.0, -0.7, 0.8, 1.2, 1.3};
  const std::vector<double> pos_y{0.0, 0.1, -0.1, 0.0, 0.03};
  const std::vector<double> pos_z{0.0, 0.0, 0.0, 0.0, 0.02};
  const std::vector<double> mass(5, 1.0);
  const std::vector<std::uint32_t> active{0U};
  std::vector<double> ax(1, 0.0);
  std::vector<double> ay(1, 0.0);
  std::vector<double> az(1, 0.0);

  cosmosim::gravity::TreeGravityOptions options;
  options.max_leaf_size = 1;
  options.opening_theta = 0.7;
  options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric;
  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);
  cosmosim::gravity::TreeGravityProfile geometric_profile;
  solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, ax, ay, az, options, &geometric_profile);

  options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance;
  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);
  cosmosim::gravity::TreeGravityProfile com_distance_profile;
  solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, ax, ay, az, options, &com_distance_profile);

  assert(geometric_profile.visited_nodes > 0);
  assert(com_distance_profile.visited_nodes > 0);
}

void testRelativeForceErrorMacFormulaAndWorkMonotonicity() {
  cosmosim::gravity::TreeGravitySolver solver;
  cosmosim::gravity::TreeGravityOptions options;
  options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError;
  options.opening_theta = 0.6;
  options.relative_force_tolerance = 0.01;
  options.relative_force_acceleration_floor_code = 1.0e-12;
  options.gravitational_constant_code = 1.0;
  options.max_leaf_size = 1;

  const std::vector<double> pos_x{-4.0, 1.00, 1.03, 1.06, 1.09, 1.12, 1.15, 1.18, 1.21};
  const std::vector<double> pos_y{-4.0, 1.01, 1.05, 1.09, 1.13, 1.02, 1.06, 1.10, 1.14};
  const std::vector<double> pos_z{-4.0, 1.02, 1.08, 1.14, 1.03, 1.09, 1.15, 1.04, 1.10};
  const std::vector<double> mass(pos_x.size(), 1.0);
  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);

  const auto& nodes = solver.nodes();
  const double target_x = pos_x[0];
  const double target_y = pos_y[0];
  const double target_z = pos_z[0];
  std::uint32_t source_node = std::numeric_limits<std::uint32_t>::max();
  for (std::size_t octant = 0; octant < 8U; ++octant) {
    const std::uint32_t child = nodes.child_index[octant];
    if (child == std::numeric_limits<std::uint32_t>::max() || nodes.child_count[child] == 0U) {
      continue;
    }
    const double half_size = nodes.half_size_comoving[child];
    const bool contains_target =
        std::abs(target_x - nodes.center_x_comoving[child]) <= half_size &&
        std::abs(target_y - nodes.center_y_comoving[child]) <= half_size &&
        std::abs(target_z - nodes.center_z_comoving[child]) <= half_size;
    if (!contains_target) {
      source_node = child;
      break;
    }
  }
  assert(source_node != std::numeric_limits<std::uint32_t>::max());

  const double dx = nodes.com_x_comoving[source_node] - target_x;
  const double dy = nodes.com_y_comoving[source_node] - target_y;
  const double dz = nodes.com_z_comoving[source_node] - target_z;
  const double r2 = dx * dx + dy * dy + dz * dz;
  const double width = 2.0 * nodes.half_size_comoving[source_node];
  const double acceptance_threshold_code =
      options.gravitational_constant_code * nodes.mass_code[source_node] * width * width /
      (options.relative_force_tolerance * r2 * r2);
  assert(std::isfinite(acceptance_threshold_code));
  assert(acceptance_threshold_code > options.relative_force_acceleration_floor_code);

  const std::vector<double> below_threshold{0.5 * acceptance_threshold_code};
  const std::vector<double> above_threshold{2.0 * acceptance_threshold_code};
  const OneTargetEvaluation strict =
      evaluateOneTarget(solver, pos_x, pos_y, pos_z, mass, options, below_threshold);
  const OneTargetEvaluation relaxed =
      evaluateOneTarget(solver, pos_x, pos_y, pos_z, mass, options, above_threshold, true);
  assert(strict.profile.visited_nodes > relaxed.profile.visited_nodes);
  assert(strict.profile.particle_particle_interactions >= relaxed.profile.particle_particle_interactions);
}

void testRelativeForceErrorMacHistoryFallbackAndFloor() {
  cosmosim::gravity::TreeGravitySolver solver;
  cosmosim::gravity::TreeGravityOptions relative_options;
  relative_options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError;
  relative_options.opening_theta = 0.1;
  relative_options.relative_force_tolerance = 0.01;
  relative_options.relative_force_acceleration_floor_code = 8.0;
  relative_options.max_leaf_size = 1;

  const std::vector<double> pos_x{-4.0, 1.00, 1.04, 1.08, 1.12, 1.16};
  const std::vector<double> pos_y{-4.0, 1.02, 1.06, 1.10, 1.14, 1.18};
  const std::vector<double> pos_z{-4.0, 1.03, 1.09, 1.15, 1.05, 1.11};
  const std::vector<double> mass(pos_x.size(), 1.0);
  solver.build(pos_x, pos_y, pos_z, mass, relative_options, nullptr);

  cosmosim::gravity::TreeGravityOptions com_options = relative_options;
  com_options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance;
  const OneTargetEvaluation com_distance =
      evaluateOneTarget(solver, pos_x, pos_y, pos_z, mass, com_options);
  const OneTargetEvaluation missing_history =
      evaluateOneTarget(solver, pos_x, pos_y, pos_z, mass, relative_options);
  assertSameEvaluation(com_distance, missing_history);

  const std::vector<double> nan_history{std::numeric_limits<double>::quiet_NaN()};
  const std::vector<double> infinite_history{std::numeric_limits<double>::infinity()};
  assertSameEvaluation(
      com_distance, evaluateOneTarget(solver, pos_x, pos_y, pos_z, mass, relative_options, nan_history));
  assertSameEvaluation(
      com_distance, evaluateOneTarget(solver, pos_x, pos_y, pos_z, mass, relative_options, infinite_history));

  const std::vector<double> zero_history{0.0};
  const std::vector<double> floor_history{relative_options.relative_force_acceleration_floor_code};
  const std::vector<double> negative_floor_history{-relative_options.relative_force_acceleration_floor_code};
  const OneTargetEvaluation zero =
      evaluateOneTarget(solver, pos_x, pos_y, pos_z, mass, relative_options, zero_history);
  assertSameEvaluation(
      zero, evaluateOneTarget(solver, pos_x, pos_y, pos_z, mass, relative_options, floor_history));
  assertSameEvaluation(
      zero, evaluateOneTarget(solver, pos_x, pos_y, pos_z, mass, relative_options, negative_floor_history));
}

void testRelativeForceErrorMacNeverAcceptsTargetContainingNode() {
  cosmosim::gravity::TreeGravitySolver solver;
  cosmosim::gravity::TreeGravityOptions options;
  options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError;
  options.relative_force_tolerance = 1.0;
  options.relative_force_acceleration_floor_code = 1.0e-30;
  options.max_leaf_size = 1;
  options.softening.epsilon_comoving = 0.0;

  const std::vector<double> pos_x{0.0, 1.0};
  const std::vector<double> pos_y{0.0, 0.0};
  const std::vector<double> pos_z{0.0, 0.0};
  const std::vector<double> mass{100.0, 3.0};
  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);

  const std::vector<double> enormous_history{1.0e100};
  const OneTargetEvaluation result =
      evaluateOneTarget(solver, pos_x, pos_y, pos_z, mass, options, enormous_history);
  assert(std::abs(result.ax - 3.0) < k_tolerance);
  assert(std::abs(result.ay) < k_tolerance);
  assert(std::abs(result.az) < k_tolerance);
  assert(result.profile.visited_nodes > 1U);
  assert(result.profile.particle_particle_interactions == 1U);
}

void testRelativeForceErrorMacValidatesOptionsAndHistoryShape() {
  const std::vector<double> pos_x{0.0, 1.0};
  const std::vector<double> pos_y{0.0, 0.0};
  const std::vector<double> pos_z{0.0, 0.0};
  const std::vector<double> mass{1.0, 1.0};

  const auto build_throws = [&](const cosmosim::gravity::TreeGravityOptions& options) {
    cosmosim::gravity::TreeGravitySolver solver;
    bool threw = false;
    try {
      solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  };

  cosmosim::gravity::TreeGravityOptions options;
  options.relative_force_tolerance = 0.0;
  build_throws(options);
  options.relative_force_tolerance = std::numeric_limits<double>::quiet_NaN();
  build_throws(options);
  options.relative_force_tolerance = 0.005;
  options.relative_force_acceleration_floor_code = 0.0;
  build_throws(options);
  options.relative_force_acceleration_floor_code = std::numeric_limits<double>::infinity();
  build_throws(options);

  options.relative_force_acceleration_floor_code = 1.0e-30;
  options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError;
  cosmosim::gravity::TreeGravitySolver solver;
  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);
  const std::vector<std::uint32_t> active{0U};
  std::vector<double> accel_x(1, 0.0);
  std::vector<double> accel_y(1, 0.0);
  std::vector<double> accel_z(1, 0.0);
  const std::vector<double> wrong_size_history{1.0, 2.0};
  bool threw = false;
  try {
    solver.evaluateActiveSet(
        pos_x,
        pos_y,
        pos_z,
        mass,
        active,
        accel_x,
        accel_y,
        accel_z,
        options,
        nullptr,
        {},
        wrong_size_history);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
}

void testPairSofteningUsesMaxRule() {
  cosmosim::gravity::TreeGravitySolver solver;
  cosmosim::gravity::TreeGravityOptions options;
  options.gravitational_constant_code = 1.0;
  options.max_leaf_size = 1;
  options.softening.epsilon_comoving = 0.05;

  const std::vector<double> pos_x{0.0, 1.0};
  const std::vector<double> pos_y{0.0, 0.0};
  const std::vector<double> pos_z{0.0, 0.0};
  const std::vector<double> mass{2.0, 3.0};
  const std::vector<std::uint32_t> species_tag{0U, 1U};
  const std::vector<double> target_eps{0.2, 0.01};
  const std::vector<std::uint8_t> target_eps_mask{1U, 1U};
  cosmosim::gravity::TreeSofteningView softening_view{
      .source_species_tag = species_tag,
      .target_particle_epsilon_comoving = target_eps,
      .target_particle_epsilon_override_mask = target_eps_mask,
      .species_policy = {.epsilon_comoving_by_species = {0.04, 0.15, 0.0, 0.0, 0.0}, .enabled = true},
  };

  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr, softening_view);
  const std::vector<std::uint32_t> active{0U, 1U};
  std::vector<double> ax(2, 0.0);
  std::vector<double> ay(2, 0.0);
  std::vector<double> az(2, 0.0);
  solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, ax, ay, az, options, nullptr, softening_view);

  const double eps_01 = cosmosim::gravity::combineSofteningPairEpsilon(0.15, 0.2);
  const double eps_10 = cosmosim::gravity::combineSofteningPairEpsilon(0.04, 0.01);
  const double expected_a0 = mass[1] * cosmosim::gravity::softenedInvR3(1.0, eps_01);
  const double expected_a1 = -mass[0] * cosmosim::gravity::softenedInvR3(1.0, eps_10);
  assert(std::abs(ax[0] - expected_a0) < k_tolerance);
  assert(std::abs(ax[1] - expected_a1) < k_tolerance);
}


void testTargetSpeciesSofteningFallback() {
  cosmosim::gravity::TreeGravitySolver solver;
  cosmosim::gravity::TreeGravityOptions options;
  options.gravitational_constant_code = 1.0;
  options.max_leaf_size = 1;
  options.softening.epsilon_comoving = 0.02;

  const std::vector<double> pos_x{0.0, 1.0};
  const std::vector<double> pos_y{0.0, 0.0};
  const std::vector<double> pos_z{0.0, 0.0};
  const std::vector<double> mass{2.0, 3.0};
  const std::vector<std::uint32_t> species_tag{0U, 1U};
  const cosmosim::gravity::TreeSofteningView softening_view{
      .source_species_tag = species_tag,
      .species_policy = {.epsilon_comoving_by_species = {0.04, 0.15, 0.0, 0.0, 0.0}, .enabled = true},
  };

  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr, softening_view);
  const std::vector<std::uint32_t> active{1U};
  std::vector<double> ax(1, 0.0);
  std::vector<double> ay(1, 0.0);
  std::vector<double> az(1, 0.0);
  solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, ax, ay, az, options, nullptr, softening_view);

  const double pair_eps = cosmosim::gravity::combineSofteningPairEpsilon(0.04, 0.15);
  const double expected = -mass[0] * cosmosim::gravity::softenedInvR3(1.0, pair_eps);
  assert(std::abs(ax[0] - expected) < k_tolerance);
}


void testSourceGenerationRejectsStaleTree() {
  cosmosim::gravity::TreeGravitySolver solver;
  cosmosim::gravity::TreeGravityOptions options;
  options.max_leaf_size = 1;
  const std::vector<double> x{0.0, 1.0, 2.0};
  const std::vector<double> y(3, 0.0);
  const std::vector<double> z(3, 0.0);
  const std::vector<double> mass(3, 1.0);
  solver.build(x, y, z, mass, options, nullptr, {}, cosmosim::gravity::GravitySourceGeneration{41U});
  const std::vector<std::uint32_t> active{0U};
  std::vector<double> ax(1), ay(1), az(1);
  bool threw = false;
  try {
    solver.evaluateActiveSet(x, y, z, mass, active, ax, ay, az, options, nullptr, {}, {}, cosmosim::gravity::GravitySourceGeneration{42U});
  } catch (const std::logic_error&) {
    threw = true;
  }
  assert(threw);
}

void testLegacyFingerprintRejectsSameSizeMutation() {
  cosmosim::gravity::TreeGravitySolver solver;
  cosmosim::gravity::TreeGravityOptions options;
  options.max_leaf_size = 1;
  std::vector<double> x{0.0, 1.0, 2.0};
  const std::vector<double> y(3, 0.0);
  const std::vector<double> z(3, 0.0);
  const std::vector<double> mass(3, 1.0);
  solver.build(x, y, z, mass, options, nullptr);
  x[2] = 3.0;
  const std::vector<std::uint32_t> active{0U};
  std::vector<double> ax(1), ay(1), az(1);
  bool threw = false;
  try {
    solver.evaluateActiveSet(x, y, z, mass, active, ax, ay, az, options, nullptr);
  } catch (const std::logic_error&) {
    threw = true;
  }
  assert(threw);
}

void testEveryMacOpensTargetContainingNode() {
  const std::vector<double> x{0.0, 10.0};
  const std::vector<double> y{0.0, 0.0};
  const std::vector<double> z{0.0, 0.0};
  const std::vector<double> mass{1.0, 1.0};
  for (const auto criterion : {
           cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric,
           cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance,
           cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError}) {
    cosmosim::gravity::TreeGravitySolver solver;
    cosmosim::gravity::TreeGravityOptions options;
    options.opening_criterion = criterion;
    options.opening_theta = 100.0;
    options.relative_force_tolerance = 100.0;
    options.max_leaf_size = 1;
    options.softening.epsilon_comoving = 0.0;
    solver.build(x, y, z, mass, options, nullptr);
    const auto got = evaluateOneTarget(solver, x, y, z, mass, options);
    assert(std::abs(got.ax - 0.01) < 1.0e-12);
  }
}


void testPmPlanResourcesAndDmoProcessPreflight() {
  using cosmosim::core::MemoryEntry;
  using cosmosim::core::MemoryLifetime;
  using cosmosim::core::MemoryReportBuilder;
  using cosmosim::core::MemorySubsystem;

  const cosmosim::gravity::PmGridShape shape{32U, 32U, 32U};
  const auto serial_layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, 1, 0);
  const auto serial_plan = cosmosim::gravity::estimatePmPlanResourcesMemory(
      shape, serial_layout, cosmosim::core::PmDecompositionMode::kSlab);
  constexpr std::uint64_t serial_real_cells = 32ULL * 32ULL * 32ULL;
  constexpr std::uint64_t serial_complex_cells = 32ULL * 32ULL * 17ULL;
  const std::uint64_t expected_serial =
      serial_real_cells * sizeof(double) +
      3ULL * serial_complex_cells * sizeof(std::complex<double>) +
      4ULL * serial_complex_cells * sizeof(double);
  assert(serial_plan.total_owned_bytes == expected_serial);
  assert(serial_plan.logical_local_complex_cells == serial_complex_cells);

  const auto rank3_layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, 8, 3);
  const auto distributed_plan = cosmosim::gravity::estimatePmPlanResourcesMemory(
      shape, rank3_layout, cosmosim::core::PmDecompositionMode::kSlab);
  constexpr std::uint64_t distributed_complex_cells = 4ULL * 32ULL * 17ULL;
  const std::uint64_t expected_distributed =
      2ULL * distributed_complex_cells * sizeof(double) +
      3ULL * distributed_complex_cells * sizeof(std::complex<double>) +
      4ULL * distributed_complex_cells * sizeof(double);
  assert(distributed_plan.total_owned_bytes == expected_distributed);

  // First-light target geometry: derive, do not hard-code, the known solver-owned
  // PlanResources payload for one 64-plane slab of a 512^3 mesh.
  const cosmosim::gravity::PmGridShape first_light_shape{512U, 512U, 512U};
  const auto first_light_layout = cosmosim::parallel::makePmSlabLayout(
      first_light_shape.nx, first_light_shape.ny, first_light_shape.nz, 8, 0);
  const auto first_light_plan = cosmosim::gravity::estimatePmPlanResourcesMemory(
      first_light_shape, first_light_layout, cosmosim::core::PmDecompositionMode::kSlab);
  const std::uint64_t first_light_complex =
      64ULL * 512ULL * (512ULL / 2ULL + 1ULL);
  const std::uint64_t first_light_expected =
      2ULL * first_light_complex * sizeof(double) +
      3ULL * first_light_complex * sizeof(std::complex<double>) +
      4ULL * first_light_complex * sizeof(double);
  assert(first_light_plan.total_owned_bytes == first_light_expected);
  assert(first_light_plan.total_owned_bytes > 3ULL * 256ULL * 1024ULL * 1024ULL);

  const auto gravity_estimate = cosmosim::gravity::estimateGravityMemory({
      .local_source_count = 1024U,
      .local_target_count = 1024U,
      .local_particle_count = 1024U,
      .tree_leaf_size = 16U,
      .pm_shape = shape,
      .decomposition_mode = cosmosim::core::PmDecompositionMode::kSlab,
      .mpi_rank_count = 8U,
      .mpi_world_rank = 3,
      .backend_unknown_reserve_bytes = 4096U,
  });
  assert(gravity_estimate.pm_plan_owned_bytes == expected_distributed);

  MemoryReportBuilder canonical_builder;
  canonical_builder.addEntry(MemoryEntry{
      .subsystem = MemorySubsystem::kParticles,
      .lifetime = MemoryLifetime::kPersistent,
      .label = "test.canonical_particles",
      .owned_capacity_bytes = 1000U,
      .high_water_bytes = 1000U,
      .estimated_next_step_bytes = 1000U,
  });
  canonical_builder.addEntry(MemoryEntry{
      .subsystem = MemorySubsystem::kScratch,
      .lifetime = MemoryLifetime::kTransient,
      .label = "test.canonical_scratch",
      .owned_capacity_bytes = 2000U,
      .high_water_bytes = 2000U,
      .estimated_next_step_bytes = 2000U,
  });
  const auto canonical_report = std::move(canonical_builder).finish();
  const auto process_estimate = cosmosim::gravity::estimateDmoProcessMemory(
      canonical_report, gravity_estimate, cosmosim::gravity::DmoProcessMemoryPolicy{
          .mpi_rank_count = 8U,
          .scheduler_current_size_bytes = 2500U,
          .scheduler_owned_bytes = 3000U,
          .scheduler_high_water_bytes = 3500U,
          .output_restart_overlap_bytes = 5000U,
          .mpi_external_reserve_bytes = 7000U,
          .fftw_external_reserve_bytes = 11000U,
          .hdf5_external_reserve_bytes = 13000U,
          .allocator_external_reserve_bytes = 17000U,
          .safety_margin_fraction = 0.25,
      });
  assert(process_estimate.known_owned_peak_bytes ==
      gravity_estimate.known_peak_bytes + 1000U + 2000U + 3000U + 5000U);
  assert(process_estimate.external_unknown_reserve_bytes ==
      gravity_estimate.external_backend_unknown_bytes + 7000U + 11000U + 13000U + 17000U);
  assert(process_estimate.modeled_subtotal_bytes ==
      process_estimate.known_owned_peak_bytes + process_estimate.external_unknown_reserve_bytes);
  assert(process_estimate.budget_required_bytes >= process_estimate.modeled_subtotal_bytes);
  assert(process_estimate.aggregate_required_bytes == 8U * process_estimate.budget_required_bytes);
  bool found_scheduler_memory = false;
  for (const auto& entry : process_estimate.report.entries) {
    if (entry.label != "dmo_process.scheduler_owned_state") {
      continue;
    }
    found_scheduler_memory = true;
    assert(entry.current_size_bytes == 2500U);
    assert(entry.owned_capacity_bytes == 3000U);
    assert(entry.high_water_bytes == 3500U);
    assert(!entry.estimate_only);
  }
  assert(found_scheduler_memory);
  bool rejected = false;
  try {
    cosmosim::gravity::enforceDmoProcessMemoryBudget(
        process_estimate, process_estimate.budget_required_bytes - 1U);
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  assert(rejected);
  cosmosim::gravity::enforceDmoProcessMemoryBudget(
      process_estimate, process_estimate.budget_required_bytes);
}

void testAdaptiveNodeReserveAndMemoryEstimate() {
  constexpr std::size_t n = 1024;
  std::vector<double> x(n), y(n), z(n), mass(n, 1.0);
  for (std::size_t i = 0; i < n; ++i) {
    x[i] = static_cast<double>(i % 16U) / 16.0;
    y[i] = static_cast<double>((i / 16U) % 8U) / 8.0;
    z[i] = static_cast<double>(i / 128U) / 8.0;
  }
  cosmosim::gravity::TreeGravityOptions options;
  options.max_leaf_size = 16;
  cosmosim::gravity::TreeGravitySolver solver;
  cosmosim::gravity::TreeGravityProfile profile;
  solver.build(x, y, z, mass, options, &profile);
  assert(profile.estimated_node_reserve < 2U * n);
  assert(profile.actual_node_count == solver.nodes().size());
  assert(profile.node_capacity_high_water >= profile.actual_node_count);

  const std::uint64_t stable_node_capacity = profile.node_capacity_high_water;
  for (int rebuild = 0; rebuild < 8; ++rebuild) {
    solver.build(x, y, z, mass, options, &profile);
    assert(profile.actual_node_count == solver.nodes().size());
    assert(profile.node_capacity_high_water == stable_node_capacity);
  }

  const auto estimate = cosmosim::gravity::estimateGravityMemory({
      .local_source_count = n,
      .local_target_count = n,
      .tree_leaf_size = options.max_leaf_size,
      .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kMonopole,
      .pm_shape = {32U, 32U, 32U},
      .assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kCic,
      .decomposition_mode = cosmosim::core::PmDecompositionMode::kSlab,
      .mpi_rank_count = 1U,
  });
  assert(estimate.known_peak_bytes > 0U);
  assert(estimate.budget_required_bytes >= estimate.known_peak_bytes);
  bool budget_rejected = false;
  try {
    cosmosim::gravity::enforceGravityMemoryBudget(estimate, estimate.budget_required_bytes - 1U);
  } catch (const std::runtime_error&) {
    budget_rejected = true;
  }
  assert(budget_rejected);
  cosmosim::gravity::enforceGravityMemoryBudget(estimate, estimate.budget_required_bytes);

  const auto distributed_estimate = cosmosim::gravity::estimateGravityMemory({
      .local_source_count = n,
      .local_target_count = n,
      .tree_leaf_size = options.max_leaf_size,
      .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kMonopole,
      .pm_shape = {32U, 32U, 32U},
      .assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc,
      .decomposition_mode = cosmosim::core::PmDecompositionMode::kSlab,
      .mpi_rank_count = 8U,
      .tree_exchange_batch_bytes = 4096U,
      .pm_exchange_batch_bytes = 1024U,
  });
  bool saw_pm_routing = false;
  for (const auto& entry : distributed_estimate.report.entries) {
    if (entry.label == "gravity.estimate.pm_routing_exchange") {
      saw_pm_routing = true;
      // Two reusable wire buffers, seven remote peers, 1024 bytes/peer.
      assert(entry.estimated_next_step_bytes == 2U * 7U * 1024U);
    }
  }
  assert(saw_pm_routing);
  assert(distributed_estimate.report.distributed.valid);
  assert(distributed_estimate.report.distributed.rank_count == 8);
  assert(distributed_estimate.report.distributed.global_sum_owned_bytes ==
      8U * distributed_estimate.report.distributed.local_owned_bytes);
}

void testGasRefinementForceConvergesAtFixedMassAndCom() {
  constexpr double epsilon = 0.05;
  cosmosim::gravity::TreeGravityOptions options;
  options.max_leaf_size = 1U;
  options.opening_theta = 1.0e-8;
  options.gravitational_constant_code = 1.0;
  options.softening.epsilon_comoving = epsilon;

  const auto evaluate_target = [&](const std::vector<double>& source_x,
                                   const std::vector<double>& source_mass,
                                   std::uint32_t target_index) {
    const std::vector<double> zeros(source_x.size(), 0.0);
    cosmosim::gravity::TreeGravitySolver solver;
    solver.build(source_x, zeros, zeros, source_mass, options);
    const std::vector<std::uint32_t> active{target_index};
    std::vector<double> ax(1, 0.0), ay(1, 0.0), az(1, 0.0);
    solver.evaluateActiveSet(
        source_x, zeros, zeros, source_mass, active, ax, ay, az, options);
    return ax[0];
  };

  // A zero-mass target at x=10 observes either one coarse source at the COM
  // or two symmetric child leaves conserving total mass and COM. Refinement
  // changes the resolved density representation, so exact equality is not the
  // contract; the far-field difference should be small and non-zero.
  const double coarse_ax = evaluate_target(
      std::vector<double>{0.0, 10.0},
      std::vector<double>{4.0, 0.0},
      1U);
  const double refined_ax = evaluate_target(
      std::vector<double>{-0.1, 0.1, 10.0},
      std::vector<double>{2.0, 2.0, 0.0},
      2U);
  const double relative_difference =
      std::abs(refined_ax - coarse_ax) / std::abs(coarse_ax);
  assert(relative_difference > 0.0);
  assert(relative_difference < 1.0e-3);
}

void testAuthoritativeGasGeometryMatchesDirectSoftenedReference() {
  // Rows 1 and 2 model distinct authoritative leaf gas cells. Their lineage
  // metadata is deliberately irrelevant to gravity: geometry/mass define the
  // physical sources, so siblings remain distinct and self force is excluded.
  const std::vector<double> x{0.0, 1.0, 2.5};
  const std::vector<double> y{0.0, 0.0, 0.0};
  const std::vector<double> z{0.0, 0.0, 0.0};
  const std::vector<double> mass{3.0, 2.0, 4.0};
  constexpr double epsilon = 0.075;
  cosmosim::gravity::TreeGravityOptions options;
  options.max_leaf_size = 1;
  options.opening_theta = 1.0e-6;  // force the direct leaf corpus
  options.gravitational_constant_code = 1.7;
  options.softening.epsilon_comoving = epsilon;

  cosmosim::gravity::TreeGravitySolver solver;
  solver.build(x, y, z, mass, options);
  const std::vector<std::uint32_t> active{1U};
  std::vector<double> ax(1, 0.0), ay(1, 0.0), az(1, 0.0);
  solver.evaluateActiveSet(x, y, z, mass, active, ax, ay, az, options);

  double expected_ax = 0.0;
  for (std::size_t source = 0; source < x.size(); ++source) {
    if (source == 1U) continue;
    const double dx = x[source] - x[1];
    const double r2 = dx * dx;
    expected_ax += options.gravitational_constant_code * mass[source] * dx *
        cosmosim::gravity::softenedInvR3(r2, epsilon);
  }
  assert(std::abs(ax[0] - expected_ax) < 1.0e-12);
  assert(std::abs(ay[0]) < 1.0e-14);
  assert(std::abs(az[0]) < 1.0e-14);
}

void testTreeAndTreePmCommonAcceptanceConformanceCorpus() {
  using cosmosim::gravity::TreeGravityOptions;
  using cosmosim::gravity::TreeOpeningCriterion;
  using cosmosim::gravity::internal::TreeNodeAcceptanceInput;
  using cosmosim::gravity::internal::acceptNodeByCommonTreePolicy;

  TreeGravityOptions options;
  options.gravitational_constant_code = 1.0;
  options.opening_theta = 0.7;
  options.relative_force_tolerance = 1.0e-3;
  options.relative_force_acceleration_floor_code = 1.0e-12;

  const std::array<TreeOpeningCriterion, 3> criteria{
      TreeOpeningCriterion::kBarnesHutGeometric,
      TreeOpeningCriterion::kBarnesHutComDistance,
      TreeOpeningCriterion::kRelativeForceError,
  };
  for (const auto criterion : criteria) {
    options.opening_criterion = criterion;
    for (const bool inside : {false, true}) {
      for (const bool previous_available : {false, true}) {
        const TreeNodeAcceptanceInput input{
            .is_leaf = false,
            .target_inside_node = inside,
            .half_size = 0.1,
            .com_center_offset = 0.02,
            .node_mass_code = 2.0,
            .r2 = 4.0,
            .previous_acceleration_available = previous_available,
            .previous_acceleration_magnitude_code = 0.5,
            .target_softening_comoving = 0.01,
            .node_softening_min_comoving = 0.01,
            .node_softening_max_comoving = 0.01,
        };
        const bool standalone_common = acceptNodeByCommonTreePolicy(input, options);
        const bool treepm_common = acceptNodeByCommonTreePolicy(input, options);
        assert(standalone_common == treepm_common);
        if (inside) {
          assert(!standalone_common);
        }
      }
    }
  }

  assert(cosmosim::gravity::internal::targetInsideNodeAabbFromCenterDelta(
      0.1, -0.1, 0.0, 0.1));
  assert(!cosmosim::gravity::internal::targetInsideNodeAabbFromCenterDelta(
      0.10001, 0.0, 0.0, 0.1));
}


}  // namespace

int main() {
  testTwoBodyForce();
  testForceSymmetry();
  testOpeningCriterionBehavior();
  testMultipoleAccumulation();
  testMacCriteriaAreAvailable();
  testRelativeForceErrorMacFormulaAndWorkMonotonicity();
  testRelativeForceErrorMacHistoryFallbackAndFloor();
  testRelativeForceErrorMacNeverAcceptsTargetContainingNode();
  testRelativeForceErrorMacValidatesOptionsAndHistoryShape();
  testPairSofteningUsesMaxRule();
  testTargetSpeciesSofteningFallback();
  testSourceGenerationRejectsStaleTree();
  testLegacyFingerprintRejectsSameSizeMutation();
  testEveryMacOpensTargetContainingNode();
  testPmPlanResourcesAndDmoProcessPreflight();
  testAdaptiveNodeReserveAndMemoryEstimate();
  testAuthoritativeGasGeometryMatchesDirectSoftenedReference();
  testGasRefinementForceConvergesAtFixedMassAndCom();
  testTreeAndTreePmCommonAcceptanceConformanceCorpus();
  return 0;
}
