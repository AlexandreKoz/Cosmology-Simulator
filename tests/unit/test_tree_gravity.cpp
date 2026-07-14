#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include "cosmosim/gravity/tree_gravity.hpp"

namespace {

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
  return 0;
}
