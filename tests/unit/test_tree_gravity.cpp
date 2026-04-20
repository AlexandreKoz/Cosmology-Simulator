#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <vector>

#include "cosmosim/gravity/tree_gravity.hpp"

namespace {

constexpr double k_tolerance = 1.0e-10;

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
  cosmosim::gravity::TreeSofteningView softening_view{
      .source_species_tag = species_tag,
      .target_particle_epsilon_comoving = target_eps,
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

}  // namespace

int main() {
  testTwoBodyForce();
  testForceSymmetry();
  testOpeningCriterionBehavior();
  testMultipoleAccumulation();
  testMacCriteriaAreAvailable();
  testPairSofteningUsesMaxRule();
  return 0;
}
