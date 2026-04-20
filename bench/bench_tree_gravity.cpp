#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <iostream>
#include <vector>

#include "cosmosim/gravity/tree_gravity.hpp"

int main() {
  constexpr std::size_t particle_count = 250000;
  constexpr std::size_t active_count = 32768;

  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 1.0);
  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod(static_cast<double>((37U * i + 13U) % 104729U) * 7.13e-5, 1.0);
    pos_y[i] = std::fmod(static_cast<double>((73U * i + 17U) % 130363U) * 5.11e-5, 1.0);
    pos_z[i] = std::fmod(static_cast<double>((97U * i + 19U) % 156007U) * 4.73e-5, 1.0);
    mass[i] = 0.5 + static_cast<double>((i % 11U)) * 0.02;
  }

  std::vector<std::uint32_t> active(active_count);
  for (std::size_t i = 0; i < active_count; ++i) {
    active[i] = static_cast<std::uint32_t>((17U * i) % particle_count);
  }

  std::vector<double> accel_x(active_count, 0.0);
  std::vector<double> accel_y(active_count, 0.0);
  std::vector<double> accel_z(active_count, 0.0);

  const std::array<cosmosim::gravity::TreeOpeningCriterion, 2> criteria = {
      cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric,
      cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance};
  const std::array<double, 3> theta_schedule = {0.45, 0.6, 0.8};

  for (const auto criterion : criteria) {
    for (const double theta : theta_schedule) {
      cosmosim::gravity::TreeGravityOptions options;
      options.opening_theta = theta;
      options.opening_criterion = criterion;
      options.multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole;
      options.gravitational_constant_code = 1.0;
      options.max_leaf_size = 16;
      options.softening.epsilon_comoving = 3.0e-4;

      cosmosim::gravity::TreeGravitySolver solver;
      cosmosim::gravity::TreeGravityProfile profile;

      const auto build_start = std::chrono::steady_clock::now();
      solver.build(pos_x, pos_y, pos_z, mass, options, &profile);
      const auto build_stop = std::chrono::steady_clock::now();

      const auto traversal_start = std::chrono::steady_clock::now();
      solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, accel_x, accel_y, accel_z, options, &profile);
      const auto traversal_stop = std::chrono::steady_clock::now();

      const double build_ms = std::chrono::duration<double, std::milli>(build_stop - build_start).count();
      const double traversal_ms = std::chrono::duration<double, std::milli>(traversal_stop - traversal_start).count();
      const double build_mpart_s = static_cast<double>(particle_count) / (build_ms * 1.0e3) * 1.0e-6;
      const double active_mpart_s = static_cast<double>(active_count) / (traversal_ms * 1.0e3) * 1.0e-6;

      std::cout << "bench_tree_gravity"
                << " criterion=" << (criterion == cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric ? "geometric" : "com_distance")
                << " theta=" << theta
                << " multipole=quadrupole"
                << " particle_count=" << particle_count
                << " active_count=" << active_count
                << " build_ms=" << build_ms
                << " multipole_ms=" << profile.multipole_ms
                << " traversal_ms=" << traversal_ms
                << " visited_nodes=" << profile.visited_nodes
                << " accepted_nodes=" << profile.accepted_nodes
                << " pp_interactions=" << profile.particle_particle_interactions
                << " avg_interactions_per_target=" << profile.average_interactions_per_target
                << " build_mpart_s=" << build_mpart_s
                << " active_mpart_s=" << active_mpart_s
                << '\n';
    }
  }

  return 0;
}
