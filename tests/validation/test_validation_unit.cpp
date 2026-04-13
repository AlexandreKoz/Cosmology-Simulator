#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/gravity/tree_gravity.hpp"
#include "validation_tolerance.hpp"

namespace {

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void testTwoBodySymmetry(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  std::vector<double> pos_x = {-0.25, 0.25};
  std::vector<double> pos_y = {0.0, 0.0};
  std::vector<double> pos_z = {0.0, 0.0};
  std::vector<double> mass = {2.0, 2.0};

  std::vector<std::uint32_t> active = {0U, 1U};
  std::vector<double> ax(2, 0.0);
  std::vector<double> ay(2, 0.0);
  std::vector<double> az(2, 0.0);

  cosmosim::gravity::TreeGravityOptions options;
  options.opening_theta = 0.2;
  options.gravitational_constant_code = 1.0;
  options.max_leaf_size = 1;
  options.softening.epsilon_comoving = 1.0e-6;

  cosmosim::gravity::TreeGravitySolver solver;
  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);
  solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, ax, ay, az, options, nullptr);

  const double dx = pos_x[1] - pos_x[0];
  const double r2 = dx * dx;
  const double expected = options.gravitational_constant_code * mass[1] * dx *
      cosmosim::gravity::softenedInvR3(r2, options.softening);

  const double rel_err_0 = std::abs(ax[0] - expected) / std::max(std::abs(expected), 1.0e-12);
  const double rel_err_1 = std::abs(ax[1] + expected) / std::max(std::abs(expected), 1.0e-12);

  const double max_rel_err = std::max(rel_err_0, rel_err_1);
  requireOrThrow(
      max_rel_err <= tolerances.require("gravity_two_body.max_relative_accel_error"),
      "gravity_two_body failed: relative acceleration error exceeds tolerance");
  requireOrThrow(std::abs(ay[0]) < 1.0e-14 && std::abs(az[0]) < 1.0e-14, "gravity_two_body transverse drift");
}

void testToleranceTableCompleteness(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  requireOrThrow(tolerances.has("gravity_two_body.max_relative_accel_error"), "missing key: gravity two-body");
  requireOrThrow(tolerances.has("hydro_sine_wave.l1_density_error_n64"), "missing key: hydro n64");
  requireOrThrow(tolerances.has("amr_reflux.mass_correction_abs"), "missing key: amr reflux");
}

}  // namespace

int main() {
  const auto tolerances = cosmosim::validation::ValidationToleranceTable::loadFromFile(
      std::string(COSMOSIM_SOURCE_DIR) + "/validation/reference/validation_tolerances_v1.txt");

  testToleranceTableCompleteness(tolerances);
  testTwoBodySymmetry(tolerances);
  return 0;
}
