#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/gravity/tree_pm_split_kernel.hpp"
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

void directSumAcceleration(
    const std::vector<double>& pos_x,
    const std::vector<double>& pos_y,
    const std::vector<double>& pos_z,
    const std::vector<double>& mass,
    const cosmosim::gravity::TreeGravityOptions& options,
    std::vector<double>& ax,
    std::vector<double>& ay,
    std::vector<double>& az) {
  for (std::size_t i = 0; i < pos_x.size(); ++i) {
    double acc_x = 0.0;
    double acc_y = 0.0;
    double acc_z = 0.0;
    for (std::size_t j = 0; j < pos_x.size(); ++j) {
      if (i == j) {
        continue;
      }
      const double dx = pos_x[j] - pos_x[i];
      const double dy = pos_y[j] - pos_y[i];
      const double dz = pos_z[j] - pos_z[i];
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

void testSmallNTreeVsDirectOpenBoundary(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  constexpr std::size_t particle_count = 24;
  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 0.0);
  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod((17.0 * static_cast<double>(i) + 3.0) * 0.021, 0.97) - 0.48;
    pos_y[i] = std::fmod((31.0 * static_cast<double>(i) + 5.0) * 0.017, 0.95) - 0.47;
    pos_z[i] = std::fmod((43.0 * static_cast<double>(i) + 7.0) * 0.013, 0.93) - 0.46;
    mass[i] = 0.8 + 0.05 * static_cast<double>(i % 7U);
  }

  cosmosim::gravity::TreeGravityOptions options;
  options.opening_theta = 0.45;
  options.gravitational_constant_code = 1.0;
  options.max_leaf_size = 4;
  options.softening.epsilon_comoving = 1.0e-3;

  std::vector<std::uint32_t> active(particle_count, 0U);
  for (std::size_t i = 0; i < particle_count; ++i) {
    active[i] = static_cast<std::uint32_t>(i);
  }

  std::vector<double> tree_ax(particle_count, 0.0);
  std::vector<double> tree_ay(particle_count, 0.0);
  std::vector<double> tree_az(particle_count, 0.0);

  cosmosim::gravity::TreeGravitySolver solver;
  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);
  solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, tree_ax, tree_ay, tree_az, options, nullptr);

  std::vector<double> ref_ax(particle_count, 0.0);
  std::vector<double> ref_ay(particle_count, 0.0);
  std::vector<double> ref_az(particle_count, 0.0);
  directSumAcceleration(pos_x, pos_y, pos_z, mass, options, ref_ax, ref_ay, ref_az);

  double err2 = 0.0;
  double ref2 = 0.0;
  for (std::size_t i = 0; i < particle_count; ++i) {
    const double dx = tree_ax[i] - ref_ax[i];
    const double dy = tree_ay[i] - ref_ay[i];
    const double dz = tree_az[i] - ref_az[i];
    err2 += dx * dx + dy * dy + dz * dz;
    ref2 += ref_ax[i] * ref_ax[i] + ref_ay[i] * ref_ay[i] + ref_az[i] * ref_az[i];
  }

  const double rel_l2 = std::sqrt(err2 / std::max(ref2, 1.0e-24));
  requireOrThrow(
      rel_l2 <= tolerances.require("gravity_tree_small_n_open_boundary.max_relative_l2_error"),
      "gravity_tree_small_n_open_boundary failed: tree-vs-direct relative L2 error above tolerance");
}

void testSplitKernelComposition(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  const auto split = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 4.5, 0.03125);
  double max_err = 0.0;
  for (double r : {0.01, 0.03, 0.1, 0.2, 0.4, 0.8}) {
    const double sr = cosmosim::gravity::treePmGaussianShortRangeForceFactor(r, split.split_scale_comoving);
    const double lr = cosmosim::gravity::treePmGaussianLongRangeForceFactor(r, split.split_scale_comoving);
    max_err = std::max(max_err, std::abs((sr + lr) - 1.0));
  }
  requireOrThrow(
      max_err <= tolerances.require("gravity_tree_pm_split.composition_abs_error"),
      "gravity_tree_pm_split failed: SR/LR composition exceeds tolerance");
}

void testToleranceTableCompleteness(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  requireOrThrow(tolerances.has("gravity_two_body.max_relative_accel_error"), "missing key: gravity two-body");
  requireOrThrow(tolerances.has("gravity_tree_small_n_open_boundary.max_relative_l2_error"), "missing key: small-N tree");
  requireOrThrow(tolerances.has("gravity_tree_pm_split.composition_abs_error"), "missing key: split composition");
  requireOrThrow(tolerances.has("gravity_pm_glass_like.max_rms_accel"), "missing key: glass-like pm uniformity");
  requireOrThrow(tolerances.has("hydro_sine_wave.l1_density_error_n64"), "missing key: hydro n64");
  requireOrThrow(tolerances.has("amr_reflux.mass_correction_abs"), "missing key: amr reflux");
}

}  // namespace

int main() {
  const auto tolerances = cosmosim::validation::ValidationToleranceTable::loadFromFile(
      std::string(COSMOSIM_SOURCE_DIR) + "/validation/reference/validation_tolerances_v1.txt");

  testToleranceTableCompleteness(tolerances);
  testTwoBodySymmetry(tolerances);
  testSmallNTreeVsDirectOpenBoundary(tolerances);
  testSplitKernelComposition(tolerances);
  return 0;
}
