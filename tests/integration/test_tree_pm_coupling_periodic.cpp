#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"

namespace {

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

[[nodiscard]] double minimumImageDelta(double delta, double box_size_comoving) {
  return delta - box_size_comoving * std::nearbyint(delta / box_size_comoving);
}

void testPeriodicTreePmAgainstDirectReference() {
  constexpr std::size_t particle_count = 48;
  constexpr double box_size_comoving = 1.0;

  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 1.0);

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod((17.0 * static_cast<double>(i) + 3.0) * 0.017, box_size_comoving);
    pos_y[i] = std::fmod((29.0 * static_cast<double>(i) + 5.0) * 0.013, box_size_comoving);
    pos_z[i] = std::fmod((41.0 * static_cast<double>(i) + 7.0) * 0.011, box_size_comoving);
    mass[i] = 0.8 + 0.01 * static_cast<double>(i % 11U);
  }

  std::vector<std::uint32_t> active(particle_count, 0U);
  for (std::size_t i = 0; i < particle_count; ++i) {
    active[i] = static_cast<std::uint32_t>(i);
  }

  std::vector<double> tree_pm_ax(particle_count, 0.0);
  std::vector<double> tree_pm_ay(particle_count, 0.0);
  std::vector<double> tree_pm_az(particle_count, 0.0);

  cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      active,
      tree_pm_ax,
      tree_pm_ay,
      tree_pm_az,
  };

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = box_size_comoving;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.enable_window_deconvolution = true;
  options.tree_options.opening_theta = 0.55;
  options.tree_options.max_leaf_size = 8;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 0.01;
  options.split_policy.split_scale_comoving = 0.08;
  options.split_policy.cutoff_radius_comoving = 0.4;

  const cosmosim::gravity::PmGridShape pm_shape =
#if COSMOSIM_ENABLE_FFTW
      {32, 32, 32};
#else
      {8, 8, 8};
#endif
  cosmosim::gravity::TreePmCoordinator coordinator(pm_shape);
  cosmosim::gravity::TreePmDiagnostics diagnostics;
  coordinator.solveActiveSet(pos_x, pos_y, pos_z, mass, accumulator, options, nullptr, &diagnostics);

  std::vector<double> ref_ax(particle_count, 0.0);
  std::vector<double> ref_ay(particle_count, 0.0);
  std::vector<double> ref_az(particle_count, 0.0);

  for (std::size_t target = 0; target < particle_count; ++target) {
    for (std::size_t source = 0; source < particle_count; ++source) {
      if (target == source) {
        continue;
      }
      const double dx = minimumImageDelta(pos_x[source] - pos_x[target], box_size_comoving);
      const double dy = minimumImageDelta(pos_y[source] - pos_y[target], box_size_comoving);
      const double dz = minimumImageDelta(pos_z[source] - pos_z[target], box_size_comoving);
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double softened = cosmosim::gravity::softenedInvR3(r2, options.tree_options.softening);
      const double factor = options.tree_options.gravitational_constant_code * mass[source] * softened;
      ref_ax[target] += factor * dx;
      ref_ay[target] += factor * dy;
      ref_az[target] += factor * dz;
    }
  }

  double ref_norm = 0.0;
  double err_norm = 0.0;
  for (std::size_t i = 0; i < particle_count; ++i) {
    const double dax = tree_pm_ax[i] - ref_ax[i];
    const double day = tree_pm_ay[i] - ref_ay[i];
    const double daz = tree_pm_az[i] - ref_az[i];
    err_norm += dax * dax + day * day + daz * daz;
    ref_norm += ref_ax[i] * ref_ax[i] + ref_ay[i] * ref_ay[i] + ref_az[i] * ref_az[i];
  }

  const double rel_l2 = std::sqrt(err_norm / std::max(ref_norm, 1.0e-30));
#if COSMOSIM_ENABLE_FFTW
  const double max_rel_l2 = 0.75;
#else
  const double max_rel_l2 = 1.8;
#endif

  std::ostringstream diag;
  diag << "TreePM periodic validation failed: build_flag.COSMOSIM_ENABLE_FFTW="
       << (COSMOSIM_ENABLE_FFTW ? "ON" : "OFF")
       << ", treePmSupportedByBuild()=" << (cosmosim::gravity::treePmSupportedByBuild() ? "true" : "false")
       << ", rel_l2=" << rel_l2 << " (required <= " << max_rel_l2 << ')'
       << ", composition_error_at_split=" << diagnostics.composition_error_at_split
       << ", max_relative_composition_error=" << diagnostics.max_relative_composition_error;

  requireOrThrow(rel_l2 < max_rel_l2, diag.str());
  requireOrThrow(diagnostics.composition_error_at_split < 1.0e-12, diag.str());
  requireOrThrow(diagnostics.max_relative_composition_error < 1.0e-12, diag.str());
}

}  // namespace

int main() {
  testPeriodicTreePmAgainstDirectReference();
  return 0;
}
