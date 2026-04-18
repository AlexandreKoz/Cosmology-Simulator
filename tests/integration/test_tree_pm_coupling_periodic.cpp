#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
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

struct ForceField {
  std::vector<double> ax;
  std::vector<double> ay;
  std::vector<double> az;
};

[[nodiscard]] ForceField computeDirectReference(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    double box_size_comoving,
    const cosmosim::gravity::TreeSofteningPolicy& softening,
    double gravitational_constant_code) {
  ForceField field{
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
  };

  for (std::size_t target = 0; target < pos_x.size(); ++target) {
    for (std::size_t source = 0; source < pos_x.size(); ++source) {
      if (target == source) {
        continue;
      }
      const double dx = minimumImageDelta(pos_x[source] - pos_x[target], box_size_comoving);
      const double dy = minimumImageDelta(pos_y[source] - pos_y[target], box_size_comoving);
      const double dz = minimumImageDelta(pos_z[source] - pos_z[target], box_size_comoving);
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double softened = cosmosim::gravity::softenedInvR3(r2, softening);
      const double factor = gravitational_constant_code * mass[source] * softened;
      field.ax[target] += factor * dx;
      field.ay[target] += factor * dy;
      field.az[target] += factor * dz;
    }
  }

  return field;
}

[[nodiscard]] double relativeL2Error(const ForceField& value, const ForceField& reference) {
  double ref_norm = 0.0;
  double err_norm = 0.0;
  for (std::size_t i = 0; i < value.ax.size(); ++i) {
    const double dax = value.ax[i] - reference.ax[i];
    const double day = value.ay[i] - reference.ay[i];
    const double daz = value.az[i] - reference.az[i];
    err_norm += dax * dax + day * day + daz * daz;
    ref_norm += reference.ax[i] * reference.ax[i] + reference.ay[i] * reference.ay[i] + reference.az[i] * reference.az[i];
  }
  return std::sqrt(err_norm / std::max(ref_norm, 1.0e-30));
}

[[nodiscard]] ForceField solveTreePm(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    const cosmosim::gravity::PmGridShape& pm_shape,
    cosmosim::gravity::TreePmOptions options,
    cosmosim::gravity::TreePmDiagnostics* diagnostics,
    cosmosim::gravity::TreePmProfileEvent* profile = nullptr) {
  std::vector<std::uint32_t> active(pos_x.size(), 0U);
  for (std::size_t i = 0; i < pos_x.size(); ++i) {
    active[i] = static_cast<std::uint32_t>(i);
  }

  ForceField field{
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
  };

  cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      active,
      field.ax,
      field.ay,
      field.az,
  };

  cosmosim::gravity::TreePmCoordinator coordinator(pm_shape);
  coordinator.solveActiveSet(pos_x, pos_y, pos_z, mass, accumulator, options, profile, diagnostics);
  return field;
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

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = box_size_comoving;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.enable_window_deconvolution = true;
  options.tree_options.opening_theta = 0.55;
  options.tree_options.max_leaf_size = 8;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 0.01;

  const cosmosim::gravity::PmGridShape pm_shape =
#if COSMOSIM_ENABLE_FFTW
      {32, 32, 32};
#else
      {8, 8, 8};
#endif
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 4.5, mesh_spacing);

  cosmosim::gravity::TreePmDiagnostics diagnostics;
  const ForceField tree_pm_force = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, options, &diagnostics);

  const ForceField reference = computeDirectReference(
      pos_x,
      pos_y,
      pos_z,
      mass,
      box_size_comoving,
      options.tree_options.softening,
      options.tree_options.gravitational_constant_code);

  const double rel_l2 = relativeL2Error(tree_pm_force, reference);
#if COSMOSIM_ENABLE_FFTW
  const double max_rel_l2 = 0.75;
#else
  const double max_rel_l2 = 1.8;
#endif

  std::ostringstream diag;
  diag << "TreePM periodic validation failed (direct reference is minimum-image periodic, not Ewald): build_flag.COSMOSIM_ENABLE_FFTW="
       << (COSMOSIM_ENABLE_FFTW ? "ON" : "OFF")
       << ", treePmSupportedByBuild()=" << (cosmosim::gravity::treePmSupportedByBuild() ? "true" : "false")
       << ", rel_l2=" << rel_l2 << " (required <= " << max_rel_l2 << ')'
       << ", pm_grid=" << pm_shape.nx
       << ", mesh_spacing=" << diagnostics.mesh_spacing_comoving
       << ", asmth_cells=" << diagnostics.asmth_cells
       << ", rcut_cells=" << diagnostics.rcut_cells
       << ", split_scale_comoving=" << diagnostics.split_scale_comoving
       << ", cutoff_radius_comoving=" << diagnostics.cutoff_radius_comoving
       << ", residual_pruned_nodes=" << diagnostics.residual_pruned_nodes
       << ", residual_pair_skips_cutoff=" << diagnostics.residual_pair_skips_cutoff
       << ", composition_error_at_split=" << diagnostics.composition_error_at_split
       << ", max_relative_composition_error=" << diagnostics.max_relative_composition_error;

  requireOrThrow(rel_l2 < max_rel_l2, diag.str());
  requireOrThrow(diagnostics.composition_error_at_split < 1.0e-12, diag.str());
  requireOrThrow(diagnostics.max_relative_composition_error < 1.0e-12, diag.str());
  requireOrThrow(diagnostics.residual_pruned_nodes > 0, diag.str());
  requireOrThrow(diagnostics.residual_pair_skips_cutoff > 0, diag.str());
}

void testResidualCutoffPrunesTreePath() {
  constexpr double box_size_comoving = 1.0;
  std::vector<double> pos_x = {0.10, 0.18, 0.62, 0.88};
  std::vector<double> pos_y = {0.10, 0.12, 0.60, 0.90};
  std::vector<double> pos_z = {0.10, 0.09, 0.58, 0.86};
  std::vector<double> mass = {1.0, 1.0, 1.0, 1.0};

  const cosmosim::gravity::PmGridShape pm_shape{32, 32, 32};
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = box_size_comoving;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.opening_theta = 0.55;
  options.tree_options.max_leaf_size = 2;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.0e-3, 3.0, mesh_spacing);

  cosmosim::gravity::TreePmDiagnostics diagnostics;
  [[maybe_unused]] const ForceField field = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, options, &diagnostics);

  std::ostringstream msg;
  msg << "Residual cutoff pruning counters invalid: pruned_nodes=" << diagnostics.residual_pruned_nodes
      << ", pair_skips=" << diagnostics.residual_pair_skips_cutoff
      << ", pair_evaluations=" << diagnostics.residual_pair_evaluations;

  requireOrThrow(diagnostics.residual_pruned_nodes > 0, msg.str());
  requireOrThrow(diagnostics.residual_pair_skips_cutoff > 0, msg.str());
  requireOrThrow(diagnostics.residual_pair_evaluations > 0, msg.str());
}

void testPmOnlyTreeOnlyAndTreePmConsistency() {
  constexpr std::size_t particle_count = 32;
  constexpr double box_size_comoving = 1.0;

  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 1.0);

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod((13.0 * static_cast<double>(i) + 1.0) * 0.023, box_size_comoving);
    pos_y[i] = std::fmod((19.0 * static_cast<double>(i) + 2.0) * 0.019, box_size_comoving);
    pos_z[i] = std::fmod((23.0 * static_cast<double>(i) + 3.0) * 0.017, box_size_comoving);
    mass[i] = 0.9 + 0.05 * static_cast<double>(i % 5U);
  }

  const cosmosim::gravity::PmGridShape pm_shape{32, 32, 32};
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);

  cosmosim::gravity::TreePmOptions base_options;
  base_options.pm_options.box_size_mpc_comoving = box_size_comoving;
  base_options.pm_options.scale_factor = 1.0;
  base_options.pm_options.gravitational_constant_code = 1.0;
  base_options.pm_options.enable_window_deconvolution = true;
  base_options.tree_options.opening_theta = 0.55;
  base_options.tree_options.max_leaf_size = 8;
  base_options.tree_options.gravitational_constant_code = 1.0;
  base_options.tree_options.softening.epsilon_comoving = 0.01;

  const ForceField direct = computeDirectReference(
      pos_x,
      pos_y,
      pos_z,
      mass,
      box_size_comoving,
      base_options.tree_options.softening,
      base_options.tree_options.gravitational_constant_code);

  cosmosim::gravity::TreePmOptions tree_only_options = base_options;
  tree_only_options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(100.0, 40.0, mesh_spacing);
  cosmosim::gravity::TreePmDiagnostics tree_only_diagnostics;
  const ForceField tree_only = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, tree_only_options, &tree_only_diagnostics);

  cosmosim::gravity::TreePmOptions pm_only_options = base_options;
  pm_only_options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(0.2, 4.5, mesh_spacing);
  cosmosim::gravity::TreePmDiagnostics pm_only_diagnostics;
  const ForceField pm_only = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, pm_only_options, &pm_only_diagnostics);

  cosmosim::gravity::TreePmOptions split_options = base_options;
  split_options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 4.5, mesh_spacing);
  cosmosim::gravity::TreePmDiagnostics split_diagnostics;
  const ForceField split = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, split_options, &split_diagnostics);

  const double tree_only_rel = relativeL2Error(tree_only, direct);
  const double pm_only_rel = relativeL2Error(pm_only, direct);
  const double split_rel = relativeL2Error(split, direct);

  std::ostringstream message;
  message << "TreePM consistency failure: tree_only_rel=" << tree_only_rel << ", split_rel=" << split_rel
          << ", pm_only_rel=" << pm_only_rel << ", split_pruned_nodes=" << split_diagnostics.residual_pruned_nodes
          << ", split_pair_skips=" << split_diagnostics.residual_pair_skips_cutoff;

  requireOrThrow(tree_only_rel < 0.35, message.str());
  requireOrThrow(split_rel < 0.90, message.str());
  requireOrThrow(pm_only_rel < 2.2, message.str());
  requireOrThrow(tree_only_rel <= split_rel + 1.0e-9, message.str());
  requireOrThrow(split_rel <= pm_only_rel + 1.0e-9, message.str());
}

}  // namespace

int main() {
  testPeriodicTreePmAgainstDirectReference();
  testResidualCutoffPrunesTreePath();
  testPmOnlyTreeOnlyAndTreePmConsistency();
  return 0;
}
