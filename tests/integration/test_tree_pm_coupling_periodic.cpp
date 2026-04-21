#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

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
  requireOrThrow(
      diagnostics.residual_pruned_nodes > 0 || diagnostics.residual_pair_skips_cutoff > 0,
      diag.str());
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

void testActiveSubsetMatchesFullSolveSingleRank() {
  constexpr double box_size_comoving = 1.0;
  const cosmosim::gravity::PmGridShape pm_shape{16, 16, 16};
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);

  const std::vector<double> pos_x = {0.05, 0.12, 0.19, 0.28, 0.34, 0.43, 0.55, 0.61, 0.74, 0.88};
  const std::vector<double> pos_y = {0.08, 0.16, 0.24, 0.31, 0.37, 0.46, 0.52, 0.67, 0.79, 0.91};
  const std::vector<double> pos_z = {0.03, 0.11, 0.21, 0.29, 0.41, 0.49, 0.58, 0.69, 0.82, 0.94};
  std::vector<double> mass(pos_x.size(), 1.0);
  for (std::size_t i = 0; i < mass.size(); ++i) {
    mass[i] = 0.95 + 0.015 * static_cast<double>(i % 5U);
  }

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = box_size_comoving;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.enable_window_deconvolution = true;
  options.tree_options.opening_theta = 0.6;
  options.tree_options.max_leaf_size = 4;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.5, 4.0, mesh_spacing);

  const ForceField full_field = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, options, nullptr);

  std::vector<std::uint32_t> active_subset = {1U, 4U, 7U};
  ForceField subset_field{
      std::vector<double>(active_subset.size(), 0.0),
      std::vector<double>(active_subset.size(), 0.0),
      std::vector<double>(active_subset.size(), 0.0),
  };
  cosmosim::gravity::TreePmForceAccumulatorView subset_accumulator{
      active_subset,
      subset_field.ax,
      subset_field.ay,
      subset_field.az,
  };
  cosmosim::gravity::TreePmCoordinator subset_coordinator(pm_shape);
  subset_coordinator.solveActiveSet(pos_x, pos_y, pos_z, mass, subset_accumulator, options, nullptr, nullptr);

  for (std::size_t i = 0; i < active_subset.size(); ++i) {
    const std::uint32_t global_index = active_subset[i];
    requireOrThrow(std::abs(subset_field.ax[i] - full_field.ax[global_index]) <= 1.0e-12, "active-subset ax mismatch");
    requireOrThrow(std::abs(subset_field.ay[i] - full_field.ay[global_index]) <= 1.0e-12, "active-subset ay mismatch");
    requireOrThrow(std::abs(subset_field.az[i] - full_field.az[global_index]) <= 1.0e-12, "active-subset az mismatch");
  }
}

#if COSMOSIM_ENABLE_MPI
void testDistributedShortRangeExportImportMatchesSingleRankReference() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    return;
  }

  constexpr double box_size_comoving = 1.0;
  const cosmosim::gravity::PmGridShape pm_shape{16, 16, 16};
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);

  const std::vector<double> global_x = {
      0.08, 0.14, 0.22, 0.30, 0.41, 0.47, 0.492, 0.498,
      0.502, 0.508, 0.53, 0.61, 0.72, 0.81, 0.90, 0.96,
  };
  const std::vector<double> global_y = {
      0.11, 0.19, 0.27, 0.33, 0.39, 0.46, 0.49, 0.52,
      0.55, 0.58, 0.63, 0.67, 0.71, 0.76, 0.84, 0.93,
  };
  const std::vector<double> global_z = {
      0.07, 0.13, 0.20, 0.28, 0.35, 0.44, 0.495, 0.499,
      0.501, 0.505, 0.57, 0.64, 0.70, 0.79, 0.87, 0.95,
  };
  std::vector<double> global_mass(global_x.size(), 1.0);
  for (std::size_t i = 0; i < global_mass.size(); ++i) {
    global_mass[i] = 0.9 + 0.02 * static_cast<double>(i % 7U);
  }

  std::vector<double> local_x;
  std::vector<double> local_y;
  std::vector<double> local_z;
  std::vector<double> local_mass;
  std::vector<std::size_t> local_to_global;
  for (std::size_t i = 0; i < global_x.size(); ++i) {
    const int owner_rank = (global_x[i] < 0.5) ? 0 : 1;
    if (owner_rank == world_rank) {
      local_x.push_back(global_x[i]);
      local_y.push_back(global_y[i]);
      local_z.push_back(global_z[i]);
      local_mass.push_back(global_mass[i]);
      local_to_global.push_back(i);
    }
  }

  std::vector<std::uint32_t> local_active(local_x.size(), 0U);
  for (std::size_t i = 0; i < local_active.size(); ++i) {
    local_active[i] = static_cast<std::uint32_t>(i);
  }
  ForceField local_field{
      std::vector<double>(local_x.size(), 0.0),
      std::vector<double>(local_x.size(), 0.0),
      std::vector<double>(local_x.size(), 0.0),
  };
  cosmosim::gravity::TreePmForceAccumulatorView local_accumulator{
      local_active,
      local_field.ax,
      local_field.ay,
      local_field.az,
  };

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = box_size_comoving;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.opening_theta = 0.55;
  options.tree_options.max_leaf_size = 4;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(2.0, 4.0, mesh_spacing);
  options.tree_exchange_batch_bytes = sizeof(double) * 3U * 3U;

  const auto local_layout =
      cosmosim::parallel::makePmSlabLayout(pm_shape.nx, pm_shape.ny, pm_shape.nz, world_size, world_rank);
  cosmosim::gravity::TreePmCoordinator distributed(pm_shape, local_layout);
  cosmosim::gravity::TreePmDiagnostics distributed_diag;
  distributed.solveActiveSet(local_x, local_y, local_z, local_mass, local_accumulator, options, nullptr, &distributed_diag);

  std::vector<std::uint32_t> ref_active(global_x.size(), 0U);
  for (std::size_t i = 0; i < ref_active.size(); ++i) {
    ref_active[i] = static_cast<std::uint32_t>(i);
  }
  ForceField reference_field{
      std::vector<double>(global_x.size(), 0.0),
      std::vector<double>(global_x.size(), 0.0),
      std::vector<double>(global_x.size(), 0.0),
  };
  cosmosim::gravity::TreePmForceAccumulatorView ref_accumulator{
      ref_active,
      reference_field.ax,
      reference_field.ay,
      reference_field.az,
  };
  cosmosim::gravity::TreePmCoordinator single_rank(pm_shape);
  cosmosim::gravity::TreePmDiagnostics reference_diag;
  single_rank.solveActiveSet(global_x, global_y, global_z, global_mass, ref_accumulator, options, nullptr, &reference_diag);

  double local_diff2 = 0.0;
  double local_ref2 = 0.0;
  for (std::size_t i = 0; i < local_x.size(); ++i) {
    const std::size_t gi = local_to_global[i];
    const double dx = local_field.ax[i] - reference_field.ax[gi];
    const double dy = local_field.ay[i] - reference_field.ay[gi];
    const double dz = local_field.az[i] - reference_field.az[gi];
    local_diff2 += dx * dx + dy * dy + dz * dz;
    local_ref2 += reference_field.ax[gi] * reference_field.ax[gi] +
        reference_field.ay[gi] * reference_field.ay[gi] +
        reference_field.az[gi] * reference_field.az[gi];
  }

  double global_diff2 = 0.0;
  double global_ref2 = 0.0;
  MPI_Allreduce(&local_diff2, &global_diff2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_ref2, &global_ref2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double rel_l2 = std::sqrt(global_diff2 / std::max(global_ref2, 1.0e-30));

  std::ostringstream msg;
  msg << "Distributed TreePM short-range export/import drift: rel_l2=" << rel_l2
      << ", world_rank=" << world_rank
      << ", residual_pruned_nodes=" << distributed_diag.residual_pruned_nodes
      << ", residual_pair_skips_cutoff=" << distributed_diag.residual_pair_skips_cutoff
      << ", reference_pair_skips_cutoff=" << reference_diag.residual_pair_skips_cutoff;
  requireOrThrow(rel_l2 <= 1.0e-9, msg.str());
  requireOrThrow(distributed_diag.residual_pair_skips_cutoff > 0, msg.str());
}

void testDistributedActiveSubsetMatchesSingleRankReference() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    return;
  }

  constexpr double box_size_comoving = 1.0;
  const cosmosim::gravity::PmGridShape pm_shape{16, 16, 16};
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);

  const std::vector<double> global_x = {
      0.06, 0.13, 0.22, 0.29, 0.38, 0.46, 0.495, 0.499,
      0.501, 0.507, 0.56, 0.62, 0.71, 0.80, 0.89, 0.97,
  };
  const std::vector<double> global_y = {
      0.09, 0.18, 0.26, 0.32, 0.40, 0.48, 0.50, 0.53,
      0.55, 0.59, 0.64, 0.68, 0.73, 0.78, 0.85, 0.92,
  };
  const std::vector<double> global_z = {
      0.04, 0.14, 0.23, 0.30, 0.36, 0.45, 0.494, 0.498,
      0.502, 0.509, 0.58, 0.65, 0.74, 0.82, 0.88, 0.96,
  };
  std::vector<double> global_mass(global_x.size(), 1.0);
  for (std::size_t i = 0; i < global_mass.size(); ++i) {
    global_mass[i] = 0.92 + 0.02 * static_cast<double>(i % 6U);
  }

  std::vector<double> local_x;
  std::vector<double> local_y;
  std::vector<double> local_z;
  std::vector<double> local_mass;
  std::vector<std::size_t> local_to_global;
  for (std::size_t i = 0; i < global_x.size(); ++i) {
    const int owner_rank = (global_x[i] < 0.5) ? 0 : 1;
    if (owner_rank == world_rank) {
      local_x.push_back(global_x[i]);
      local_y.push_back(global_y[i]);
      local_z.push_back(global_z[i]);
      local_mass.push_back(global_mass[i]);
      local_to_global.push_back(i);
    }
  }

  std::vector<std::uint32_t> local_active;
  std::vector<std::size_t> local_active_to_global;
  for (std::size_t li = 0; li < local_to_global.size(); ++li) {
    if ((local_to_global[li] % 2U) == 0U) {
      local_active.push_back(static_cast<std::uint32_t>(li));
      local_active_to_global.push_back(local_to_global[li]);
    }
  }

  ForceField local_field{
      std::vector<double>(local_active.size(), 0.0),
      std::vector<double>(local_active.size(), 0.0),
      std::vector<double>(local_active.size(), 0.0),
  };
  cosmosim::gravity::TreePmForceAccumulatorView local_accumulator{
      local_active,
      local_field.ax,
      local_field.ay,
      local_field.az,
  };

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = box_size_comoving;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.opening_theta = 0.55;
  options.tree_options.max_leaf_size = 4;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(2.0, 4.0, mesh_spacing);
  options.tree_exchange_batch_bytes = sizeof(double) * 3U * 3U;

  const auto local_layout =
      cosmosim::parallel::makePmSlabLayout(pm_shape.nx, pm_shape.ny, pm_shape.nz, world_size, world_rank);
  cosmosim::gravity::TreePmCoordinator distributed(pm_shape, local_layout);
  distributed.solveActiveSet(local_x, local_y, local_z, local_mass, local_accumulator, options, nullptr, nullptr);

  std::vector<std::uint32_t> ref_active(global_x.size(), 0U);
  for (std::size_t i = 0; i < ref_active.size(); ++i) {
    ref_active[i] = static_cast<std::uint32_t>(i);
  }
  ForceField reference_field{
      std::vector<double>(global_x.size(), 0.0),
      std::vector<double>(global_x.size(), 0.0),
      std::vector<double>(global_x.size(), 0.0),
  };
  cosmosim::gravity::TreePmForceAccumulatorView ref_accumulator{
      ref_active,
      reference_field.ax,
      reference_field.ay,
      reference_field.az,
  };
  cosmosim::gravity::TreePmCoordinator single_rank(pm_shape);
  single_rank.solveActiveSet(global_x, global_y, global_z, global_mass, ref_accumulator, options, nullptr, nullptr);

  double local_diff2 = 0.0;
  double local_ref2 = 0.0;
  for (std::size_t i = 0; i < local_active_to_global.size(); ++i) {
    const std::size_t gi = local_active_to_global[i];
    const double dx = local_field.ax[i] - reference_field.ax[gi];
    const double dy = local_field.ay[i] - reference_field.ay[gi];
    const double dz = local_field.az[i] - reference_field.az[gi];
    local_diff2 += dx * dx + dy * dy + dz * dz;
    local_ref2 += reference_field.ax[gi] * reference_field.ax[gi] +
        reference_field.ay[gi] * reference_field.ay[gi] +
        reference_field.az[gi] * reference_field.az[gi];
  }
  double global_diff2 = 0.0;
  double global_ref2 = 0.0;
  MPI_Allreduce(&local_diff2, &global_diff2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_ref2, &global_ref2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double rel_l2 = std::sqrt(global_diff2 / std::max(global_ref2, 1.0e-30));
  requireOrThrow(rel_l2 <= 1.0e-9, "distributed active-subset vs single-rank mismatch");
}
#endif


void testZoomMembershipPeriodicWrapAffectsCorrection() {
  constexpr double box_size_comoving = 1.0;
  std::vector<double> pos_x = {0.98, 0.02, 0.50};
  std::vector<double> pos_y = {0.50, 0.50, 0.50};
  std::vector<double> pos_z = {0.50, 0.50, 0.50};
  std::vector<double> mass = {1.0, 1.0, 2.0};
  std::vector<std::uint8_t> source_is_high_res = {1, 1, 0};
  std::vector<std::uint8_t> active_is_high_res = source_is_high_res;

  cosmosim::gravity::TreePmOptions opts;
  opts.pm_options.box_size_mpc_comoving = box_size_comoving;
  opts.pm_options.box_size_x_mpc_comoving = box_size_comoving;
  opts.pm_options.box_size_y_mpc_comoving = box_size_comoving;
  opts.pm_options.box_size_z_mpc_comoving = box_size_comoving;
  opts.pm_options.scale_factor = 1.0;
  opts.pm_options.gravitational_constant_code = 1.0;
  opts.tree_options.opening_theta = 0.55;
  opts.tree_options.max_leaf_size = 4;
  opts.tree_options.gravitational_constant_code = 1.0;
  opts.tree_options.softening.epsilon_comoving = 1.0e-3;
  const cosmosim::gravity::PmGridShape coarse_shape{16,16,16};
  const double mesh_spacing = box_size_comoving / static_cast<double>(coarse_shape.nx);
  opts.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 4.5, mesh_spacing);
  opts.enable_zoom_long_range_correction = true;
  opts.zoom_focused_pm_shape = cosmosim::gravity::PmGridShape{24,24,24};
  opts.source_is_high_res = source_is_high_res;
  opts.active_is_high_res = active_is_high_res;
  opts.zoom_region_center_x_comoving = 0.0;
  opts.zoom_region_center_y_comoving = 0.5;
  opts.zoom_region_center_z_comoving = 0.5;
  opts.zoom_region_radius_comoving = 0.06;
  opts.zoom_contamination_radius_comoving = 0.1;
  cosmosim::gravity::TreePmDiagnostics diag;
  const ForceField force = solveTreePm(pos_x, pos_y, pos_z, mass, coarse_shape, opts, &diag);
  requireOrThrow(diag.force_l2_pm_zoom_correction > 0.0, "wrapped zoom pair should produce focused correction");
  requireOrThrow(std::abs(force.ax[0]) > 0.0 || std::abs(force.ax[1]) > 0.0, "wrapped high-res targets should feel non-zero force");
}

void testZoomFocusedPmCorrectionAndContaminationDiagnostics() {
  constexpr double box_size_comoving = 1.0;
  std::vector<double> pos_x = {0.49, 0.52, 0.51, 0.60, 0.80};
  std::vector<double> pos_y = {0.50, 0.50, 0.53, 0.50, 0.80};
  std::vector<double> pos_z = {0.50, 0.50, 0.47, 0.50, 0.80};
  std::vector<double> mass = {1.0, 1.1, 0.9, 2.0, 2.5};
  std::vector<std::uint8_t> source_is_high_res = {1, 1, 1, 0, 0};
  std::vector<std::uint8_t> active_is_high_res = source_is_high_res;

  cosmosim::gravity::TreePmOptions no_zoom;
  no_zoom.pm_options.box_size_mpc_comoving = box_size_comoving;
  no_zoom.pm_options.scale_factor = 1.0;
  no_zoom.pm_options.gravitational_constant_code = 1.0;
  no_zoom.tree_options.opening_theta = 0.55;
  no_zoom.tree_options.max_leaf_size = 4;
  no_zoom.tree_options.gravitational_constant_code = 1.0;
  no_zoom.tree_options.softening.epsilon_comoving = 1.0e-3;
  const cosmosim::gravity::PmGridShape coarse_shape{16, 16, 16};
  const double mesh_spacing = box_size_comoving / static_cast<double>(coarse_shape.nx);
  no_zoom.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 4.5, mesh_spacing);
  cosmosim::gravity::TreePmDiagnostics no_zoom_diag;
  const ForceField no_zoom_force =
      solveTreePm(pos_x, pos_y, pos_z, mass, coarse_shape, no_zoom, &no_zoom_diag);

  cosmosim::gravity::TreePmOptions with_zoom = no_zoom;
  with_zoom.enable_zoom_long_range_correction = true;
  with_zoom.zoom_focused_pm_shape = cosmosim::gravity::PmGridShape{32, 32, 32};
  with_zoom.source_is_high_res = source_is_high_res;
  with_zoom.active_is_high_res = active_is_high_res;
  with_zoom.zoom_region_center_x_comoving = 0.5;
  with_zoom.zoom_region_center_y_comoving = 0.5;
  with_zoom.zoom_region_center_z_comoving = 0.5;
  with_zoom.zoom_region_radius_comoving = 0.08;
  with_zoom.zoom_contamination_radius_comoving = 0.12;
  cosmosim::gravity::TreePmDiagnostics with_zoom_diag;
  const ForceField with_zoom_force =
      solveTreePm(pos_x, pos_y, pos_z, mass, coarse_shape, with_zoom, &with_zoom_diag);

  const double highres_delta = std::abs(with_zoom_force.ax[0] - no_zoom_force.ax[0]) +
      std::abs(with_zoom_force.ay[1] - no_zoom_force.ay[1]) +
      std::abs(with_zoom_force.az[2] - no_zoom_force.az[2]);
  const double lowres_delta = std::abs(with_zoom_force.ax[3] - no_zoom_force.ax[3]) +
      std::abs(with_zoom_force.ay[4] - no_zoom_force.ay[4]);
  requireOrThrow(highres_delta > 1.0e-6, "zoom PM correction should change high-res target forces");
  requireOrThrow(lowres_delta < 1.0e-12, "zoom PM correction must not modify low-res target forces");
  requireOrThrow(with_zoom_diag.zoom_high_res_source_count == 3, "unexpected zoom high-res source count");
  requireOrThrow(with_zoom_diag.zoom_low_res_source_count == 2, "unexpected zoom low-res source count");
  requireOrThrow(with_zoom_diag.zoom_low_res_contamination_count == 1, "expected exactly one low-res contaminant");
  requireOrThrow(with_zoom_diag.force_l2_pm_zoom_correction > 0.0, "zoom PM correction norm should be positive");
  requireOrThrow(std::isfinite(with_zoom_diag.force_l2_pm_global), "pm-global force norm must be finite");
  requireOrThrow(std::isfinite(with_zoom_diag.force_l2_tree_short_range), "tree-short force norm must be finite");
  requireOrThrow(std::isfinite(with_zoom_diag.force_l2_total), "total force norm must be finite");
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
#endif
  testPeriodicTreePmAgainstDirectReference();
  testResidualCutoffPrunesTreePath();
  testPmOnlyTreeOnlyAndTreePmConsistency();
  testZoomMembershipPeriodicWrapAffectsCorrection();
  testZoomFocusedPmCorrectionAndContaminationDiagnostics();
  testActiveSubsetMatchesFullSolveSingleRank();
#if COSMOSIM_ENABLE_MPI
  testDistributedShortRangeExportImportMatchesSingleRankReference();
  testDistributedActiveSubsetMatchesSingleRankReference();
  MPI_Finalize();
#endif
  return 0;
}
