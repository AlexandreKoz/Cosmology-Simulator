#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
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

[[nodiscard]] constexpr cosmosim::gravity::PmGridShape pmShapeForAvailableBackend(
    std::size_t fftw_n,
    std::size_t fallback_n) noexcept {
#if COSMOSIM_ENABLE_FFTW
  return {fftw_n, fftw_n, fftw_n};
#else
  return {fallback_n, fallback_n, fallback_n};
#endif
}

[[nodiscard]] constexpr double rcutCellsForAvailableBackend(
    double fftw_rcut_cells,
    double fallback_rcut_cells) noexcept {
  // Naive-DFT fixtures intentionally use tiny meshes. Keep their diagnostic
  // split inside the production minimum-image cutoff domain without raising
  // the O(N_mesh^2) fallback cost; FFTW fixtures retain their named profile.
#if COSMOSIM_ENABLE_FFTW
  return fftw_rcut_cells;
#else
  return fallback_rcut_cells;
#endif
}

[[nodiscard]] double minimumImageDelta(double delta, double box_size_comoving) {
  return delta - box_size_comoving * std::nearbyint(delta / box_size_comoving);
}

struct ForceField {
  std::vector<double> ax;
  std::vector<double> ay;
  std::vector<double> az;
};

struct AxisBoxLengths {
  double x = 1.0;
  double y = 1.0;
  double z = 1.0;
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
    cosmosim::gravity::TreePmProfileEvent* profile = nullptr,
    std::span<const double> previous_acceleration_magnitude_code = {}) {
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
      previous_acceleration_magnitude_code,
  };

  cosmosim::gravity::TreePmCoordinator coordinator(pm_shape);
  coordinator.solveActiveSet(pos_x, pos_y, pos_z, mass, accumulator, options, profile, diagnostics);
  return field;
}

[[nodiscard]] ForceField solveTreePmActiveSubset(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    std::span<const std::uint32_t> active,
    const cosmosim::gravity::PmGridShape& pm_shape,
    cosmosim::gravity::TreePmOptions options,
    cosmosim::gravity::TreePmDiagnostics* diagnostics,
    cosmosim::gravity::TreePmProfileEvent* profile = nullptr,
    std::span<const double> previous_acceleration_magnitude_code = {}) {
  ForceField field{
      std::vector<double>(active.size(), 0.0),
      std::vector<double>(active.size(), 0.0),
      std::vector<double>(active.size(), 0.0),
  };
  cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      active, field.ax, field.ay, field.az, previous_acceleration_magnitude_code};
  cosmosim::gravity::TreePmCoordinator coordinator(pm_shape);
  coordinator.solveActiveSet(pos_x, pos_y, pos_z, mass, accumulator, options, profile, diagnostics);
  return field;
}

void requireEquivalentForces(
    const ForceField& lhs,
    const ForceField& rhs,
    double tolerance,
    const std::string& context) {
  requireOrThrow(lhs.ax.size() == rhs.ax.size(), context + ": force extent mismatch");
  for (std::size_t i = 0; i < lhs.ax.size(); ++i) {
    const std::array lhs_components{lhs.ax[i], lhs.ay[i], lhs.az[i]};
    const std::array rhs_components{rhs.ax[i], rhs.ay[i], rhs.az[i]};
    for (std::size_t axis = 0; axis < lhs_components.size(); ++axis) {
      const double scale = std::max({1.0, std::abs(lhs_components[axis]), std::abs(rhs_components[axis])});
      if (std::abs(lhs_components[axis] - rhs_components[axis]) > tolerance * scale) {
        std::ostringstream message;
        message << context << ": force mismatch at target=" << i << ", axis=" << axis
                << ", lhs=" << lhs_components[axis] << ", rhs=" << rhs_components[axis]
                << ", tolerance=" << tolerance * scale;
        throw std::runtime_error(message.str());
      }
    }
  }
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

  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(32, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, rcutCellsForAvailableBackend(4.5, 3.9), mesh_spacing);

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

void testPeriodicTreeGeometryIsSeamSafeAndTranslationInvariant() {
  constexpr std::uint8_t x_seam = 1U << 0U;
  constexpr std::uint8_t y_seam = 1U << 1U;
  constexpr std::uint8_t z_seam = 1U << 2U;
  const AxisBoxLengths box{.x = 1.0, .y = 1.5, .z = 2.0};
  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};

  enum class TargetPlacement {
    kAllClusterParticles,
    kInteriorTarget,
    kSeamTargetWithInteriorSources,
    kSeamTargetWithSeamSources,
  };
  struct Case {
    std::string label;
    std::uint8_t source_seams = 0U;
    TargetPlacement target_placement = TargetPlacement::kAllClusterParticles;
    cosmosim::gravity::TreeMultipoleOrder multipole_order = cosmosim::gravity::TreeMultipoleOrder::kMonopole;
    cosmosim::gravity::TreeOpeningCriterion opening_criterion =
        cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric;
    bool direct_leaf = false;
    bool require_compact_root = false;
    bool require_internal_acceptance = false;
  };

  std::vector<Case> cases;
  for (const auto order : {cosmosim::gravity::TreeMultipoleOrder::kMonopole,
                           cosmosim::gravity::TreeMultipoleOrder::kQuadrupole}) {
    for (const auto criterion : {cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric,
                                 cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance}) {
      for (const bool direct_leaf : {true, false}) {
        cases.push_back(Case{
            .label = "x_seam_full_numerics_matrix",
            .source_seams = x_seam,
            .target_placement = TargetPlacement::kAllClusterParticles,
            .multipole_order = order,
            .opening_criterion = criterion,
            .direct_leaf = direct_leaf,
            .require_compact_root = true,
        });
      }
    }
  }
  cases.push_back(Case{.label = "y_seam", .source_seams = y_seam, .require_compact_root = true});
  cases.push_back(Case{
      .label = "z_seam_quadrupole",
      .source_seams = z_seam,
      .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole,
      .require_compact_root = true});
  cases.push_back(Case{
      .label = "xy_edge_com_mac",
      .source_seams = static_cast<std::uint8_t>(x_seam | y_seam),
      .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance,
      .require_compact_root = true});
  cases.push_back(Case{
      .label = "xyz_corner_quadrupole",
      .source_seams = static_cast<std::uint8_t>(x_seam | y_seam | z_seam),
      .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole,
      .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance,
      .require_compact_root = true});
  cases.push_back(Case{
      .label = "interior_target_seam_sources",
      .source_seams = static_cast<std::uint8_t>(x_seam | y_seam),
      .target_placement = TargetPlacement::kInteriorTarget,
      .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole,
      .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance,
      .require_internal_acceptance = true});
  cases.push_back(Case{
      .label = "seam_target_interior_sources",
      .source_seams = 0U,
      .target_placement = TargetPlacement::kSeamTargetWithInteriorSources,
      .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole});
  cases.push_back(Case{
      .label = "seam_target_seam_sources",
      .source_seams = static_cast<std::uint8_t>(x_seam | y_seam | z_seam),
      .target_placement = TargetPlacement::kSeamTargetWithSeamSources,
      .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance});

  for (const Case& test_case : cases) {
    constexpr std::size_t source_count = 12U;
    std::vector<double> pos_x(source_count, 0.0);
    std::vector<double> pos_y(source_count, 0.0);
    std::vector<double> pos_z(source_count, 0.0);
    std::vector<double> mass(source_count, 1.0);
    const auto coordinate = [](std::size_t i, double length, bool crosses_seam) {
      const double offset = (0.006 + 0.0015 * static_cast<double>(i % 4U)) * length;
      if (crosses_seam) {
        return (i % 2U) == 0U ? offset : length - offset;
      }
      return (0.5 + 0.004 * static_cast<double>(static_cast<int>(i % 5U) - 2)) * length;
    };
    for (std::size_t i = 0; i < source_count; ++i) {
      pos_x[i] = coordinate(i, box.x, (test_case.source_seams & x_seam) != 0U);
      pos_y[i] = coordinate(i + 1U, box.y, (test_case.source_seams & y_seam) != 0U);
      pos_z[i] = coordinate(i + 2U, box.z, (test_case.source_seams & z_seam) != 0U);
      mass[i] = 0.7 + 0.07 * static_cast<double>(i % 6U);
    }

    std::vector<std::uint32_t> active;
    if (test_case.target_placement == TargetPlacement::kAllClusterParticles) {
      active.resize(source_count);
      for (std::size_t i = 0; i < source_count; ++i) {
        active[i] = static_cast<std::uint32_t>(i);
      }
    } else {
      const bool target_on_seam =
          test_case.target_placement == TargetPlacement::kSeamTargetWithInteriorSources ||
          test_case.target_placement == TargetPlacement::kSeamTargetWithSeamSources;
      pos_x.push_back(target_on_seam ? 0.998 * box.x : 0.47 * box.x);
      pos_y.push_back(target_on_seam ? 0.002 * box.y : 0.44 * box.y);
      pos_z.push_back(target_on_seam ? 0.999 * box.z : 0.53 * box.z);
      mass.push_back(0.0);
      active.push_back(static_cast<std::uint32_t>(source_count));
    }

    std::vector<double> lifted_x = pos_x;
    std::vector<double> lifted_y = pos_y;
    std::vector<double> lifted_z = pos_z;
    for (std::size_t i = 0; i < pos_x.size(); ++i) {
      lifted_x[i] += static_cast<double>(static_cast<int>(i % 3U) - 1) * box.x;
      lifted_y[i] += static_cast<double>(static_cast<int>((i + 1U) % 3U) - 1) * box.y;
      lifted_z[i] += static_cast<double>(static_cast<int>((i + 2U) % 3U) - 1) * box.z;
    }

    cosmosim::gravity::TreePmOptions options;
    options.pm_options.box_size_x_mpc_comoving = box.x;
    options.pm_options.box_size_y_mpc_comoving = box.y;
    options.pm_options.box_size_z_mpc_comoving = box.z;
    options.pm_options.scale_factor = 1.0;
    options.pm_options.gravitational_constant_code = 1.0;
    options.tree_options.opening_theta = 0.8;
    options.tree_options.opening_criterion = test_case.opening_criterion;
    options.tree_options.multipole_order = test_case.multipole_order;
    options.tree_options.max_leaf_size = test_case.direct_leaf ? 64U : 1U;
    options.tree_options.gravitational_constant_code = 1.0;
    options.tree_options.softening.epsilon_comoving = 2.0e-3;
    options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
        1.25, 3.9, box.x / static_cast<double>(pm_shape.nx));

    cosmosim::gravity::TreePmDiagnostics wrapped_diagnostics;
    cosmosim::gravity::TreePmDiagnostics lifted_diagnostics;
    cosmosim::gravity::TreePmProfileEvent wrapped_profile;
    const ForceField wrapped_force = solveTreePmActiveSubset(
        pos_x, pos_y, pos_z, mass, active, pm_shape, options, &wrapped_diagnostics, &wrapped_profile);
    const ForceField lifted_force = solveTreePmActiveSubset(
        lifted_x, lifted_y, lifted_z, mass, active, pm_shape, options, &lifted_diagnostics);
    requireEquivalentForces(wrapped_force, lifted_force, 5.0e-11, test_case.label);

    const double geometry_scale = std::max({1.0,
                                            wrapped_diagnostics.tree_root_half_size_comoving,
                                            lifted_diagnostics.tree_root_half_size_comoving});
    requireOrThrow(
        std::abs(wrapped_diagnostics.tree_root_half_size_comoving -
                 lifted_diagnostics.tree_root_half_size_comoving) <= 5.0e-12 * geometry_scale,
        test_case.label + ": periodic translation changed root extent");
    if (test_case.require_compact_root) {
      requireOrThrow(
          wrapped_diagnostics.tree_root_half_size_comoving < 0.03 * box.z,
          test_case.label + ": seam cluster produced a noncompact tree root");
    }
    if (test_case.require_internal_acceptance) {
      requireOrThrow(
          wrapped_profile.tree_profile.particle_particle_interactions < source_count,
          test_case.label + ": distant compact seam cluster was not accepted as an internal node");
    }
  }
}

void testRelativeForceErrorMacUsesOwnerHistoryAcrossPeriodicSeam() {
  const AxisBoxLengths box{.x = 1.0, .y = 1.4, .z = 1.8};
  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};
  std::vector<double> pos_x{0.006, 0.994, 0.009, 0.991, 0.012, 0.988, 0.015, 0.985, 0.35};
  std::vector<double> pos_y{0.69, 0.71, 0.68, 0.72, 0.695, 0.705, 0.685, 0.715, 0.70};
  std::vector<double> pos_z{0.89, 0.91, 0.88, 0.92, 0.895, 0.905, 0.885, 0.915, 0.90};
  std::vector<double> mass(pos_x.size(), 1.0);
  mass.back() = 0.0;
  const std::vector<std::uint32_t> active{static_cast<std::uint32_t>(pos_x.size() - 1U)};

  cosmosim::gravity::TreePmOptions relative_options;
  relative_options.pm_options.box_size_x_mpc_comoving = box.x;
  relative_options.pm_options.box_size_y_mpc_comoving = box.y;
  relative_options.pm_options.box_size_z_mpc_comoving = box.z;
  relative_options.pm_options.gravitational_constant_code = 1.0;
  relative_options.tree_options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError;
  relative_options.tree_options.opening_theta = 0.03;
  relative_options.tree_options.multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole;
  relative_options.tree_options.max_leaf_size = 1U;
  relative_options.tree_options.relative_force_tolerance = 0.005;
  relative_options.tree_options.relative_force_acceleration_floor_code = 1.0e-12;
  relative_options.tree_options.softening.epsilon_comoving = 1.0e-3;
  relative_options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, 3.9, box.x / static_cast<double>(pm_shape.nx));

  cosmosim::gravity::TreePmOptions com_options = relative_options;
  com_options.tree_options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance;
  cosmosim::gravity::TreePmProfileEvent missing_profile;
  cosmosim::gravity::TreePmProfileEvent com_profile;
  const ForceField missing_history = solveTreePmActiveSubset(
      pos_x, pos_y, pos_z, mass, active, pm_shape, relative_options, nullptr, &missing_profile);
  const ForceField com_fallback = solveTreePmActiveSubset(
      pos_x, pos_y, pos_z, mass, active, pm_shape, com_options, nullptr, &com_profile);
  requireEquivalentForces(missing_history, com_fallback, 1.0e-13, "relative MAC missing-history fallback");
  requireOrThrow(
      missing_profile.tree_profile.visited_nodes == com_profile.tree_profile.visited_nodes,
      "relative MAC missing history must reproduce COM-distance tree work");

  const std::vector<double> previous_acceleration{1.0e8};
  cosmosim::gravity::TreePmProfileEvent history_profile;
  const ForceField with_history = solveTreePmActiveSubset(
      pos_x,
      pos_y,
      pos_z,
      mass,
      active,
      pm_shape,
      relative_options,
      nullptr,
      &history_profile,
      previous_acceleration);
  requireOrThrow(
      history_profile.tree_profile.particle_particle_interactions <
          missing_profile.tree_profile.particle_particle_interactions,
      "relative MAC owner history should permit less tree work for a distant compact node");

  std::vector<double> lifted_x = pos_x;
  std::vector<double> lifted_y = pos_y;
  std::vector<double> lifted_z = pos_z;
  for (std::size_t i = 0; i < pos_x.size(); ++i) {
    lifted_x[i] += (i % 2U == 0U ? box.x : -box.x);
    lifted_y[i] += (i % 3U == 0U ? box.y : -box.y);
    lifted_z[i] += (i % 2U == 0U ? -box.z : box.z);
  }
  const ForceField lifted_history = solveTreePmActiveSubset(
      lifted_x,
      lifted_y,
      lifted_z,
      mass,
      active,
      pm_shape,
      relative_options,
      nullptr,
      nullptr,
      previous_acceleration);
  requireEquivalentForces(with_history, lifted_history, 5.0e-11, "relative MAC periodic translation");
}

void testIndependentTargetAndEmptySourceContractSingleRank() {
  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 3.9, 0.125);

  const std::vector<double> empty;
  const std::vector<std::uint32_t> target_identity{std::numeric_limits<std::uint32_t>::max()};
  const std::vector<double> target_x{0.01};
  const std::vector<double> target_y{0.25};
  const std::vector<double> target_z{0.75};
  ForceField empty_source_force{{0.0}, {0.0}, {0.0}};
  const cosmosim::gravity::TreePmForceAccumulatorView empty_source_target{
      .active_particle_index = target_identity,
      .accel_x_comoving = empty_source_force.ax,
      .accel_y_comoving = empty_source_force.ay,
      .accel_z_comoving = empty_source_force.az,
      .target_pos_x_comoving = target_x,
      .target_pos_y_comoving = target_y,
      .target_pos_z_comoving = target_z,
  };
  cosmosim::gravity::TreePmCoordinator empty_source_coordinator(pm_shape);
  empty_source_coordinator.solveActiveSet(
      empty, empty, empty, empty, empty_source_target, options, nullptr, nullptr);
  requireOrThrow(
      empty_source_force.ax[0] == 0.0 && empty_source_force.ay[0] == 0.0 && empty_source_force.az[0] == 0.0,
      "an independent target with no sources must receive exactly zero force");

  const std::vector<std::uint32_t> no_targets;
  ForceField all_empty_force;
  const cosmosim::gravity::TreePmForceAccumulatorView all_empty{
      .active_particle_index = no_targets,
      .accel_x_comoving = all_empty_force.ax,
      .accel_y_comoving = all_empty_force.ay,
      .accel_z_comoving = all_empty_force.az,
  };
  cosmosim::gravity::TreePmCoordinator all_empty_coordinator(pm_shape);
  all_empty_coordinator.solveActiveSet(empty, empty, empty, empty, all_empty, options, nullptr, nullptr);
}

void testResidualCutoffPrunesTreePath() {
  constexpr double box_size_comoving = 1.0;
  std::vector<double> pos_x = {0.10, 0.18, 0.62, 0.88};
  std::vector<double> pos_y = {0.10, 0.12, 0.60, 0.90};
  std::vector<double> pos_z = {0.10, 0.09, 0.58, 0.86};
  std::vector<double> mass = {1.0, 1.0, 1.0, 1.0};

  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(32, 8);
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

void testPmOnlyAndTreePmConsistency() {
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

  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(32, 8);
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

  cosmosim::gravity::TreePmOptions pm_only_options = base_options;
  pm_only_options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      0.2, rcutCellsForAvailableBackend(4.5, 3.9), mesh_spacing);
  cosmosim::gravity::TreePmDiagnostics pm_only_diagnostics;
  const ForceField pm_only = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, pm_only_options, &pm_only_diagnostics);

  cosmosim::gravity::TreePmOptions split_options = base_options;
  split_options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, rcutCellsForAvailableBackend(4.5, 3.9), mesh_spacing);
  cosmosim::gravity::TreePmDiagnostics split_diagnostics;
  const ForceField split = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, split_options, &split_diagnostics);

  const double pm_only_rel = relativeL2Error(pm_only, direct);
  const double split_rel = relativeL2Error(split, direct);

  std::ostringstream message;
  message << "TreePM consistency failure: split_rel=" << split_rel
          << ", pm_only_rel=" << pm_only_rel << ", split_pruned_nodes=" << split_diagnostics.residual_pruned_nodes
          << ", split_pair_skips=" << split_diagnostics.residual_pair_skips_cutoff;

  requireOrThrow(split_rel < 0.90, message.str());
  requireOrThrow(pm_only_rel < 2.2, message.str());
  requireOrThrow(split_rel <= pm_only_rel + 1.0e-9, message.str());
}

void testCosmologicalTreePmNormalizationAndSplitComplementarity() {
  constexpr double box_size_comoving = 1.0;
  constexpr double non_unit_scale_factor = 0.5;
  const std::vector<double> pos_x = {
      0.030, 0.052, 0.105, 0.188, 0.327, 0.469,
      0.531, 0.674, 0.812, 0.895, 0.948, 0.970,
  };
  const std::vector<double> pos_y = {
      0.120, 0.137, 0.490, 0.715, 0.238, 0.881,
      0.862, 0.401, 0.643, 0.292, 0.108, 0.091,
  };
  const std::vector<double> pos_z = {
      0.760, 0.744, 0.205, 0.349, 0.917, 0.558,
      0.579, 0.083, 0.426, 0.681, 0.816, 0.834,
  };
  const std::vector<double> mass = {
      0.80, 1.15, 0.95, 1.30, 0.72, 1.08,
      0.91, 1.22, 0.86, 1.11, 0.77, 1.04,
  };

  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(32, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);
  cosmosim::gravity::TreePmOptions unit_options;
  unit_options.pm_options.box_size_mpc_comoving = box_size_comoving;
  unit_options.pm_options.scale_factor = 1.0;
  unit_options.pm_options.gravitational_constant_code = 1.0;
  unit_options.pm_options.enable_window_deconvolution = true;
  unit_options.tree_options.opening_theta = 0.55;
  unit_options.tree_options.max_leaf_size = 2U;
  unit_options.tree_options.gravitational_constant_code = 1.0;
  unit_options.tree_options.softening.epsilon_comoving = 2.0e-3;
  unit_options.split_policy =
      cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
          1.25, rcutCellsForAvailableBackend(6.25, 3.9), mesh_spacing);

  cosmosim::gravity::TreePmDiagnostics unit_diagnostics;
  const ForceField unit_force =
      solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, unit_options, &unit_diagnostics);

  cosmosim::gravity::TreePmOptions scaled_options = unit_options;
  scaled_options.pm_options.scale_factor = non_unit_scale_factor;
  cosmosim::gravity::TreePmDiagnostics scaled_diagnostics;
  const ForceField scaled_force =
      solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, scaled_options, &scaled_diagnostics);

  const double expected_force_scale = 1.0;
  ForceField expected_scaled_force = unit_force;
  for (std::size_t i = 0; i < expected_scaled_force.ax.size(); ++i) {
    expected_scaled_force.ax[i] *= expected_force_scale;
    expected_scaled_force.ay[i] *= expected_force_scale;
    expected_scaled_force.az[i] *= expected_force_scale;
  }
  const double total_scaling_rel_l2 = relativeL2Error(scaled_force, expected_scaled_force);
  const auto component_scaling_error = [expected_force_scale](double scaled, double unit) {
    const double expected = expected_force_scale * unit;
    return std::abs(scaled - expected) / std::max(std::abs(expected), 1.0e-30);
  };
  const auto norm_scaling_error = [expected_force_scale](double scaled, double unit) {
    return std::abs(scaled - expected_force_scale * unit) /
        std::max(expected_force_scale * unit, 1.0e-30);
  };
  const double pm_norm_scaling_error = norm_scaling_error(
      scaled_diagnostics.force_l2_pm_global, unit_diagnostics.force_l2_pm_global);
  const double tree_norm_scaling_error = norm_scaling_error(
      scaled_diagnostics.force_l2_tree_short_range,
      unit_diagnostics.force_l2_tree_short_range);
  const double total_norm_scaling_error = norm_scaling_error(
      scaled_diagnostics.force_l2_total, unit_diagnostics.force_l2_total);
  double max_component_scaling_error = 0.0;
  for (std::size_t i = 0; i < unit_force.ax.size(); ++i) {
    for (const auto [scaled, unit] : std::array{
             std::pair{scaled_force.ax[i], unit_force.ax[i]},
             std::pair{scaled_force.ay[i], unit_force.ay[i]},
             std::pair{scaled_force.az[i], unit_force.az[i]}}) {
      if (std::abs(expected_force_scale * unit) > 1.0e-12 * unit_diagnostics.force_l2_total) {
        max_component_scaling_error =
            std::max(max_component_scaling_error, component_scaling_error(scaled, unit));
      }
    }
  }

  const double pm_fraction = unit_diagnostics.force_l2_pm_global /
      (unit_diagnostics.force_l2_pm_global + unit_diagnostics.force_l2_tree_short_range);
  const double tree_fraction = unit_diagnostics.force_l2_tree_short_range /
      (unit_diagnostics.force_l2_pm_global + unit_diagnostics.force_l2_tree_short_range);
  std::ostringstream message;
  message << "cosmological TreePM normalization/split regression failed: a="
          << non_unit_scale_factor
          << ", expected_force_scale=" << expected_force_scale
          << ", total_scaling_rel_l2=" << total_scaling_rel_l2
          << ", max_component_scaling_error=" << max_component_scaling_error
          << ", pm_norm_scaling_error=" << pm_norm_scaling_error
          << ", tree_norm_scaling_error=" << tree_norm_scaling_error
          << ", total_norm_scaling_error=" << total_norm_scaling_error
          << ", unit_pm_l2=" << unit_diagnostics.force_l2_pm_global
          << ", unit_tree_l2=" << unit_diagnostics.force_l2_tree_short_range
          << ", pm_fraction=" << pm_fraction
          << ", tree_fraction=" << tree_fraction
          << ", unit_composition_error=" << unit_diagnostics.max_relative_composition_error
          << ", scaled_composition_error=" << scaled_diagnostics.max_relative_composition_error;

  requireOrThrow(unit_diagnostics.force_l2_total > 0.0, message.str());
  requireOrThrow(pm_fraction > 1.0e-2 && tree_fraction > 1.0e-2, message.str());
  requireOrThrow(total_scaling_rel_l2 <= 2.0e-12, message.str());
  requireOrThrow(max_component_scaling_error <= 2.0e-11, message.str());
  requireOrThrow(pm_norm_scaling_error <= 2.0e-12, message.str());
  requireOrThrow(tree_norm_scaling_error <= 2.0e-12, message.str());
  requireOrThrow(total_norm_scaling_error <= 2.0e-12, message.str());
  requireOrThrow(unit_diagnostics.max_relative_composition_error <= 1.0e-12, message.str());
  requireOrThrow(scaled_diagnostics.max_relative_composition_error <= 1.0e-12, message.str());

  bool emit_metrics = true;
#if COSMOSIM_ENABLE_MPI
  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  emit_metrics = world_rank == 0;
#endif
  if (emit_metrics) {
    std::cout << "TREEPM_COSMOLOGICAL_NORMALIZATION_PASS"
              << " a=" << non_unit_scale_factor
              << " expected_force_scale=" << expected_force_scale
              << " total_scaling_rel_l2=" << total_scaling_rel_l2
              << " max_component_scaling_error=" << max_component_scaling_error
              << " pm_norm_scaling_error=" << pm_norm_scaling_error
              << " tree_norm_scaling_error=" << tree_norm_scaling_error
              << " total_norm_scaling_error=" << total_norm_scaling_error
              << " pm_fraction=" << pm_fraction
              << " tree_fraction=" << tree_fraction
              << '\n';
  }
}

void testTreePmRejectsMismatchedGravitationalConstants() {
  const std::vector<double> pos_x = {0.20, 0.35};
  const std::vector<double> pos_y = {0.30, 0.42};
  const std::vector<double> pos_z = {0.40, 0.53};
  const std::vector<double> mass = {1.0, 1.0};
  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.gravitational_constant_code = 1.25;
  options.split_policy =
      cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 3.9, 0.125);

  bool rejected = false;
  try {
    [[maybe_unused]] const ForceField field =
        solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, options, nullptr);
  } catch (const std::invalid_argument& error) {
    rejected = std::string(error.what()).find("gravitational_constant_code values must match") !=
        std::string::npos;
  }
  requireOrThrow(
      rejected,
      "TreePM must reject mismatched PM/tree gravitational constants before composing the force split");
}

void testPeriodicTreePmRejectsAmbiguousMinimumImageCutoff() {
  const std::vector<double> pos_x = {0.20, 0.35};
  const std::vector<double> pos_y = {0.30, 0.42};
  const std::vector<double> pos_z = {0.40, 0.53};
  const std::vector<double> mass = {1.0, 1.0};
  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.gravitational_constant_code = 1.0;
  options.split_policy =
      cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 4.0, 0.125);

  bool rejected = false;
  try {
    [[maybe_unused]] const ForceField field =
        solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, options, nullptr);
  } catch (const std::invalid_argument& error) {
    rejected = std::string(error.what()).find("< half the shortest box axis") !=
        std::string::npos;
  }
  requireOrThrow(
      rejected,
      "periodic TreePM must reject r_cut=L_min/2 before choosing one of two equidistant images");
}

void testActiveSubsetMatchesFullSolveSingleRank() {
  constexpr double box_size_comoving = 1.0;
  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(16, 8);
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
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.5, rcutCellsForAvailableBackend(4.0, 3.9), mesh_spacing);

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

void testLongRangeCadenceCacheIsFailClosedAndOwnershipCompatible() {
#if COSMOSIM_ENABLE_MPI
  int world_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size != 1) {
    return;
  }
#endif
  const cosmosim::gravity::PmGridShape shape =
      pmShapeForAvailableBackend(8U, 8U);
  const std::vector<double> pos_x{0.15, 0.62};
  const std::vector<double> pos_y{0.28, 0.71};
  const std::vector<double> pos_z{0.41, 0.84};
  const std::vector<double> mass{1.0, 1.5};
  const std::vector<std::uint32_t> active{0U, 1U};
  ForceField force{
      std::vector<double>(active.size(), 0.0),
      std::vector<double>(active.size(), 0.0),
      std::vector<double>(active.size(), 0.0)};
  const cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      .active_particle_index = active,
      .accel_x_comoving = force.ax,
      .accel_y_comoving = force.ay,
      .accel_z_comoving = force.az,
  };
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.box_size_x_mpc_comoving = 1.0;
  options.pm_options.box_size_y_mpc_comoving = 1.0;
  options.pm_options.box_size_z_mpc_comoving = 1.0;
  options.pm_options.scale_factor = 0.5;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, 3.9, 1.0 / static_cast<double>(shape.nx));
  options.decomposition_epoch = 3U;
  options.force_epoch = 7U;

  cosmosim::gravity::TreePmCoordinator coordinator(shape);
  cosmosim::gravity::TreePmDiagnostics diagnostics;
  coordinator.solveActiveSetWithPmCadence(
      pos_x, pos_y, pos_z, mass, accumulator, options, true, nullptr,
      &diagnostics);
  requireOrThrow(diagnostics.pm_solve_count == 1U,
                 "authorized cadence refresh did not solve PM");

  force.ax.assign(active.size(), 0.0);
  force.ay.assign(active.size(), 0.0);
  force.az.assign(active.size(), 0.0);
  options.decomposition_epoch = 4U;
  coordinator.solveActiveSetWithPmCadence(
      pos_x, pos_y, pos_z, mass, accumulator, options, false, nullptr,
      &diagnostics);
  requireOrThrow(diagnostics.pm_reuse_count == 1U,
                 "particle-ownership epoch incorrectly invalidated fixed-slab PM cache");

  options.force_epoch = 8U;
  bool incompatible_reuse_threw = false;
  try {
    coordinator.solveActiveSetWithPmCadence(
        pos_x, pos_y, pos_z, mass, accumulator, options, false, nullptr,
        &diagnostics);
  } catch (const std::runtime_error& ex) {
    incompatible_reuse_threw =
        std::string(ex.what()).find("without a compatible PM field") !=
        std::string::npos;
  }
  requireOrThrow(incompatible_reuse_threw,
                 "incompatible PM cache reuse did not fail closed");
}

#if COSMOSIM_ENABLE_MPI
void testDivergentTreePmLayoutFailsCoordinately() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size < 2 || world_size > 4) {
    return;
  }

  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, 3.9, 1.0 / static_cast<double>(pm_shape.nx));

  const std::vector<double> empty;
  const std::vector<std::uint32_t> no_targets;
  ForceField no_force;
  const cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      .active_particle_index = no_targets,
      .accel_x_comoving = no_force.ax,
      .accel_y_comoving = no_force.ay,
      .accel_z_comoving = no_force.az,
  };

  bool rejected = false;
  std::string diagnostic;
  try {
    if (world_rank == 0) {
      // Deliberately advertise a rank-local serial layout while every peer
      // advertises the communicator-wide slab layout. All ranks must vote the
      // mismatch before any rank enters a layout-selected PM collective.
      cosmosim::gravity::TreePmCoordinator coordinator(pm_shape);
      coordinator.solveActiveSet(
          empty, empty, empty, empty, accumulator, options, nullptr, nullptr);
    } else {
      const auto layout = cosmosim::parallel::makePmSlabLayout(
          pm_shape.nx, pm_shape.ny, pm_shape.nz, world_size, world_rank);
      cosmosim::gravity::TreePmCoordinator coordinator(pm_shape, layout);
      coordinator.solveActiveSet(
          empty, empty, empty, empty, accumulator, options, nullptr, nullptr);
    }
  } catch (const std::exception& ex) {
    rejected = true;
    diagnostic = ex.what();
  }

  const std::uint64_t local_rejection = rejected ? 1U : 0U;
  std::uint64_t rejection_count = 0U;
  MPI_Allreduce(
      &local_rejection,
      &rejection_count,
      1,
      MPI_UINT64_T,
      MPI_SUM,
      MPI_COMM_WORLD);
  requireOrThrow(
      rejection_count == static_cast<std::uint64_t>(world_size),
      "divergent TreePM layouts were not rejected on every rank");
  requireOrThrow(
      diagnostic.find("world metadata") != std::string::npos ||
          diagnostic.find("disagree") != std::string::npos ||
          diagnostic.find("peer rank rejected") != std::string::npos,
      "divergent TreePM layout rejection omitted the coordinated preflight diagnostic");

  const cosmosim::gravity::PmGridShape rank_local_shape =
      world_rank == 0 ? cosmosim::gravity::PmGridShape{8U, 8U, 8U}
                      : cosmosim::gravity::PmGridShape{10U, 8U, 8U};
  const auto rank_local_layout = cosmosim::parallel::makePmSlabLayout(
      rank_local_shape.nx,
      rank_local_shape.ny,
      rank_local_shape.nz,
      world_size,
      world_rank);
  bool shape_rejected = false;
  std::string shape_diagnostic;
  try {
    cosmosim::gravity::TreePmCoordinator coordinator(
        rank_local_shape, rank_local_layout);
    coordinator.solveActiveSet(
        empty, empty, empty, empty, accumulator, options, nullptr, nullptr);
  } catch (const std::exception& ex) {
    shape_rejected = true;
    shape_diagnostic = ex.what();
  }
  const std::uint64_t local_shape_rejection = shape_rejected ? 1U : 0U;
  std::uint64_t shape_rejection_count = 0U;
  MPI_Allreduce(
      &local_shape_rejection,
      &shape_rejection_count,
      1,
      MPI_UINT64_T,
      MPI_SUM,
      MPI_COMM_WORLD);
  requireOrThrow(
      shape_rejection_count == static_cast<std::uint64_t>(world_size),
      "divergent TreePM mesh shapes were not rejected on every rank");
  requireOrThrow(
      shape_diagnostic.find("mesh shape") != std::string::npos ||
          shape_diagnostic.find("diverged") != std::string::npos,
      "divergent TreePM mesh shape rejection omitted the protocol diagnostic");
}

void testDistributedShortRangeExportImportMatchesSingleRankReference() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    return;
  }

  constexpr double box_size_comoving = 1.0;
  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(16, 8);
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
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      2.0, rcutCellsForAvailableBackend(4.0, 3.9), mesh_spacing);
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
  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(16, 8);
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
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      2.0, rcutCellsForAvailableBackend(4.0, 3.9), mesh_spacing);
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

void testDistributedEmptyRanksIndependentTargetsAndMigration() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size < 2 || world_size > 4) {
    return;
  }

  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc;
  options.pm_options.enable_window_deconvolution = true;
  options.tree_options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError;
  options.tree_options.relative_force_tolerance = 0.005;
  options.tree_options.relative_force_acceleration_floor_code = 1.0e-12;
  options.tree_options.max_leaf_size = 1U;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 3.9, 0.125);
  options.decomposition_epoch = 7U;
  options.force_epoch = 11U;
  options.tree_exchange_batch_bytes = 256U;

  const std::size_t global_source_count = world_size == 2 ? 1U : 2U;
  std::vector<double> global_x(global_source_count, 0.0);
  std::vector<double> global_y(global_source_count, 0.0);
  std::vector<double> global_z(global_source_count, 0.0);
  std::vector<double> global_mass(global_source_count, 0.0);
  for (std::size_t i = 0; i < global_source_count; ++i) {
    global_x[i] = 0.992 - 0.014 * static_cast<double>(i);
    global_y[i] = 0.30 + 0.025 * static_cast<double>(i);
    global_z[i] = 0.70 - 0.018 * static_cast<double>(i);
    global_mass[i] = 0.9 + 0.2 * static_cast<double>(i);
  }

  std::vector<double> local_x;
  std::vector<double> local_y;
  std::vector<double> local_z;
  std::vector<double> local_mass;
  if (world_rank == world_size - 1) {
    local_x = global_x;
    local_y = global_y;
    local_z = global_z;
    local_mass = global_mass;
  }
  const std::vector<std::uint32_t> independent_target = world_rank == 0
      ? std::vector<std::uint32_t>{std::numeric_limits<std::uint32_t>::max()}
      : std::vector<std::uint32_t>{};
  const std::vector<double> target_x = world_rank == 0 ? std::vector<double>{0.008} : std::vector<double>{};
  const std::vector<double> target_y = world_rank == 0 ? std::vector<double>{0.31} : std::vector<double>{};
  const std::vector<double> target_z = world_rank == 0 ? std::vector<double>{0.69} : std::vector<double>{};
  const std::vector<double> previous_acceleration = world_rank == 0
      ? std::vector<double>{1.0e6}
      : std::vector<double>{};
  ForceField distributed_target_force{
      std::vector<double>(independent_target.size(), 0.0),
      std::vector<double>(independent_target.size(), 0.0),
      std::vector<double>(independent_target.size(), 0.0),
  };
  const cosmosim::gravity::TreePmForceAccumulatorView target_accumulator{
      .active_particle_index = independent_target,
      .accel_x_comoving = distributed_target_force.ax,
      .accel_y_comoving = distributed_target_force.ay,
      .accel_z_comoving = distributed_target_force.az,
      .previous_acceleration_magnitude_code = previous_acceleration,
      .target_pos_x_comoving = target_x,
      .target_pos_y_comoving = target_y,
      .target_pos_z_comoving = target_z,
  };
  const auto layout =
      cosmosim::parallel::makePmSlabLayout(pm_shape.nx, pm_shape.ny, pm_shape.nz, world_size, world_rank);
  cosmosim::gravity::TreePmCoordinator distributed(pm_shape, layout);
  cosmosim::gravity::TreePmDiagnostics distributed_diagnostics;
  distributed.solveActiveSet(
      local_x,
      local_y,
      local_z,
      local_mass,
      target_accumulator,
      options,
      nullptr,
      &distributed_diagnostics);

  const std::vector<std::uint32_t> reference_target{std::numeric_limits<std::uint32_t>::max()};
  const std::vector<double> reference_target_x{0.008};
  const std::vector<double> reference_target_y{0.31};
  const std::vector<double> reference_target_z{0.69};
  const std::vector<double> reference_previous{1.0e6};
  ForceField reference_force{{0.0}, {0.0}, {0.0}};
  const cosmosim::gravity::TreePmForceAccumulatorView reference_accumulator{
      .active_particle_index = reference_target,
      .accel_x_comoving = reference_force.ax,
      .accel_y_comoving = reference_force.ay,
      .accel_z_comoving = reference_force.az,
      .previous_acceleration_magnitude_code = reference_previous,
      .target_pos_x_comoving = reference_target_x,
      .target_pos_y_comoving = reference_target_y,
      .target_pos_z_comoving = reference_target_z,
  };
  cosmosim::gravity::TreePmCoordinator reference(pm_shape);
  reference.solveActiveSet(
      global_x, global_y, global_z, global_mass, reference_accumulator, options, nullptr, nullptr);

  double local_max_error = 0.0;
  if (world_rank == 0) {
    local_max_error = std::max({
        std::abs(distributed_target_force.ax[0] - reference_force.ax[0]),
        std::abs(distributed_target_force.ay[0] - reference_force.ay[0]),
        std::abs(distributed_target_force.az[0] - reference_force.az[0]),
    });
    requireOrThrow(
        distributed_diagnostics.residual_remote_request_packets > 0U,
        "empty-source rank with an independent target did not issue a remote tree request");
  } else {
    requireOrThrow(independent_target.empty(), "non-target rank unexpectedly owns an active target");
  }
  double global_max_error = 0.0;
  MPI_Allreduce(&local_max_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  requireOrThrow(global_max_error <= 2.0e-9, "empty-rank independent target differs from single-rank TreePM");

  // Exercise an all-empty collective round and distinguish explicit zero-source
  // descriptors from missing rank participation.
  const std::vector<double> empty;
  const std::vector<std::uint32_t> no_targets;
  ForceField no_force;
  const cosmosim::gravity::TreePmForceAccumulatorView empty_accumulator{
      .active_particle_index = no_targets,
      .accel_x_comoving = no_force.ax,
      .accel_y_comoving = no_force.ay,
      .accel_z_comoving = no_force.az,
  };
  cosmosim::gravity::TreePmCoordinator all_empty(pm_shape, layout);
  options.force_epoch = 12U;
  all_empty.solveActiveSet(empty, empty, empty, empty, empty_accumulator, options, nullptr, nullptr);

  // Reuse one coordinator across a decomposition transition in which every
  // nonzero rank becomes source-empty after the first force epoch.
  cosmosim::gravity::TreePmCoordinator migration(pm_shape, layout);
  std::vector<double> before_x{(static_cast<double>(world_rank) + 0.5) / static_cast<double>(world_size)};
  std::vector<double> before_y{0.2 + 0.1 * static_cast<double>(world_rank)};
  std::vector<double> before_z{0.4};
  std::vector<double> before_mass{1.0};
  options.decomposition_epoch = 20U;
  options.force_epoch = 20U;
  migration.solveActiveSet(
      before_x, before_y, before_z, before_mass, empty_accumulator, options, nullptr, nullptr);
  std::vector<double> after_x;
  std::vector<double> after_y;
  std::vector<double> after_z;
  std::vector<double> after_mass;
  if (world_rank == 0) {
    for (int rank = 0; rank < world_size; ++rank) {
      after_x.push_back((static_cast<double>(rank) + 0.5) / static_cast<double>(world_size));
      after_y.push_back(0.2 + 0.1 * static_cast<double>(rank));
      after_z.push_back(0.4);
      after_mass.push_back(1.0);
    }
  }
  options.decomposition_epoch = 21U;
  options.force_epoch = 21U;
  migration.solveActiveSet(
      after_x, after_y, after_z, after_mass, empty_accumulator, options, nullptr, nullptr);
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
  const cosmosim::gravity::PmGridShape coarse_shape = pmShapeForAvailableBackend(16, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(coarse_shape.nx);
  opts.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, rcutCellsForAvailableBackend(4.5, 3.9), mesh_spacing);
  opts.enable_zoom_long_range_correction = true;
  opts.zoom_focused_pm_shape = pmShapeForAvailableBackend(24, 8);
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
  std::vector<double> pos_x = {0.49, 0.52, 0.51, 0.20, 0.80};
  std::vector<double> pos_y = {0.50, 0.50, 0.53, 0.20, 0.80};
  std::vector<double> pos_z = {0.50, 0.50, 0.47, 0.20, 0.80};
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
  const cosmosim::gravity::PmGridShape coarse_shape = pmShapeForAvailableBackend(16, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(coarse_shape.nx);
  no_zoom.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, rcutCellsForAvailableBackend(4.5, 3.9), mesh_spacing);
  cosmosim::gravity::TreePmDiagnostics no_zoom_diag;
  const ForceField no_zoom_force =
      solveTreePm(pos_x, pos_y, pos_z, mass, coarse_shape, no_zoom, &no_zoom_diag);

  cosmosim::gravity::TreePmOptions with_zoom = no_zoom;
  with_zoom.enable_zoom_long_range_correction = true;
  with_zoom.zoom_focused_pm_shape = pmShapeForAvailableBackend(32, 8);
  with_zoom.source_is_high_res = source_is_high_res;
  with_zoom.active_is_high_res = active_is_high_res;
  with_zoom.zoom_region_center_x_comoving = 0.5;
  with_zoom.zoom_region_center_y_comoving = 0.5;
  with_zoom.zoom_region_center_z_comoving = 0.5;
  with_zoom.zoom_region_radius_comoving = 0.08;
  with_zoom.zoom_contamination_radius_comoving = 0.55;
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
  requireOrThrow(with_zoom_diag.zoom_low_res_contamination_count >= 1, "expected at least one low-res contaminant");
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
  testPeriodicTreeGeometryIsSeamSafeAndTranslationInvariant();
  testRelativeForceErrorMacUsesOwnerHistoryAcrossPeriodicSeam();
  testIndependentTargetAndEmptySourceContractSingleRank();
  testResidualCutoffPrunesTreePath();
  testPmOnlyAndTreePmConsistency();
  testCosmologicalTreePmNormalizationAndSplitComplementarity();
  testTreePmRejectsMismatchedGravitationalConstants();
  testPeriodicTreePmRejectsAmbiguousMinimumImageCutoff();
  testZoomMembershipPeriodicWrapAffectsCorrection();
  testZoomFocusedPmCorrectionAndContaminationDiagnostics();
  testActiveSubsetMatchesFullSolveSingleRank();
  testLongRangeCadenceCacheIsFailClosedAndOwnershipCompatible();
#if COSMOSIM_ENABLE_MPI
  testDivergentTreePmLayoutFailsCoordinately();
  testDistributedShortRangeExportImportMatchesSingleRankReference();
  testDistributedActiveSubsetMatchesSingleRankReference();
  testDistributedEmptyRanksIndependentTargetsAndMigration();
  MPI_Finalize();
#endif
  return 0;
}
