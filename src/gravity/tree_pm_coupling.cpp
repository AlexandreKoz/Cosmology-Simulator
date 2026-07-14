#include "cosmosim/gravity/tree_pm_coupling.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::gravity {
namespace {

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
bool queryActiveMpiWorld(int& world_size, int& world_rank) noexcept {
  world_size = 1;
  world_rank = 0;
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized == 0) {
    return false;
  }
  int finalized = 0;
  MPI_Finalized(&finalized);
  if (finalized != 0) {
    return false;
  }
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  return true;
}
#endif

[[nodiscard]] double minimumImageDelta(double delta, double box_size_comoving) {
  if (box_size_comoving <= 0.0) {
    return delta;
  }
  return delta - box_size_comoving * std::nearbyint(delta / box_size_comoving);
}

struct PeriodicBoxLengths {
  double lx = 0.0;
  double ly = 0.0;
  double lz = 0.0;
};

[[nodiscard]] PeriodicBoxLengths effectivePeriodicBoxLengths(const PmSolveOptions& options) {
  const double scalar = options.box_size_mpc_comoving;
  return PeriodicBoxLengths{
      .lx = options.box_size_x_mpc_comoving > 0.0 ? options.box_size_x_mpc_comoving : scalar,
      .ly = options.box_size_y_mpc_comoving > 0.0 ? options.box_size_y_mpc_comoving : scalar,
      .lz = options.box_size_z_mpc_comoving > 0.0 ? options.box_size_z_mpc_comoving : scalar};
}

// Map one periodic coordinate lane into the shortest contiguous interval that
// contains all sources. The interval begins immediately after the largest
// circular gap. This gives topology, COMs, quadrupoles, MAC geometry, and
// exported hierarchy bounds one coherent unwrapped frame.
void unwrapPeriodicAxis(
    std::span<const double> input,
    double box_size_comoving,
    std::vector<double>& output) {
  if (!std::isfinite(box_size_comoving) || box_size_comoving <= 0.0) {
    throw std::invalid_argument("Periodic TreePM tree geometry requires finite positive axis lengths");
  }
  output.resize(input.size());
  if (input.empty()) {
    return;
  }

  std::vector<double> wrapped(input.size(), 0.0);
  for (std::size_t i = 0; i < input.size(); ++i) {
    if (!std::isfinite(input[i])) {
      throw std::invalid_argument("Periodic TreePM tree geometry requires finite source coordinates");
    }
    double value = input[i] - box_size_comoving * std::floor(input[i] / box_size_comoving);
    if (value >= box_size_comoving) {
      value = 0.0;
    } else if (value < 0.0) {
      value += box_size_comoving;
    }
    wrapped[i] = value;
  }

  std::vector<double> ordered = wrapped;
  std::sort(ordered.begin(), ordered.end());
  double anchor = ordered.front();
  if (ordered.size() > 1U) {
    double largest_gap = -1.0;
    double best_anchor = ordered.front();
    for (std::size_t i = 0; i < ordered.size(); ++i) {
      const std::size_t next_index = (i + 1U) % ordered.size();
      const double next_value = next_index == 0U ? ordered.front() + box_size_comoving : ordered[next_index];
      const double gap = next_value - ordered[i];
      const double candidate_anchor = ordered[next_index];
      if (gap > largest_gap || (gap == largest_gap && candidate_anchor < best_anchor)) {
        largest_gap = gap;
        best_anchor = candidate_anchor;
      }
    }
    anchor = best_anchor;
  }

  for (std::size_t i = 0; i < wrapped.size(); ++i) {
    output[i] = wrapped[i] < anchor ? wrapped[i] + box_size_comoving : wrapped[i];
  }
}

[[nodiscard]] double norm3(double x, double y, double z) {
  return std::sqrt(x * x + y * y + z * z);
}

[[nodiscard]] bool forceAccumulatorShapeValid(const TreePmForceAccumulatorView& accumulator) {
  const std::size_t target_count = accumulator.active_particle_index.size();
  const bool explicit_targets_absent = accumulator.target_pos_x_comoving.empty() &&
      accumulator.target_pos_y_comoving.empty() && accumulator.target_pos_z_comoving.empty();
  const bool explicit_targets_complete = accumulator.target_pos_x_comoving.size() == target_count &&
      accumulator.target_pos_y_comoving.size() == target_count &&
      accumulator.target_pos_z_comoving.size() == target_count;
  return target_count == accumulator.accel_x_comoving.size() &&
      accumulator.active_particle_index.size() == accumulator.accel_y_comoving.size() &&
      accumulator.active_particle_index.size() == accumulator.accel_z_comoving.size() &&
      (accumulator.previous_acceleration_magnitude_code.empty() ||
       accumulator.active_particle_index.size() == accumulator.previous_acceleration_magnitude_code.size()) &&
      (explicit_targets_absent || explicit_targets_complete);
}

void validateTreePmPreflight(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    const TreePmForceAccumulatorView& accumulator,
    const TreePmOptions& options,
    const TreeSofteningView& softening_view) {
  if (!forceAccumulatorShapeValid(accumulator)) {
    throw std::invalid_argument(
        "TreePM force accumulator spans must have matching active-set extent");
  }
  if (pos_x_comoving.size() != pos_y_comoving.size() ||
      pos_x_comoving.size() != pos_z_comoving.size() ||
      pos_x_comoving.size() != mass_code.size()) {
    throw std::invalid_argument(
        "TreePM requires position and mass spans with equal extent");
  }
  validateTreePmSplitPolicy(options.split_policy);
  if (options.tree_exchange_batch_bytes == 0U) {
    throw std::invalid_argument("TreePM tree_exchange_batch_bytes must be > 0");
  }

  const PmSolveOptions& pm = options.pm_options;
  const PeriodicBoxLengths box = effectivePeriodicBoxLengths(pm);
  if (!std::isfinite(box.lx) || !std::isfinite(box.ly) ||
      !std::isfinite(box.lz) || box.lx <= 0.0 || box.ly <= 0.0 ||
      box.lz <= 0.0 || !std::isfinite(pm.scale_factor) ||
      pm.scale_factor <= 0.0 ||
      !std::isfinite(pm.gravitational_constant_code) ||
      pm.gravitational_constant_code <= 0.0 ||
      !std::isfinite(options.tree_options.gravitational_constant_code) ||
      options.tree_options.gravitational_constant_code <= 0.0) {
    throw std::invalid_argument(
        "TreePM requires finite positive box axes, scale factor, and gravitational constants");
  }
  if (options.tree_options.gravitational_constant_code !=
      pm.gravitational_constant_code) {
    throw std::invalid_argument(
        "TreePM PM and tree gravitational_constant_code values must match");
  }
  if (!std::isfinite(options.tree_options.opening_theta) ||
      options.tree_options.opening_theta <= 0.0 ||
      !std::isfinite(options.tree_options.relative_force_tolerance) ||
      options.tree_options.relative_force_tolerance <= 0.0 ||
      !std::isfinite(
          options.tree_options.relative_force_acceleration_floor_code) ||
      options.tree_options.relative_force_acceleration_floor_code <= 0.0 ||
      options.tree_options.max_leaf_size == 0U ||
      !std::isfinite(options.tree_options.softening.epsilon_comoving) ||
      options.tree_options.softening.epsilon_comoving < 0.0 ||
      options.tree_options.softening.kernel != TreeSofteningKernel::kPlummer) {
    throw std::invalid_argument("TreePM tree options are invalid");
  }
  switch (options.tree_options.opening_criterion) {
    case TreeOpeningCriterion::kBarnesHutGeometric:
    case TreeOpeningCriterion::kBarnesHutComDistance:
    case TreeOpeningCriterion::kRelativeForceError:
      break;
    default:
      throw std::invalid_argument("TreePM tree opening criterion is invalid");
  }
  switch (options.tree_options.multipole_order) {
    case TreeMultipoleOrder::kMonopole:
    case TreeMultipoleOrder::kQuadrupole:
      break;
    default:
      throw std::invalid_argument("TreePM tree multipole order is invalid");
  }
  switch (pm.assignment_scheme) {
    case PmAssignmentScheme::kCic:
    case PmAssignmentScheme::kTsc:
      break;
    default:
      throw std::invalid_argument("TreePM PM assignment scheme is invalid");
  }
  switch (pm.boundary_condition) {
    case PmBoundaryCondition::kPeriodic:
    case PmBoundaryCondition::kIsolatedOpen:
      break;
    default:
      throw std::invalid_argument("TreePM PM boundary condition is invalid");
  }
  const double half_shortest_periodic_axis =
      0.5 * std::min({box.lx, box.ly, box.lz});
  if (pm.boundary_condition == PmBoundaryCondition::kPeriodic &&
      options.split_policy.cutoff_radius_comoving >=
          std::nextafter(half_shortest_periodic_axis, 0.0)) {
    throw std::invalid_argument(
        "Periodic TreePM requires cutoff_radius_comoving < half the shortest box axis; "
        "cutoff_radius_comoving=" +
        std::to_string(options.split_policy.cutoff_radius_comoving) +
        ", half_shortest_axis=" +
        std::to_string(half_shortest_periodic_axis) +
        "; the short-range tree evaluates one minimum image per source");
  }

  for (std::size_t i = 0; i < mass_code.size(); ++i) {
    if (!std::isfinite(pos_x_comoving[i]) ||
        !std::isfinite(pos_y_comoving[i]) ||
        !std::isfinite(pos_z_comoving[i]) || !std::isfinite(mass_code[i]) ||
        mass_code[i] < 0.0) {
      throw std::invalid_argument(
          "TreePM sources require finite coordinates and finite non-negative masses");
    }
  }

  const bool explicit_targets = !accumulator.target_pos_x_comoving.empty();
  for (std::size_t i = 0; i < accumulator.active_particle_index.size(); ++i) {
    const std::uint32_t source_index = accumulator.active_particle_index[i];
    if (!explicit_targets && source_index >= pos_x_comoving.size()) {
      throw std::out_of_range("TreePM active index exceeds particle count");
    }
    if (explicit_targets && source_index >= pos_x_comoving.size() &&
        source_index != std::numeric_limits<std::uint32_t>::max()) {
      throw std::out_of_range(
          "TreePM explicit target source identity must be local or UINT32_MAX");
    }
    const double x = explicit_targets ? accumulator.target_pos_x_comoving[i]
                                      : pos_x_comoving[source_index];
    const double y = explicit_targets ? accumulator.target_pos_y_comoving[i]
                                      : pos_y_comoving[source_index];
    const double z = explicit_targets ? accumulator.target_pos_z_comoving[i]
                                      : pos_z_comoving[source_index];
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
      throw std::invalid_argument("TreePM target positions must be finite");
    }
  }

  const std::size_t source_count = pos_x_comoving.size();
  const std::size_t target_count = accumulator.active_particle_index.size();
  if ((!softening_view.source_particle_epsilon_comoving.empty() &&
       softening_view.source_particle_epsilon_comoving.size() != source_count) ||
      (!softening_view.source_particle_epsilon_override_mask.empty() &&
       softening_view.source_particle_epsilon_override_mask.size() != source_count) ||
      (!softening_view.source_species_tag.empty() &&
       softening_view.source_species_tag.size() != source_count) ||
      (!softening_view.target_particle_epsilon_comoving.empty() &&
       softening_view.target_particle_epsilon_comoving.size() != target_count) ||
      (!softening_view.target_particle_epsilon_override_mask.empty() &&
       softening_view.target_particle_epsilon_override_mask.size() != target_count) ||
      (!softening_view.target_species_tag.empty() &&
       softening_view.target_species_tag.size() != target_count)) {
    throw std::invalid_argument("TreePM softening/species sidecar extent mismatch");
  }
  const auto validate_epsilon_lane = [](std::span<const double> lane) {
    for (const double epsilon : lane) {
      if (!std::isfinite(epsilon) || epsilon < 0.0) {
        throw std::invalid_argument(
            "TreePM softening sidecars require finite non-negative values");
      }
    }
  };
  validate_epsilon_lane(softening_view.source_particle_epsilon_comoving);
  validate_epsilon_lane(softening_view.target_particle_epsilon_comoving);
  if (softening_view.species_policy.enabled) {
    for (const double epsilon :
         softening_view.species_policy.epsilon_comoving_by_species) {
      if (!std::isfinite(epsilon) || epsilon < 0.0) {
        throw std::invalid_argument(
            "TreePM species softening requires finite non-negative values");
      }
    }
  }
  if (explicit_targets &&
      ((!softening_view.source_particle_epsilon_comoving.empty() &&
        softening_view.target_particle_epsilon_comoving.empty()) ||
       (!softening_view.source_species_tag.empty() &&
        softening_view.target_species_tag.empty()))) {
    for (const std::uint32_t source_index : accumulator.active_particle_index) {
      if (source_index == std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument(
            "TreePM independent targets require target-owned softening/species sidecars");
      }
    }
  }

  if (options.enable_zoom_long_range_correction) {
    if (!options.zoom_focused_pm_shape.isValid() ||
        options.source_is_high_res.size() != source_count ||
        options.active_is_high_res.size() != target_count ||
        !std::isfinite(options.zoom_region_center_x_comoving) ||
        !std::isfinite(options.zoom_region_center_y_comoving) ||
        !std::isfinite(options.zoom_region_center_z_comoving) ||
        !std::isfinite(options.zoom_region_radius_comoving) ||
        !std::isfinite(options.zoom_contamination_radius_comoving) ||
        options.zoom_region_radius_comoving <= 0.0 ||
        options.zoom_contamination_radius_comoving < 0.0 ||
        options.zoom_high_res_allgather_limit_bytes == 0U) {
      throw std::invalid_argument("TreePM zoom correction options are invalid");
    }
  }
}

[[nodiscard]] std::uint64_t treePmOptionFingerprint(
    const PmGridShape& coordinator_pm_shape,
    const TreePmOptions& options,
    const PmSolveOptions& effective_pm_options,
    const TreeSofteningView& softening_view) {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&hash](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  const auto mix_double = [&mix](double value) {
    mix(std::bit_cast<std::uint64_t>(value));
  };
  mix(static_cast<std::uint64_t>(coordinator_pm_shape.nx));
  mix(static_cast<std::uint64_t>(coordinator_pm_shape.ny));
  mix(static_cast<std::uint64_t>(coordinator_pm_shape.nz));
  const PeriodicBoxLengths box = effectivePeriodicBoxLengths(effective_pm_options);
  mix_double(box.lx);
  mix_double(box.ly);
  mix_double(box.lz);
  mix_double(effective_pm_options.scale_factor);
  mix_double(effective_pm_options.gravitational_constant_code);
  mix_double(effective_pm_options.tree_pm_split_scale_comoving);
  mix(static_cast<std::uint64_t>(effective_pm_options.assignment_scheme));
  mix(effective_pm_options.enable_window_deconvolution ? 1U : 0U);
  mix(static_cast<std::uint64_t>(effective_pm_options.decomposition_mode));
  mix(static_cast<std::uint64_t>(effective_pm_options.boundary_condition));
  mix(static_cast<std::uint64_t>(effective_pm_options.execution_policy));
  mix(static_cast<std::uint64_t>(effective_pm_options.data_residency));
  mix(effective_pm_options.isolated_open_root_workspace_limit_bytes);
  mix(static_cast<std::uint64_t>(options.tree_options.opening_criterion));
  mix(static_cast<std::uint64_t>(options.tree_options.multipole_order));
  mix_double(options.tree_options.opening_theta);
  mix_double(options.tree_options.relative_force_tolerance);
  mix_double(options.tree_options.relative_force_acceleration_floor_code);
  mix_double(options.tree_options.gravitational_constant_code);
  mix(static_cast<std::uint64_t>(options.tree_options.max_leaf_size));
  mix(static_cast<std::uint64_t>(options.tree_options.softening.kernel));
  mix_double(options.tree_options.softening.epsilon_comoving);
  mix_double(options.split_policy.mesh_spacing_comoving);
  mix_double(options.split_policy.asmth_cells);
  mix_double(options.split_policy.rcut_cells);
  mix_double(options.split_policy.split_scale_comoving);
  mix_double(options.split_policy.cutoff_radius_comoving);
  mix(options.decomposition_epoch);
  mix(options.force_epoch);
  mix(options.tree_exchange_batch_bytes);
  mix(options.enable_zoom_long_range_correction ? 1U : 0U);
  mix(static_cast<std::uint64_t>(options.zoom_focused_pm_shape.nx));
  mix(static_cast<std::uint64_t>(options.zoom_focused_pm_shape.ny));
  mix(static_cast<std::uint64_t>(options.zoom_focused_pm_shape.nz));
  mix_double(options.zoom_region_center_x_comoving);
  mix_double(options.zoom_region_center_y_comoving);
  mix_double(options.zoom_region_center_z_comoving);
  mix_double(options.zoom_region_radius_comoving);
  mix_double(options.zoom_contamination_radius_comoving);
  mix(options.zoom_high_res_allgather_limit_bytes);
  mix(softening_view.species_policy.enabled ? 1U : 0U);
  for (const double epsilon :
       softening_view.species_policy.epsilon_comoving_by_species) {
    mix_double(epsilon);
  }
  return hash;
}

[[nodiscard]] bool acceptNodeByMac(
    bool is_leaf,
    bool target_inside_node,
    double half_size,
    double com_center_offset,
    double node_mass_code,
    double r2,
    bool previous_acceleration_available,
    double previous_acceleration_magnitude_code,
    const TreeGravityOptions& options) {
  if (is_leaf) {
    return true;
  }
  const double r = std::sqrt(r2 + 1.0e-30);
  const double width = 2.0 * half_size;
  if (options.opening_criterion == TreeOpeningCriterion::kBarnesHutGeometric) {
    return (width / r) < options.opening_theta;
  }
  if (options.opening_criterion == TreeOpeningCriterion::kBarnesHutComDistance) {
    const double effective_size = width + com_center_offset;
    return (effective_size / r) < options.opening_theta;
  }
  if (target_inside_node) {
    return false;
  }
  if (!previous_acceleration_available) {
    const double effective_size = width + com_center_offset;
    return (effective_size / r) < options.opening_theta;
  }
  const double acceleration_scale_code = std::max(
      std::abs(previous_acceleration_magnitude_code), options.relative_force_acceleration_floor_code);
  const double estimated_error_scale =
      options.gravitational_constant_code * node_mass_code * width * width;
  const double allowed_error_scale =
      options.relative_force_tolerance * acceleration_scale_code * r2 * r2;
  return estimated_error_scale <= allowed_error_scale;
}


[[nodiscard]] bool passesSofteningEnvelopeGuard(
    bool is_leaf,
    double half_size,
    double r,
    double target_softening_comoving,
    double node_softening_min_comoving,
    double node_softening_max_comoving) {
  if (is_leaf) {
    return true;
  }
  const double heterogeneity = node_softening_max_comoving - node_softening_min_comoving;
  if (heterogeneity <= 1.0e-12) {
    return true;
  }
  const double pair_softening_max =
      combineSofteningPairEpsilon(node_softening_max_comoving, target_softening_comoving);
  const double envelope_radius = 2.0 * half_size + 2.0 * pair_softening_max;
  return r > envelope_radius;
}

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
[[nodiscard]] std::uint64_t checkedZoomGatherBytes(
    std::size_t particle_count,
    std::string_view context) {
  constexpr std::uint64_t k_fields = 4U;
  if (particle_count > std::numeric_limits<std::uint64_t>::max() / (k_fields * sizeof(double))) {
    throw std::overflow_error(std::string(context) + " byte estimate overflows uint64_t");
  }
  return static_cast<std::uint64_t>(particle_count) * k_fields * sizeof(double);
}

struct GatheredParticleField {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> m;
};

[[nodiscard]] GatheredParticleField allGatherParticleField(
    std::span<const double> x,
    std::span<const double> y,
    std::span<const double> z,
    std::span<const double> m,
    std::uint64_t gather_limit_bytes,
    std::uint64_t* gathered_bytes = nullptr) {
  if (x.size() != y.size() || x.size() != z.size() || x.size() != m.size()) {
    throw std::invalid_argument("allGatherParticleField requires equal local span lengths");
  }
  int world_size = 1;
  int world_rank = 0;
  if (!queryActiveMpiWorld(world_size, world_rank)) {
    throw std::runtime_error("TreePM zoom high-res all-gather requires an active MPI runtime");
  }
  if (x.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("TreePM zoom high-res all-gather local particle count exceeds MPI int limit");
  }
  const int local_count = static_cast<int>(x.size());
  std::vector<int> counts(static_cast<std::size_t>(world_size), 0);
  MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  std::vector<int> displs(static_cast<std::size_t>(world_size), 0);
  int total = 0;
  for (int i = 0; i < world_size; ++i) {
    if (counts[static_cast<std::size_t>(i)] < 0 ||
        counts[static_cast<std::size_t>(i)] > std::numeric_limits<int>::max() - total) {
      throw std::overflow_error("TreePM zoom high-res all-gather count/displacement exceeds MPI int limit");
    }
    displs[static_cast<std::size_t>(i)] = total;
    total += counts[static_cast<std::size_t>(i)];
  }
  const std::uint64_t total_bytes =
      checkedZoomGatherBytes(static_cast<std::size_t>(total), "TreePM zoom high-res all-gather");
  if (gathered_bytes != nullptr) {
    *gathered_bytes = total_bytes;
  }
  if (gather_limit_bytes == 0U || total_bytes > gather_limit_bytes) {
    throw std::runtime_error(
        "TreePM zoom high-res source all-gather guard exceeded on rank " + std::to_string(world_rank) +
        ": gathered_high_res_sources=" + std::to_string(total) +
        " ranks=" + std::to_string(world_size) +
        " estimated_bytes=" + std::to_string(total_bytes) +
        " configured_limit_bytes=" + std::to_string(gather_limit_bytes) +
        " route=allgather policy=bounded_zoom_correction_only");
  }
  GatheredParticleField out;
  out.x.resize(static_cast<std::size_t>(total));
  out.y.resize(static_cast<std::size_t>(total));
  out.z.resize(static_cast<std::size_t>(total));
  out.m.resize(static_cast<std::size_t>(total));
  MPI_Allgatherv(x.data(), local_count, MPI_DOUBLE, out.x.data(), counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgatherv(y.data(), local_count, MPI_DOUBLE, out.y.data(), counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgatherv(z.data(), local_count, MPI_DOUBLE, out.z.data(), counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgatherv(m.data(), local_count, MPI_DOUBLE, out.m.data(), counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  return out;
}
#endif

[[nodiscard]] std::array<double, 3> monopolePlusQuadrupoleAccelPeriodic(
    const TreeNodeSoa& nodes,
    std::uint32_t node_index,
    double dx,
    double dy,
    double dz,
    double target_softening_comoving,
    const TreeGravityOptions& options,
    double split_scale_comoving) {
  const double r2 = dx * dx + dy * dy + dz * dz;
  const double r = std::sqrt(std::max(r2, 1.0e-30));
  const double pair_epsilon =
      combineSofteningPairEpsilon(nodes.softening_max_comoving[node_index], target_softening_comoving);
  const double eps2 = pair_epsilon * pair_epsilon;
  const double denom = std::max(r2 + eps2, 1.0e-30);
  const double softened_inv_r3 = 1.0 / (denom * std::sqrt(denom));
  const double split_factor = treePmGaussianShortRangeForceFactor(r, split_scale_comoving);
  const double prefactor = options.gravitational_constant_code;

  double ax = prefactor * nodes.mass_code[node_index] * split_factor * softened_inv_r3 * dx;
  double ay = prefactor * nodes.mass_code[node_index] * split_factor * softened_inv_r3 * dy;
  double az = prefactor * nodes.mass_code[node_index] * split_factor * softened_inv_r3 * dz;
  if (options.multipole_order != TreeMultipoleOrder::kQuadrupole) {
    return {ax, ay, az};
  }

  const double qxx = nodes.quad_xx[node_index];
  const double qxy = nodes.quad_xy[node_index];
  const double qxz = nodes.quad_xz[node_index];
  const double qyy = nodes.quad_yy[node_index];
  const double qyz = nodes.quad_yz[node_index];
  const double qzz = nodes.quad_zz[node_index];
  const double qrx = qxx * dx + qxy * dy + qxz * dz;
  const double qry = qxy * dx + qyy * dy + qyz * dz;
  const double qrz = qxz * dx + qyz * dy + qzz * dz;
  const double rqr = dx * qrx + dy * qry + dz * qrz;
  (void)qrx;
  (void)qry;
  (void)qrz;
  (void)rqr;

  // A screened/softened radial kernel is not harmonic, so its second-order
  // expansion needs both the traceless quadrupole and the trace of the raw
  // central second moment.  Multiplying the Newtonian quadrupole by the split
  // factor would omit derivatives of that factor and is not a valid TreePM
  // multipole expansion near the split scale.
  const double moment_trace = nodes.second_moment_trace[node_index];
  const double mxx = (qxx + moment_trace) / 3.0;
  const double mxy = qxy / 3.0;
  const double mxz = qxz / 3.0;
  const double myy = (qyy + moment_trace) / 3.0;
  const double myz = qyz / 3.0;
  const double mzz = (qzz + moment_trace) / 3.0;
  const double mdx = mxx * dx + mxy * dy + mxz * dz;
  const double mdy = mxy * dx + myy * dy + myz * dz;
  const double mdz = mxz * dx + myz * dy + mzz * dz;
  const double dmd = dx * mdx + dy * mdy + dz * mdz;

  constexpr double inv_sqrt_pi = 0.564189583547756286948079451560772586;
  const double split_scale2 = split_scale_comoving * split_scale_comoving;
  const double exponential = std::exp(-r2 / (4.0 * split_scale2));
  const double split_first =
      -0.5 * inv_sqrt_pi * r2 * exponential /
      (split_scale_comoving * split_scale2);
  const double split_second = inv_sqrt_pi * exponential *
      (-r / (split_scale_comoving * split_scale2) +
       0.25 * r * r2 / (split_scale_comoving * split_scale2 * split_scale2));
  const double softened_first = -3.0 * r / (denom * denom * std::sqrt(denom));
  const double softened_second =
      -3.0 / (denom * denom * std::sqrt(denom)) +
      15.0 * r2 / (denom * denom * denom * std::sqrt(denom));
  const double radial_first = split_first * softened_inv_r3 + split_factor * softened_first;
  const double radial_second = split_second * softened_inv_r3 +
      2.0 * split_first * softened_first + split_factor * softened_second;
  const double first_over_r = radial_first / r;
  const double contracted_radial = radial_second / r2 - radial_first / (r2 * r);
  ax += prefactor *
      (mdx * first_over_r + 0.5 * moment_trace * dx * first_over_r +
       0.5 * dx * dmd * contracted_radial);
  ay += prefactor *
      (mdy * first_over_r + 0.5 * moment_trace * dy * first_over_r +
       0.5 * dy * dmd * contracted_radial);
  az += prefactor *
      (mdz * first_over_r + 0.5 * moment_trace * dz * first_over_r +
       0.5 * dz * dmd * contracted_radial);
  return {ax, ay, az};
}

void pushChildrenNearFirstPeriodic(
    const TreeNodeSoa& nodes,
    std::uint32_t node_index,
    double px,
    double py,
    double pz,
    const PeriodicBoxLengths& box_lengths,
    std::vector<std::uint32_t>& stack) {
  std::array<std::pair<double, std::uint32_t>, 8> child_dist2{};
  std::size_t count = 0;
  const std::size_t child_offset = static_cast<std::size_t>(node_index) * 8U;
  for (std::uint8_t octant = 0; octant < 8U; ++octant) {
    const std::uint32_t child = nodes.child_index[child_offset + octant];
    if (child == std::numeric_limits<std::uint32_t>::max()) {
      continue;
    }
    const double dx = minimumImageDelta(nodes.center_x_comoving[child] - px, box_lengths.lx);
    const double dy = minimumImageDelta(nodes.center_y_comoving[child] - py, box_lengths.ly);
    const double dz = minimumImageDelta(nodes.center_z_comoving[child] - pz, box_lengths.lz);
    child_dist2[count++] = {dx * dx + dy * dy + dz * dz, child};
  }
  std::sort(child_dist2.begin(), child_dist2.begin() + static_cast<std::ptrdiff_t>(count),
      [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
  for (std::size_t i = count; i > 0; --i) {
    stack.push_back(child_dist2[i - 1].second);
  }
}

void resizeCompactSidecars(std::vector<double>& first, std::vector<double>& second, std::vector<double>& third, std::size_t size) {
  first.resize(size);
  second.resize(size);
  third.resize(size);
}

[[nodiscard]] double l2NormFromComponents(
    std::span<const double> ax,
    std::span<const double> ay,
    std::span<const double> az) {
  double sum = 0.0;
  for (std::size_t i = 0; i < ax.size(); ++i) {
    sum += ax[i] * ax[i] + ay[i] * ay[i] + az[i] * az[i];
  }
  return std::sqrt(sum);
}

void appendWireU32(std::vector<std::uint8_t>& bytes, std::uint32_t value) {
  for (unsigned shift = 0; shift < 32U; shift += 8U) {
    bytes.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

void appendWireU64(std::vector<std::uint8_t>& bytes, std::uint64_t value) {
  for (unsigned shift = 0; shift < 64U; shift += 8U) {
    bytes.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

[[nodiscard]] std::uint32_t readWireU32(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  if (bytes.size() - std::min(bytes.size(), offset) < sizeof(std::uint32_t)) {
    throw std::runtime_error("TreePM short-range wire record is truncated");
  }
  std::uint32_t value = 0U;
  for (unsigned shift = 0; shift < 32U; shift += 8U) {
    value |= static_cast<std::uint32_t>(bytes[offset++]) << shift;
  }
  return value;
}

[[nodiscard]] std::uint64_t readWireU64(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  if (bytes.size() - std::min(bytes.size(), offset) < sizeof(std::uint64_t)) {
    throw std::runtime_error("TreePM short-range wire record is truncated");
  }
  std::uint64_t value = 0U;
  for (unsigned shift = 0; shift < 64U; shift += 8U) {
    value |= static_cast<std::uint64_t>(bytes[offset++]) << shift;
  }
  return value;
}

void appendWireDouble(std::vector<std::uint8_t>& bytes, double value) {
  appendWireU64(bytes, std::bit_cast<std::uint64_t>(value));
}

[[nodiscard]] double readWireDouble(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  return std::bit_cast<double>(readWireU64(bytes, offset));
}

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
[[nodiscard]] int checkedMpiByteCount(
    std::size_t record_count,
    std::size_t record_bytes,
    std::string_view context) {
  if (record_bytes != 0U &&
      record_count > static_cast<std::size_t>(std::numeric_limits<int>::max()) / record_bytes) {
    throw std::overflow_error(std::string(context) + " exceeds MPI int byte-count capacity");
  }
  return static_cast<int>(record_count * record_bytes);
}

[[nodiscard]] int populateMpiByteDisplacements(
    std::span<const int> counts,
    std::span<int> displacements,
    std::string_view context) {
  if (counts.size() != displacements.size()) {
    throw std::invalid_argument(std::string(context) + " count/displacement vector sizes differ");
  }
  int total_bytes = 0;
  for (std::size_t peer = 0; peer < counts.size(); ++peer) {
    const int count = counts[peer];
    if (count < 0) {
      throw std::runtime_error(std::string(context) + " contains a negative byte count");
    }
    displacements[peer] = total_bytes;
    if (count > std::numeric_limits<int>::max() - total_bytes) {
      throw std::overflow_error(std::string(context) + " total byte count exceeds MPI int capacity");
    }
    total_bytes += count;
  }
  return total_bytes;
}
#endif

[[nodiscard]] double minimumDistanceToNodeAabb(
    double px,
    double py,
    double pz,
    double cx,
    double cy,
    double cz,
    double half_size_comoving,
    const PeriodicBoxLengths& box_lengths) {
  const double dx_abs = std::abs(minimumImageDelta(cx - px, box_lengths.lx));
  const double dy_abs = std::abs(minimumImageDelta(cy - py, box_lengths.ly));
  const double dz_abs = std::abs(minimumImageDelta(cz - pz, box_lengths.lz));
  const double ex = std::max(0.0, dx_abs - half_size_comoving);
  const double ey = std::max(0.0, dy_abs - half_size_comoving);
  const double ez = std::max(0.0, dz_abs - half_size_comoving);
  return std::sqrt(ex * ex + ey * ey + ez * ez);
}

[[nodiscard]] double maximumDistanceToNodeAabb(
    double px,
    double py,
    double pz,
    double cx,
    double cy,
    double cz,
    double half_size_comoving,
    const PeriodicBoxLengths& box_lengths) {
  const double dx_abs = std::abs(minimumImageDelta(cx - px, box_lengths.lx));
  const double dy_abs = std::abs(minimumImageDelta(cy - py, box_lengths.ly));
  const double dz_abs = std::abs(minimumImageDelta(cz - pz, box_lengths.lz));
  const double max_x = dx_abs + half_size_comoving;
  const double max_y = dy_abs + half_size_comoving;
  const double max_z = dz_abs + half_size_comoving;
  return std::sqrt(max_x * max_x + max_y * max_y + max_z * max_z);
}

struct ShortRangeTargetRequestPacket {
  std::uint32_t wire_version = 1U;
  int origin_rank = -1;
  int destination_rank = -1;
  std::uint32_t flags = 0U;
  std::uint32_t batch_token = 0;
  std::uint32_t request_id = 0;
  std::uint64_t exchange_sequence = 0U;
  std::uint64_t decomposition_epoch = 0U;
  std::uint64_t force_epoch = 0U;
  std::uint64_t target_identity = 0U;
  double target_x_comoving = 0.0;
  double target_y_comoving = 0.0;
  double target_z_comoving = 0.0;
  double target_softening_epsilon_comoving = 0.0;
  double previous_acceleration_magnitude = 0.0;
};

struct ShortRangeTargetResponsePacket {
  std::uint32_t wire_version = 1U;
  int target_owner_rank = -1;
  int source_rank = -1;
  std::uint32_t flags = 0U;
  std::uint32_t batch_token = 0;
  std::uint32_t request_id = 0;
  std::uint64_t exchange_sequence = 0U;
  std::uint64_t decomposition_epoch = 0U;
  std::uint64_t force_epoch = 0U;
  std::uint64_t target_identity = 0U;
  double accel_x_comoving = 0.0;
  double accel_y_comoving = 0.0;
  double accel_z_comoving = 0.0;
};

constexpr std::uint32_t k_short_range_wire_version = 1U;
constexpr std::uint32_t k_short_range_flag_previous_acceleration = 1U << 0U;
constexpr std::size_t k_short_range_request_wire_bytes = 96U;
constexpr std::size_t k_short_range_response_wire_bytes = 80U;

[[nodiscard]] std::vector<std::uint8_t> encodeShortRangeRequests(
    std::span<const ShortRangeTargetRequestPacket> records) {
  if (records.size() > std::numeric_limits<std::size_t>::max() / k_short_range_request_wire_bytes) {
    throw std::overflow_error("TreePM short-range request serialization size overflows size_t");
  }
  std::vector<std::uint8_t> bytes;
  bytes.reserve(records.size() * k_short_range_request_wire_bytes);
  for (const ShortRangeTargetRequestPacket& record : records) {
    appendWireU32(bytes, record.wire_version);
    appendWireU32(bytes, static_cast<std::uint32_t>(record.origin_rank));
    appendWireU32(bytes, static_cast<std::uint32_t>(record.destination_rank));
    appendWireU32(bytes, record.flags);
    appendWireU32(bytes, record.batch_token);
    appendWireU32(bytes, record.request_id);
    appendWireU64(bytes, record.exchange_sequence);
    appendWireU64(bytes, record.decomposition_epoch);
    appendWireU64(bytes, record.force_epoch);
    appendWireU64(bytes, record.target_identity);
    appendWireDouble(bytes, record.target_x_comoving);
    appendWireDouble(bytes, record.target_y_comoving);
    appendWireDouble(bytes, record.target_z_comoving);
    appendWireDouble(bytes, record.target_softening_epsilon_comoving);
    appendWireDouble(bytes, record.previous_acceleration_magnitude);
  }
  return bytes;
}

[[nodiscard]] std::vector<ShortRangeTargetRequestPacket> decodeShortRangeRequests(
    std::span<const std::uint8_t> bytes) {
  if (bytes.size() % k_short_range_request_wire_bytes != 0U) {
    throw std::runtime_error("TreePM short-range request wire payload is misaligned");
  }
  std::vector<ShortRangeTargetRequestPacket> records;
  records.reserve(bytes.size() / k_short_range_request_wire_bytes);
  std::size_t offset = 0U;
  while (offset < bytes.size()) {
    ShortRangeTargetRequestPacket record;
    record.wire_version = readWireU32(bytes, offset);
    record.origin_rank = static_cast<int>(readWireU32(bytes, offset));
    record.destination_rank = static_cast<int>(readWireU32(bytes, offset));
    record.flags = readWireU32(bytes, offset);
    record.batch_token = readWireU32(bytes, offset);
    record.request_id = readWireU32(bytes, offset);
    record.exchange_sequence = readWireU64(bytes, offset);
    record.decomposition_epoch = readWireU64(bytes, offset);
    record.force_epoch = readWireU64(bytes, offset);
    record.target_identity = readWireU64(bytes, offset);
    record.target_x_comoving = readWireDouble(bytes, offset);
    record.target_y_comoving = readWireDouble(bytes, offset);
    record.target_z_comoving = readWireDouble(bytes, offset);
    record.target_softening_epsilon_comoving = readWireDouble(bytes, offset);
    record.previous_acceleration_magnitude = readWireDouble(bytes, offset);
    records.push_back(record);
  }
  return records;
}

[[nodiscard]] std::vector<std::uint8_t> encodeShortRangeResponses(
    std::span<const ShortRangeTargetResponsePacket> records) {
  if (records.size() > std::numeric_limits<std::size_t>::max() / k_short_range_response_wire_bytes) {
    throw std::overflow_error("TreePM short-range response serialization size overflows size_t");
  }
  std::vector<std::uint8_t> bytes;
  bytes.reserve(records.size() * k_short_range_response_wire_bytes);
  for (const ShortRangeTargetResponsePacket& record : records) {
    appendWireU32(bytes, record.wire_version);
    appendWireU32(bytes, static_cast<std::uint32_t>(record.target_owner_rank));
    appendWireU32(bytes, static_cast<std::uint32_t>(record.source_rank));
    appendWireU32(bytes, record.flags);
    appendWireU32(bytes, record.batch_token);
    appendWireU32(bytes, record.request_id);
    appendWireU64(bytes, record.exchange_sequence);
    appendWireU64(bytes, record.decomposition_epoch);
    appendWireU64(bytes, record.force_epoch);
    appendWireU64(bytes, record.target_identity);
    appendWireDouble(bytes, record.accel_x_comoving);
    appendWireDouble(bytes, record.accel_y_comoving);
    appendWireDouble(bytes, record.accel_z_comoving);
  }
  return bytes;
}

[[nodiscard]] std::vector<ShortRangeTargetResponsePacket> decodeShortRangeResponses(
    std::span<const std::uint8_t> bytes) {
  if (bytes.size() % k_short_range_response_wire_bytes != 0U) {
    throw std::runtime_error("TreePM short-range response wire payload is misaligned");
  }
  std::vector<ShortRangeTargetResponsePacket> records;
  records.reserve(bytes.size() / k_short_range_response_wire_bytes);
  std::size_t offset = 0U;
  while (offset < bytes.size()) {
    ShortRangeTargetResponsePacket record;
    record.wire_version = readWireU32(bytes, offset);
    record.target_owner_rank = static_cast<int>(readWireU32(bytes, offset));
    record.source_rank = static_cast<int>(readWireU32(bytes, offset));
    record.flags = readWireU32(bytes, offset);
    record.batch_token = readWireU32(bytes, offset);
    record.request_id = readWireU32(bytes, offset);
    record.exchange_sequence = readWireU64(bytes, offset);
    record.decomposition_epoch = readWireU64(bytes, offset);
    record.force_epoch = readWireU64(bytes, offset);
    record.target_identity = readWireU64(bytes, offset);
    record.accel_x_comoving = readWireDouble(bytes, offset);
    record.accel_y_comoving = readWireDouble(bytes, offset);
    record.accel_z_comoving = readWireDouble(bytes, offset);
    records.push_back(record);
  }
  return records;
}

struct SourceDomainBoundsPacket {
  double min_x_comoving = 0.0;
  double max_x_comoving = 0.0;
  double min_y_comoving = 0.0;
  double max_y_comoving = 0.0;
  double min_z_comoving = 0.0;
  double max_z_comoving = 0.0;
  std::uint8_t wraps_x = 0;
  std::uint8_t wraps_y = 0;
  std::uint8_t wraps_z = 0;
  std::uint8_t reserved0 = 0;
  std::uint64_t source_particle_count = 0;
};

static_assert(std::is_trivially_copyable_v<SourceDomainBoundsPacket>);

struct PeriodicInterval {
  double min = 0.0;
  double max = 0.0;
  bool wraps = false;
};

[[nodiscard]] PeriodicInterval tightPeriodicInterval(std::span<const double> values, double box_size_comoving) {
  PeriodicInterval interval{};
  if (values.empty()) {
    return interval;
  }
  if (box_size_comoving <= 0.0) {
    interval.min = *std::min_element(values.begin(), values.end());
    interval.max = *std::max_element(values.begin(), values.end());
    interval.wraps = false;
    return interval;
  }
  std::vector<double> wrapped(values.begin(), values.end());
  for (double& value : wrapped) {
    value = std::fmod(value, box_size_comoving);
    if (value < 0.0) {
      value += box_size_comoving;
    }
  }
  std::sort(wrapped.begin(), wrapped.end());
  if (wrapped.size() == 1) {
    interval.min = wrapped.front();
    interval.max = wrapped.front();
    interval.wraps = false;
    return interval;
  }
  std::size_t largest_gap_start = 0;
  double largest_gap = -1.0;
  for (std::size_t i = 0; i < wrapped.size(); ++i) {
    const std::size_t j = (i + 1U) % wrapped.size();
    const double next = (j == 0) ? (wrapped[j] + box_size_comoving) : wrapped[j];
    const double gap = next - wrapped[i];
    if (gap > largest_gap) {
      largest_gap = gap;
      largest_gap_start = i;
    }
  }
  const std::size_t interval_begin_index = (largest_gap_start + 1U) % wrapped.size();
  const std::size_t interval_end_index = largest_gap_start;
  interval.min = wrapped[interval_begin_index];
  interval.max = wrapped[interval_end_index];
  interval.wraps = interval.min > interval.max;
  return interval;
}

[[nodiscard]] double minimumDistanceToPeriodicInterval(
    double coordinate,
    double interval_min,
    double interval_max,
    bool interval_wraps,
    double box_size_comoving) {
  if (interval_max < interval_min && !interval_wraps) {
    return 0.0;
  }
  auto interval_distance = [](double value, double lower, double upper) {
    if (value < lower) {
      return lower - value;
    }
    if (value > upper) {
      return value - upper;
    }
    return 0.0;
  };
  if (interval_wraps) {
    if (box_size_comoving <= 0.0) {
      return 0.0;
    }
    const double d0 = interval_distance(coordinate, interval_min, box_size_comoving);
    const double d1 = interval_distance(coordinate, 0.0, interval_max);
    return std::min(d0, d1);
  }
  if (box_size_comoving <= 0.0) {
    return interval_distance(coordinate, interval_min, interval_max);
  }
  const double width = interval_max - interval_min;
  if (width >= box_size_comoving) {
    return 0.0;
  }
  const double center = 0.5 * (interval_min + interval_max);
  const double center_distance = std::abs(minimumImageDelta(coordinate - center, box_size_comoving));
  return std::max(0.0, center_distance - 0.5 * width);
}

[[nodiscard]] double minimumDistanceToPeriodicBounds(
    double px,
    double py,
    double pz,
    const SourceDomainBoundsPacket& bounds,
    const PeriodicBoxLengths& box_lengths) {
  if (bounds.source_particle_count == 0) {
    return std::numeric_limits<double>::infinity();
  }
  const double dx = minimumDistanceToPeriodicInterval(
      px,
      bounds.min_x_comoving,
      bounds.max_x_comoving,
      bounds.wraps_x != 0U,
      box_lengths.lx);
  const double dy = minimumDistanceToPeriodicInterval(
      py,
      bounds.min_y_comoving,
      bounds.max_y_comoving,
      bounds.wraps_y != 0U,
      box_lengths.ly);
  const double dz = minimumDistanceToPeriodicInterval(
      pz,
      bounds.min_z_comoving,
      bounds.max_z_comoving,
      bounds.wraps_z != 0U,
      box_lengths.lz);
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

[[nodiscard]] SourceDomainBoundsPacket computeLocalSourceBounds(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    const PeriodicBoxLengths& box_lengths) {
  SourceDomainBoundsPacket bounds;
  bounds.source_particle_count = static_cast<std::uint64_t>(pos_x_comoving.size());
  if (pos_x_comoving.empty()) {
    return bounds;
  }
  const PeriodicInterval interval_x = tightPeriodicInterval(pos_x_comoving, box_lengths.lx);
  const PeriodicInterval interval_y = tightPeriodicInterval(pos_y_comoving, box_lengths.ly);
  const PeriodicInterval interval_z = tightPeriodicInterval(pos_z_comoving, box_lengths.lz);
  bounds.min_x_comoving = interval_x.min;
  bounds.max_x_comoving = interval_x.max;
  bounds.min_y_comoving = interval_y.min;
  bounds.max_y_comoving = interval_y.max;
  bounds.min_z_comoving = interval_z.min;
  bounds.max_z_comoving = interval_z.max;
  bounds.wraps_x = interval_x.wraps ? 1U : 0U;
  bounds.wraps_y = interval_y.wraps ? 1U : 0U;
  bounds.wraps_z = interval_z.wraps ? 1U : 0U;
  return bounds;
}

[[nodiscard]] parallel::TreePseudoParticlePacket makeLocalTreePseudoParticlePacket(
    int world_rank,
    std::uint64_t decomposition_epoch,
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    const PeriodicBoxLengths& box_lengths) {
  const SourceDomainBoundsPacket bounds = computeLocalSourceBounds(
      pos_x_comoving, pos_y_comoving, pos_z_comoving, box_lengths);
  double mass_sum = 0.0;
  double mass_x = 0.0;
  double mass_y = 0.0;
  double mass_z = 0.0;
  for (std::size_t i = 0; i < mass_code.size(); ++i) {
    const double mass = std::max(0.0, mass_code[i]);
    mass_sum += mass;
    mass_x += mass * pos_x_comoving[i];
    mass_y += mass * pos_y_comoving[i];
    mass_z += mass * pos_z_comoving[i];
  }
  const double inv_mass = (mass_sum > 0.0) ? (1.0 / mass_sum) : 0.0;
  parallel::TreePseudoParticlePacket packet;
  packet.descriptor = parallel::TreePseudoParticleDescriptor{
      .pseudo_particle_id = (static_cast<std::uint64_t>(std::max(world_rank, 0)) << 32U) ^ decomposition_epoch,
      .source_rank = world_rank,
      .decomposition_epoch = decomposition_epoch,
      .derived_not_authoritative = true,
  };
  packet.mass_code = mass_sum;
  packet.center_x_comoving = mass_x * inv_mass;
  packet.center_y_comoving = mass_y * inv_mass;
  packet.center_z_comoving = mass_z * inv_mass;
  packet.min_x_comoving = bounds.min_x_comoving;
  packet.max_x_comoving = bounds.max_x_comoving;
  packet.min_y_comoving = bounds.min_y_comoving;
  packet.max_y_comoving = bounds.max_y_comoving;
  packet.min_z_comoving = bounds.min_z_comoving;
  packet.max_z_comoving = bounds.max_z_comoving;
  packet.source_count = bounds.source_particle_count;
  parallel::validateTreePseudoParticlePacket(packet);
  return packet;
}

[[nodiscard]] SourceDomainBoundsPacket boundsFromTreePseudoParticlePacket(
    const parallel::TreePseudoParticlePacket& packet) {
  parallel::validateTreePseudoParticlePacket(packet);
  SourceDomainBoundsPacket bounds;
  bounds.min_x_comoving = packet.min_x_comoving;
  bounds.max_x_comoving = packet.max_x_comoving;
  bounds.min_y_comoving = packet.min_y_comoving;
  bounds.max_y_comoving = packet.max_y_comoving;
  bounds.min_z_comoving = packet.min_z_comoving;
  bounds.max_z_comoving = packet.max_z_comoving;
  bounds.source_particle_count = packet.source_count;
  return bounds;
}

[[nodiscard]] std::vector<parallel::TreePseudoParticlePacket> makeLocalTreePseudoParticleHierarchyPackets(
    int world_rank,
    std::uint64_t decomposition_epoch,
    std::uint64_t force_epoch,
    std::uint64_t exchange_sequence,
    bool periodic_unwrapped_geometry,
    const TreeNodeSoa& nodes,
  std::size_t max_packets = 128) {
  std::vector<parallel::TreePseudoParticlePacket> packets;
  if (max_packets == 0) {
    return packets;
  }
  if (nodes.size() == 0) {
    parallel::TreePseudoParticlePacket empty_packet;
    empty_packet.descriptor = parallel::TreePseudoParticleDescriptor{
        .wire_version = 1U,
        .pseudo_particle_id = (static_cast<std::uint64_t>(std::max(world_rank, 0)) << 48U) ^ decomposition_epoch,
        .source_rank = world_rank,
        .decomposition_epoch = decomposition_epoch,
        .force_epoch = force_epoch,
        .exchange_sequence = exchange_sequence,
        .derived_not_authoritative = true,
    };
    empty_packet.source_count = 0U;
    empty_packet.hierarchy_level = 0U;
    empty_packet.local_node_index = 0U;
    empty_packet.child_count = 0U;
    empty_packet.is_leaf = 1U;
    empty_packet.geometry_frame = periodic_unwrapped_geometry ? 1U : 0U;
    parallel::validateTreePseudoParticlePacket(empty_packet);
    packets.push_back(empty_packet);
    return packets;
  }
  struct QueueEntry {
    std::uint32_t node_index = 0;
    std::uint32_t level = 0;
  };
  std::vector<QueueEntry> queue;
  queue.reserve(std::min<std::size_t>(nodes.size(), max_packets));
  queue.push_back(QueueEntry{.node_index = 0U, .level = 0U});
  for (std::size_t cursor = 0; cursor < queue.size() && packets.size() < max_packets; ++cursor) {
    const QueueEntry entry = queue[cursor];
    if (entry.node_index >= nodes.size()) {
      continue;
    }
    const double half_size = nodes.half_size_comoving[entry.node_index];
    parallel::TreePseudoParticlePacket packet;
    packet.descriptor = parallel::TreePseudoParticleDescriptor{
        .wire_version = 1U,
        .pseudo_particle_id = (static_cast<std::uint64_t>(std::max(world_rank, 0)) << 48U) ^
            (static_cast<std::uint64_t>(entry.level) << 40U) ^
            static_cast<std::uint64_t>(entry.node_index) ^ decomposition_epoch,
        .source_rank = world_rank,
        .decomposition_epoch = decomposition_epoch,
        .force_epoch = force_epoch,
        .exchange_sequence = exchange_sequence,
        .derived_not_authoritative = true,
    };
    packet.mass_code = nodes.mass_code[entry.node_index];
    packet.center_x_comoving = nodes.com_x_comoving[entry.node_index];
    packet.center_y_comoving = nodes.com_y_comoving[entry.node_index];
    packet.center_z_comoving = nodes.com_z_comoving[entry.node_index];
    packet.min_x_comoving = nodes.center_x_comoving[entry.node_index] - half_size;
    packet.max_x_comoving = nodes.center_x_comoving[entry.node_index] + half_size;
    packet.min_y_comoving = nodes.center_y_comoving[entry.node_index] - half_size;
    packet.max_y_comoving = nodes.center_y_comoving[entry.node_index] + half_size;
    packet.min_z_comoving = nodes.center_z_comoving[entry.node_index] - half_size;
    packet.max_z_comoving = nodes.center_z_comoving[entry.node_index] + half_size;
    packet.source_count = nodes.particle_count[entry.node_index];
    packet.hierarchy_level = entry.level;
    packet.local_node_index = entry.node_index;
    packet.child_count = nodes.child_count[entry.node_index];
    packet.is_leaf = nodes.child_count[entry.node_index] == 0U ? 1U : 0U;
    packet.geometry_frame = periodic_unwrapped_geometry ? 1U : 0U;
    parallel::validateTreePseudoParticlePacket(packet);
    packets.push_back(packet);

    if (nodes.child_count[entry.node_index] == 0U || queue.size() >= max_packets) {
      continue;
    }
    const std::size_t child_slot_offset = static_cast<std::size_t>(entry.node_index) * 8U;
    for (std::uint8_t octant = 0; octant < 8U && queue.size() < max_packets; ++octant) {
      const std::uint32_t child = nodes.child_index[child_slot_offset + octant];
      if (child == std::numeric_limits<std::uint32_t>::max()) {
        continue;
      }
      queue.push_back(QueueEntry{.node_index = child, .level = entry.level + 1U});
    }
  }
  return packets;
}

[[nodiscard]] bool remoteTreeHierarchyIntersectsCutoff(
    double px,
    double py,
    double pz,
    std::span<const parallel::TreePseudoParticlePacket> packets,
    const PeriodicBoxLengths& box_lengths,
  double cutoff_radius_comoving) {
  for (const parallel::TreePseudoParticlePacket& packet : packets) {
    if (packet.source_count == 0U || packet.mass_code == 0.0) {
      continue;
    }
    if (minimumDistanceToPeriodicBounds(px, py, pz, boundsFromTreePseudoParticlePacket(packet), box_lengths) <=
        cutoff_radius_comoving) {
      return true;
    }
  }
  return false;
}

}  // namespace

void TreePmForceAccumulatorView::reset() const {
  std::fill(accel_x_comoving.begin(), accel_x_comoving.end(), 0.0);
  std::fill(accel_y_comoving.begin(), accel_y_comoving.end(), 0.0);
  std::fill(accel_z_comoving.begin(), accel_z_comoving.end(), 0.0);
}

void TreePmForceAccumulatorView::addToActiveSlot(
    std::size_t active_slot,
    double ax_comoving,
    double ay_comoving,
    double az_comoving) const {
  accel_x_comoving[active_slot] += ax_comoving;
  accel_y_comoving[active_slot] += ay_comoving;
  accel_z_comoving[active_slot] += az_comoving;
}

TreePmCoordinator::TreePmCoordinator(PmGridShape pm_shape)
    : m_shape(pm_shape),
      m_mpi_context(),
      m_grid(pm_shape),
      m_pm_solver(pm_shape),
      m_tree_solver() {}

TreePmCoordinator::TreePmCoordinator(PmGridShape pm_shape, parallel::PmSlabLayout pm_layout)
    : TreePmCoordinator(pm_shape, std::move(pm_layout), parallel::MpiContext()) {}

TreePmCoordinator::TreePmCoordinator(
    PmGridShape pm_shape,
    parallel::PmSlabLayout pm_layout,
    parallel::MpiContext mpi_context)
    : m_shape(pm_shape),
      m_mpi_context(std::move(mpi_context)),
      m_grid(pm_shape, std::move(pm_layout)),
      m_pm_solver(pm_shape),
      m_tree_solver() {}

const parallel::PmSlabLayout& TreePmCoordinator::slabLayout() const noexcept {
  return m_grid.slabLayout();
}

bool TreePmCoordinator::ownsFullPmDomain() const noexcept {
  return m_grid.ownsFullDomain();
}

const parallel::PmSlabHaloExchangeResult& TreePmCoordinator::lastPmSlabHaloExchange() const noexcept {
  return m_last_pm_slab_halo_exchange;
}

core::MemoryReport TreePmCoordinator::memoryReport() const {
  core::MemoryReportBuilder builder;
  m_grid.appendMemoryReport(builder);
  m_tree_solver.nodes().appendMemoryReport(builder);
  const auto add_active = [&builder](std::string label, const auto& container) {
    const std::uint64_t bytes = core::ownedCapacityBytesForContainer(container);
    builder.addEntry(core::MemoryEntry{.subsystem = core::MemorySubsystem::kActiveSets,
                                       .lifetime = core::MemoryLifetime::kTransient,
                                       .label = std::move(label),
                                       .owned_capacity_bytes = bytes,
                                       .high_water_bytes = bytes});
  };
  add_active("treepm.active_pos_x_comoving", m_active_pos_x_comoving);
  add_active("treepm.active_pos_y_comoving", m_active_pos_y_comoving);
  add_active("treepm.active_pos_z_comoving", m_active_pos_z_comoving);
  add_active("treepm.active_pm_ax_comoving", m_active_pm_ax_comoving);
  add_active("treepm.active_pm_ay_comoving", m_active_pm_ay_comoving);
  add_active("treepm.active_pm_az_comoving", m_active_pm_az_comoving);
  add_active("treepm.active_zoom_corr_ax_comoving", m_active_zoom_corr_ax_comoving);
  add_active("treepm.active_zoom_corr_ay_comoving", m_active_zoom_corr_ay_comoving);
  add_active("treepm.active_zoom_corr_az_comoving", m_active_zoom_corr_az_comoving);

  const auto add_tree_scratch = [&builder](std::string label, const auto& container) {
    const std::uint64_t bytes = core::ownedCapacityBytesForContainer(container);
    builder.addEntry(core::MemoryEntry{.subsystem = core::MemorySubsystem::kTree,
                                       .lifetime = core::MemoryLifetime::kTransient,
                                       .label = std::move(label),
                                       .owned_capacity_bytes = bytes,
                                       .high_water_bytes = bytes});
  };
  add_tree_scratch("treepm.periodic_tree_source_x_comoving", m_tree_source_x_comoving);
  add_tree_scratch("treepm.periodic_tree_source_y_comoving", m_tree_source_y_comoving);
  add_tree_scratch("treepm.periodic_tree_source_z_comoving", m_tree_source_z_comoving);

  const auto add_mpi = [&builder](std::string label, const auto& container) {
    const std::uint64_t bytes = core::ownedCapacityBytesForContainer(container);
    builder.addEntry(core::MemoryEntry{.subsystem = core::MemorySubsystem::kMpiBuffers,
                                       .lifetime = core::MemoryLifetime::kTransient,
                                       .label = std::move(label),
                                       .owned_capacity_bytes = bytes,
                                       .high_water_bytes = bytes});
  };
  add_mpi("treepm.exchange.send_counts", m_tree_exchange_workspace.send_counts);
  add_mpi("treepm.exchange.recv_counts", m_tree_exchange_workspace.recv_counts);
  add_mpi("treepm.exchange.send_displs", m_tree_exchange_workspace.send_displs);
  add_mpi("treepm.exchange.recv_displs", m_tree_exchange_workspace.recv_displs);
  add_mpi("treepm.exchange.response_send_counts", m_tree_exchange_workspace.response_send_counts);
  add_mpi("treepm.exchange.response_recv_counts", m_tree_exchange_workspace.response_recv_counts);
  add_mpi("treepm.exchange.response_send_displs", m_tree_exchange_workspace.response_send_displs);
  add_mpi("treepm.exchange.response_recv_displs", m_tree_exchange_workspace.response_recv_displs);
  add_mpi("treepm.exchange.send_payload", m_tree_exchange_workspace.send_payload);
  add_mpi("treepm.exchange.recv_payload", m_tree_exchange_workspace.recv_payload);
  add_mpi("treepm.exchange.response_send_payload", m_tree_exchange_workspace.response_send_payload);
  add_mpi("treepm.exchange.response_recv_payload", m_tree_exchange_workspace.response_recv_payload);
  add_mpi("treepm.exchange.remote_batch_ax", m_tree_exchange_workspace.remote_batch_ax);
  add_mpi("treepm.exchange.remote_batch_ay", m_tree_exchange_workspace.remote_batch_ay);
  add_mpi("treepm.exchange.remote_batch_az", m_tree_exchange_workspace.remote_batch_az);
  add_mpi("treepm.exchange.expected_response_count", m_tree_exchange_workspace.expected_response_count);
  add_mpi("treepm.exchange.received_response_count", m_tree_exchange_workspace.received_response_count);

  builder.addEntry(core::MemoryEntry{.subsystem = core::MemorySubsystem::kPmMesh,
                                     .lifetime = core::MemoryLifetime::kUnknown,
                                     .label = "pm_solver.external_fftw_or_cuda_plan_cache",
                                     .estimate_only = true,
                                     .uncertainty_note = "FFTW/cuFFT plan internals are backend-owned; cached plan count is reported elsewhere."});
  core::MemoryReport report = std::move(builder).finish();
  report.notes.push_back("TreePM report covers owned host mesh, tree-node, active-set, and exchange-buffer capacities.");
  report.notes.push_back("FFTW/GPU library internals remain unknown unless a backend exposes exact allocation hooks.");
  return report;
}

void TreePmCoordinator::solveActiveSet(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    const TreePmForceAccumulatorView& accumulator,
    const TreePmOptions& options,
    TreePmProfileEvent* profile,
    TreePmDiagnostics* diagnostics,
    const TreeSofteningView& softening_view) {
  solveActiveSetWithPmCadence(
      pos_x_comoving,
      pos_y_comoving,
      pos_z_comoving,
      mass_code,
      accumulator,
      options,
      true,
      profile,
      diagnostics,
      softening_view);
}

void TreePmCoordinator::solveActiveSetWithPmCadence(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    const TreePmForceAccumulatorView& accumulator,
    const TreePmOptions& options,
    bool refresh_long_range_field,
    TreePmProfileEvent* profile,
    TreePmDiagnostics* diagnostics,
    const TreeSofteningView& softening_view) {
  int tree_pm_entry_world_size = 1;
  bool tree_pm_rank_local_serial_mode = true;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  int tree_pm_entry_world_rank = 0;
  bool tree_pm_layout_mode_diverged = false;
  queryActiveMpiWorld(tree_pm_entry_world_size, tree_pm_entry_world_rank);
  const int local_rank_serial_mode =
      m_grid.slabLayout().world_size == 1 && m_grid.ownsFullDomain() ? 1 : 0;
  int minimum_rank_serial_mode = local_rank_serial_mode;
  int maximum_rank_serial_mode = local_rank_serial_mode;
  if (tree_pm_entry_world_size > 1) {
    MPI_Allreduce(
        &local_rank_serial_mode,
        &minimum_rank_serial_mode,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        &local_rank_serial_mode,
        &maximum_rank_serial_mode,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD);
  }
  tree_pm_layout_mode_diverged =
      minimum_rank_serial_mode != maximum_rank_serial_mode;
  tree_pm_rank_local_serial_mode = minimum_rank_serial_mode == 1;
#endif
  std::exception_ptr local_preflight_failure;
  try {
    validateTreePmPreflight(
        pos_x_comoving,
        pos_y_comoving,
        pos_z_comoving,
        mass_code,
        accumulator,
        options,
        softening_view);
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    const parallel::PmSlabLayout& layout = m_grid.slabLayout();
    if (tree_pm_layout_mode_diverged) {
      throw std::invalid_argument(
          "TreePM ranks disagree on rank-local-serial versus distributed PM layout mode");
    }
    const bool layout_matches_selected_mode = tree_pm_rank_local_serial_mode
        ? layout.world_size == 1 && layout.world_rank == 0 &&
              m_grid.ownsFullDomain()
        : layout.world_size == tree_pm_entry_world_size &&
              layout.world_rank == tree_pm_entry_world_rank;
    if (!layout_matches_selected_mode ||
        m_mpi_context.worldSize() != tree_pm_entry_world_size ||
        m_mpi_context.worldRank() != tree_pm_entry_world_rank ||
        (tree_pm_entry_world_size > 1 && !m_mpi_context.isEnabled())) {
      throw std::invalid_argument(
          "TreePM PM layout and MPI context world metadata must match MPI_COMM_WORLD");
    }
#else
    if (m_grid.slabLayout().world_size != 1 ||
        m_grid.slabLayout().world_rank != 0) {
      throw std::invalid_argument(
          "TreePM distributed PM layout requires COSMOSIM_ENABLE_MPI=ON");
    }
#endif
  } catch (...) {
    local_preflight_failure = std::current_exception();
  }
  std::uint64_t preflight_failure_count = local_preflight_failure ? 1U : 0U;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (tree_pm_entry_world_size > 1) {
    std::uint64_t global_preflight_failure_count = 0U;
    MPI_Allreduce(
        &preflight_failure_count,
        &global_preflight_failure_count,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    preflight_failure_count = global_preflight_failure_count;
  }
#endif
  if (preflight_failure_count != 0U) {
    if (local_preflight_failure) {
      std::rethrow_exception(local_preflight_failure);
    }
    throw std::runtime_error(
        "TreePM peer rank rejected its local input during coordinated preflight");
  }
  const auto coordinate_tree_pm_failure = [&](std::exception_ptr local_failure,
                                               std::string_view phase) {
    std::uint64_t failure_count = local_failure ? 1U : 0U;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    if (tree_pm_entry_world_size > 1) {
      std::uint64_t global_failure_count = 0U;
      MPI_Allreduce(
          &failure_count,
          &global_failure_count,
          1,
          MPI_UINT64_T,
          MPI_SUM,
          MPI_COMM_WORLD);
      failure_count = global_failure_count;
    }
#endif
    if (failure_count == 0U) {
      return;
    }
    if (local_failure) {
      std::rethrow_exception(local_failure);
    }
    throw std::runtime_error(
        "TreePM peer rank rejected " + std::string(phase));
  };
  const auto start = std::chrono::steady_clock::now();
  accumulator.reset();

  // PM owns long-range force via explicit Gaussian Fourier filter. Cadence-aware callers
  // may choose to reuse the previously solved PM mesh field.
  PmSolveOptions pm_options = options.pm_options;
  pm_options.tree_pm_split_scale_comoving = options.split_policy.split_scale_comoving;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (tree_pm_entry_world_size > 1) {
    const std::uint64_t local_fingerprint =
        treePmOptionFingerprint(m_shape, options, pm_options, softening_view);
    std::uint64_t minimum_fingerprint = 0U;
    std::uint64_t maximum_fingerprint = 0U;
    MPI_Allreduce(
        &local_fingerprint,
        &minimum_fingerprint,
        1,
        MPI_UINT64_T,
        MPI_MIN,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        &local_fingerprint,
        &maximum_fingerprint,
        1,
        MPI_UINT64_T,
        MPI_MAX,
        MPI_COMM_WORLD);
    if (minimum_fingerprint != maximum_fingerprint) {
      throw std::runtime_error(
          "TreePM mesh shape, physical options, or protocol epochs diverged across ranks before PM collectives");
    }
  }
#endif
  const double comoving_gravity_prefactor =
      pm_options.gravitational_constant_code;
  if (!std::isfinite(comoving_gravity_prefactor) ||
      comoving_gravity_prefactor <= 0.0) {
    throw std::overflow_error(
        "TreePM code-unit gravitational prefactor is non-finite or non-positive");
  }

  // TreePM returns the scale-free comoving Newtonian kernel A. The KDK
  // integrator applies its cosmological 1/a^2 kick factor to physical peculiar
  // velocity, so neither PM nor the complementary residual may insert another
  // scale-factor power here. Keep the caller-owned standalone-tree options
  // unchanged and derive the effective TreePM short-range options locally.
  TreePmOptions short_range_options = options;
  short_range_options.pm_options = pm_options;
  short_range_options.tree_options.gravitational_constant_code =
      comoving_gravity_prefactor;
  const PeriodicBoxLengths cache_box_lengths = effectivePeriodicBoxLengths(pm_options);
  // Particle ownership changes advance the tree protocol epoch, but do not
  // invalidate a mesh field whose FFT slab ownership, physical force epoch,
  // scale factor, and mesh options are unchanged. Interpolation is routed to
  // the fixed slab owners for the caller's current target decomposition.
  const auto cache_matches_request = [&]() {
    return m_long_range_field_validity.valid &&
        m_long_range_field_validity.force_epoch == options.force_epoch &&
        m_long_range_field_validity.scale_factor == pm_options.scale_factor &&
        m_long_range_field_validity.gravitational_constant_code == pm_options.gravitational_constant_code &&
        m_long_range_field_validity.split_scale_comoving == pm_options.tree_pm_split_scale_comoving &&
        m_long_range_field_validity.box_size_x_comoving == cache_box_lengths.lx &&
        m_long_range_field_validity.box_size_y_comoving == cache_box_lengths.ly &&
        m_long_range_field_validity.box_size_z_comoving == cache_box_lengths.lz &&
        m_long_range_field_validity.assignment_scheme == pm_options.assignment_scheme &&
        m_long_range_field_validity.boundary_condition == pm_options.boundary_condition &&
        m_long_range_field_validity.decomposition_mode == pm_options.decomposition_mode &&
        m_long_range_field_validity.window_deconvolution == pm_options.enable_window_deconvolution;
  };

  bool perform_long_range_refresh = refresh_long_range_field;
  if (tree_pm_entry_world_size > 1) {
    const std::uint64_t requested_refresh_votes =
        m_mpi_context.allreduceSumUint64(refresh_long_range_field ? 1U : 0U);
    const std::uint64_t world_size =
        static_cast<std::uint64_t>(tree_pm_entry_world_size);
    if (requested_refresh_votes != 0U && requested_refresh_votes != world_size) {
      throw std::runtime_error(
          "TreePM long-range refresh/reuse request diverged across ranks before PM collectives");
    }
    const std::uint64_t compatible_cache_votes =
        m_mpi_context.allreduceSumUint64(cache_matches_request() ? 1U : 0U);
    if (requested_refresh_votes == 0U && compatible_cache_votes != world_size) {
      throw std::runtime_error(
          "TreePM long-range reuse requested without a compatible PM field on every rank");
    }
    perform_long_range_refresh = requested_refresh_votes == world_size;
  } else if (!refresh_long_range_field && !cache_matches_request()) {
    throw std::runtime_error(
        "TreePM long-range reuse requested without a compatible PM field");
  }

  if (perform_long_range_refresh) {
    m_pm_solver.assignDensity(
        m_grid,
        pos_x_comoving,
        pos_y_comoving,
        pos_z_comoving,
        mass_code,
        pm_options,
        profile != nullptr ? &profile->pm_profile : nullptr);
    m_last_pm_slab_halo_exchange = {};
    m_grid.clearForceHaloCache();
    if (pm_options.boundary_condition == PmBoundaryCondition::kPeriodic) {
      m_pm_solver.solvePoissonPeriodic(m_grid, pm_options, profile != nullptr ? &profile->pm_profile : nullptr);
    } else {
      m_pm_solver.solvePoissonIsolatedOpen(m_grid, pm_options, profile != nullptr ? &profile->pm_profile : nullptr);
    }
    if (!m_grid.ownsFullDomain() && m_grid.slabLayout().world_size > 1) {
      const std::uint64_t exchange_sequence = ++m_pm_halo_exchange_sequence;
      const parallel::PmSlabHaloExchangeResult force_x_halo = parallel::executeBlockingPmSlabHaloExchange(
          m_mpi_context,
          m_grid.slabLayout(),
          m_grid.force_x(),
          /*halo_depth_x=*/1,
          pm_options.boundary_condition == PmBoundaryCondition::kPeriodic,
          exchange_sequence * 3U + 0U);
      const parallel::PmSlabHaloExchangeResult force_y_halo = parallel::executeBlockingPmSlabHaloExchange(
          m_mpi_context,
          m_grid.slabLayout(),
          m_grid.force_y(),
          /*halo_depth_x=*/1,
          pm_options.boundary_condition == PmBoundaryCondition::kPeriodic,
          exchange_sequence * 3U + 1U);
      const parallel::PmSlabHaloExchangeResult force_z_halo = parallel::executeBlockingPmSlabHaloExchange(
          m_mpi_context,
          m_grid.slabLayout(),
          m_grid.force_z(),
          /*halo_depth_x=*/1,
          pm_options.boundary_condition == PmBoundaryCondition::kPeriodic,
          exchange_sequence * 3U + 2U);
      std::exception_ptr halo_cache_commit_failure;
      try {
        m_grid.setForceHaloCache(force_x_halo, force_y_halo, force_z_halo, exchange_sequence);
        m_last_pm_slab_halo_exchange = force_x_halo;
        m_last_pm_slab_halo_exchange.sent_bytes += force_y_halo.sent_bytes + force_z_halo.sent_bytes;
        m_last_pm_slab_halo_exchange.received_bytes += force_y_halo.received_bytes + force_z_halo.received_bytes;
        if (profile != nullptr) {
          profile->pm_profile.bytes_moved +=
              m_last_pm_slab_halo_exchange.sent_bytes + m_last_pm_slab_halo_exchange.received_bytes;
        }
      } catch (...) {
        halo_cache_commit_failure = std::current_exception();
      }
      coordinate_tree_pm_failure(
          halo_cache_commit_failure, "PM halo-cache commit");
    }
    m_long_range_field_validity = LongRangeFieldValidity{
        .valid = true,
        .decomposition_epoch = options.decomposition_epoch,
        .force_epoch = options.force_epoch,
        .scale_factor = pm_options.scale_factor,
        .gravitational_constant_code = pm_options.gravitational_constant_code,
        .split_scale_comoving = pm_options.tree_pm_split_scale_comoving,
        .box_size_x_comoving = cache_box_lengths.lx,
        .box_size_y_comoving = cache_box_lengths.ly,
        .box_size_z_comoving = cache_box_lengths.lz,
        .assignment_scheme = pm_options.assignment_scheme,
        .boundary_condition = pm_options.boundary_condition,
        .decomposition_mode = pm_options.decomposition_mode,
        .window_deconvolution = pm_options.enable_window_deconvolution,
    };
  }
  if (!m_long_range_field_validity.valid) {
    throw std::runtime_error("TreePM long-range mesh field is unavailable for reuse");
  }

  const std::size_t active_count = accumulator.active_particle_index.size();
  std::uint64_t zoom_high_res_allgather_bytes = 0;
  const bool has_explicit_target_positions = !accumulator.target_pos_x_comoving.empty();
  std::exception_ptr compact_target_preparation_failure;
  try {
    resizeCompactSidecars(m_active_pos_x_comoving, m_active_pos_y_comoving, m_active_pos_z_comoving, active_count);
    resizeCompactSidecars(m_active_pm_ax_comoving, m_active_pm_ay_comoving, m_active_pm_az_comoving, active_count);
    resizeCompactSidecars(
        m_active_zoom_corr_ax_comoving,
        m_active_zoom_corr_ay_comoving,
        m_active_zoom_corr_az_comoving,
        active_count);
    std::fill(m_active_zoom_corr_ax_comoving.begin(), m_active_zoom_corr_ax_comoving.end(), 0.0);
    std::fill(m_active_zoom_corr_ay_comoving.begin(), m_active_zoom_corr_ay_comoving.end(), 0.0);
    std::fill(m_active_zoom_corr_az_comoving.begin(), m_active_zoom_corr_az_comoving.end(), 0.0);

    for (std::size_t i = 0; i < accumulator.active_particle_index.size(); ++i) {
      const std::uint32_t particle_index = accumulator.active_particle_index[i];
      if (!has_explicit_target_positions && particle_index >= pos_x_comoving.size()) {
        throw std::out_of_range("TreePM active index exceeds particle count");
      }
      m_active_pos_x_comoving[i] = has_explicit_target_positions
          ? accumulator.target_pos_x_comoving[i]
          : pos_x_comoving[particle_index];
      m_active_pos_y_comoving[i] = has_explicit_target_positions
          ? accumulator.target_pos_y_comoving[i]
          : pos_y_comoving[particle_index];
      m_active_pos_z_comoving[i] = has_explicit_target_positions
          ? accumulator.target_pos_z_comoving[i]
          : pos_z_comoving[particle_index];
      if (!std::isfinite(m_active_pos_x_comoving[i]) || !std::isfinite(m_active_pos_y_comoving[i]) ||
          !std::isfinite(m_active_pos_z_comoving[i])) {
        throw std::invalid_argument("TreePM target positions must be finite");
      }
    }
  } catch (...) {
    compact_target_preparation_failure = std::current_exception();
  }
  coordinate_tree_pm_failure(
      compact_target_preparation_failure, "compact PM target preparation");

  m_pm_solver.interpolateForces(
      m_grid,
      m_active_pos_x_comoving,
      m_active_pos_y_comoving,
      m_active_pos_z_comoving,
      m_active_pm_ax_comoving,
      m_active_pm_ay_comoving,
      m_active_pm_az_comoving,
      pm_options,
      profile != nullptr ? &profile->pm_profile : nullptr);

  for (std::size_t i = 0; i < accumulator.active_particle_index.size(); ++i) {
    accumulator.addToActiveSlot(i, m_active_pm_ax_comoving[i], m_active_pm_ay_comoving[i], m_active_pm_az_comoving[i]);
  }

  if (options.enable_zoom_long_range_correction) {
    if (!options.zoom_focused_pm_shape.isValid()) {
      throw std::invalid_argument("zoom_focused_pm_shape must be valid when zoom correction is enabled");
    }
    if (options.source_is_high_res.size() != pos_x_comoving.size()) {
      throw std::invalid_argument("zoom source high-res mask size must match source particle count");
    }
    if (options.active_is_high_res.size() != accumulator.active_particle_index.size()) {
      throw std::invalid_argument("zoom active high-res mask size must match active set");
    }
    std::vector<double> high_res_source_x;
    std::vector<double> high_res_source_y;
    std::vector<double> high_res_source_z;
    std::vector<double> high_res_source_mass;
    high_res_source_x.reserve(pos_x_comoving.size());
    high_res_source_y.reserve(pos_y_comoving.size());
    high_res_source_z.reserve(pos_z_comoving.size());
    high_res_source_mass.reserve(mass_code.size());
    for (std::size_t i = 0; i < pos_x_comoving.size(); ++i) {
      if (options.source_is_high_res[i] == 0U) {
        continue;
      }
      high_res_source_x.push_back(pos_x_comoving[i]);
      high_res_source_y.push_back(pos_y_comoving[i]);
      high_res_source_z.push_back(pos_z_comoving[i]);
      high_res_source_mass.push_back(mass_code[i]);
    }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    if (m_grid.slabLayout().world_size > 1) {
      const auto gathered = allGatherParticleField(
          high_res_source_x,
          high_res_source_y,
          high_res_source_z,
          high_res_source_mass,
          options.zoom_high_res_allgather_limit_bytes,
          &zoom_high_res_allgather_bytes);
      high_res_source_x = std::move(gathered.x);
      high_res_source_y = std::move(gathered.y);
      high_res_source_z = std::move(gathered.z);
      high_res_source_mass = std::move(gathered.m);
      if (profile != nullptr) {
        profile->pm_profile.bytes_moved += zoom_high_res_allgather_bytes;
      }
    }
#endif

    if (!high_res_source_x.empty()) {
      PmGridStorage high_res_coarse_grid(m_shape);
      m_pm_solver.assignDensity(
          high_res_coarse_grid,
          high_res_source_x,
          high_res_source_y,
          high_res_source_z,
          high_res_source_mass,
          pm_options,
          profile != nullptr ? &profile->pm_profile : nullptr);
      if (pm_options.boundary_condition == PmBoundaryCondition::kPeriodic) {
        m_pm_solver.solvePoissonPeriodic(
            high_res_coarse_grid, pm_options, profile != nullptr ? &profile->pm_profile : nullptr);
      } else {
        m_pm_solver.solvePoissonIsolatedOpen(
            high_res_coarse_grid, pm_options, profile != nullptr ? &profile->pm_profile : nullptr);
      }

      const PeriodicBoxLengths box_lengths = effectivePeriodicBoxLengths(pm_options);
      const double focused_half_extent = std::max({
          options.zoom_region_radius_comoving,
          options.zoom_contamination_radius_comoving,
          options.split_policy.cutoff_radius_comoving});
      PmSolveOptions zoom_pm_options = pm_options;
      zoom_pm_options.boundary_condition = PmBoundaryCondition::kIsolatedOpen;
      // The focused correction is an isolated/open convolution. Assignment-
      // window deconvolution is currently certified only for the periodic
      // spectral operator, so this branch must select its documented
      // non-deconvolved isolated policy explicitly rather than inheriting the
      // global periodic PM setting.
      zoom_pm_options.enable_window_deconvolution = false;
      zoom_pm_options.box_size_x_mpc_comoving = 2.0 * focused_half_extent;
      zoom_pm_options.box_size_y_mpc_comoving = 2.0 * focused_half_extent;
      zoom_pm_options.box_size_z_mpc_comoving = 2.0 * focused_half_extent;
      zoom_pm_options.box_size_mpc_comoving = 2.0 * focused_half_extent;
      PmSolver zoom_solver(options.zoom_focused_pm_shape);
      PmGridStorage high_res_focused_grid(options.zoom_focused_pm_shape);
      std::vector<double> focused_source_x;
      std::vector<double> focused_source_y;
      std::vector<double> focused_source_z;
      std::vector<double> focused_source_mass;
      focused_source_x.reserve(high_res_source_x.size());
      focused_source_y.reserve(high_res_source_y.size());
      focused_source_z.reserve(high_res_source_z.size());
      focused_source_mass.reserve(high_res_source_mass.size());
      for (std::size_t i = 0; i < high_res_source_x.size(); ++i) {
        const double dx_local = (pm_options.boundary_condition == PmBoundaryCondition::kPeriodic)
            ? minimumImageDelta(high_res_source_x[i] - options.zoom_region_center_x_comoving, box_lengths.lx)
            : (high_res_source_x[i] - options.zoom_region_center_x_comoving);
        const double dy_local = (pm_options.boundary_condition == PmBoundaryCondition::kPeriodic)
            ? minimumImageDelta(high_res_source_y[i] - options.zoom_region_center_y_comoving, box_lengths.ly)
            : (high_res_source_y[i] - options.zoom_region_center_y_comoving);
        const double dz_local = (pm_options.boundary_condition == PmBoundaryCondition::kPeriodic)
            ? minimumImageDelta(high_res_source_z[i] - options.zoom_region_center_z_comoving, box_lengths.lz)
            : (high_res_source_z[i] - options.zoom_region_center_z_comoving);
        if (std::abs(dx_local) > focused_half_extent ||
            std::abs(dy_local) > focused_half_extent ||
            std::abs(dz_local) > focused_half_extent) {
          continue;
        }
        focused_source_x.push_back(dx_local + focused_half_extent);
        focused_source_y.push_back(dy_local + focused_half_extent);
        focused_source_z.push_back(dz_local + focused_half_extent);
        focused_source_mass.push_back(high_res_source_mass[i]);
      }
      zoom_solver.assignDensity(
          high_res_focused_grid,
          focused_source_x,
          focused_source_y,
          focused_source_z,
          focused_source_mass,
          zoom_pm_options,
          profile != nullptr ? &profile->pm_profile : nullptr);
      zoom_solver.solvePoissonIsolatedOpen(
          high_res_focused_grid, zoom_pm_options, profile != nullptr ? &profile->pm_profile : nullptr);

      std::vector<double> high_res_pm_coarse_ax(active_count, 0.0);
      std::vector<double> high_res_pm_coarse_ay(active_count, 0.0);
      std::vector<double> high_res_pm_coarse_az(active_count, 0.0);
      std::vector<double> high_res_pm_focused_ax(active_count, 0.0);
      std::vector<double> high_res_pm_focused_ay(active_count, 0.0);
      std::vector<double> high_res_pm_focused_az(active_count, 0.0);
      m_pm_solver.interpolateForces(
          high_res_coarse_grid,
          m_active_pos_x_comoving,
          m_active_pos_y_comoving,
          m_active_pos_z_comoving,
          high_res_pm_coarse_ax,
          high_res_pm_coarse_ay,
          high_res_pm_coarse_az,
          pm_options,
          profile != nullptr ? &profile->pm_profile : nullptr);
      std::vector<double> focused_active_x(active_count, focused_half_extent);
      std::vector<double> focused_active_y(active_count, focused_half_extent);
      std::vector<double> focused_active_z(active_count, focused_half_extent);
      for (std::size_t i = 0; i < active_count; ++i) {
        const double dx_local = (pm_options.boundary_condition == PmBoundaryCondition::kPeriodic)
            ? minimumImageDelta(m_active_pos_x_comoving[i] - options.zoom_region_center_x_comoving, box_lengths.lx)
            : (m_active_pos_x_comoving[i] - options.zoom_region_center_x_comoving);
        const double dy_local = (pm_options.boundary_condition == PmBoundaryCondition::kPeriodic)
            ? minimumImageDelta(m_active_pos_y_comoving[i] - options.zoom_region_center_y_comoving, box_lengths.ly)
            : (m_active_pos_y_comoving[i] - options.zoom_region_center_y_comoving);
        const double dz_local = (pm_options.boundary_condition == PmBoundaryCondition::kPeriodic)
            ? minimumImageDelta(m_active_pos_z_comoving[i] - options.zoom_region_center_z_comoving, box_lengths.lz)
            : (m_active_pos_z_comoving[i] - options.zoom_region_center_z_comoving);
        focused_active_x[i] = dx_local + focused_half_extent;
        focused_active_y[i] = dy_local + focused_half_extent;
        focused_active_z[i] = dz_local + focused_half_extent;
      }
      zoom_solver.interpolateForces(
          high_res_focused_grid,
          focused_active_x,
          focused_active_y,
          focused_active_z,
          high_res_pm_focused_ax,
          high_res_pm_focused_ay,
          high_res_pm_focused_az,
          zoom_pm_options,
          profile != nullptr ? &profile->pm_profile : nullptr);

      for (std::size_t i = 0; i < active_count; ++i) {
        if (options.active_is_high_res[i] == 0U) {
          continue;
        }
        const double corr_x = high_res_pm_focused_ax[i] - high_res_pm_coarse_ax[i];
        const double corr_y = high_res_pm_focused_ay[i] - high_res_pm_coarse_ay[i];
        const double corr_z = high_res_pm_focused_az[i] - high_res_pm_coarse_az[i];
        m_active_zoom_corr_ax_comoving[i] = corr_x;
        m_active_zoom_corr_ay_comoving[i] = corr_y;
        m_active_zoom_corr_az_comoving[i] = corr_z;
        accumulator.addToActiveSlot(i, corr_x, corr_y, corr_z);
      }
    }
  }

  // Tree owns short-range residual with the complementary real-space kernel.
  const auto tree_start = std::chrono::steady_clock::now();
  std::span<const double> tree_source_x = pos_x_comoving;
  std::span<const double> tree_source_y = pos_y_comoving;
  std::span<const double> tree_source_z = pos_z_comoving;
  std::exception_ptr local_tree_build_failure;
  try {
    if (options.pm_options.boundary_condition == PmBoundaryCondition::kPeriodic) {
      const PeriodicBoxLengths box_lengths = effectivePeriodicBoxLengths(options.pm_options);
      unwrapPeriodicAxis(pos_x_comoving, box_lengths.lx, m_tree_source_x_comoving);
      unwrapPeriodicAxis(pos_y_comoving, box_lengths.ly, m_tree_source_y_comoving);
      unwrapPeriodicAxis(pos_z_comoving, box_lengths.lz, m_tree_source_z_comoving);
      tree_source_x = m_tree_source_x_comoving;
      tree_source_y = m_tree_source_y_comoving;
      tree_source_z = m_tree_source_z_comoving;
    }
    m_tree_solver.build(
        tree_source_x,
        tree_source_y,
        tree_source_z,
        mass_code,
        short_range_options.tree_options,
        profile != nullptr ? &profile->tree_profile : nullptr,
        softening_view);
  } catch (...) {
    local_tree_build_failure = std::current_exception();
  }
  std::uint64_t tree_build_failure_count = local_tree_build_failure ? 1U : 0U;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (tree_pm_entry_world_size > 1) {
    std::uint64_t global_tree_build_failure_count = 0U;
    MPI_Allreduce(
        &tree_build_failure_count,
        &global_tree_build_failure_count,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    tree_build_failure_count = global_tree_build_failure_count;
  }
#endif
  if (tree_build_failure_count != 0U) {
    if (local_tree_build_failure) {
      std::rethrow_exception(local_tree_build_failure);
    }
    throw std::runtime_error(
        "TreePM peer rank rejected periodic unwrapping or local tree build");
  }
  evaluateShortRangeResidual(
      pos_x_comoving,
      pos_y_comoving,
      pos_z_comoving,
      mass_code,
      accumulator,
      short_range_options,
      softening_view,
      tree_pm_rank_local_serial_mode,
      profile != nullptr ? &profile->tree_profile : nullptr);
  const auto tree_stop = std::chrono::steady_clock::now();

  std::uint64_t empty_source_rank_count = pos_x_comoving.empty() ? 1U : 0U;
  std::uint64_t empty_target_rank_count = active_count == 0U ? 1U : 0U;
  if (m_grid.slabLayout().world_size > 1) {
    empty_source_rank_count = m_mpi_context.allreduceSumUint64(empty_source_rank_count);
    empty_target_rank_count = m_mpi_context.allreduceSumUint64(empty_target_rank_count);
  }

  if (diagnostics != nullptr) {
    const PeriodicBoxLengths diagnostic_box_lengths = effectivePeriodicBoxLengths(options.pm_options);
    *diagnostics = computeTreePmDiagnostics(options.split_policy);
    diagnostics->local_source_count = static_cast<std::uint64_t>(pos_x_comoving.size());
    diagnostics->local_active_target_count = static_cast<std::uint64_t>(active_count);
    diagnostics->local_tree_node_count = static_cast<std::uint64_t>(m_tree_solver.nodes().size());
    diagnostics->pm_solve_count = perform_long_range_refresh ? 1U : 0U;
    diagnostics->pm_reuse_count = perform_long_range_refresh ? 0U : 1U;
    diagnostics->pm_halo_value_count = static_cast<std::uint64_t>(
        m_last_pm_slab_halo_exchange.left_halo.size() + m_last_pm_slab_halo_exchange.right_halo.size());
    diagnostics->pm_local_nx = static_cast<std::uint64_t>(m_grid.slabLayout().local_nx());
    diagnostics->pm_local_ny = static_cast<std::uint64_t>(m_shape.ny);
    diagnostics->pm_local_nz = static_cast<std::uint64_t>(m_shape.nz);
    diagnostics->empty_source_rank_count = empty_source_rank_count;
    diagnostics->empty_target_rank_count = empty_target_rank_count;
    diagnostics->zoom_high_res_allgather_bytes = zoom_high_res_allgather_bytes;
    diagnostics->zoom_high_res_allgather_limit_bytes = options.zoom_high_res_allgather_limit_bytes;
    diagnostics->residual_pruned_nodes = m_last_residual_stats.pruned_nodes;
    diagnostics->residual_pair_skips_cutoff = m_last_residual_stats.pair_skips_cutoff;
    diagnostics->residual_pair_evaluations = m_last_residual_stats.pair_evaluations;
    diagnostics->residual_remote_request_packets = m_last_residual_stats.remote_request_packets;
    diagnostics->residual_remote_response_packets = m_last_residual_stats.remote_response_packets;
    diagnostics->residual_remote_request_bytes = m_last_residual_stats.remote_request_bytes;
    diagnostics->residual_remote_response_bytes = m_last_residual_stats.remote_response_bytes;
    diagnostics->residual_remote_request_batches = m_last_residual_stats.remote_request_batches;
    diagnostics->residual_remote_peer_participations = m_last_residual_stats.remote_peer_participations;
    diagnostics->residual_remote_targets_with_requests = m_last_residual_stats.remote_targets_with_requests;
    diagnostics->residual_remote_targets_without_requests = m_last_residual_stats.remote_targets_without_requests;
    diagnostics->residual_remote_pairs_pruned_by_bounds = m_last_residual_stats.remote_pairs_pruned_by_bounds;
    diagnostics->residual_remote_request_packets_max_peer = m_last_residual_stats.remote_request_packets_max_peer;
    diagnostics->residual_remote_response_packets_max_peer = m_last_residual_stats.remote_response_packets_max_peer;
    diagnostics->residual_remote_request_packet_imbalance_ratio =
        m_last_residual_stats.remote_request_packet_imbalance_ratio;
    diagnostics->remote_hierarchy_packet_count = m_last_residual_stats.remote_hierarchy_packets;
    diagnostics->communicating_peer_count = m_last_residual_stats.communicating_peer_count;
    diagnostics->force_l2_pm_global = l2NormFromComponents(
        m_active_pm_ax_comoving, m_active_pm_ay_comoving, m_active_pm_az_comoving);
    diagnostics->force_l2_pm_zoom_correction = l2NormFromComponents(
        m_active_zoom_corr_ax_comoving,
        m_active_zoom_corr_ay_comoving,
        m_active_zoom_corr_az_comoving);
    double tree_sq = 0.0;
    double total_sq = 0.0;
    for (std::size_t i = 0; i < active_count; ++i) {
      const double total_x = accumulator.accel_x_comoving[i];
      const double total_y = accumulator.accel_y_comoving[i];
      const double total_z = accumulator.accel_z_comoving[i];
      const double tree_x =
          total_x - m_active_pm_ax_comoving[i] - m_active_zoom_corr_ax_comoving[i];
      const double tree_y =
          total_y - m_active_pm_ay_comoving[i] - m_active_zoom_corr_ay_comoving[i];
      const double tree_z =
          total_z - m_active_pm_az_comoving[i] - m_active_zoom_corr_az_comoving[i];
      tree_sq += tree_x * tree_x + tree_y * tree_y + tree_z * tree_z;
      total_sq += total_x * total_x + total_y * total_y + total_z * total_z;
    }
    diagnostics->force_l2_tree_short_range = std::sqrt(tree_sq);
    diagnostics->force_l2_tree_short_range_local = std::sqrt(m_last_residual_stats.local_short_range_sum_sq);
    diagnostics->force_l2_tree_short_range_remote = std::sqrt(m_last_residual_stats.remote_short_range_sum_sq);
    diagnostics->force_l2_total = std::sqrt(total_sq);
    if (m_tree_solver.nodes().size() > 0U) {
      diagnostics->tree_root_half_size_comoving = m_tree_solver.nodes().half_size_comoving.front();
      diagnostics->tree_root_com_x_comoving = m_tree_solver.nodes().com_x_comoving.front();
      diagnostics->tree_root_com_y_comoving = m_tree_solver.nodes().com_y_comoving.front();
      diagnostics->tree_root_com_z_comoving = m_tree_solver.nodes().com_z_comoving.front();
    }
    for (std::size_t i = 0; i < options.source_is_high_res.size(); ++i) {
      if (options.source_is_high_res[i] != 0U) {
        ++diagnostics->zoom_high_res_source_count;
      } else {
        ++diagnostics->zoom_low_res_source_count;
        const double dx = minimumImageDelta(
            pos_x_comoving[i] - options.zoom_region_center_x_comoving,
            diagnostic_box_lengths.lx);
        const double dy = minimumImageDelta(
            pos_y_comoving[i] - options.zoom_region_center_y_comoving,
            diagnostic_box_lengths.ly);
        const double dz = minimumImageDelta(
            pos_z_comoving[i] - options.zoom_region_center_z_comoving,
            diagnostic_box_lengths.lz);
        const double r = norm3(dx, dy, dz);
        if (r <= options.zoom_contamination_radius_comoving) {
          ++diagnostics->zoom_low_res_contamination_count;
          diagnostics->zoom_low_res_contamination_mass_code += mass_code[i];
        }
      }
    }
  }

  if (profile != nullptr) {
    profile->tree_short_range_ms += std::chrono::duration<double, std::milli>(tree_stop - tree_start).count();
    profile->coupling_overhead_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
  }
}

void TreePmCoordinator::evaluateShortRangeResidual(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    const TreePmForceAccumulatorView& accumulator,
    const TreePmOptions& options,
    const TreeSofteningView& softening_view,
    bool rank_local_serial_mode,
    TreeGravityProfile* tree_profile) {
  std::uint64_t visited_nodes = 0;
  std::uint64_t accepted_nodes = 0;
  std::uint64_t opened_nodes = 0;
  std::uint64_t pp_interactions = 0;
  std::uint64_t cutoff_pruned_nodes = 0;
  std::uint64_t cutoff_pair_skips = 0;
  m_last_residual_stats = {};

  const auto traversal_start = std::chrono::steady_clock::now();
  const TreeNodeSoa& nodes = m_tree_solver.nodes();
  const TreeMortonOrdering& ordering = m_tree_solver.ordering();
  const PeriodicBoxLengths box_lengths = options.pm_options.boundary_condition == PmBoundaryCondition::kPeriodic
      ? effectivePeriodicBoxLengths(options.pm_options)
      : PeriodicBoxLengths{};
  const double cutoff_radius_comoving = options.split_policy.cutoff_radius_comoving;
  const double cutoff_radius2_comoving = cutoff_radius_comoving * cutoff_radius_comoving;
  if (!softening_view.source_particle_epsilon_comoving.empty() &&
      softening_view.source_particle_epsilon_comoving.size() != pos_x_comoving.size()) {
    throw std::invalid_argument("TreePM source softening sidecar size must match source particle count");
  }
  if (!softening_view.source_particle_epsilon_override_mask.empty() &&
      softening_view.source_particle_epsilon_override_mask.size() != pos_x_comoving.size()) {
    throw std::invalid_argument("TreePM source softening override mask size must match source particle count");
  }
  if (!softening_view.source_species_tag.empty() && softening_view.source_species_tag.size() != pos_x_comoving.size()) {
    throw std::invalid_argument("TreePM source species sidecar size must match source particle count");
  }
  if (!softening_view.target_particle_epsilon_comoving.empty() &&
      softening_view.target_particle_epsilon_comoving.size() != accumulator.active_particle_index.size()) {
    throw std::invalid_argument("TreePM target softening sidecar size must match active-set size");
  }
  if (!softening_view.target_particle_epsilon_override_mask.empty() &&
      softening_view.target_particle_epsilon_override_mask.size() != accumulator.active_particle_index.size()) {
    throw std::invalid_argument("TreePM target softening override mask size must match active-set size");
  }
  if (!softening_view.target_species_tag.empty() &&
      softening_view.target_species_tag.size() != accumulator.active_particle_index.size()) {
    throw std::invalid_argument("TreePM target species sidecar size must match active-set size");
  }
  const auto resolve_target_softening = [&](std::size_t active_slot, std::uint32_t source_index) {
    const bool has_local_source_identity = source_index < pos_x_comoving.size();
    if (!has_local_source_identity &&
        ((softening_view.target_particle_epsilon_comoving.empty() &&
          !softening_view.source_particle_epsilon_comoving.empty()) ||
         (softening_view.target_species_tag.empty() && !softening_view.source_species_tag.empty()))) {
      throw std::invalid_argument(
          "TreePM independent targets require target-owned softening/species sidecars when source sidecars are present");
    }
    return resolveTargetSofteningEpsilon(
        active_slot,
        has_local_source_identity ? source_index : 0U,
        options.tree_options.softening,
        softening_view);
  };
  int tree_mpi_world_size = 1;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  int tree_mpi_world_rank = 0;
  queryActiveMpiWorld(tree_mpi_world_size, tree_mpi_world_rank);
#endif
  std::vector<std::uint32_t> stack;
  std::exception_ptr traversal_workspace_failure;
  try {
    stack.reserve(nodes.size());
  } catch (...) {
    traversal_workspace_failure = std::current_exception();
  }
  std::uint64_t traversal_workspace_failure_count =
      traversal_workspace_failure ? 1U : 0U;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (tree_mpi_world_size > 1) {
    std::uint64_t global_traversal_workspace_failure_count = 0U;
    MPI_Allreduce(
        &traversal_workspace_failure_count,
        &global_traversal_workspace_failure_count,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    traversal_workspace_failure_count =
        global_traversal_workspace_failure_count;
  }
#endif
  if (traversal_workspace_failure_count != 0U) {
    if (traversal_workspace_failure) {
      std::rethrow_exception(traversal_workspace_failure);
    }
    throw std::runtime_error(
        "TreePM peer rank rejected residual traversal workspace preparation");
  }

  auto evaluateTargetAgainstLocalTree = [&](
                                         double px,
                                         double py,
                                         double pz,
                                         std::uint32_t self_index,
                                         bool skip_self,
                                         double target_softening_comoving,
                                         bool previous_acceleration_available,
                                         double previous_acceleration_magnitude_code) {
    double ax = 0.0;
    if (nodes.size() == 0) {
      return std::array<double, 3>{0.0, 0.0, 0.0};
    }
    double ay = 0.0;
    double az = 0.0;
    stack.clear();
    stack.push_back(0U);

    while (!stack.empty()) {
      const std::uint32_t node_index = stack.back();
      stack.pop_back();
      ++visited_nodes;

      const double half_size = nodes.half_size_comoving[node_index];
      const double min_node_distance = minimumDistanceToNodeAabb(
          px,
          py,
          pz,
          nodes.center_x_comoving[node_index],
          nodes.center_y_comoving[node_index],
          nodes.center_z_comoving[node_index],
          half_size,
          box_lengths);
      if (min_node_distance > cutoff_radius_comoving) {
        ++cutoff_pruned_nodes;
        if (!skip_self) {
          ++m_last_residual_stats.remote_pairs_pruned_by_bounds;
        }
        continue;
      }

      const double dx = minimumImageDelta(nodes.com_x_comoving[node_index] - px, box_lengths.lx);
      const double dy = minimumImageDelta(nodes.com_y_comoving[node_index] - py, box_lengths.ly);
      const double dz = minimumImageDelta(nodes.com_z_comoving[node_index] - pz, box_lengths.lz);
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double r = std::sqrt(std::max(r2, 1.0e-30));

      const bool is_leaf = nodes.child_count[node_index] == 0;
      const double center_dx = nodes.center_x_comoving[node_index] - nodes.com_x_comoving[node_index];
      const double center_dy = nodes.center_y_comoving[node_index] - nodes.com_y_comoving[node_index];
      const double center_dz = nodes.center_z_comoving[node_index] - nodes.com_z_comoving[node_index];
      const double com_offset = std::sqrt(center_dx * center_dx + center_dy * center_dy + center_dz * center_dz);
      const bool target_inside_node = skip_self && !is_leaf &&
          std::abs(minimumImageDelta(nodes.center_x_comoving[node_index] - px, box_lengths.lx)) <= half_size &&
          std::abs(minimumImageDelta(nodes.center_y_comoving[node_index] - py, box_lengths.ly)) <= half_size &&
          std::abs(minimumImageDelta(nodes.center_z_comoving[node_index] - pz, box_lengths.lz)) <= half_size;
      const bool mac_accept = acceptNodeByMac(
          is_leaf,
          target_inside_node,
          half_size,
          com_offset,
          nodes.mass_code[node_index],
          r2,
          previous_acceleration_available,
          previous_acceleration_magnitude_code,
          options.tree_options);
      const bool node_within_cutoff = is_leaf || (maximumDistanceToNodeAabb(
          px,
          py,
          pz,
          nodes.center_x_comoving[node_index],
          nodes.center_y_comoving[node_index],
          nodes.center_z_comoving[node_index],
          half_size,
          box_lengths) <= cutoff_radius_comoving);
      const bool softening_accept = passesSofteningEnvelopeGuard(
          is_leaf,
          half_size,
          r,
          target_softening_comoving,
          nodes.softening_min_comoving[node_index],
          nodes.softening_max_comoving[node_index]);
      // A rank-local forest does not have the same topology as the serial
      // tree.  Bound the residual multipole truncation independently of that
      // topology so changing rank ownership cannot amplify the configured
      // MAC's approximation error past the distributed-equivalence floor.
      // Monopoles recurse to exact leaves; screened quadrupoles retain a small
      // decomposition-independent geometric envelope.
      constexpr double k_screened_quadrupole_width_over_distance_limit = 0.08;
      const bool decomposition_stable_accept = is_leaf ||
          (options.tree_options.multipole_order == TreeMultipoleOrder::kQuadrupole &&
           (2.0 * half_size / r) < k_screened_quadrupole_width_over_distance_limit);
      const bool accept =
          mac_accept && node_within_cutoff && softening_accept && decomposition_stable_accept;

      if (accept) {
        ++accepted_nodes;
        if (is_leaf) {
          const std::uint32_t begin = nodes.particle_begin[node_index];
          const std::uint32_t end = begin + nodes.particle_count[node_index];
          for (std::uint32_t sorted_i = begin; sorted_i < end; ++sorted_i) {
            const std::uint32_t source_index = ordering.sorted_particle_index[sorted_i];
            if (skip_self && source_index == self_index) {
              continue;
            }
            const double sx = minimumImageDelta(pos_x_comoving[source_index] - px, box_lengths.lx);
            const double sy = minimumImageDelta(pos_y_comoving[source_index] - py, box_lengths.ly);
            const double sz = minimumImageDelta(pos_z_comoving[source_index] - pz, box_lengths.lz);
            const double sr2 = sx * sx + sy * sy + sz * sz;
            if (sr2 > cutoff_radius2_comoving) {
              ++cutoff_pair_skips;
              continue;
            }
            const double sr = std::sqrt(std::max(sr2, 1.0e-30));
            const double split_factor = treePmGaussianShortRangeForceFactor(sr, options.split_policy.split_scale_comoving);
            const double source_softening =
                resolveSourceSofteningEpsilon(source_index, options.tree_options.softening, softening_view);
            const double pair_epsilon = combineSofteningPairEpsilon(source_softening, target_softening_comoving);
            // Contract: short-range residual is the softened tree force multiplied by the
            // Gaussian real-space residual factor so that tree+PM composes to the unsplit
            // softened force before explicit r_cut truncation.
            const double softened_factor =
                softenedInvR3(sr2, pair_epsilon) * split_factor * options.tree_options.gravitational_constant_code;
            ax += softened_factor * mass_code[source_index] * sx;
            ay += softened_factor * mass_code[source_index] * sy;
            az += softened_factor * mass_code[source_index] * sz;
            ++pp_interactions;
          }
        } else {
          // Same softened-residual contract as the leaf pair path, applied to accepted nodes.
          const auto contrib =
              monopolePlusQuadrupoleAccelPeriodic(
                  nodes,
                  node_index,
                  dx,
                  dy,
                  dz,
                  target_softening_comoving,
                  options.tree_options,
                  options.split_policy.split_scale_comoving);
          ax += contrib[0];
          ay += contrib[1];
          az += contrib[2];
        }
      } else {
        ++opened_nodes;
        pushChildrenNearFirstPeriodic(nodes, node_index, px, py, pz, box_lengths, stack);
      }
    }
    return std::array<double, 3>{ax, ay, az};
  };

  bool distributed_short_range = false;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  distributed_short_range =
      tree_mpi_world_size > 1 && !rank_local_serial_mode;
#endif
  if (!distributed_short_range) {
    for (std::size_t active_i = 0; active_i < accumulator.active_particle_index.size(); ++active_i) {
      const std::uint32_t particle_index = accumulator.active_particle_index[active_i];
      const bool has_local_source_identity = particle_index < pos_x_comoving.size();
      const double px = m_active_pos_x_comoving[active_i];
      const double py = m_active_pos_y_comoving[active_i];
      const double pz = m_active_pos_z_comoving[active_i];
      const double target_softening =
          resolve_target_softening(active_i, particle_index);
      const bool previous_acceleration_available =
          !accumulator.previous_acceleration_magnitude_code.empty() &&
          std::isfinite(accumulator.previous_acceleration_magnitude_code[active_i]);
      const double previous_acceleration_magnitude_code = previous_acceleration_available
          ? accumulator.previous_acceleration_magnitude_code[active_i]
          : 0.0;
      const auto local_accel = evaluateTargetAgainstLocalTree(
          px,
          py,
          pz,
          has_local_source_identity ? particle_index : 0U,
          has_local_source_identity,
          target_softening,
          previous_acceleration_available,
          previous_acceleration_magnitude_code);
      accumulator.addToActiveSlot(active_i, local_accel[0], local_accel[1], local_accel[2]);
      m_last_residual_stats.local_short_range_sum_sq +=
          local_accel[0] * local_accel[0] + local_accel[1] * local_accel[1] + local_accel[2] * local_accel[2];
    }
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const parallel::PmSlabLayout& layout = m_grid.slabLayout();
  if (distributed_short_range) {
    const int mpi_world_size = tree_mpi_world_size;
    const int mpi_world_rank = tree_mpi_world_rank;
    const auto coordinate_protocol_failure = [&](std::exception_ptr local_failure,
                                                  std::string_view phase) {
      const std::uint64_t local_vote = local_failure ? 1U : 0U;
      std::uint64_t failure_count = 0U;
      MPI_Allreduce(
          &local_vote,
          &failure_count,
          1,
          MPI_UINT64_T,
          MPI_SUM,
          MPI_COMM_WORLD);
      if (failure_count == 0U) {
        return;
      }
      if (local_failure) {
        std::rethrow_exception(local_failure);
      }
      throw std::runtime_error(
          "TreePM peer rejected the " + std::string(phase) +
          " protocol phase");
    };

    std::size_t max_requests_per_peer = 0U;
    std::uint64_t exchange_sequence = 0U;
    std::exception_ptr exchange_preflight_failure;
    try {
      if (mpi_world_size != layout.world_size || mpi_world_rank != layout.world_rank) {
        throw std::invalid_argument("TreePM short-range exchange layout world metadata must match MPI_COMM_WORLD");
      }
      max_requests_per_peer = std::max<std::size_t>(
          1U,
          static_cast<std::size_t>(options.tree_exchange_batch_bytes / k_short_range_request_wire_bytes));
      if (m_tree_exchange_sequence == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("TreePM short-range exchange sequence exhausted");
      }
      exchange_sequence = m_tree_exchange_sequence + 1U;

      if (m_tree_exchange_workspace.world_size != mpi_world_size) {
        m_tree_exchange_workspace.send_counts.assign(static_cast<std::size_t>(mpi_world_size), 0);
        m_tree_exchange_workspace.recv_counts.assign(static_cast<std::size_t>(mpi_world_size), 0);
        m_tree_exchange_workspace.send_displs.assign(static_cast<std::size_t>(mpi_world_size), 0);
        m_tree_exchange_workspace.recv_displs.assign(static_cast<std::size_t>(mpi_world_size), 0);
        m_tree_exchange_workspace.response_send_counts.assign(static_cast<std::size_t>(mpi_world_size), 0);
        m_tree_exchange_workspace.response_recv_counts.assign(static_cast<std::size_t>(mpi_world_size), 0);
        m_tree_exchange_workspace.response_send_displs.assign(static_cast<std::size_t>(mpi_world_size), 0);
        m_tree_exchange_workspace.response_recv_displs.assign(static_cast<std::size_t>(mpi_world_size), 0);
        m_tree_exchange_workspace.world_size = mpi_world_size;
      }
    } catch (...) {
      exchange_preflight_failure = std::current_exception();
    }
    coordinate_protocol_failure(exchange_preflight_failure, "exchange preflight");

    const std::array<std::uint64_t, 4> local_protocol_identity{
        exchange_sequence, options.decomposition_epoch, options.force_epoch, options.tree_exchange_batch_bytes};
    std::array<std::uint64_t, 4> min_protocol_identity{};
    std::array<std::uint64_t, 4> max_protocol_identity{};
    MPI_Allreduce(
        local_protocol_identity.data(),
        min_protocol_identity.data(),
        static_cast<int>(local_protocol_identity.size()),
        MPI_UINT64_T,
        MPI_MIN,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        local_protocol_identity.data(),
        max_protocol_identity.data(),
        static_cast<int>(local_protocol_identity.size()),
        MPI_UINT64_T,
        MPI_MAX,
        MPI_COMM_WORLD);
    if (min_protocol_identity != max_protocol_identity) {
      throw std::runtime_error("TreePM short-range ranks disagree on exchange sequence, runtime epochs, or batch policy");
    }
    m_tree_exchange_sequence = exchange_sequence;

    auto& send_counts = m_tree_exchange_workspace.send_counts;
    auto& recv_counts = m_tree_exchange_workspace.recv_counts;
    auto& send_displs = m_tree_exchange_workspace.send_displs;
    auto& recv_displs = m_tree_exchange_workspace.recv_displs;
    auto& response_send_counts = m_tree_exchange_workspace.response_send_counts;
    auto& response_recv_counts = m_tree_exchange_workspace.response_recv_counts;
    auto& response_send_displs = m_tree_exchange_workspace.response_send_displs;
    auto& response_recv_displs = m_tree_exchange_workspace.response_recv_displs;
    auto& send_payload = m_tree_exchange_workspace.send_payload;
    auto& recv_payload = m_tree_exchange_workspace.recv_payload;
    auto& response_send_payload = m_tree_exchange_workspace.response_send_payload;
    auto& response_recv_payload = m_tree_exchange_workspace.response_recv_payload;
    auto& remote_batch_ax = m_tree_exchange_workspace.remote_batch_ax;
    auto& remote_batch_ay = m_tree_exchange_workspace.remote_batch_ay;
    auto& remote_batch_az = m_tree_exchange_workspace.remote_batch_az;
    auto& expected_response_count = m_tree_exchange_workspace.expected_response_count;
    auto& received_response_count = m_tree_exchange_workspace.received_response_count;

    std::vector<parallel::TreePseudoParticlePacket> local_pseudo_hierarchy;
    std::exception_ptr local_hierarchy_failure;
    try {
      local_pseudo_hierarchy = makeLocalTreePseudoParticleHierarchyPackets(
          mpi_world_rank,
          options.decomposition_epoch,
          options.force_epoch,
          exchange_sequence,
          options.pm_options.boundary_condition == PmBoundaryCondition::kPeriodic,
          m_tree_solver.nodes(),
          /*max_packets=*/std::max<std::size_t>(32U, static_cast<std::size_t>(mpi_world_size) * 16U));
    } catch (...) {
      local_hierarchy_failure = std::current_exception();
    }
    coordinate_protocol_failure(local_hierarchy_failure, "local hierarchy preparation");

    const std::vector<parallel::TreePseudoParticlePacket> peer_pseudo_packets =
        parallel::executeBlockingTreePseudoParticleHierarchyExchange(
            m_mpi_context, local_pseudo_hierarchy, exchange_sequence);
    std::vector<std::vector<parallel::TreePseudoParticlePacket>> peer_hierarchy_by_rank;
    std::vector<std::vector<ShortRangeTargetRequestPacket>> requests_by_peer;
    std::vector<std::vector<std::uint8_t>> response_expected_by_peer;
    std::vector<std::vector<std::uint8_t>> response_seen_by_peer;
    std::vector<std::uint8_t> communicated_with_peer;
    std::exception_ptr peer_hierarchy_failure;
    try {
      if (peer_pseudo_packets.empty()) {
        throw std::runtime_error("TreePM pseudo-particle hierarchy exchange returned no derived nodes");
      }
      peer_hierarchy_by_rank.resize(static_cast<std::size_t>(mpi_world_size));
      for (const parallel::TreePseudoParticlePacket& packet : peer_pseudo_packets) {
        if (packet.descriptor.source_rank < 0 || packet.descriptor.source_rank >= mpi_world_size) {
          throw std::runtime_error("TreePM pseudo-particle hierarchy packet has invalid source rank");
        }
        peer_hierarchy_by_rank[static_cast<std::size_t>(packet.descriptor.source_rank)].push_back(packet);
      }
      for (int peer = 0; peer < mpi_world_size; ++peer) {
        if (peer_hierarchy_by_rank[static_cast<std::size_t>(peer)].empty()) {
          throw std::runtime_error("TreePM pseudo-particle hierarchy exchange returned incomplete rank coverage");
        }
        if (peer != mpi_world_rank) {
          m_last_residual_stats.remote_hierarchy_packets += static_cast<std::uint64_t>(
              peer_hierarchy_by_rank[static_cast<std::size_t>(peer)].size());
        }
      }

      requests_by_peer.resize(static_cast<std::size_t>(mpi_world_size));
      response_expected_by_peer.resize(static_cast<std::size_t>(mpi_world_size));
      response_seen_by_peer.resize(static_cast<std::size_t>(mpi_world_size));
      communicated_with_peer.assign(static_cast<std::size_t>(mpi_world_size), 0U);
    } catch (...) {
      peer_hierarchy_failure = std::current_exception();
    }
    coordinate_protocol_failure(peer_hierarchy_failure, "received hierarchy validation");

    std::uint64_t local_active_count_u64 =
        static_cast<std::uint64_t>(accumulator.active_particle_index.size());
    std::uint64_t global_active_count_max_u64 = 0;
    MPI_Allreduce(
        &local_active_count_u64,
        &global_active_count_max_u64,
        1,
        MPI_UINT64_T,
        MPI_MAX,
        MPI_COMM_WORLD);
    const std::size_t global_active_count_max =
        static_cast<std::size_t>(global_active_count_max_u64);
    for (std::size_t batch_begin = 0; batch_begin < global_active_count_max;) {
      const std::size_t local_remaining =
          batch_begin < accumulator.active_particle_index.size()
          ? accumulator.active_particle_index.size() - batch_begin
          : 0U;
      const std::size_t batch_size = std::min(
          max_requests_per_peer,
          local_remaining);
      std::uint32_t batch_token = 0U;
      int total_send_bytes = 0;
      std::exception_ptr request_preparation_failure;
      try {
        if (batch_begin > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
            batch_size > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
          throw std::overflow_error("TreePM short-range batch token or request ID exceeds packet capacity");
        }
        batch_token = static_cast<std::uint32_t>(batch_begin);

        for (auto& peer_requests : requests_by_peer) {
          peer_requests.clear();
        }
        for (int peer = 0; peer < mpi_world_size; ++peer) {
          response_expected_by_peer[static_cast<std::size_t>(peer)].assign(batch_size, 0U);
          response_seen_by_peer[static_cast<std::size_t>(peer)].assign(batch_size, 0U);
          if (peer != mpi_world_rank) {
            requests_by_peer[static_cast<std::size_t>(peer)].reserve(batch_size);
          }
        }

        expected_response_count.assign(batch_size, 0U);
        received_response_count.assign(batch_size, 0U);
        for (std::size_t batch_slot = 0; batch_slot < batch_size; ++batch_slot) {
        const std::uint32_t particle_index = accumulator.active_particle_index[batch_begin + batch_slot];
        const double px = m_active_pos_x_comoving[batch_begin + batch_slot];
        const double py = m_active_pos_y_comoving[batch_begin + batch_slot];
        const double pz = m_active_pos_z_comoving[batch_begin + batch_slot];
        const double target_softening =
            resolve_target_softening(batch_begin + batch_slot, particle_index);
        const bool previous_acceleration_available =
            !accumulator.previous_acceleration_magnitude_code.empty() &&
            std::isfinite(accumulator.previous_acceleration_magnitude_code[batch_begin + batch_slot]);
        const double previous_acceleration_magnitude_code = previous_acceleration_available
            ? accumulator.previous_acceleration_magnitude_code[batch_begin + batch_slot]
            : 0.0;
        for (int peer = 0; peer < mpi_world_size; ++peer) {
          if (peer == mpi_world_rank) {
            continue;
          }
          const auto& peer_hierarchy = peer_hierarchy_by_rank[static_cast<std::size_t>(peer)];
          if (!remoteTreeHierarchyIntersectsCutoff(
                  px, py, pz, peer_hierarchy, box_lengths, cutoff_radius_comoving)) {
            m_last_residual_stats.remote_pairs_pruned_by_bounds +=
                static_cast<std::uint64_t>(std::count_if(
                    peer_hierarchy.begin(),
                    peer_hierarchy.end(),
                    [](const parallel::TreePseudoParticlePacket& packet) {
                      return packet.source_count > 0U && packet.mass_code > 0.0;
                    }));
            continue;
          }
          requests_by_peer[static_cast<std::size_t>(peer)].push_back(ShortRangeTargetRequestPacket{
              .wire_version = k_short_range_wire_version,
              .origin_rank = mpi_world_rank,
              .destination_rank = peer,
              .flags = previous_acceleration_available ? k_short_range_flag_previous_acceleration : 0U,
              .batch_token = batch_token,
              .request_id = static_cast<std::uint32_t>(batch_slot),
              .exchange_sequence = exchange_sequence,
              .decomposition_epoch = options.decomposition_epoch,
              .force_epoch = options.force_epoch,
              .target_identity = static_cast<std::uint64_t>(batch_begin + batch_slot),
              .target_x_comoving = px,
              .target_y_comoving = py,
              .target_z_comoving = pz,
              .target_softening_epsilon_comoving = target_softening,
              .previous_acceleration_magnitude = previous_acceleration_magnitude_code,
          });
          response_expected_by_peer[static_cast<std::size_t>(peer)][batch_slot] = 1U;
          if (expected_response_count[batch_slot] == std::numeric_limits<std::uint32_t>::max()) {
            throw std::overflow_error("TreePM short-range expected-response counter exceeds packet coverage capacity");
          }
          ++expected_response_count[batch_slot];
        }
        if (expected_response_count[batch_slot] > 0U) {
          ++m_last_residual_stats.remote_targets_with_requests;
        } else {
          ++m_last_residual_stats.remote_targets_without_requests;
        }
        }

        std::fill(send_counts.begin(), send_counts.end(), 0);
        std::fill(recv_counts.begin(), recv_counts.end(), 0);
        std::fill(send_displs.begin(), send_displs.end(), 0);
        std::fill(recv_displs.begin(), recv_displs.end(), 0);
        send_payload.clear();
        for (int peer = 0; peer < mpi_world_size; ++peer) {
        if (peer == mpi_world_rank) {
          continue;
        }
        send_counts[static_cast<std::size_t>(peer)] = checkedMpiByteCount(
            requests_by_peer[static_cast<std::size_t>(peer)].size(),
            k_short_range_request_wire_bytes,
            "TreePM short-range request payload");
        }
        total_send_bytes = populateMpiByteDisplacements(
            send_counts, send_displs, "TreePM short-range request send layout");
        send_payload.resize(static_cast<std::size_t>(total_send_bytes), 0U);
        for (int peer = 0; peer < mpi_world_size; ++peer) {
        if (peer == mpi_world_rank || requests_by_peer[static_cast<std::size_t>(peer)].empty()) {
          continue;
        }
        const auto request_bytes = encodeShortRangeRequests(std::span<const ShortRangeTargetRequestPacket>(
            requests_by_peer[static_cast<std::size_t>(peer)].data(),
            requests_by_peer[static_cast<std::size_t>(peer)].size()));
        std::copy(
            request_bytes.begin(),
            request_bytes.end(),
            send_payload.begin() + send_displs[static_cast<std::size_t>(peer)]);
        }
      } catch (...) {
        request_preparation_failure = std::current_exception();
      }
      coordinate_protocol_failure(
          request_preparation_failure, "outbound request preparation");

      MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
      for (int peer = 0; peer < mpi_world_size; ++peer) {
        if (peer != mpi_world_rank &&
            (send_counts[static_cast<std::size_t>(peer)] > 0 || recv_counts[static_cast<std::size_t>(peer)] > 0)) {
          communicated_with_peer[static_cast<std::size_t>(peer)] = 1U;
        }
      }
      int total_recv_bytes = 0;
      std::exception_ptr request_transport_preparation_failure;
      try {
        total_recv_bytes = populateMpiByteDisplacements(
            recv_counts, recv_displs, "TreePM short-range request receive layout");
        recv_payload.resize(static_cast<std::size_t>(total_recv_bytes), 0U);

        for (std::size_t batch_slot = 0; batch_slot < batch_size; ++batch_slot) {
        const std::uint32_t particle_index = accumulator.active_particle_index[batch_begin + batch_slot];
        const bool has_local_source_identity = particle_index < pos_x_comoving.size();
        const double target_softening =
            resolve_target_softening(batch_begin + batch_slot, particle_index);
        const bool previous_acceleration_available =
            !accumulator.previous_acceleration_magnitude_code.empty() &&
            std::isfinite(accumulator.previous_acceleration_magnitude_code[batch_begin + batch_slot]);
        const double previous_acceleration_magnitude_code = previous_acceleration_available
            ? accumulator.previous_acceleration_magnitude_code[batch_begin + batch_slot]
            : 0.0;
        const auto local_accel = evaluateTargetAgainstLocalTree(
            m_active_pos_x_comoving[batch_begin + batch_slot],
            m_active_pos_y_comoving[batch_begin + batch_slot],
            m_active_pos_z_comoving[batch_begin + batch_slot],
            has_local_source_identity ? particle_index : 0U,
            has_local_source_identity,
            target_softening,
            previous_acceleration_available,
            previous_acceleration_magnitude_code);
        accumulator.addToActiveSlot(batch_begin + batch_slot, local_accel[0], local_accel[1], local_accel[2]);
        m_last_residual_stats.local_short_range_sum_sq +=
            local_accel[0] * local_accel[0] + local_accel[1] * local_accel[1] + local_accel[2] * local_accel[2];
        }
      } catch (...) {
        request_transport_preparation_failure = std::current_exception();
      }
      coordinate_protocol_failure(
          request_transport_preparation_failure,
          "request receive/local-force preparation");

      // All ranks enter both transport collectives for every globally coordinated batch.
      // A rank may have no local payload while peers exchange requests and responses.
      std::uint8_t empty_request_payload = 0U;
      const std::uint8_t* request_send_buffer =
          send_payload.empty() ? &empty_request_payload : send_payload.data();
      std::uint8_t* request_receive_buffer =
          recv_payload.empty() ? &empty_request_payload : recv_payload.data();
      MPI_Alltoallv(
          request_send_buffer,
          send_counts.data(),
          send_displs.data(),
          MPI_BYTE,
          request_receive_buffer,
          recv_counts.data(),
          recv_displs.data(),
          MPI_BYTE,
          MPI_COMM_WORLD);

      m_last_residual_stats.remote_request_batches += 1;
      m_last_residual_stats.remote_request_packets +=
          static_cast<std::uint64_t>(total_send_bytes / static_cast<int>(k_short_range_request_wire_bytes));
      m_last_residual_stats.remote_request_bytes += static_cast<std::uint64_t>(total_send_bytes);
      m_last_residual_stats.remote_peer_participations += static_cast<std::uint64_t>(std::count_if(send_counts.begin(), send_counts.end(), [](int bytes) { return bytes > 0; }));
      m_last_residual_stats.remote_peer_participations += static_cast<std::uint64_t>(std::count_if(recv_counts.begin(), recv_counts.end(), [](int bytes) { return bytes > 0; }));
      {
        std::uint64_t packet_sum = 0;
        std::uint64_t packet_max = 0;
        for (int peer = 0; peer < mpi_world_size; ++peer) {
          if (peer == mpi_world_rank) {
            continue;
          }
          const std::uint64_t packets = static_cast<std::uint64_t>(
              send_counts[static_cast<std::size_t>(peer)] / static_cast<int>(k_short_range_request_wire_bytes));
          packet_sum += packets;
          packet_max = std::max(packet_max, packets);
        }
        m_last_residual_stats.remote_request_packets_max_peer =
            std::max(m_last_residual_stats.remote_request_packets_max_peer, packet_max);
        const double mean_packets = (mpi_world_size > 1)
            ? static_cast<double>(packet_sum) / static_cast<double>(mpi_world_size - 1)
            : 0.0;
        const double imbalance = (mean_packets > 0.0) ? static_cast<double>(packet_max) / mean_packets : 0.0;
        m_last_residual_stats.remote_request_packet_imbalance_ratio =
            std::max(m_last_residual_stats.remote_request_packet_imbalance_ratio, imbalance);
      }

      int total_response_send_bytes = 0;
      int total_response_recv_bytes = 0;
      std::exception_ptr request_validation_failure;
      try {
      std::fill(response_send_counts.begin(), response_send_counts.end(), 0);
      std::fill(response_send_displs.begin(), response_send_displs.end(), 0);
      std::fill(response_recv_counts.begin(), response_recv_counts.end(), 0);
      std::fill(response_recv_displs.begin(), response_recv_displs.end(), 0);
      for (int peer = 0; peer < mpi_world_size; ++peer) {
        if (recv_counts[static_cast<std::size_t>(peer)] % static_cast<int>(k_short_range_request_wire_bytes) != 0) {
          throw std::runtime_error("TreePM short-range request payload size is not request-record aligned");
        }
        if (send_counts[static_cast<std::size_t>(peer)] % static_cast<int>(k_short_range_request_wire_bytes) != 0) {
          throw std::runtime_error("TreePM short-range outbound request payload size is not request-record aligned");
        }
        const int inbound_request_records =
            recv_counts[static_cast<std::size_t>(peer)] / static_cast<int>(k_short_range_request_wire_bytes);
        const int outbound_request_records =
            send_counts[static_cast<std::size_t>(peer)] / static_cast<int>(k_short_range_request_wire_bytes);
        response_send_counts[static_cast<std::size_t>(peer)] = checkedMpiByteCount(
            static_cast<std::size_t>(inbound_request_records),
            k_short_range_response_wire_bytes,
            "TreePM short-range response send payload");
        response_recv_counts[static_cast<std::size_t>(peer)] = checkedMpiByteCount(
            static_cast<std::size_t>(outbound_request_records),
            k_short_range_response_wire_bytes,
            "TreePM short-range response receive payload");
      }
      total_response_send_bytes = populateMpiByteDisplacements(
          response_send_counts, response_send_displs, "TreePM short-range response send layout");
      total_response_recv_bytes = populateMpiByteDisplacements(
          response_recv_counts, response_recv_displs, "TreePM short-range response receive layout");

      response_send_payload.assign(static_cast<std::size_t>(total_response_send_bytes), 0U);
      for (int peer = 0; peer < mpi_world_size; ++peer) {
        const int peer_recv_bytes = recv_counts[static_cast<std::size_t>(peer)];
        if (peer_recv_bytes == 0) {
          continue;
        }
        const std::span<const std::uint8_t> peer_request_bytes(
            recv_payload.data() + recv_displs[static_cast<std::size_t>(peer)],
            static_cast<std::size_t>(peer_recv_bytes));
        const std::vector<ShortRangeTargetRequestPacket> peer_requests =
            decodeShortRangeRequests(peer_request_bytes);

        std::vector<ShortRangeTargetResponsePacket> peer_responses;
        peer_responses.reserve(peer_requests.size());
        std::unordered_set<std::uint64_t> peer_target_identities;
        peer_target_identities.reserve(peer_requests.size());
        for (const ShortRangeTargetRequestPacket& request : peer_requests) {
          if (request.wire_version != k_short_range_wire_version || request.origin_rank != peer ||
              request.destination_rank != mpi_world_rank || request.exchange_sequence != exchange_sequence ||
              request.decomposition_epoch != options.decomposition_epoch || request.force_epoch != options.force_epoch ||
              request.batch_token != batch_token || request.request_id >= max_requests_per_peer ||
              request.target_identity != static_cast<std::uint64_t>(request.batch_token) + request.request_id ||
              (request.flags & ~k_short_range_flag_previous_acceleration) != 0U) {
            throw std::runtime_error("TreePM short-range request protocol identity mismatch");
          }
          if (!peer_target_identities.insert(request.target_identity).second) {
            throw std::runtime_error("TreePM short-range request contains a duplicate target identity");
          }
          if (!std::isfinite(request.target_x_comoving) || !std::isfinite(request.target_y_comoving) ||
              !std::isfinite(request.target_z_comoving) ||
              !std::isfinite(request.target_softening_epsilon_comoving) ||
              request.target_softening_epsilon_comoving < 0.0 ||
              !std::isfinite(request.previous_acceleration_magnitude)) {
            throw std::runtime_error("TreePM short-range request contains invalid numeric data");
          }
          const auto remote_accel = evaluateTargetAgainstLocalTree(
              request.target_x_comoving,
              request.target_y_comoving,
              request.target_z_comoving,
              0U,
              false,
              request.target_softening_epsilon_comoving,
              (request.flags & k_short_range_flag_previous_acceleration) != 0U,
              request.previous_acceleration_magnitude);
          peer_responses.push_back(ShortRangeTargetResponsePacket{
              .wire_version = k_short_range_wire_version,
              .target_owner_rank = request.origin_rank,
              .source_rank = mpi_world_rank,
              .flags = request.flags,
              .batch_token = request.batch_token,
              .request_id = request.request_id,
              .exchange_sequence = request.exchange_sequence,
              .decomposition_epoch = request.decomposition_epoch,
              .force_epoch = request.force_epoch,
              .target_identity = request.target_identity,
              .accel_x_comoving = remote_accel[0],
              .accel_y_comoving = remote_accel[1],
              .accel_z_comoving = remote_accel[2],
          });
        }
        const auto encoded_responses =
            encodeShortRangeResponses(
                std::span<const ShortRangeTargetResponsePacket>(peer_responses.data(), peer_responses.size()));
        if (encoded_responses.size() != static_cast<std::size_t>(response_send_counts[static_cast<std::size_t>(peer)])) {
          throw std::runtime_error("TreePM short-range response payload size mismatch");
        }
        std::copy(
            encoded_responses.begin(),
            encoded_responses.end(),
            response_send_payload.begin() + response_send_displs[static_cast<std::size_t>(peer)]);
      }
      response_recv_payload.resize(static_cast<std::size_t>(total_response_recv_bytes), 0U);
      } catch (...) {
        request_validation_failure = std::current_exception();
      }
      coordinate_protocol_failure(
          request_validation_failure, "received request");

      std::uint8_t empty_response_payload = 0U;
      const std::uint8_t* response_send_buffer =
          response_send_payload.empty() ? &empty_response_payload : response_send_payload.data();
      std::uint8_t* response_receive_buffer =
          response_recv_payload.empty() ? &empty_response_payload : response_recv_payload.data();
      MPI_Alltoallv(
          response_send_buffer,
          response_send_counts.data(),
          response_send_displs.data(),
          MPI_BYTE,
          response_receive_buffer,
          response_recv_counts.data(),
          response_recv_displs.data(),
          MPI_BYTE,
          MPI_COMM_WORLD);

      std::exception_ptr response_validation_failure;
      try {
      remote_batch_ax.assign(batch_size, 0.0);
      remote_batch_ay.assign(batch_size, 0.0);
      remote_batch_az.assign(batch_size, 0.0);
      m_last_residual_stats.remote_response_bytes += static_cast<std::uint64_t>(response_recv_payload.size());
      {
        std::uint64_t response_max = 0;
        for (int peer = 0; peer < mpi_world_size; ++peer) {
          if (peer == mpi_world_rank) {
            continue;
          }
          const std::uint64_t packets = static_cast<std::uint64_t>(
              response_recv_counts[static_cast<std::size_t>(peer)] /
              static_cast<int>(k_short_range_response_wire_bytes));
          response_max = std::max(response_max, packets);
        }
        m_last_residual_stats.remote_response_packets_max_peer =
            std::max(m_last_residual_stats.remote_response_packets_max_peer, response_max);
      }
      for (int peer = 0; peer < mpi_world_size; ++peer) {
        if (peer == mpi_world_rank) {
          continue;
        }
        const int expected_bytes = response_recv_counts[static_cast<std::size_t>(peer)];
        if (expected_bytes == 0) {
          continue;
        }
        const std::span<const std::uint8_t> peer_response_bytes(
            response_recv_payload.data() + response_recv_displs[static_cast<std::size_t>(peer)],
            static_cast<std::size_t>(expected_bytes));
        const std::vector<ShortRangeTargetResponsePacket> responses =
            decodeShortRangeResponses(peer_response_bytes);
        for (const ShortRangeTargetResponsePacket& response : responses) {
          const std::uint32_t expected_flags =
              response.request_id < batch_size &&
                  !accumulator.previous_acceleration_magnitude_code.empty() &&
                  std::isfinite(accumulator.previous_acceleration_magnitude_code[batch_begin + response.request_id])
              ? k_short_range_flag_previous_acceleration
              : 0U;
          if (response.wire_version != k_short_range_wire_version || response.target_owner_rank != mpi_world_rank ||
              response.source_rank != peer || response.flags != expected_flags || response.batch_token != batch_token ||
              response.request_id >= batch_size || response.exchange_sequence != exchange_sequence ||
              response.decomposition_epoch != options.decomposition_epoch || response.force_epoch != options.force_epoch ||
              response.target_identity != static_cast<std::uint64_t>(batch_begin + response.request_id)) {
            throw std::runtime_error("TreePM short-range response protocol identity mismatch");
          }
          if (response_expected_by_peer[static_cast<std::size_t>(peer)][response.request_id] == 0U) {
            throw std::runtime_error("TreePM short-range response was not requested from its source peer");
          }
          if (response_seen_by_peer[static_cast<std::size_t>(peer)][response.request_id] != 0U) {
            throw std::runtime_error("TreePM short-range response duplicates a peer/request identity");
          }
          if (!std::isfinite(response.accel_x_comoving) || !std::isfinite(response.accel_y_comoving) ||
              !std::isfinite(response.accel_z_comoving)) {
            throw std::runtime_error("TreePM short-range response contains non-finite acceleration");
          }
          response_seen_by_peer[static_cast<std::size_t>(peer)][response.request_id] = 1U;
          remote_batch_ax[response.request_id] += response.accel_x_comoving;
          remote_batch_ay[response.request_id] += response.accel_y_comoving;
          remote_batch_az[response.request_id] += response.accel_z_comoving;
          if (received_response_count[response.request_id] == std::numeric_limits<std::uint32_t>::max()) {
            throw std::overflow_error("TreePM short-range received-response counter overflow");
          }
          ++received_response_count[response.request_id];
          ++m_last_residual_stats.remote_response_packets;
        }
      }

      for (int peer = 0; peer < mpi_world_size; ++peer) {
        if (peer == mpi_world_rank) {
          continue;
        }
        for (std::size_t batch_slot = 0; batch_slot < batch_size; ++batch_slot) {
          if (response_expected_by_peer[static_cast<std::size_t>(peer)][batch_slot] !=
              response_seen_by_peer[static_cast<std::size_t>(peer)][batch_slot]) {
            throw std::runtime_error(
                "TreePM short-range response missing for expected peer/request identity: peer=" +
                std::to_string(peer) + ", batch=" + std::to_string(batch_token) +
                ", request=" + std::to_string(batch_slot));
          }
        }
      }

      for (std::size_t batch_slot = 0; batch_slot < batch_size; ++batch_slot) {
        if (received_response_count[batch_slot] != expected_response_count[batch_slot]) {
          throw std::runtime_error(
              "TreePM short-range response coverage mismatch for batch token " + std::to_string(batch_token) +
              ", slot=" + std::to_string(batch_slot) +
              ": expected=" + std::to_string(expected_response_count[batch_slot]) +
              ", received=" + std::to_string(received_response_count[batch_slot]));
        }
        if (received_response_count[batch_slot] == 0U) {
          continue;
        }
        accumulator.addToActiveSlot(
            batch_begin + batch_slot,
            remote_batch_ax[batch_slot],
            remote_batch_ay[batch_slot],
            remote_batch_az[batch_slot]);
        m_last_residual_stats.remote_short_range_sum_sq +=
            remote_batch_ax[batch_slot] * remote_batch_ax[batch_slot] +
            remote_batch_ay[batch_slot] * remote_batch_ay[batch_slot] +
            remote_batch_az[batch_slot] * remote_batch_az[batch_slot];
      }
      } catch (...) {
        response_validation_failure = std::current_exception();
      }
      coordinate_protocol_failure(
          response_validation_failure, "received response");
      batch_begin += max_requests_per_peer;
    }
    m_last_residual_stats.communicating_peer_count = static_cast<std::uint64_t>(std::count(
        communicated_with_peer.begin(), communicated_with_peer.end(), static_cast<std::uint8_t>(1U)));
  }
#else
  if (m_grid.slabLayout().world_size > 1) {
    throw std::invalid_argument("TreePM short-range distributed exchange requires COSMOSIM_ENABLE_MPI=ON");
  }
#endif

  if (tree_profile != nullptr) {
    tree_profile->visited_nodes += visited_nodes;
    tree_profile->accepted_nodes += accepted_nodes;
    tree_profile->opened_nodes += opened_nodes;
    tree_profile->particle_particle_interactions += pp_interactions;
    tree_profile->traversal_ms +=
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - traversal_start).count();
    if (!accumulator.active_particle_index.empty()) {
      tree_profile->average_interactions_per_target = static_cast<double>(tree_profile->particle_particle_interactions) /
          static_cast<double>(accumulator.active_particle_index.size());
    }
  }
  m_last_residual_stats.pruned_nodes = cutoff_pruned_nodes;
  m_last_residual_stats.pair_skips_cutoff = cutoff_pair_skips;
  m_last_residual_stats.pair_evaluations = pp_interactions;
}

TreePmDiagnostics computeTreePmDiagnostics(const TreePmSplitPolicy& split_policy) {
  validateTreePmSplitPolicy(split_policy);

  TreePmDiagnostics diagnostics;
  diagnostics.mesh_spacing_comoving = split_policy.mesh_spacing_comoving;
  diagnostics.asmth_cells = split_policy.asmth_cells;
  diagnostics.rcut_cells = split_policy.rcut_cells;
  diagnostics.split_scale_comoving = split_policy.split_scale_comoving;
  diagnostics.cutoff_radius_comoving = split_policy.cutoff_radius_comoving;
  diagnostics.short_range_factor_at_split =
      treePmGaussianShortRangeForceFactor(split_policy.split_scale_comoving, split_policy.split_scale_comoving);
  diagnostics.long_range_factor_at_split =
      treePmGaussianLongRangeForceFactor(split_policy.split_scale_comoving, split_policy.split_scale_comoving);
  diagnostics.short_range_factor_at_cutoff =
      treePmGaussianShortRangeForceFactor(split_policy.cutoff_radius_comoving, split_policy.split_scale_comoving);
  diagnostics.long_range_factor_at_cutoff =
      treePmGaussianLongRangeForceFactor(split_policy.cutoff_radius_comoving, split_policy.split_scale_comoving);
  diagnostics.composition_error_at_split = std::abs(
      diagnostics.short_range_factor_at_split + diagnostics.long_range_factor_at_split - 1.0);

  const double radii[] = {
      0.25 * split_policy.split_scale_comoving,
      0.5 * split_policy.split_scale_comoving,
      split_policy.split_scale_comoving,
      2.0 * split_policy.split_scale_comoving,
      4.0 * split_policy.split_scale_comoving,
  };
  for (const double radius_comoving : radii) {
    const double composed = treePmGaussianShortRangeForceFactor(radius_comoving, split_policy.split_scale_comoving) +
        treePmGaussianLongRangeForceFactor(radius_comoving, split_policy.split_scale_comoving);
    diagnostics.max_relative_composition_error = std::max(
        diagnostics.max_relative_composition_error,
        std::abs(composed - 1.0) / std::max(1.0, std::abs(composed)));
  }

  return diagnostics;
}

}  // namespace cosmosim::gravity
