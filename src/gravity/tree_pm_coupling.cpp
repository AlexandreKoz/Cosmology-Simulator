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

// shared hot tree interaction invariants
#include "internal/tree_interaction_common.hpp"

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

// Wrapped periodic coordinates are finite and non-negative, so IEEE-754 bit
// order equals numerical order. Eight stable byte passes remove the three
// O(N log N) comparison sorts from every periodic tree rebuild while reusing
// buffers already owned by the coordinator.
void radixSortNonNegativeDoubles(
    std::vector<double>& values,
    std::vector<double>& scratch) {
  if (values.size() < 2U) {
    return;
  }
  scratch.resize(values.size());
  for (unsigned pass = 0; pass < 8U; ++pass) {
    std::array<std::size_t, 256> counts{};
    const auto& source = (pass % 2U == 0U) ? values : scratch;
    auto& destination = (pass % 2U == 0U) ? scratch : values;
    const unsigned shift = pass * 8U;
    for (const double value : source) {
      const std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
      ++counts[(bits >> shift) & 0xffU];
    }
    std::array<std::size_t, 256> offsets{};
    for (std::size_t bucket = 1U; bucket < offsets.size(); ++bucket) {
      offsets[bucket] = offsets[bucket - 1U] + counts[bucket - 1U];
    }
    for (const double value : source) {
      const std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
      const std::size_t bucket = static_cast<std::size_t>((bits >> shift) & 0xffU);
      destination[offsets[bucket]++] = value;
    }
  }
}

// Map one periodic coordinate lane into the shortest contiguous interval that
// contains all sources. The interval begins immediately after the largest
// circular gap. This gives topology, COMs, quadrupoles, MAC geometry, and
// exported hierarchy bounds one coherent unwrapped frame.
void unwrapPeriodicAxis(
    std::span<const double> input,
    double box_size_comoving,
    std::vector<double>& output,
    std::vector<double>& wrapped,
    std::vector<double>& ordered) {
  if (!std::isfinite(box_size_comoving) || box_size_comoving <= 0.0) {
    throw std::invalid_argument("Periodic TreePM tree geometry requires finite positive axis lengths");
  }
  output.resize(input.size());
  if (input.empty()) {
    return;
  }

  wrapped.assign(input.size(), 0.0);
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

  ordered.assign(wrapped.begin(), wrapped.end());
  // output is not yet authoritative here, so reuse it as radix scratch.
  radixSortNonNegativeDoubles(ordered, output);
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
    const TreeLocalIndex source_index = accumulator.active_particle_index[i];
    if (!explicit_targets && source_index >= pos_x_comoving.size()) {
      throw std::out_of_range("TreePM active index exceeds particle count");
    }
    if (explicit_targets && source_index >= pos_x_comoving.size() &&
        source_index != kInvalidTreeLocalIndex) {
      throw std::out_of_range(
          "TreePM explicit target source identity must be local or the invalid local-index sentinel");
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
  const auto validate_species_lane = [](std::span<const std::uint32_t> lane) {
    for (const std::uint32_t species_tag : lane) {
      if (!core::isValidParticleSpeciesTag(species_tag)) {
        throw std::invalid_argument("TreePM species sidecars contain an invalid species tag");
      }
    }
  };
  validate_species_lane(softening_view.source_species_tag);
  validate_species_lane(softening_view.target_species_tag);
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
    for (const TreeLocalIndex source_index : accumulator.active_particle_index) {
      if (source_index == kInvalidTreeLocalIndex) {
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

[[nodiscard, maybe_unused]] std::uint64_t treePmOptionFingerprint(
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
  mix(options.decomposition_epoch.value);
  mix(options.source_generation.value);
  mix(options.force_epoch.sequence);
  mix_double(options.force_epoch.scale_factor);
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
    TreeLocalIndex node_index,
    double dx,
    double dy,
    double dz,
    double target_softening_comoving,
    const TreeGravityOptions& options,
    double split_scale_comoving) {
  const double r2 = dx * dx + dy * dy + dz * dz;
  const double r = std::sqrt(std::max(r2, 1.0e-30));
  const double pair_epsilon =
      combineSofteningPairEpsilonUnchecked(nodes.softening_max_comoving[node_index], target_softening_comoving);
  const double eps2 = pair_epsilon * pair_epsilon;
  const double denom = std::max(r2 + eps2, 1.0e-30);
  const double softened_inv_r3 = 1.0 / (denom * std::sqrt(denom));
  const double split_factor = treePmGaussianShortRangeForceFactorUnchecked(r, split_scale_comoving);
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
    TreeLocalIndex node_index,
    double px,
    double py,
    double pz,
    const PeriodicBoxLengths& box_lengths,
    std::vector<TreeLocalIndex>& stack) {
  std::array<std::pair<double, TreeLocalIndex>, 8> child_dist2{};
  std::size_t count = 0;
  const std::size_t child_offset = static_cast<std::size_t>(node_index) * 8U;
  for (std::uint8_t octant = 0; octant < 8U; ++octant) {
    const TreeLocalIndex child = nodes.child_index[child_offset + octant];
    if (child == kInvalidTreeLocalIndex) {
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

struct SparsePeerGraph {
  MPI_Comm communicator = MPI_COMM_NULL;
  std::vector<int> incoming_peers;
  std::vector<int> outgoing_peers;

  static void releaseCommunicator(MPI_Comm& communicator) noexcept {
    if (communicator == MPI_COMM_NULL) {
      return;
    }
    int initialized = 0;
    int finalized = 0;
    (void)MPI_Initialized(&initialized);
    if (initialized != 0) {
      (void)MPI_Finalized(&finalized);
    }
    if (initialized != 0 && finalized == 0) {
      (void)MPI_Comm_free(&communicator);
    } else {
      // A cached communicator may outlive MPI finalization during process
      // teardown. MPI objects cannot be freed after MPI_Finalize; dropping the
      // handle here is the only legal teardown and process exit reclaims it.
      communicator = MPI_COMM_NULL;
    }
  }

  SparsePeerGraph() = default;
  SparsePeerGraph(const SparsePeerGraph&) = delete;
  SparsePeerGraph& operator=(const SparsePeerGraph&) = delete;
  SparsePeerGraph(SparsePeerGraph&& other) noexcept
      : communicator(std::exchange(other.communicator, MPI_COMM_NULL)),
        incoming_peers(std::move(other.incoming_peers)),
        outgoing_peers(std::move(other.outgoing_peers)) {}
  SparsePeerGraph& operator=(SparsePeerGraph&& other) noexcept {
    if (this != &other) {
      releaseCommunicator(communicator);
      communicator = std::exchange(other.communicator, MPI_COMM_NULL);
      incoming_peers = std::move(other.incoming_peers);
      outgoing_peers = std::move(other.outgoing_peers);
    }
    return *this;
  }
  ~SparsePeerGraph() { releaseCommunicator(communicator); }
};

[[nodiscard]] SparsePeerGraph makeSymmetricSparsePeerGraph(
    int world_rank,
    std::span<const int> requested_outgoing_peers) {
  MPI_Comm directed = MPI_COMM_NULL;
  const int source = world_rank;
  if (requested_outgoing_peers.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("TreePM sparse peer degree exceeds MPI int capacity");
  }
  const int degree = static_cast<int>(requested_outgoing_peers.size());
  const int* destinations = requested_outgoing_peers.empty()
      ? nullptr
      : requested_outgoing_peers.data();
  MPI_Dist_graph_create(
      MPI_COMM_WORLD,
      1,
      &source,
      &degree,
      destinations,
      MPI_UNWEIGHTED,
      MPI_INFO_NULL,
      0,
      &directed);
  if (directed == MPI_COMM_NULL) {
    throw std::runtime_error("TreePM failed to create directed sparse peer graph");
  }

  int indegree = 0;
  int outdegree = 0;
  int weighted = 0;
  MPI_Dist_graph_neighbors_count(directed, &indegree, &outdegree, &weighted);
  std::vector<int> incoming(static_cast<std::size_t>(std::max(indegree, 0)));
  std::vector<int> outgoing(static_cast<std::size_t>(std::max(outdegree, 0)));
  MPI_Dist_graph_neighbors(
      directed,
      indegree,
      incoming.empty() ? nullptr : incoming.data(),
      MPI_UNWEIGHTED,
      outdegree,
      outgoing.empty() ? nullptr : outgoing.data(),
      MPI_UNWEIGHTED);
  MPI_Comm_free(&directed);

  std::vector<int> symmetric_peers = incoming;
  symmetric_peers.insert(symmetric_peers.end(), outgoing.begin(), outgoing.end());
  std::sort(symmetric_peers.begin(), symmetric_peers.end());
  symmetric_peers.erase(
      std::unique(symmetric_peers.begin(), symmetric_peers.end()),
      symmetric_peers.end());
  symmetric_peers.erase(
      std::remove(symmetric_peers.begin(), symmetric_peers.end(), world_rank),
      symmetric_peers.end());
  if (symmetric_peers.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("TreePM symmetric sparse peer degree exceeds MPI int capacity");
  }

  SparsePeerGraph graph;
  const int symmetric_degree = static_cast<int>(symmetric_peers.size());
  const int* symmetric_data = symmetric_peers.empty() ? nullptr : symmetric_peers.data();
  MPI_Dist_graph_create_adjacent(
      MPI_COMM_WORLD,
      symmetric_degree,
      symmetric_data,
      MPI_UNWEIGHTED,
      symmetric_degree,
      symmetric_data,
      MPI_UNWEIGHTED,
      MPI_INFO_NULL,
      0,
      &graph.communicator);
  if (graph.communicator == MPI_COMM_NULL) {
    throw std::runtime_error("TreePM failed to create symmetric sparse peer graph");
  }
  MPI_Dist_graph_neighbors_count(
      graph.communicator, &indegree, &outdegree, &weighted);
  graph.incoming_peers.resize(static_cast<std::size_t>(std::max(indegree, 0)));
  graph.outgoing_peers.resize(static_cast<std::size_t>(std::max(outdegree, 0)));
  MPI_Dist_graph_neighbors(
      graph.communicator,
      indegree,
      graph.incoming_peers.empty() ? nullptr : graph.incoming_peers.data(),
      MPI_UNWEIGHTED,
      outdegree,
      graph.outgoing_peers.empty() ? nullptr : graph.outgoing_peers.data(),
      MPI_UNWEIGHTED);
  return graph;
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

[[nodiscard, maybe_unused]] std::vector<std::uint8_t> encodeShortRangeRequests(
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

[[nodiscard, maybe_unused]] std::vector<ShortRangeTargetRequestPacket> decodeShortRangeRequests(
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

[[nodiscard, maybe_unused]] std::vector<std::uint8_t> encodeShortRangeResponses(
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

[[nodiscard, maybe_unused]] std::vector<ShortRangeTargetResponsePacket> decodeShortRangeResponses(
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

[[nodiscard, maybe_unused]] parallel::TreePseudoParticlePacket makeLocalTreePseudoParticlePacket(
    int world_rank,
    std::uint64_t decomposition_epoch,
    std::uint64_t force_epoch,
    std::uint64_t exchange_sequence,
    const TreeNodeSoa& nodes,
    std::size_t source_count) {
  parallel::TreePseudoParticlePacket packet;
  packet.descriptor = parallel::TreePseudoParticleDescriptor{
      .pseudo_particle_id = (static_cast<std::uint64_t>(std::max(world_rank, 0)) << 32U) ^ decomposition_epoch,
      .source_rank = world_rank,
      .decomposition_epoch = decomposition_epoch,
      .force_epoch = force_epoch,
      .exchange_sequence = exchange_sequence,
      .derived_not_authoritative = true,
  };
  packet.source_count = static_cast<std::uint64_t>(source_count);
  packet.geometry_frame = 1U;
  packet.is_leaf = 1U;
  if (source_count == 0U || nodes.size() == 0U) {
    parallel::validateTreePseudoParticlePacket(packet);
    return packet;
  }
  const std::size_t root = 0U;
  const double half_size = nodes.half_size_comoving[root];
  packet.mass_code = nodes.mass_code[root];
  packet.center_x_comoving = nodes.com_x_comoving[root];
  packet.center_y_comoving = nodes.com_y_comoving[root];
  packet.center_z_comoving = nodes.com_z_comoving[root];
  packet.min_x_comoving = nodes.center_x_comoving[root] - half_size;
  packet.max_x_comoving = nodes.center_x_comoving[root] + half_size;
  packet.min_y_comoving = nodes.center_y_comoving[root] - half_size;
  packet.max_y_comoving = nodes.center_y_comoving[root] + half_size;
  packet.min_z_comoving = nodes.center_z_comoving[root] - half_size;
  packet.max_z_comoving = nodes.center_z_comoving[root] + half_size;
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

[[nodiscard, maybe_unused]] bool remoteTreeHierarchyIntersectsCutoff(
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


class TopLevelDomainHierarchy {
 public:
  explicit TopLevelDomainHierarchy(
      std::span<const parallel::TreePseudoParticlePacket> leaves)
      : m_leaves(leaves.begin(), leaves.end()) {
    m_order.resize(m_leaves.size());
    std::iota(m_order.begin(), m_order.end(), 0U);
    if (!m_order.empty()) {
      m_nodes.reserve(m_order.size() * 2U);
      (void)buildNode(0U, m_order.size());
    }
  }

  [[nodiscard]] std::size_t nodeCount() const noexcept { return m_nodes.size(); }

  [[nodiscard]] std::vector<int> ownersWithinCutoff(
      double px, double py, double pz,
      double cutoff_radius_comoving,
      const PeriodicBoxLengths& box_lengths,
      int excluded_rank) const {
    std::vector<int> owners;
    if (m_nodes.empty()) {
      return owners;
    }
    std::vector<std::size_t> stack{0U};
    while (!stack.empty()) {
      const std::size_t node_index = stack.back();
      stack.pop_back();
      const Node& node = m_nodes[node_index];
      if (minimumDistanceToPeriodicBounds(
              px, py, pz, node.bounds, box_lengths) > cutoff_radius_comoving) {
        continue;
      }
      if (node.leaf_count != 0U) {
        for (std::size_t slot = node.leaf_begin;
             slot < node.leaf_begin + node.leaf_count; ++slot) {
          const auto& leaf = m_leaves[m_order[slot]];
          if (leaf.descriptor.source_rank == excluded_rank ||
              leaf.source_count == 0U || leaf.mass_code == 0.0) {
            continue;
          }
          if (minimumDistanceToPeriodicBounds(
                  px, py, pz, boundsFromTreePseudoParticlePacket(leaf),
                  box_lengths) <= cutoff_radius_comoving) {
            owners.push_back(leaf.descriptor.source_rank);
          }
        }
      } else {
        stack.push_back(node.left);
        stack.push_back(node.right);
      }
    }
    std::sort(owners.begin(), owners.end());
    owners.erase(std::unique(owners.begin(), owners.end()), owners.end());
    return owners;
  }

 private:
  struct Node {
    SourceDomainBoundsPacket bounds{};
    std::size_t left = 0U;
    std::size_t right = 0U;
    std::size_t leaf_begin = 0U;
    std::size_t leaf_count = 0U;
  };

  [[nodiscard]] std::size_t buildNode(std::size_t begin, std::size_t end) {
    Node node;
    const auto first_bounds = boundsFromTreePseudoParticlePacket(
        m_leaves[m_order[begin]]);
    node.bounds = first_bounds;
    for (std::size_t slot = begin + 1U; slot < end; ++slot) {
      const auto b = boundsFromTreePseudoParticlePacket(m_leaves[m_order[slot]]);
      node.bounds.min_x_comoving = std::min(node.bounds.min_x_comoving, b.min_x_comoving);
      node.bounds.max_x_comoving = std::max(node.bounds.max_x_comoving, b.max_x_comoving);
      node.bounds.min_y_comoving = std::min(node.bounds.min_y_comoving, b.min_y_comoving);
      node.bounds.max_y_comoving = std::max(node.bounds.max_y_comoving, b.max_y_comoving);
      node.bounds.min_z_comoving = std::min(node.bounds.min_z_comoving, b.min_z_comoving);
      node.bounds.max_z_comoving = std::max(node.bounds.max_z_comoving, b.max_z_comoving);
      node.bounds.source_particle_count += b.source_particle_count;
    }
    const std::size_t node_index = m_nodes.size();
    m_nodes.push_back(node);
    const std::size_t count = end - begin;
    if (count <= 2U) {
      m_nodes[node_index].leaf_begin = begin;
      m_nodes[node_index].leaf_count = count;
      return node_index;
    }
    const double ex = node.bounds.max_x_comoving - node.bounds.min_x_comoving;
    const double ey = node.bounds.max_y_comoving - node.bounds.min_y_comoving;
    const double ez = node.bounds.max_z_comoving - node.bounds.min_z_comoving;
    const int axis = ey > ex && ey >= ez ? 1 : (ez > ex && ez > ey ? 2 : 0);
    const std::size_t middle = begin + count / 2U;
    std::nth_element(
        m_order.begin() + static_cast<std::ptrdiff_t>(begin),
        m_order.begin() + static_cast<std::ptrdiff_t>(middle),
        m_order.begin() + static_cast<std::ptrdiff_t>(end),
        [&](std::size_t lhs, std::size_t rhs) {
          const auto center = [&](const parallel::TreePseudoParticlePacket& p) {
            if (axis == 1) return p.center_y_comoving;
            if (axis == 2) return p.center_z_comoving;
            return p.center_x_comoving;
          };
          return center(m_leaves[lhs]) < center(m_leaves[rhs]);
        });
    m_nodes[node_index].left = buildNode(begin, middle);
    m_nodes[node_index].right = buildNode(middle, end);
    return node_index;
  }

  std::vector<parallel::TreePseudoParticlePacket> m_leaves;
  std::vector<std::size_t> m_order;
  std::vector<Node> m_nodes;
};

[[nodiscard]] std::uint64_t domainGeometryFingerprint(
    std::span<const parallel::TreePseudoParticlePacket> leaves) {
  std::uint64_t hash = 1469598103934665603ULL;
  auto mix = [&](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  for (const auto& leaf : leaves) {
    mix(static_cast<std::uint64_t>(static_cast<std::uint32_t>(leaf.descriptor.source_rank)));
    mix(leaf.descriptor.decomposition_epoch);
    mix(std::bit_cast<std::uint64_t>(leaf.min_x_comoving));
    mix(std::bit_cast<std::uint64_t>(leaf.max_x_comoving));
    mix(std::bit_cast<std::uint64_t>(leaf.min_y_comoving));
    mix(std::bit_cast<std::uint64_t>(leaf.max_y_comoving));
    mix(std::bit_cast<std::uint64_t>(leaf.min_z_comoving));
    mix(std::bit_cast<std::uint64_t>(leaf.max_z_comoving));
  }
  return hash;
}

}  // namespace

struct TreePmCoordinator::SparsePeerGraphCacheOpaque {
  bool valid = false;
  DecompositionEpoch decomposition_epoch{};
  int world_size = 1;
  std::vector<int> requested_outgoing_peers;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  SparsePeerGraph graph;
#endif
};

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
      m_tree_solver(),
      m_sparse_peer_graph_cache(std::make_unique<SparsePeerGraphCacheOpaque>()) {}

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
      m_tree_solver(),
      m_sparse_peer_graph_cache(std::make_unique<SparsePeerGraphCacheOpaque>()) {}

TreePmCoordinator::~TreePmCoordinator() = default;

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
  m_tree_solver.appendMemoryReport(builder);
  const auto add_active = [&builder](
                              std::string label,
                              const auto& container,
                              std::uint64_t high_water_bytes) {
    const std::uint64_t bytes = core::ownedCapacityBytesForContainer(container);
    builder.addEntry(core::MemoryEntry{.subsystem = core::MemorySubsystem::kActiveSets,
                                       .lifetime = core::MemoryLifetime::kTransient,
                                       .label = std::move(label),
                                       .current_size_bytes = core::currentSizeBytesForContainer(container),
                                       .owned_capacity_bytes = bytes,
                                       .high_water_bytes = high_water_bytes,
                                       .estimated_next_step_bytes = 0U,
                                       .uncertainty_note = "historical high-water tracked; next-step zoom requirement is profile-dependent"});
  };
  add_active("treepm.active_zoom_corr_ax_comoving", m_active_zoom_corr_ax_comoving, m_zoom_corr_high_water_bytes[0]);
  add_active("treepm.active_zoom_corr_ay_comoving", m_active_zoom_corr_ay_comoving, m_zoom_corr_high_water_bytes[1]);
  add_active("treepm.active_zoom_corr_az_comoving", m_active_zoom_corr_az_comoving, m_zoom_corr_high_water_bytes[2]);

  const auto add_tree_scratch = [&builder](std::string label, const auto& container) {
    const std::uint64_t bytes = core::ownedCapacityBytesForContainer(container);
    builder.addEntry(core::MemoryEntry{.subsystem = core::MemorySubsystem::kTree,
                                       .lifetime = core::MemoryLifetime::kTransient,
                                       .label = std::move(label),
                                       .current_size_bytes = core::currentSizeBytesForContainer(container),
                                       .owned_capacity_bytes = bytes,
                                       .high_water_bytes = bytes,
                                       .estimated_next_step_bytes = 0U,
                                       .uncertainty_note = "next-step requirement not predicted from retained capacity"});
  };
  add_tree_scratch("treepm.periodic_tree_source_x_comoving", m_tree_source_x_comoving);
  add_tree_scratch("treepm.periodic_tree_source_y_comoving", m_tree_source_y_comoving);
  add_tree_scratch("treepm.periodic_tree_source_z_comoving", m_tree_source_z_comoving);

  const auto add_mpi = [&builder](std::string label, const auto& container) {
    const std::uint64_t bytes = core::ownedCapacityBytesForContainer(container);
    builder.addEntry(core::MemoryEntry{.subsystem = core::MemorySubsystem::kMpiBuffers,
                                       .lifetime = core::MemoryLifetime::kTransient,
                                       .label = std::move(label),
                                       .current_size_bytes = core::currentSizeBytesForContainer(container),
                                       .owned_capacity_bytes = bytes,
                                       .high_water_bytes = bytes,
                                       .estimated_next_step_bytes = 0U,
                                       .uncertainty_note = "next-step requirement not predicted from retained capacity"});
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
        m_long_range_field_validity.source_generation == options.source_generation &&
        m_long_range_field_validity.last_force_epoch == options.force_epoch &&
        m_long_range_field_validity.pm_field_version == options.pm_field_version &&
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
        .source_generation = options.source_generation,
        .pm_field_version = options.pm_field_version,
        .last_force_epoch = options.force_epoch,
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
    if (has_explicit_target_positions) {
      if (accumulator.target_pos_x_comoving.size() != active_count ||
          accumulator.target_pos_y_comoving.size() != active_count ||
          accumulator.target_pos_z_comoving.size() != active_count) {
        throw std::invalid_argument(
            "TreePM independent target coordinate extents must match active count");
      }
    } else if (!accumulator.target_pos_y_comoving.empty() ||
               !accumulator.target_pos_z_comoving.empty()) {
      throw std::invalid_argument(
          "TreePM independent target coordinate triplet must be entirely present or absent");
    }
    for (std::size_t i = 0; i < active_count; ++i) {
      const TreeLocalIndex source_index = accumulator.active_particle_index[i];
      if (!has_explicit_target_positions &&
          static_cast<std::size_t>(source_index) >= pos_x_comoving.size()) {
        throw std::out_of_range("TreePM active source index exceeds source count");
      }
      const double px = has_explicit_target_positions
          ? accumulator.target_pos_x_comoving[i]
          : pos_x_comoving[source_index];
      const double py = has_explicit_target_positions
          ? accumulator.target_pos_y_comoving[i]
          : pos_y_comoving[source_index];
      const double pz = has_explicit_target_positions
          ? accumulator.target_pos_z_comoving[i]
          : pos_z_comoving[source_index];
      if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz)) {
        throw std::invalid_argument("TreePM target positions must be finite");
      }
    }
    if (options.enable_zoom_long_range_correction) {
      resizeCompactSidecars(
          m_active_zoom_corr_ax_comoving,
          m_active_zoom_corr_ay_comoving,
          m_active_zoom_corr_az_comoving,
          active_count);
      m_zoom_corr_high_water_bytes[0] = std::max(
          m_zoom_corr_high_water_bytes[0],
          core::ownedCapacityBytesForContainer(m_active_zoom_corr_ax_comoving));
      m_zoom_corr_high_water_bytes[1] = std::max(
          m_zoom_corr_high_water_bytes[1],
          core::ownedCapacityBytesForContainer(m_active_zoom_corr_ay_comoving));
      m_zoom_corr_high_water_bytes[2] = std::max(
          m_zoom_corr_high_water_bytes[2],
          core::ownedCapacityBytesForContainer(m_active_zoom_corr_az_comoving));
      std::fill(m_active_zoom_corr_ax_comoving.begin(), m_active_zoom_corr_ax_comoving.end(), 0.0);
      std::fill(m_active_zoom_corr_ay_comoving.begin(), m_active_zoom_corr_ay_comoving.end(), 0.0);
      std::fill(m_active_zoom_corr_az_comoving.begin(), m_active_zoom_corr_az_comoving.end(), 0.0);
    } else {
      std::vector<double>().swap(m_active_zoom_corr_ax_comoving);
      std::vector<double>().swap(m_active_zoom_corr_ay_comoving);
      std::vector<double>().swap(m_active_zoom_corr_az_comoving);
    }
  } catch (...) {
    compact_target_preparation_failure = std::current_exception();
  }
  coordinate_tree_pm_failure(
      compact_target_preparation_failure, "compact PM target preparation");

  const auto target_x = [&](std::size_t active_i) -> double {
    return has_explicit_target_positions
        ? accumulator.target_pos_x_comoving[active_i]
        : pos_x_comoving[accumulator.active_particle_index[active_i]];
  };
  const auto target_y = [&](std::size_t active_i) -> double {
    return has_explicit_target_positions
        ? accumulator.target_pos_y_comoving[active_i]
        : pos_y_comoving[accumulator.active_particle_index[active_i]];
  };
  const auto target_z = [&](std::size_t active_i) -> double {
    return has_explicit_target_positions
        ? accumulator.target_pos_z_comoving[active_i]
        : pos_z_comoving[accumulator.active_particle_index[active_i]];
  };

  m_pm_solver.interpolateForces(
      m_grid,
      PmSolver::PmForceTargetView{
          .active_particle_index = accumulator.active_particle_index,
          .pos_x_comoving = has_explicit_target_positions
              ? accumulator.target_pos_x_comoving : pos_x_comoving,
          .pos_y_comoving = has_explicit_target_positions
              ? accumulator.target_pos_y_comoving : pos_y_comoving,
          .pos_z_comoving = has_explicit_target_positions
              ? accumulator.target_pos_z_comoving : pos_z_comoving,
          .accel_x_comoving = accumulator.accel_x_comoving,
          .accel_y_comoving = accumulator.accel_y_comoving,
          .accel_z_comoving = accumulator.accel_z_comoving,
          .output_layout = PmSolver::PmForceOutputLayout::kCompactActive,
          .coordinate_source_index = has_explicit_target_positions
              ? std::span<const TreeLocalIndex>{}
              : accumulator.active_particle_index,
      },
      pm_options,
      profile != nullptr ? &profile->pm_profile : nullptr);
  const double pm_force_l2_global = l2NormFromComponents(
      accumulator.accel_x_comoving,
      accumulator.accel_y_comoving,
      accumulator.accel_z_comoving);

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
      std::vector<double> high_res_pm_coarse_ax(active_count, 0.0);
      std::vector<double> high_res_pm_coarse_ay(active_count, 0.0);
      std::vector<double> high_res_pm_coarse_az(active_count, 0.0);

      // The coarse and focused PM grids can each be a dominant temporary.
      // Interpolate the coarse contribution while its grid is live and release
      // it before allocating the focused grid so peak memory is bounded by one
      // correction mesh rather than their sum.
      {
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
        m_pm_solver.interpolateForces(
            high_res_coarse_grid,
            PmSolver::PmForceTargetView{
                .active_particle_index = accumulator.active_particle_index,
                .pos_x_comoving = has_explicit_target_positions
                    ? accumulator.target_pos_x_comoving : pos_x_comoving,
                .pos_y_comoving = has_explicit_target_positions
                    ? accumulator.target_pos_y_comoving : pos_y_comoving,
                .pos_z_comoving = has_explicit_target_positions
                    ? accumulator.target_pos_z_comoving : pos_z_comoving,
                .accel_x_comoving = high_res_pm_coarse_ax,
                .accel_y_comoving = high_res_pm_coarse_ay,
                .accel_z_comoving = high_res_pm_coarse_az,
                .output_layout = PmSolver::PmForceOutputLayout::kCompactActive,
                .coordinate_source_index = has_explicit_target_positions
                    ? std::span<const TreeLocalIndex>{}
                    : accumulator.active_particle_index,
            },
            pm_options,
            profile != nullptr ? &profile->pm_profile : nullptr);
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

      std::vector<double> focused_active_x(active_count, focused_half_extent);
      std::vector<double> focused_active_y(active_count, focused_half_extent);
      std::vector<double> focused_active_z(active_count, focused_half_extent);
      for (std::size_t i = 0; i < active_count; ++i) {
        const double dx_local = (pm_options.boundary_condition == PmBoundaryCondition::kPeriodic)
            ? minimumImageDelta(target_x(i) - options.zoom_region_center_x_comoving, box_lengths.lx)
            : (target_x(i) - options.zoom_region_center_x_comoving);
        const double dy_local = (pm_options.boundary_condition == PmBoundaryCondition::kPeriodic)
            ? minimumImageDelta(target_y(i) - options.zoom_region_center_y_comoving, box_lengths.ly)
            : (target_y(i) - options.zoom_region_center_y_comoving);
        const double dz_local = (pm_options.boundary_condition == PmBoundaryCondition::kPeriodic)
            ? minimumImageDelta(target_z(i) - options.zoom_region_center_z_comoving, box_lengths.lz)
            : (target_z(i) - options.zoom_region_center_z_comoving);
        focused_active_x[i] = dx_local + focused_half_extent;
        focused_active_y[i] = dy_local + focused_half_extent;
        focused_active_z[i] = dz_local + focused_half_extent;
      }
      zoom_solver.interpolateForces(
          high_res_focused_grid,
          focused_active_x,
          focused_active_y,
          focused_active_z,
          m_active_zoom_corr_ax_comoving,
          m_active_zoom_corr_ay_comoving,
          m_active_zoom_corr_az_comoving,
          zoom_pm_options,
          profile != nullptr ? &profile->pm_profile : nullptr);

      for (std::size_t i = 0; i < active_count; ++i) {
        if (options.active_is_high_res[i] == 0U) {
          m_active_zoom_corr_ax_comoving[i] = 0.0;
          m_active_zoom_corr_ay_comoving[i] = 0.0;
          m_active_zoom_corr_az_comoving[i] = 0.0;
          continue;
        }
        const double corr_x = m_active_zoom_corr_ax_comoving[i] - high_res_pm_coarse_ax[i];
        const double corr_y = m_active_zoom_corr_ay_comoving[i] - high_res_pm_coarse_ay[i];
        const double corr_z = m_active_zoom_corr_az_comoving[i] - high_res_pm_coarse_az[i];
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
    const auto source_preprocess_start = std::chrono::steady_clock::now();
    if (options.pm_options.boundary_condition == PmBoundaryCondition::kPeriodic) {
      const PeriodicBoxLengths box_lengths = effectivePeriodicBoxLengths(options.pm_options);
      unwrapPeriodicAxis(
          pos_x_comoving, box_lengths.lx, m_tree_source_x_comoving,
          m_periodic_wrapped_axis_scratch, m_periodic_ordered_axis_scratch);
      unwrapPeriodicAxis(
          pos_y_comoving, box_lengths.ly, m_tree_source_y_comoving,
          m_periodic_wrapped_axis_scratch, m_periodic_ordered_axis_scratch);
      unwrapPeriodicAxis(
          pos_z_comoving, box_lengths.lz, m_tree_source_z_comoving,
          m_periodic_wrapped_axis_scratch, m_periodic_ordered_axis_scratch);
      tree_source_x = m_tree_source_x_comoving;
      tree_source_y = m_tree_source_y_comoving;
      tree_source_z = m_tree_source_z_comoving;
    }
    if (profile != nullptr) {
      profile->source_preprocess_ms += std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - source_preprocess_start).count();
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
    diagnostics->top_level_domain_leaf_count = m_last_residual_stats.top_level_domain_leaf_count;
    diagnostics->let_candidate_peer_count = m_last_residual_stats.let_candidate_peer_count;
    diagnostics->let_exported_target_count = m_last_residual_stats.let_exported_target_count;
    diagnostics->let_imported_target_count = m_last_residual_stats.let_imported_target_count;
    diagnostics->let_wire_bytes_sent = m_last_residual_stats.let_wire_bytes_sent;
    diagnostics->let_wire_bytes_received = m_last_residual_stats.let_wire_bytes_received;
    diagnostics->let_high_water_bytes = m_last_residual_stats.let_high_water_bytes;
    diagnostics->let_discovery_ms = m_last_residual_stats.let_discovery_ms;
    diagnostics->let_graph_setup_ms = m_last_residual_stats.let_graph_setup_ms;
    diagnostics->let_communication_ms = m_last_residual_stats.let_communication_ms;
    diagnostics->let_overlap_local_work_ms = m_last_residual_stats.let_overlap_local_work_ms;
    diagnostics->let_communication_wait_ms = m_last_residual_stats.let_communication_wait_ms;
    diagnostics->let_overlap_efficiency = m_last_residual_stats.let_overlap_efficiency;
    diagnostics->let_remote_traversal_ms = m_last_residual_stats.let_remote_traversal_ms;
    diagnostics->force_l2_pm_global = pm_force_l2_global;
    diagnostics->force_l2_pm_zoom_correction = l2NormFromComponents(
        m_active_zoom_corr_ax_comoving,
        m_active_zoom_corr_ay_comoving,
        m_active_zoom_corr_az_comoving);
    double total_sq = 0.0;
    for (std::size_t i = 0; i < active_count; ++i) {
      const double total_x = accumulator.accel_x_comoving[i];
      const double total_y = accumulator.accel_y_comoving[i];
      const double total_z = accumulator.accel_z_comoving[i];
      total_sq += total_x * total_x + total_y * total_y + total_z * total_z;
    }
    diagnostics->force_l2_tree_short_range = std::sqrt(
        m_last_residual_stats.local_short_range_sum_sq +
        m_last_residual_stats.remote_short_range_sum_sq);
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
    profile->let_discovery_ms += m_last_residual_stats.let_discovery_ms;
    profile->let_graph_setup_ms += m_last_residual_stats.let_graph_setup_ms;
    profile->let_communication_ms += m_last_residual_stats.let_communication_ms;
    profile->let_overlap_local_work_ms += m_last_residual_stats.let_overlap_local_work_ms;
    profile->let_communication_wait_ms += m_last_residual_stats.let_communication_wait_ms;
    profile->let_overlap_efficiency = m_last_residual_stats.let_overlap_efficiency;
    profile->remote_traversal_ms += m_last_residual_stats.let_remote_traversal_ms;
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

  const bool has_explicit_target_positions = !accumulator.target_pos_x_comoving.empty();
  const auto target_x = [&](std::size_t active_i) -> double {
    return has_explicit_target_positions
        ? accumulator.target_pos_x_comoving[active_i]
        : pos_x_comoving[accumulator.active_particle_index[active_i]];
  };
  const auto target_y = [&](std::size_t active_i) -> double {
    return has_explicit_target_positions
        ? accumulator.target_pos_y_comoving[active_i]
        : pos_y_comoving[accumulator.active_particle_index[active_i]];
  };
  const auto target_z = [&](std::size_t active_i) -> double {
    return has_explicit_target_positions
        ? accumulator.target_pos_z_comoving[active_i]
        : pos_z_comoving[accumulator.active_particle_index[active_i]];
  };

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
  const auto resolve_target_softening = [&](std::size_t active_slot, TreeLocalIndex source_index) {
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
#else
  static_cast<void>(rank_local_serial_mode);
  static_cast<void>(tree_mpi_world_size);
#endif
  std::vector<TreeLocalIndex> stack;
  std::exception_ptr traversal_workspace_failure;
  try {
    // DFS storage scales with traversal frontier/depth, not total node count.
    // Start small and let vector growth reflect actual traversal complexity.
    stack.reserve(std::min<std::size_t>(nodes.size(), 256U));
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
                                         TreeLocalIndex self_index,
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
      const TreeLocalIndex node_index = stack.back();
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
      const bool mac_accept = internal::acceptNodeByMac(
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
      const bool softening_accept = internal::passesSofteningEnvelopeGuard(
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
          const TreeLocalIndex begin = nodes.particle_begin[node_index];
          const TreeLocalIndex end = begin + nodes.particle_count[node_index];
          for (TreeLocalIndex sorted_i = begin; sorted_i < end; ++sorted_i) {
            const TreeLocalIndex source_index = ordering.sorted_particle_index[sorted_i];
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
            const double split_factor = treePmGaussianShortRangeForceFactorUnchecked(sr, options.split_policy.split_scale_comoving);
            const double source_softening =
                resolveSourceSofteningEpsilon(source_index, options.tree_options.softening, softening_view);
            const double pair_epsilon = combineSofteningPairEpsilonUnchecked(source_softening, target_softening_comoving);
            // Contract: short-range residual is the softened tree force multiplied by the
            // Gaussian real-space residual factor so that tree+PM composes to the unsplit
            // softened force before explicit r_cut truncation.
            const double softened_factor =
                softenedInvR3Unchecked(sr2, pair_epsilon) * split_factor * options.tree_options.gravitational_constant_code;
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
      const TreeLocalIndex particle_index = accumulator.active_particle_index[active_i];
      const bool has_local_source_identity = particle_index < pos_x_comoving.size();
      const double px = target_x(active_i);
      const double py = target_y(active_i);
      const double pz = target_z(active_i);
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
        exchange_sequence, options.decomposition_epoch.value, options.force_epoch.sequence, options.tree_exchange_batch_bytes};
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

    // Every rank holds only one compact top-level source-domain leaf per rank.
    // Detailed remote tree nodes are no longer globally replicated. Traversal
    // discovers target-export work against these ownership leaves, and the
    // remote owner evaluates the request against its authoritative local tree.
    std::vector<parallel::TreePseudoParticlePacket> local_top_level_domain;
    std::exception_ptr local_domain_failure;
    try {
      local_top_level_domain.push_back(makeLocalTreePseudoParticlePacket(
          mpi_world_rank,
          options.decomposition_epoch.value,
          options.force_epoch.sequence,
          exchange_sequence,
          m_tree_solver.nodes(),
          pos_x_comoving.size()));
    } catch (...) {
      local_domain_failure = std::current_exception();
    }
    coordinate_protocol_failure(local_domain_failure, "top-level domain preparation");

    // The current parallel decomposition API exposes one compact gravity-derived
    // source envelope per rank rather than authoritative spatial top leaves.
    // These envelopes may change whenever source positions change, so do not key
    // them to local tree-build generation or reuse them across force evaluations.
    // Keeping topology validity separate from numerical tree validity avoids the
    // former guaranteed cache invalidation on every rebuild and leaves a clean
    // upgrade point for authoritative decomposition leaves later.
    std::vector<parallel::TreePseudoParticlePacket> peer_pseudo_packets =
        parallel::executeBlockingTreePseudoParticleHierarchyExchange(
            m_mpi_context, local_top_level_domain, exchange_sequence);
    const std::uint64_t domain_geometry_fingerprint =
        domainGeometryFingerprint(peer_pseudo_packets);
    const bool authoritative_geometry = std::all_of(
        peer_pseudo_packets.begin(),
        peer_pseudo_packets.end(),
        [](const parallel::TreePseudoParticlePacket& packet) {
          return !packet.descriptor.derived_not_authoritative;
        });
    m_let_domain_cache.valid = true;
    m_let_domain_cache.decomposition_epoch = options.decomposition_epoch;
    m_let_domain_cache.geometry_fingerprint = domain_geometry_fingerprint;
    m_let_domain_cache.authoritative_geometry = authoritative_geometry;
    m_let_domain_cache.top_level_domain_leaves = peer_pseudo_packets;

    std::vector<std::vector<parallel::TreePseudoParticlePacket>> peer_hierarchy_by_rank;
    std::vector<std::vector<ShortRangeTargetRequestPacket>> requests_by_peer;
    std::vector<std::vector<std::uint8_t>> response_expected_by_peer;
    std::vector<std::vector<std::uint8_t>> response_seen_by_peer;
    std::vector<std::uint8_t> communicated_with_peer;
    std::vector<int> requested_outgoing_peers;
    TopLevelDomainHierarchy top_level_domain_hierarchy(peer_pseudo_packets);
    SparsePeerGraph* sparse_peer_graph = nullptr;
    std::exception_ptr peer_domain_failure;
    const auto let_discovery_start = std::chrono::steady_clock::now();
    try {
      if (peer_pseudo_packets.size() != static_cast<std::size_t>(mpi_world_size)) {
        throw std::runtime_error("TreePM top-level domain exchange returned incomplete rank coverage");
      }
      peer_hierarchy_by_rank.resize(static_cast<std::size_t>(mpi_world_size));
      for (const parallel::TreePseudoParticlePacket& packet : peer_pseudo_packets) {
        if (packet.descriptor.source_rank < 0 || packet.descriptor.source_rank >= mpi_world_size) {
          throw std::runtime_error("TreePM top-level domain packet has invalid source rank");
        }
        auto& rank_packets = peer_hierarchy_by_rank[static_cast<std::size_t>(packet.descriptor.source_rank)];
        if (!rank_packets.empty()) {
          throw std::runtime_error("TreePM top-level domain exchange produced duplicate owner leaves");
        }
        rank_packets.push_back(packet);
      }
      for (int peer = 0; peer < mpi_world_size; ++peer) {
        if (peer_hierarchy_by_rank[static_cast<std::size_t>(peer)].size() != 1U) {
          throw std::runtime_error("TreePM top-level domain exchange did not produce exactly one leaf per rank");
        }
      }

      // Route target discovery through a compact top-domain BVH rather than
      // scanning every remote rank envelope for every local target. The current
      // leaves are still derived rank envelopes; when src/parallel exposes
      // authoritative decomposition top leaves the same hierarchy can consume
      // them without changing request transport.
      std::vector<std::uint8_t> peer_needed(static_cast<std::size_t>(mpi_world_size), 0U);
      for (std::size_t active_i = 0; active_i < accumulator.active_particle_index.size(); ++active_i) {
        const std::vector<int> owners = top_level_domain_hierarchy.ownersWithinCutoff(
            target_x(active_i),
            target_y(active_i),
            target_z(active_i),
            cutoff_radius_comoving,
            box_lengths,
            mpi_world_rank);
        for (const int peer : owners) {
          peer_needed[static_cast<std::size_t>(peer)] = 1U;
        }
      }
      for (int peer = 0; peer < mpi_world_size; ++peer) {
        if (peer != mpi_world_rank && peer_needed[static_cast<std::size_t>(peer)] != 0U) {
          requested_outgoing_peers.push_back(peer);
        }
      }

      m_last_residual_stats.top_level_domain_leaf_count =
          static_cast<std::uint64_t>(peer_pseudo_packets.size());
      m_last_residual_stats.remote_hierarchy_packets =
          mpi_world_size > 0 ? static_cast<std::uint64_t>(mpi_world_size - 1) : 0U;
      requests_by_peer.resize(static_cast<std::size_t>(mpi_world_size));
      response_expected_by_peer.resize(static_cast<std::size_t>(mpi_world_size));
      response_seen_by_peer.resize(static_cast<std::size_t>(mpi_world_size));
      communicated_with_peer.assign(static_cast<std::size_t>(mpi_world_size), 0U);
    } catch (...) {
      peer_domain_failure = std::current_exception();
    }
    coordinate_protocol_failure(peer_domain_failure, "top-level domain and sparse-peer discovery");
    m_last_residual_stats.let_discovery_ms += std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - let_discovery_start).count();

    // Graph communicators are expensive MPI objects. Reuse the cached graph while
    // decomposition epoch and this rank's adjacency are unchanged. The all-rank
    // agreement makes rebuild/reuse a collective decision, preventing one rank
    // from entering graph construction while a peer reuses an old communicator.
    const auto graph_setup_start = std::chrono::steady_clock::now();
    auto& graph_cache = *m_sparse_peer_graph_cache;
    const bool local_graph_cache_match =
        graph_cache.valid &&
        graph_cache.decomposition_epoch == options.decomposition_epoch &&
        graph_cache.world_size == mpi_world_size &&
        graph_cache.requested_outgoing_peers == requested_outgoing_peers;
    int local_graph_cache_match_int = local_graph_cache_match ? 1 : 0;
    int all_graph_cache_match = 0;
    MPI_Allreduce(
        &local_graph_cache_match_int,
        &all_graph_cache_match,
        1,
        MPI_INT,
        MPI_MIN,
        MPI_COMM_WORLD);
    if (all_graph_cache_match == 0) {
      graph_cache.graph = makeSymmetricSparsePeerGraph(
          mpi_world_rank,
          std::span<const int>(requested_outgoing_peers.data(), requested_outgoing_peers.size()));
      graph_cache.valid = true;
      graph_cache.decomposition_epoch = options.decomposition_epoch;
      graph_cache.world_size = mpi_world_size;
      graph_cache.requested_outgoing_peers = requested_outgoing_peers;
    }
    sparse_peer_graph = &graph_cache.graph;
    m_last_residual_stats.let_candidate_peer_count =
        static_cast<std::uint64_t>(sparse_peer_graph->outgoing_peers.size());
    m_last_residual_stats.let_graph_setup_ms += std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - graph_setup_start).count();

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
        const TreeLocalIndex particle_index = accumulator.active_particle_index[batch_begin + batch_slot];
        const double px = target_x(batch_begin + batch_slot);
        const double py = target_y(batch_begin + batch_slot);
        const double pz = target_z(batch_begin + batch_slot);
        const double target_softening =
            resolve_target_softening(batch_begin + batch_slot, particle_index);
        const bool previous_acceleration_available =
            !accumulator.previous_acceleration_magnitude_code.empty() &&
            std::isfinite(accumulator.previous_acceleration_magnitude_code[batch_begin + batch_slot]);
        const double previous_acceleration_magnitude_code = previous_acceleration_available
            ? accumulator.previous_acceleration_magnitude_code[batch_begin + batch_slot]
            : 0.0;
        const std::vector<int> target_peers = top_level_domain_hierarchy.ownersWithinCutoff(
            px, py, pz, cutoff_radius_comoving, box_lengths, mpi_world_rank);
        if (sparse_peer_graph->outgoing_peers.size() > target_peers.size()) {
          m_last_residual_stats.remote_pairs_pruned_by_bounds +=
              static_cast<std::uint64_t>(
                  sparse_peer_graph->outgoing_peers.size() - target_peers.size());
        }
        for (const int peer : target_peers) {
          requests_by_peer[static_cast<std::size_t>(peer)].push_back(ShortRangeTargetRequestPacket{
              .wire_version = k_short_range_wire_version,
              .origin_rank = mpi_world_rank,
              .destination_rank = peer,
              .flags = previous_acceleration_available ? k_short_range_flag_previous_acceleration : 0U,
              .batch_token = batch_token,
              .request_id = static_cast<std::uint32_t>(batch_slot),
              .exchange_sequence = exchange_sequence,
              .decomposition_epoch = options.decomposition_epoch.value,
              .force_epoch = options.force_epoch.sequence,
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

      // Exchange request counts only across the locality-derived graph. The
      // graph is symmetric so source-only ranks can receive requests from
      // target-only ranks without communicator-wide all-to-all participation.
      std::vector<int> neighbor_request_send_counts(sparse_peer_graph->outgoing_peers.size(), 0);
      std::vector<int> neighbor_request_recv_counts(sparse_peer_graph->incoming_peers.size(), 0);
      for (std::size_t i = 0; i < sparse_peer_graph->outgoing_peers.size(); ++i) {
        neighbor_request_send_counts[i] =
            send_counts[static_cast<std::size_t>(sparse_peer_graph->outgoing_peers[i])];
      }
      const auto request_count_communication_start = std::chrono::steady_clock::now();
      MPI_Neighbor_alltoall(
          neighbor_request_send_counts.empty() ? nullptr : neighbor_request_send_counts.data(),
          1,
          MPI_INT,
          neighbor_request_recv_counts.empty() ? nullptr : neighbor_request_recv_counts.data(),
          1,
          MPI_INT,
          sparse_peer_graph->communicator);
      m_last_residual_stats.let_communication_ms += std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - request_count_communication_start).count();
      for (std::size_t i = 0; i < sparse_peer_graph->incoming_peers.size(); ++i) {
        recv_counts[static_cast<std::size_t>(sparse_peer_graph->incoming_peers[i])] =
            neighbor_request_recv_counts[i];
      }
      for (const int peer : sparse_peer_graph->outgoing_peers) {
        if (send_counts[static_cast<std::size_t>(peer)] > 0 ||
            recv_counts[static_cast<std::size_t>(peer)] > 0) {
          communicated_with_peer[static_cast<std::size_t>(peer)] = 1U;
        }
      }
      int total_recv_bytes = 0;
      std::exception_ptr request_transport_preparation_failure;
      MPI_Request request_payload_exchange = MPI_REQUEST_NULL;
      try {
        total_recv_bytes = populateMpiByteDisplacements(
            recv_counts, recv_displs, "TreePM short-range request receive layout");
        recv_payload.resize(static_cast<std::size_t>(total_recv_bytes), 0U);

        std::vector<int> neighbor_request_send_displs(sparse_peer_graph->outgoing_peers.size(), 0);
        std::vector<int> neighbor_request_recv_displs(sparse_peer_graph->incoming_peers.size(), 0);
        for (std::size_t i = 0; i < sparse_peer_graph->outgoing_peers.size(); ++i) {
          neighbor_request_send_displs[i] =
              send_displs[static_cast<std::size_t>(sparse_peer_graph->outgoing_peers[i])];
        }
        for (std::size_t i = 0; i < sparse_peer_graph->incoming_peers.size(); ++i) {
          neighbor_request_recv_displs[i] =
              recv_displs[static_cast<std::size_t>(sparse_peer_graph->incoming_peers[i])];
        }

        std::uint8_t empty_request_payload = 0U;
        const std::uint8_t* request_send_buffer =
            send_payload.empty() ? &empty_request_payload : send_payload.data();
        std::uint8_t* request_receive_buffer =
            recv_payload.empty() ? &empty_request_payload : recv_payload.data();
        MPI_Ineighbor_alltoallv(
            request_send_buffer,
            neighbor_request_send_counts.empty() ? nullptr : neighbor_request_send_counts.data(),
            neighbor_request_send_displs.empty() ? nullptr : neighbor_request_send_displs.data(),
            MPI_BYTE,
            request_receive_buffer,
            neighbor_request_recv_counts.empty() ? nullptr : neighbor_request_recv_counts.data(),
            neighbor_request_recv_displs.empty() ? nullptr : neighbor_request_recv_displs.data(),
            MPI_BYTE,
            sparse_peer_graph->communicator,
            &request_payload_exchange);

        // Useful local traversal overlaps the sparse remote request transfer.
        const auto overlap_local_work_start = std::chrono::steady_clock::now();
        for (std::size_t batch_slot = 0; batch_slot < batch_size; ++batch_slot) {
        const TreeLocalIndex particle_index = accumulator.active_particle_index[batch_begin + batch_slot];
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
            target_x(batch_begin + batch_slot),
            target_y(batch_begin + batch_slot),
            target_z(batch_begin + batch_slot),
            has_local_source_identity ? particle_index : 0U,
            has_local_source_identity,
            target_softening,
            previous_acceleration_available,
            previous_acceleration_magnitude_code);
        accumulator.addToActiveSlot(batch_begin + batch_slot, local_accel[0], local_accel[1], local_accel[2]);
        m_last_residual_stats.local_short_range_sum_sq +=
            local_accel[0] * local_accel[0] + local_accel[1] * local_accel[1] + local_accel[2] * local_accel[2];
        }
        const double overlap_local_work_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - overlap_local_work_start).count();
        m_last_residual_stats.let_overlap_local_work_ms += overlap_local_work_ms;
        const auto request_wait_start = std::chrono::steady_clock::now();
        MPI_Wait(&request_payload_exchange, MPI_STATUS_IGNORE);
        const double request_wait_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - request_wait_start).count();
        m_last_residual_stats.let_communication_wait_ms += request_wait_ms;
        m_last_residual_stats.let_communication_ms += request_wait_ms;
      } catch (...) {
        if (request_payload_exchange != MPI_REQUEST_NULL) {
          MPI_Wait(&request_payload_exchange, MPI_STATUS_IGNORE);
        }
        request_transport_preparation_failure = std::current_exception();
      }
      coordinate_protocol_failure(
          request_transport_preparation_failure,
          "request receive/local-force preparation");
      const double overlap_denominator_ms =
          m_last_residual_stats.let_overlap_local_work_ms +
          m_last_residual_stats.let_communication_wait_ms;
      m_last_residual_stats.let_overlap_efficiency = overlap_denominator_ms > 0.0
          ? m_last_residual_stats.let_overlap_local_work_ms / overlap_denominator_ms
          : 0.0;

      m_last_residual_stats.remote_request_batches += 1;
      m_last_residual_stats.remote_request_packets +=
          static_cast<std::uint64_t>(total_send_bytes / static_cast<int>(k_short_range_request_wire_bytes));
      m_last_residual_stats.remote_request_bytes += static_cast<std::uint64_t>(total_send_bytes);
      m_last_residual_stats.let_exported_target_count +=
          static_cast<std::uint64_t>(total_send_bytes / static_cast<int>(k_short_range_request_wire_bytes));
      m_last_residual_stats.let_imported_target_count +=
          static_cast<std::uint64_t>(total_recv_bytes / static_cast<int>(k_short_range_request_wire_bytes));
      m_last_residual_stats.let_wire_bytes_sent += static_cast<std::uint64_t>(total_send_bytes);
      m_last_residual_stats.let_wire_bytes_received += static_cast<std::uint64_t>(total_recv_bytes);
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
      const auto remote_traversal_start = std::chrono::steady_clock::now();
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
              request.decomposition_epoch != options.decomposition_epoch.value || request.force_epoch != options.force_epoch.sequence ||
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
      m_last_residual_stats.let_remote_traversal_ms += std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - remote_traversal_start).count();

      std::vector<int> neighbor_response_send_counts(sparse_peer_graph->outgoing_peers.size(), 0);
      std::vector<int> neighbor_response_send_displs(sparse_peer_graph->outgoing_peers.size(), 0);
      std::vector<int> neighbor_response_recv_counts(sparse_peer_graph->incoming_peers.size(), 0);
      std::vector<int> neighbor_response_recv_displs(sparse_peer_graph->incoming_peers.size(), 0);
      for (std::size_t i = 0; i < sparse_peer_graph->outgoing_peers.size(); ++i) {
        const int peer = sparse_peer_graph->outgoing_peers[i];
        neighbor_response_send_counts[i] = response_send_counts[static_cast<std::size_t>(peer)];
        neighbor_response_send_displs[i] = response_send_displs[static_cast<std::size_t>(peer)];
      }
      for (std::size_t i = 0; i < sparse_peer_graph->incoming_peers.size(); ++i) {
        const int peer = sparse_peer_graph->incoming_peers[i];
        neighbor_response_recv_counts[i] = response_recv_counts[static_cast<std::size_t>(peer)];
        neighbor_response_recv_displs[i] = response_recv_displs[static_cast<std::size_t>(peer)];
      }
      std::uint8_t empty_response_payload = 0U;
      const std::uint8_t* response_send_buffer =
          response_send_payload.empty() ? &empty_response_payload : response_send_payload.data();
      std::uint8_t* response_receive_buffer =
          response_recv_payload.empty() ? &empty_response_payload : response_recv_payload.data();
      const auto response_communication_start = std::chrono::steady_clock::now();
      MPI_Neighbor_alltoallv(
          response_send_buffer,
          neighbor_response_send_counts.empty() ? nullptr : neighbor_response_send_counts.data(),
          neighbor_response_send_displs.empty() ? nullptr : neighbor_response_send_displs.data(),
          MPI_BYTE,
          response_receive_buffer,
          neighbor_response_recv_counts.empty() ? nullptr : neighbor_response_recv_counts.data(),
          neighbor_response_recv_displs.empty() ? nullptr : neighbor_response_recv_displs.data(),
          MPI_BYTE,
          sparse_peer_graph->communicator);
      m_last_residual_stats.let_communication_ms += std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - response_communication_start).count();

      std::exception_ptr response_validation_failure;
      try {
      remote_batch_ax.assign(batch_size, 0.0);
      remote_batch_ay.assign(batch_size, 0.0);
      remote_batch_az.assign(batch_size, 0.0);
      m_last_residual_stats.remote_response_bytes += static_cast<std::uint64_t>(response_recv_payload.size());
      m_last_residual_stats.let_wire_bytes_sent += static_cast<std::uint64_t>(total_response_send_bytes);
      m_last_residual_stats.let_wire_bytes_received += static_cast<std::uint64_t>(total_response_recv_bytes);
      const std::uint64_t let_workspace_bytes =
          static_cast<std::uint64_t>(send_payload.capacity()) +
          static_cast<std::uint64_t>(recv_payload.capacity()) +
          static_cast<std::uint64_t>(response_send_payload.capacity()) +
          static_cast<std::uint64_t>(response_recv_payload.capacity());
      m_last_residual_stats.let_high_water_bytes =
          std::max(m_last_residual_stats.let_high_water_bytes, let_workspace_bytes);
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
              response.decomposition_epoch != options.decomposition_epoch.value || response.force_epoch != options.force_epoch.sequence ||
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
    tree_profile->traversed_target_count +=
        static_cast<std::uint64_t>(accumulator.active_particle_index.size());
    if (tree_profile->traversed_target_count > 0U) {
      tree_profile->average_interactions_per_target =
          static_cast<double>(tree_profile->particle_particle_interactions) /
          static_cast<double>(tree_profile->traversed_target_count);
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
      treePmGaussianShortRangeForceFactorUnchecked(split_policy.split_scale_comoving, split_policy.split_scale_comoving);
  diagnostics.long_range_factor_at_split =
      treePmGaussianLongRangeForceFactorUnchecked(split_policy.split_scale_comoving, split_policy.split_scale_comoving);
  diagnostics.short_range_factor_at_cutoff =
      treePmGaussianShortRangeForceFactorUnchecked(split_policy.cutoff_radius_comoving, split_policy.split_scale_comoving);
  diagnostics.long_range_factor_at_cutoff =
      treePmGaussianLongRangeForceFactorUnchecked(split_policy.cutoff_radius_comoving, split_policy.split_scale_comoving);
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
    const double composed = treePmGaussianShortRangeForceFactorUnchecked(radius_comoving, split_policy.split_scale_comoving) +
        treePmGaussianLongRangeForceFactorUnchecked(radius_comoving, split_policy.split_scale_comoving);
    diagnostics.max_relative_composition_error = std::max(
        diagnostics.max_relative_composition_error,
        std::abs(composed - 1.0) / std::max(1.0, std::abs(composed)));
  }

  return diagnostics;
}

}  // namespace cosmosim::gravity
