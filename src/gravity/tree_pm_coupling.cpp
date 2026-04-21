#include "cosmosim/gravity/tree_pm_coupling.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::gravity {
namespace {

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

[[nodiscard]] double norm3(double x, double y, double z) {
  return std::sqrt(x * x + y * y + z * z);
}

[[nodiscard]] bool forceAccumulatorShapeValid(const TreePmForceAccumulatorView& accumulator) {
  return accumulator.active_particle_index.size() == accumulator.accel_x_comoving.size() &&
      accumulator.active_particle_index.size() == accumulator.accel_y_comoving.size() &&
      accumulator.active_particle_index.size() == accumulator.accel_z_comoving.size();
}

[[nodiscard]] bool acceptNodeByMac(
    TreeOpeningCriterion criterion,
    bool is_leaf,
    double theta,
    double half_size,
    double com_center_offset,
    double r2) {
  if (is_leaf) {
    return true;
  }
  const double r = std::sqrt(r2 + 1.0e-30);
  const double width = 2.0 * half_size;
  if (criterion == TreeOpeningCriterion::kBarnesHutGeometric) {
    return (width / r) < theta;
  }
  const double effective_size = width + com_center_offset;
  return (effective_size / r) < theta;
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
    std::span<const double> m) {
  if (x.size() != y.size() || x.size() != z.size() || x.size() != m.size()) {
    throw std::invalid_argument("allGatherParticleField requires equal local span lengths");
  }
  int world_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  const int local_count = static_cast<int>(x.size());
  std::vector<int> counts(static_cast<std::size_t>(world_size), 0);
  MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  std::vector<int> displs(static_cast<std::size_t>(world_size), 0);
  int total = 0;
  for (int i = 0; i < world_size; ++i) {
    displs[static_cast<std::size_t>(i)] = total;
    total += counts[static_cast<std::size_t>(i)];
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
    double split_factor) {
  const double r2 = dx * dx + dy * dy + dz * dz;
  const double pair_epsilon =
      combineSofteningPairEpsilon(nodes.softening_max_comoving[node_index], target_softening_comoving);
  const double eps2 = pair_epsilon * pair_epsilon;
  const double denom = std::max(r2 + eps2, 1.0e-30);
  const double inv_r3 = 1.0 / (denom * std::sqrt(denom));
  const double inv_r5 = inv_r3 / denom;
  const double inv_r7 = inv_r5 / denom;
  const double prefactor = options.gravitational_constant_code * split_factor;

  double ax = prefactor * nodes.mass_code[node_index] * inv_r3 * dx;
  double ay = prefactor * nodes.mass_code[node_index] * inv_r3 * dy;
  double az = prefactor * nodes.mass_code[node_index] * inv_r3 * dz;
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
  const double quad_prefactor = 0.5 * prefactor;
  ax += quad_prefactor * (2.0 * qrx * inv_r5 - 5.0 * rqr * dx * inv_r7);
  ay += quad_prefactor * (2.0 * qry * inv_r5 - 5.0 * rqr * dy * inv_r7);
  az += quad_prefactor * (2.0 * qrz * inv_r5 - 5.0 * rqr * dz * inv_r7);
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

template <typename T>
[[nodiscard]] std::vector<std::uint8_t> asBytes(std::span<const T> records) {
  static_assert(std::is_trivially_copyable_v<T>);
  std::vector<std::uint8_t> bytes(records.size_bytes(), 0U);
  if (!records.empty()) {
    std::memcpy(bytes.data(), records.data(), records.size_bytes());
  }
  return bytes;
}

template <typename T>
[[nodiscard]] std::vector<T> decodeRecords(std::span<const std::uint8_t> bytes) {
  static_assert(std::is_trivially_copyable_v<T>);
  if (bytes.size() % sizeof(T) != 0) {
    throw std::runtime_error("TreePM short-range exchange payload is misaligned");
  }
  std::vector<T> records(bytes.size() / sizeof(T));
  if (!records.empty()) {
    std::memcpy(records.data(), bytes.data(), bytes.size());
  }
  return records;
}

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
  std::uint32_t batch_token = 0;
  std::uint32_t request_id = 0;
  double target_x_comoving = 0.0;
  double target_y_comoving = 0.0;
  double target_z_comoving = 0.0;
  double target_softening_epsilon_comoving = 0.0;
};

struct ShortRangeTargetResponsePacket {
  std::uint32_t batch_token = 0;
  std::uint32_t request_id = 0;
  double accel_x_comoving = 0.0;
  double accel_y_comoving = 0.0;
  double accel_z_comoving = 0.0;
};

static_assert(std::is_trivially_copyable_v<ShortRangeTargetRequestPacket>);
static_assert(std::is_trivially_copyable_v<ShortRangeTargetResponsePacket>);

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
  double best = std::numeric_limits<double>::max();
  for (const double shift : {-box_size_comoving, 0.0, box_size_comoving}) {
    best = std::min(best, interval_distance(coordinate, interval_min + shift, interval_max + shift));
  }
  return best;
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
    : m_shape(pm_shape), m_grid(pm_shape), m_pm_solver(pm_shape), m_tree_solver() {}

TreePmCoordinator::TreePmCoordinator(PmGridShape pm_shape, parallel::PmSlabLayout pm_layout)
    : m_shape(pm_shape), m_grid(pm_shape, std::move(pm_layout)), m_pm_solver(pm_shape), m_tree_solver() {}

const parallel::PmSlabLayout& TreePmCoordinator::slabLayout() const noexcept {
  return m_grid.slabLayout();
}

bool TreePmCoordinator::ownsFullPmDomain() const noexcept {
  return m_grid.ownsFullDomain();
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
  if (!forceAccumulatorShapeValid(accumulator)) {
    throw std::invalid_argument("TreePM force accumulator spans must have matching active-set extent");
  }
  if (pos_x_comoving.size() != pos_y_comoving.size() || pos_x_comoving.size() != pos_z_comoving.size() ||
      pos_x_comoving.size() != mass_code.size()) {
    throw std::invalid_argument("TreePM requires position and mass spans with equal extent");
  }
  validateTreePmSplitPolicy(options.split_policy);
  if (options.tree_exchange_batch_bytes == 0) {
    throw std::invalid_argument("TreePM tree_exchange_batch_bytes must be > 0");
  }
  if (options.pm_options.boundary_condition == PmBoundaryCondition::kIsolatedOpen &&
      m_grid.slabLayout().world_size > 1) {
    throw std::invalid_argument("TreePM isolated PM path currently supports only single-rank execution");
  }

  const auto start = std::chrono::steady_clock::now();
  accumulator.reset();

  // PM owns long-range force via explicit Gaussian Fourier filter. Cadence-aware callers
  // may choose to reuse the previously solved PM mesh field.
  PmSolveOptions pm_options = options.pm_options;
  pm_options.tree_pm_split_scale_comoving = options.split_policy.split_scale_comoving;
  if (refresh_long_range_field || !m_has_cached_long_range_field) {
    m_pm_solver.assignDensity(
        m_grid,
        pos_x_comoving,
        pos_y_comoving,
        pos_z_comoving,
        mass_code,
        pm_options,
        profile != nullptr ? &profile->pm_profile : nullptr);
    if (pm_options.boundary_condition == PmBoundaryCondition::kPeriodic) {
      m_pm_solver.solvePoissonPeriodic(m_grid, pm_options, profile != nullptr ? &profile->pm_profile : nullptr);
    } else {
      m_pm_solver.solvePoissonIsolatedOpen(m_grid, pm_options, profile != nullptr ? &profile->pm_profile : nullptr);
    }
    m_has_cached_long_range_field = true;
  }
  if (!m_has_cached_long_range_field) {
    throw std::runtime_error("TreePM long-range mesh field is unavailable for reuse");
  }

  const std::size_t active_count = accumulator.active_particle_index.size();
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
    if (particle_index >= pos_x_comoving.size()) {
      throw std::out_of_range("TreePM active index exceeds particle count");
    }
    m_active_pos_x_comoving[i] = pos_x_comoving[particle_index];
    m_active_pos_y_comoving[i] = pos_y_comoving[particle_index];
    m_active_pos_z_comoving[i] = pos_z_comoving[particle_index];
  }

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
      const auto gathered = allGatherParticleField(high_res_source_x, high_res_source_y, high_res_source_z, high_res_source_mass);
      high_res_source_x = std::move(gathered.x);
      high_res_source_y = std::move(gathered.y);
      high_res_source_z = std::move(gathered.z);
      high_res_source_mass = std::move(gathered.m);
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
  m_tree_solver.build(
      pos_x_comoving,
      pos_y_comoving,
      pos_z_comoving,
      mass_code,
      options.tree_options,
      profile != nullptr ? &profile->tree_profile : nullptr,
      softening_view);
  evaluateShortRangeResidual(
      pos_x_comoving,
      pos_y_comoving,
      pos_z_comoving,
      mass_code,
      accumulator,
      options,
      softening_view,
      profile != nullptr ? &profile->tree_profile : nullptr);
  const auto tree_stop = std::chrono::steady_clock::now();

  if (diagnostics != nullptr) {
    const PeriodicBoxLengths diagnostic_box_lengths = effectivePeriodicBoxLengths(options.pm_options);
    *diagnostics = computeTreePmDiagnostics(options.split_policy);
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
    diagnostics->force_l2_total = std::sqrt(total_sq);
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
    TreeGravityProfile* tree_profile) {
  std::uint64_t visited_nodes = 0;
  std::uint64_t accepted_nodes = 0;
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
  if (!softening_view.source_species_tag.empty() && softening_view.source_species_tag.size() != pos_x_comoving.size()) {
    throw std::invalid_argument("TreePM source species sidecar size must match source particle count");
  }
  if (!softening_view.target_particle_epsilon_comoving.empty() &&
      softening_view.target_particle_epsilon_comoving.size() != accumulator.active_particle_index.size()) {
    throw std::invalid_argument("TreePM target softening sidecar size must match active-set size");
  }
  std::vector<std::uint32_t> stack;
  stack.reserve(256);

  auto evaluateTargetAgainstLocalTree = [&](
                                         double px,
                                         double py,
                                         double pz,
                                         std::uint32_t self_index,
                                         bool skip_self,
                                         double target_softening_comoving) {
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
      const bool mac_accept = acceptNodeByMac(
          options.tree_options.opening_criterion,
          is_leaf,
          options.tree_options.opening_theta,
          half_size,
          com_offset,
          r2);
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
      const bool accept = mac_accept && node_within_cutoff && softening_accept;

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
          const double split_factor = treePmGaussianShortRangeForceFactor(r, options.split_policy.split_scale_comoving);
          // Same softened-residual contract as the leaf pair path, applied to accepted nodes.
          const auto contrib =
              monopolePlusQuadrupoleAccelPeriodic(
                  nodes, node_index, dx, dy, dz, target_softening_comoving, options.tree_options, split_factor);
          ax += contrib[0];
          ay += contrib[1];
          az += contrib[2];
        }
      } else {
        pushChildrenNearFirstPeriodic(nodes, node_index, px, py, pz, box_lengths, stack);
      }
    }
    return std::array<double, 3>{ax, ay, az};
  };

  bool distributed_short_range = false;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  distributed_short_range = m_grid.slabLayout().world_size > 1;
#endif
  if (!distributed_short_range) {
    for (std::size_t active_i = 0; active_i < accumulator.active_particle_index.size(); ++active_i) {
      const std::uint32_t particle_index = accumulator.active_particle_index[active_i];
      const double px = pos_x_comoving[particle_index];
      const double py = pos_y_comoving[particle_index];
      const double pz = pos_z_comoving[particle_index];
      const double target_softening = resolveTargetSofteningEpsilon(active_i, options.tree_options.softening, softening_view);
      const auto local_accel = evaluateTargetAgainstLocalTree(px, py, pz, particle_index, true, target_softening);
      accumulator.addToActiveSlot(active_i, local_accel[0], local_accel[1], local_accel[2]);
    }
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const parallel::PmSlabLayout& layout = m_grid.slabLayout();
  if (layout.world_size > 1) {
    int mpi_world_size = 1;
    int mpi_world_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_world_rank);
    if (mpi_world_size != layout.world_size || mpi_world_rank != layout.world_rank) {
      throw std::invalid_argument("TreePM short-range exchange layout world metadata must match MPI_COMM_WORLD");
    }

    const std::size_t max_requests_per_peer = std::max<std::size_t>(
        1U,
        static_cast<std::size_t>(options.tree_exchange_batch_bytes / sizeof(ShortRangeTargetRequestPacket)));

    if (m_tree_exchange_workspace.world_size != mpi_world_size) {
      m_tree_exchange_workspace.world_size = mpi_world_size;
      m_tree_exchange_workspace.send_counts.assign(static_cast<std::size_t>(mpi_world_size), 0);
      m_tree_exchange_workspace.recv_counts.assign(static_cast<std::size_t>(mpi_world_size), 0);
      m_tree_exchange_workspace.send_displs.assign(static_cast<std::size_t>(mpi_world_size), 0);
      m_tree_exchange_workspace.recv_displs.assign(static_cast<std::size_t>(mpi_world_size), 0);
    }
    auto& send_counts = m_tree_exchange_workspace.send_counts;
    auto& recv_counts = m_tree_exchange_workspace.recv_counts;
    auto& send_displs = m_tree_exchange_workspace.send_displs;
    auto& recv_displs = m_tree_exchange_workspace.recv_displs;
    auto& send_payload = m_tree_exchange_workspace.send_payload;
    auto& recv_payload = m_tree_exchange_workspace.recv_payload;
    auto& response_send_payload = m_tree_exchange_workspace.response_send_payload;
    auto& response_recv_payload = m_tree_exchange_workspace.response_recv_payload;
    auto& remote_batch_ax = m_tree_exchange_workspace.remote_batch_ax;
    auto& remote_batch_ay = m_tree_exchange_workspace.remote_batch_ay;
    auto& remote_batch_az = m_tree_exchange_workspace.remote_batch_az;
    auto& expected_response_count = m_tree_exchange_workspace.expected_response_count;
    auto& received_response_count = m_tree_exchange_workspace.received_response_count;

    SourceDomainBoundsPacket local_bounds =
        computeLocalSourceBounds(pos_x_comoving, pos_y_comoving, pos_z_comoving, box_lengths);
    std::vector<SourceDomainBoundsPacket> peer_bounds(static_cast<std::size_t>(mpi_world_size));
    MPI_Allgather(
        &local_bounds,
        static_cast<int>(sizeof(SourceDomainBoundsPacket)),
        MPI_BYTE,
        peer_bounds.data(),
        static_cast<int>(sizeof(SourceDomainBoundsPacket)),
        MPI_BYTE,
        MPI_COMM_WORLD);

    std::vector<std::vector<ShortRangeTargetRequestPacket>> requests_by_peer(static_cast<std::size_t>(mpi_world_size));

    for (std::size_t batch_begin = 0; batch_begin < accumulator.active_particle_index.size();) {
      const std::size_t batch_size = std::min(
          max_requests_per_peer,
          accumulator.active_particle_index.size() - batch_begin);
      const std::uint32_t batch_token = static_cast<std::uint32_t>(batch_begin);

      for (auto& peer_requests : requests_by_peer) {
        peer_requests.clear();
      }
      for (int peer = 0; peer < mpi_world_size; ++peer) {
        if (peer != mpi_world_rank) {
          requests_by_peer[static_cast<std::size_t>(peer)].reserve(batch_size);
        }
      }

      expected_response_count.assign(batch_size, 0U);
      received_response_count.assign(batch_size, 0U);
      for (std::size_t batch_slot = 0; batch_slot < batch_size; ++batch_slot) {
        const std::uint32_t particle_index = accumulator.active_particle_index[batch_begin + batch_slot];
        const double px = pos_x_comoving[particle_index];
        const double py = pos_y_comoving[particle_index];
        const double pz = pos_z_comoving[particle_index];
        const double target_softening =
            resolveTargetSofteningEpsilon(batch_begin + batch_slot, options.tree_options.softening, softening_view);
        for (int peer = 0; peer < mpi_world_size; ++peer) {
          if (peer == mpi_world_rank) {
            continue;
          }
          if (minimumDistanceToPeriodicBounds(px, py, pz, peer_bounds[static_cast<std::size_t>(peer)], box_lengths) >
              cutoff_radius_comoving) {
            ++m_last_residual_stats.remote_pairs_pruned_by_bounds;
            continue;
          }
          requests_by_peer[static_cast<std::size_t>(peer)].push_back(ShortRangeTargetRequestPacket{
              .batch_token = batch_token,
              .request_id = static_cast<std::uint32_t>(batch_slot),
              .target_x_comoving = px,
              .target_y_comoving = py,
              .target_z_comoving = pz,
              .target_softening_epsilon_comoving = target_softening,
          });
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
        send_counts[static_cast<std::size_t>(peer)] = static_cast<int>(
            requests_by_peer[static_cast<std::size_t>(peer)].size() * sizeof(ShortRangeTargetRequestPacket));
      }
      std::partial_sum(send_counts.begin(), send_counts.end() - 1, send_displs.begin() + 1);
      const int total_send_bytes = std::accumulate(send_counts.begin(), send_counts.end(), 0);
      send_payload.resize(static_cast<std::size_t>(total_send_bytes), 0U);
      for (int peer = 0; peer < mpi_world_size; ++peer) {
        if (peer == mpi_world_rank || requests_by_peer[static_cast<std::size_t>(peer)].empty()) {
          continue;
        }
        const auto request_bytes = asBytes(std::span<const ShortRangeTargetRequestPacket>(
            requests_by_peer[static_cast<std::size_t>(peer)].data(),
            requests_by_peer[static_cast<std::size_t>(peer)].size()));
        std::copy(
            request_bytes.begin(),
            request_bytes.end(),
            send_payload.begin() + send_displs[static_cast<std::size_t>(peer)]);
      }

      MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
      std::partial_sum(recv_counts.begin(), recv_counts.end() - 1, recv_displs.begin() + 1);
      const int total_recv_bytes = std::accumulate(recv_counts.begin(), recv_counts.end(), 0);
      recv_payload.resize(static_cast<std::size_t>(total_recv_bytes), 0U);

      for (std::size_t batch_slot = 0; batch_slot < batch_size; ++batch_slot) {
        const std::uint32_t particle_index = accumulator.active_particle_index[batch_begin + batch_slot];
        const double target_softening =
            resolveTargetSofteningEpsilon(batch_begin + batch_slot, options.tree_options.softening, softening_view);
        const auto local_accel = evaluateTargetAgainstLocalTree(
            pos_x_comoving[particle_index],
            pos_y_comoving[particle_index],
            pos_z_comoving[particle_index],
            particle_index,
            true,
            target_softening);
        accumulator.addToActiveSlot(batch_begin + batch_slot, local_accel[0], local_accel[1], local_accel[2]);
      }

      if (total_send_bytes == 0 && total_recv_bytes == 0) {
        batch_begin += batch_size;
        continue;
      }

      MPI_Request request_exchange = MPI_REQUEST_NULL;
      MPI_Ialltoallv(
          send_payload.data(),
          send_counts.data(),
          send_displs.data(),
          MPI_BYTE,
          recv_payload.data(),
          recv_counts.data(),
          recv_displs.data(),
          MPI_BYTE,
          MPI_COMM_WORLD,
          &request_exchange);
      MPI_Wait(&request_exchange, MPI_STATUS_IGNORE);

      m_last_residual_stats.remote_request_batches += 1;
      m_last_residual_stats.remote_request_packets += static_cast<std::uint64_t>(total_send_bytes / sizeof(ShortRangeTargetRequestPacket));
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
              send_counts[static_cast<std::size_t>(peer)] / static_cast<int>(sizeof(ShortRangeTargetRequestPacket)));
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

      response_send_payload.assign(recv_payload.size(), 0U);
      for (int peer = 0; peer < mpi_world_size; ++peer) {
        const int peer_recv_bytes = recv_counts[static_cast<std::size_t>(peer)];
        if (peer_recv_bytes == 0) {
          continue;
        }
        const std::span<const std::uint8_t> peer_request_bytes(
            recv_payload.data() + recv_displs[static_cast<std::size_t>(peer)],
            static_cast<std::size_t>(peer_recv_bytes));
        const std::vector<ShortRangeTargetRequestPacket> peer_requests =
            decodeRecords<ShortRangeTargetRequestPacket>(peer_request_bytes);

        std::vector<ShortRangeTargetResponsePacket> peer_responses;
        peer_responses.reserve(peer_requests.size());
        for (const ShortRangeTargetRequestPacket& request : peer_requests) {
          const auto remote_accel = evaluateTargetAgainstLocalTree(
              request.target_x_comoving,
              request.target_y_comoving,
              request.target_z_comoving,
              0U,
              false,
              request.target_softening_epsilon_comoving);
          peer_responses.push_back(ShortRangeTargetResponsePacket{
              .batch_token = request.batch_token,
              .request_id = request.request_id,
              .accel_x_comoving = remote_accel[0],
              .accel_y_comoving = remote_accel[1],
              .accel_z_comoving = remote_accel[2],
          });
        }
        const auto encoded_responses =
            asBytes(std::span<const ShortRangeTargetResponsePacket>(peer_responses.data(), peer_responses.size()));
        if (encoded_responses.size() != static_cast<std::size_t>(peer_recv_bytes)) {
          throw std::runtime_error("TreePM short-range response payload size mismatch");
        }
        std::copy(
            encoded_responses.begin(),
            encoded_responses.end(),
            response_send_payload.begin() + recv_displs[static_cast<std::size_t>(peer)]);
      }

      response_recv_payload.resize(send_payload.size(), 0U);
      MPI_Alltoallv(
          response_send_payload.data(),
          recv_counts.data(),
          recv_displs.data(),
          MPI_BYTE,
          response_recv_payload.data(),
          send_counts.data(),
          send_displs.data(),
          MPI_BYTE,
          MPI_COMM_WORLD);

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
              send_counts[static_cast<std::size_t>(peer)] / static_cast<int>(sizeof(ShortRangeTargetResponsePacket)));
          response_max = std::max(response_max, packets);
        }
        m_last_residual_stats.remote_response_packets_max_peer =
            std::max(m_last_residual_stats.remote_response_packets_max_peer, response_max);
      }
      for (int peer = 0; peer < mpi_world_size; ++peer) {
        if (peer == mpi_world_rank) {
          continue;
        }
        const int expected_bytes = send_counts[static_cast<std::size_t>(peer)];
        if (expected_bytes == 0) {
          continue;
        }
        const std::span<const std::uint8_t> peer_response_bytes(
            response_recv_payload.data() + send_displs[static_cast<std::size_t>(peer)],
            static_cast<std::size_t>(expected_bytes));
        const std::vector<ShortRangeTargetResponsePacket> responses =
            decodeRecords<ShortRangeTargetResponsePacket>(peer_response_bytes);
        for (const ShortRangeTargetResponsePacket& response : responses) {
          if (response.batch_token != batch_token || response.request_id >= batch_size) {
            throw std::runtime_error("TreePM short-range response packet ID mismatch");
          }
          remote_batch_ax[response.request_id] += response.accel_x_comoving;
          remote_batch_ay[response.request_id] += response.accel_y_comoving;
          remote_batch_az[response.request_id] += response.accel_z_comoving;
          ++received_response_count[response.request_id];
          ++m_last_residual_stats.remote_response_packets;
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
      }
      batch_begin += batch_size;
    }
  }
#else
  if (m_grid.slabLayout().world_size > 1) {
    throw std::invalid_argument("TreePM short-range distributed exchange requires COSMOSIM_ENABLE_MPI=ON");
  }
#endif

  if (tree_profile != nullptr) {
    tree_profile->visited_nodes += visited_nodes;
    tree_profile->accepted_nodes += accepted_nodes;
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
