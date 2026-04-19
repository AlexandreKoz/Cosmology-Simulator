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

[[nodiscard]] bool forceAccumulatorShapeValid(const TreePmForceAccumulatorView& accumulator) {
  return accumulator.active_particle_index.size() == accumulator.accel_x_comoving.size() &&
      accumulator.active_particle_index.size() == accumulator.accel_y_comoving.size() &&
      accumulator.active_particle_index.size() == accumulator.accel_z_comoving.size();
}

void resizeCompactSidecars(std::vector<double>& first, std::vector<double>& second, std::vector<double>& third, std::size_t size) {
  first.resize(size);
  second.resize(size);
  third.resize(size);
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
    double box_size_comoving) {
  const double dx_abs = std::abs(minimumImageDelta(cx - px, box_size_comoving));
  const double dy_abs = std::abs(minimumImageDelta(cy - py, box_size_comoving));
  const double dz_abs = std::abs(minimumImageDelta(cz - pz, box_size_comoving));
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
    double box_size_comoving) {
  const double dx_abs = std::abs(minimumImageDelta(cx - px, box_size_comoving));
  const double dy_abs = std::abs(minimumImageDelta(cy - py, box_size_comoving));
  const double dz_abs = std::abs(minimumImageDelta(cz - pz, box_size_comoving));
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
    TreePmDiagnostics* diagnostics) {
  solveActiveSetWithPmCadence(
      pos_x_comoving,
      pos_y_comoving,
      pos_z_comoving,
      mass_code,
      accumulator,
      options,
      true,
      profile,
      diagnostics);
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
    TreePmDiagnostics* diagnostics) {
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
    m_pm_solver.solvePoissonPeriodic(m_grid, pm_options, profile != nullptr ? &profile->pm_profile : nullptr);
    m_has_cached_long_range_field = true;
  }
  if (!m_has_cached_long_range_field) {
    throw std::runtime_error("TreePM long-range mesh field is unavailable for reuse");
  }

  const std::size_t active_count = accumulator.active_particle_index.size();
  resizeCompactSidecars(m_active_pos_x_comoving, m_active_pos_y_comoving, m_active_pos_z_comoving, active_count);
  resizeCompactSidecars(m_active_pm_ax_comoving, m_active_pm_ay_comoving, m_active_pm_az_comoving, active_count);

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

  // Tree owns short-range residual with the complementary real-space kernel.
  const auto tree_start = std::chrono::steady_clock::now();
  m_tree_solver.build(pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code, options.tree_options,
      profile != nullptr ? &profile->tree_profile : nullptr);
  evaluateShortRangeResidual(
      pos_x_comoving,
      pos_y_comoving,
      pos_z_comoving,
      mass_code,
      accumulator,
      options,
      profile != nullptr ? &profile->tree_profile : nullptr);
  const auto tree_stop = std::chrono::steady_clock::now();

  if (diagnostics != nullptr) {
    *diagnostics = computeTreePmDiagnostics(options.split_policy);
    diagnostics->residual_pruned_nodes = m_last_residual_stats.pruned_nodes;
    diagnostics->residual_pair_skips_cutoff = m_last_residual_stats.pair_skips_cutoff;
    diagnostics->residual_pair_evaluations = m_last_residual_stats.pair_evaluations;
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
    TreeGravityProfile* tree_profile) {
  std::uint64_t visited_nodes = 0;
  std::uint64_t accepted_nodes = 0;
  std::uint64_t pp_interactions = 0;
  std::uint64_t cutoff_pruned_nodes = 0;
  std::uint64_t cutoff_pair_skips = 0;

  const auto traversal_start = std::chrono::steady_clock::now();
  const TreeNodeSoa& nodes = m_tree_solver.nodes();
  const TreeMortonOrdering& ordering = m_tree_solver.ordering();
  const double box_size_comoving = options.pm_options.box_size_mpc_comoving;
  const double cutoff_radius_comoving = options.split_policy.cutoff_radius_comoving;
  const double cutoff_radius2_comoving = cutoff_radius_comoving * cutoff_radius_comoving;
  std::vector<std::uint32_t> stack;
  stack.reserve(256);

  auto evaluateTargetAgainstLocalTree = [&](double px, double py, double pz, std::uint32_t self_index, bool skip_self) {
    double ax = 0.0;
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
          box_size_comoving);
      if (min_node_distance > cutoff_radius_comoving) {
        ++cutoff_pruned_nodes;
        continue;
      }

      const double dx = minimumImageDelta(nodes.com_x_comoving[node_index] - px, box_size_comoving);
      const double dy = minimumImageDelta(nodes.com_y_comoving[node_index] - py, box_size_comoving);
      const double dz = minimumImageDelta(nodes.com_z_comoving[node_index] - pz, box_size_comoving);
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double r = std::sqrt(std::max(r2, 1.0e-30));

      const double l_over_r = (2.0 * half_size) / r;
      const bool is_leaf = nodes.child_count[node_index] == 0;
      const bool geometric_accept = is_leaf || (l_over_r < options.tree_options.opening_theta);
      const bool node_within_cutoff = is_leaf || (maximumDistanceToNodeAabb(
          px,
          py,
          pz,
          nodes.center_x_comoving[node_index],
          nodes.center_y_comoving[node_index],
          nodes.center_z_comoving[node_index],
          half_size,
          box_size_comoving) <= cutoff_radius_comoving);
      const bool accept = geometric_accept && node_within_cutoff;

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
            const double sx = minimumImageDelta(pos_x_comoving[source_index] - px, box_size_comoving);
            const double sy = minimumImageDelta(pos_y_comoving[source_index] - py, box_size_comoving);
            const double sz = minimumImageDelta(pos_z_comoving[source_index] - pz, box_size_comoving);
            const double sr2 = sx * sx + sy * sy + sz * sz;
            if (sr2 > cutoff_radius2_comoving) {
              ++cutoff_pair_skips;
              continue;
            }
            const double sr = std::sqrt(std::max(sr2, 1.0e-30));
            const double split_factor = treePmGaussianShortRangeForceFactor(sr, options.split_policy.split_scale_comoving);
            // Contract: short-range residual is the softened tree force multiplied by the
            // Gaussian real-space residual factor so that tree+PM composes to the unsplit
            // softened force before explicit r_cut truncation.
            const double softened_factor =
                softenedInvR3(sr2, options.tree_options.softening) * split_factor * options.tree_options.gravitational_constant_code;
            ax += softened_factor * mass_code[source_index] * sx;
            ay += softened_factor * mass_code[source_index] * sy;
            az += softened_factor * mass_code[source_index] * sz;
            ++pp_interactions;
          }
        } else {
          const double split_factor = treePmGaussianShortRangeForceFactor(r, options.split_policy.split_scale_comoving);
          // Same softened-residual contract as the leaf pair path, applied to accepted nodes.
          const double softened_factor =
              softenedInvR3(r2, options.tree_options.softening) * split_factor * options.tree_options.gravitational_constant_code;
          ax += softened_factor * nodes.mass_code[node_index] * dx;
          ay += softened_factor * nodes.mass_code[node_index] * dy;
          az += softened_factor * nodes.mass_code[node_index] * dz;
        }
      } else {
        const std::size_t child_offset = static_cast<std::size_t>(node_index) * 8U;
        for (std::uint8_t octant = 0; octant < 8U; ++octant) {
          const std::uint32_t child = nodes.child_index[child_offset + octant];
          if (child != std::numeric_limits<std::uint32_t>::max()) {
            stack.push_back(child);
          }
        }
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
      const auto local_accel = evaluateTargetAgainstLocalTree(px, py, pz, particle_index, true);
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
    auto& seen_response = m_tree_exchange_workspace.seen_response;

    for (std::size_t batch_begin = 0; batch_begin < accumulator.active_particle_index.size();) {
      const std::size_t batch_size = std::min(
          max_requests_per_peer,
          accumulator.active_particle_index.size() - batch_begin);
      const std::uint32_t batch_token = static_cast<std::uint32_t>(batch_begin);

      std::vector<ShortRangeTargetRequestPacket> local_requests;
      local_requests.reserve(batch_size);
      for (std::size_t batch_slot = 0; batch_slot < batch_size; ++batch_slot) {
        const std::uint32_t particle_index = accumulator.active_particle_index[batch_begin + batch_slot];
        local_requests.push_back(ShortRangeTargetRequestPacket{
            .batch_token = batch_token,
            .request_id = static_cast<std::uint32_t>(batch_slot),
            .target_x_comoving = pos_x_comoving[particle_index],
            .target_y_comoving = pos_y_comoving[particle_index],
            .target_z_comoving = pos_z_comoving[particle_index],
        });
      }

      std::fill(send_counts.begin(), send_counts.end(), 0);
      std::fill(recv_counts.begin(), recv_counts.end(), 0);
      std::fill(send_displs.begin(), send_displs.end(), 0);
      std::fill(recv_displs.begin(), recv_displs.end(), 0);
      send_payload.clear();
      for (int rank = 0; rank < mpi_world_size; ++rank) {
        if (rank == mpi_world_rank) {
          send_counts[static_cast<std::size_t>(rank)] = 0;
        } else {
          send_counts[static_cast<std::size_t>(rank)] = static_cast<int>(local_requests.size() * sizeof(ShortRangeTargetRequestPacket));
        }
      }
      std::partial_sum(send_counts.begin(), send_counts.end() - 1, send_displs.begin() + 1);
      send_payload.resize(static_cast<std::size_t>(std::accumulate(send_counts.begin(), send_counts.end(), 0)), 0U);
      for (int rank = 0; rank < mpi_world_size; ++rank) {
        if (rank == mpi_world_rank) {
          continue;
        }
        const auto request_bytes = asBytes(std::span<const ShortRangeTargetRequestPacket>(local_requests.data(), local_requests.size()));
        std::copy(
            request_bytes.begin(),
            request_bytes.end(),
            send_payload.begin() + send_displs[static_cast<std::size_t>(rank)]);
      }

      MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
      std::partial_sum(recv_counts.begin(), recv_counts.end() - 1, recv_displs.begin() + 1);
      recv_payload.resize(static_cast<std::size_t>(std::accumulate(recv_counts.begin(), recv_counts.end(), 0)), 0U);
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

      for (std::size_t batch_slot = 0; batch_slot < batch_size; ++batch_slot) {
        const std::uint32_t particle_index = accumulator.active_particle_index[batch_begin + batch_slot];
        const auto local_accel = evaluateTargetAgainstLocalTree(
            pos_x_comoving[particle_index], pos_y_comoving[particle_index], pos_z_comoving[particle_index], particle_index, true);
        accumulator.addToActiveSlot(batch_begin + batch_slot, local_accel[0], local_accel[1], local_accel[2]);
      }
      MPI_Wait(&request_exchange, MPI_STATUS_IGNORE);

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
              false);
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
      for (int peer = 0; peer < mpi_world_size; ++peer) {
        const int expected_bytes = send_counts[static_cast<std::size_t>(peer)];
        if (peer == mpi_world_rank) {
          continue;
        }
        if (expected_bytes == 0) {
          continue;
        }
        const std::span<const std::uint8_t> peer_response_bytes(
            response_recv_payload.data() + send_displs[static_cast<std::size_t>(peer)],
            static_cast<std::size_t>(expected_bytes));
        const std::vector<ShortRangeTargetResponsePacket> responses =
            decodeRecords<ShortRangeTargetResponsePacket>(peer_response_bytes);
        if (responses.size() != batch_size) {
          throw std::runtime_error("TreePM short-range response count mismatch");
        }
        seen_response.assign(batch_size, 0U);
        for (const ShortRangeTargetResponsePacket& response : responses) {
          if (response.batch_token != batch_token || response.request_id >= batch_size) {
            throw std::runtime_error("TreePM short-range response packet ID mismatch");
          }
          if (seen_response[response.request_id] != 0U) {
            throw std::runtime_error("TreePM short-range duplicate response packet detected");
          }
          seen_response[response.request_id] = 1U;
          remote_batch_ax[response.request_id] += response.accel_x_comoving;
          remote_batch_ay[response.request_id] += response.accel_y_comoving;
          remote_batch_az[response.request_id] += response.accel_z_comoving;
        }
        if (std::any_of(seen_response.begin(), seen_response.end(), [](std::uint8_t mark) { return mark == 0U; })) {
          throw std::runtime_error("TreePM short-range missing response packet detected");
        }
      }

      for (std::size_t batch_slot = 0; batch_slot < batch_size; ++batch_slot) {
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
