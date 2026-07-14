#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/cosmosim.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

[[nodiscard]] std::uint64_t xorRangeOneToN(std::uint64_t n) {
  switch (n & 3ULL) {
    case 0ULL: return n;
    case 1ULL: return 1ULL;
    case 2ULL: return n + 1ULL;
    default: return 0ULL;
  }
}

struct VectorErrorNorm {
  double diff_l2 = 0.0;
  double ref_l2 = 0.0;
  double max_rel = 0.0;
};

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

[[nodiscard]] bool hasFastSpectralPmBackend() {
  return cosmosim::gravity::PmSolver::fftBackendAvailable();
}

void buildDeterministicParticles(
    std::size_t particle_count,
    double box_size,
    std::vector<double>& pos_x,
    std::vector<double>& pos_y,
    std::vector<double>& pos_z,
    std::vector<double>& mass) {
  pos_x.resize(particle_count);
  pos_y.resize(particle_count);
  pos_z.resize(particle_count);
  mass.resize(particle_count);
  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod((19.0 * static_cast<double>(i) + 3.0) * 0.017, box_size);
    pos_y[i] = std::fmod((29.0 * static_cast<double>(i) + 5.0) * 0.013, box_size);
    pos_z[i] = std::fmod((37.0 * static_cast<double>(i) + 7.0) * 0.011, box_size);
    mass[i] = 0.8 + 0.04 * static_cast<double>(i % 11U);
  }
}

void buildClusteredParticles(
    std::size_t particle_count,
    double box_size,
    std::vector<double>& pos_x,
    std::vector<double>& pos_y,
    std::vector<double>& pos_z,
    std::vector<double>& mass) {
  pos_x.resize(particle_count);
  pos_y.resize(particle_count);
  pos_z.resize(particle_count);
  mass.resize(particle_count);
  for (std::size_t i = 0; i < particle_count; ++i) {
    const bool first_cluster = i < (particle_count / 2U);
    const double cx = first_cluster ? 0.08 : 0.92;
    const double cy = first_cluster ? 0.18 : 0.81;
    const double cz = first_cluster ? 0.26 : 0.74;
    const double phase = static_cast<double>(i + 1U);
    const auto wrap = [box_size](double value) {
      double out = std::fmod(value, box_size);
      if (out < 0.0) {
        out += box_size;
      }
      return out;
    };
    pos_x[i] = wrap(cx + 0.03 * std::sin(0.37 * phase));
    pos_y[i] = wrap(cy + 0.028 * std::cos(0.31 * phase));
    pos_z[i] = wrap(cz + 0.024 * std::sin(0.29 * phase + (first_cluster ? 0.3 : 1.1)));
    mass[i] = first_cluster ? 1.1 : 0.75;
  }
}

void buildPeriodicSeamCluster(
    std::size_t particle_count,
    double box_size,
    std::vector<double>& pos_x,
    std::vector<double>& pos_y,
    std::vector<double>& pos_z,
    std::vector<double>& mass) {
  pos_x.resize(particle_count);
  pos_y.resize(particle_count);
  pos_z.resize(particle_count);
  mass.resize(particle_count);
  const auto wrap = [box_size](double value) {
    double out = std::fmod(value, box_size);
    if (out < 0.0) {
      out += box_size;
    }
    return out;
  };
  for (std::size_t i = 0; i < particle_count; ++i) {
    const double phase = static_cast<double>(i + 1U);
    const double seam_x = (i & 1U) == 0U ? 0.992 : 0.008;
    const double seam_y = (i % 3U) == 0U ? 0.988 : 0.012;
    const double seam_z = (i % 5U) < 2U ? 0.985 : 0.015;
    pos_x[i] = wrap(seam_x + 0.0041 * std::sin(0.73 * phase));
    pos_y[i] = wrap(seam_y + 0.0037 * std::cos(0.61 * phase + 0.2));
    pos_z[i] = wrap(seam_z + 0.0033 * std::sin(0.47 * phase + 0.7));
    mass[i] = 0.65 + 0.035 * static_cast<double>((7U * i + 3U) % 13U);
  }
}

enum class ParticlePattern {
  kUniform,
  kClustered,
  kPeriodicSeamCluster,
};

struct TreePmEquivalenceCase {
  std::string_view name;
  ParticlePattern particle_pattern = ParticlePattern::kUniform;
  bool disjoint_uneven_ownership = false;
  cosmosim::gravity::PmAssignmentScheme assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc;
  cosmosim::gravity::TreeMultipoleOrder multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole;
  cosmosim::gravity::TreeOpeningCriterion opening_criterion =
      cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance;
  std::size_t fast_particle_count = 0;
  std::size_t fallback_particle_count = 0;
  int solve_count = 1;
  std::uint64_t exchange_batch_bytes = 4ULL * 1024ULL * 1024ULL;
  bool require_multiple_batches = false;
  bool require_cutoff_evidence = false;
};

struct ExchangeMetricTotals {
  std::uint64_t request_packets = 0;
  std::uint64_t response_packets = 0;
  std::uint64_t request_bytes = 0;
  std::uint64_t response_bytes = 0;
  std::uint64_t batches = 0;
  std::uint64_t peer_participations = 0;
  std::uint64_t pruned_nodes = 0;
  std::uint64_t cutoff_pair_skips = 0;
  std::uint64_t remote_pairs_pruned_by_bounds = 0;
};

void accumulateExchangeMetrics(
    const cosmosim::gravity::TreePmDiagnostics& diagnostics,
    ExchangeMetricTotals& totals) {
  totals.request_packets += diagnostics.residual_remote_request_packets;
  totals.response_packets += diagnostics.residual_remote_response_packets;
  totals.request_bytes += diagnostics.residual_remote_request_bytes;
  totals.response_bytes += diagnostics.residual_remote_response_bytes;
  totals.batches += diagnostics.residual_remote_request_batches;
  totals.peer_participations += diagnostics.residual_remote_peer_participations;
  totals.pruned_nodes += diagnostics.residual_pruned_nodes;
  totals.cutoff_pair_skips += diagnostics.residual_pair_skips_cutoff;
  totals.remote_pairs_pruned_by_bounds += diagnostics.residual_remote_pairs_pruned_by_bounds;
}

[[nodiscard]] std::uint64_t globalSumU64(std::uint64_t local_value) {
#if COSMOSIM_ENABLE_MPI
  std::uint64_t global_value = 0;
  MPI_Allreduce(&local_value, &global_value, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  return global_value;
#else
  return local_value;
#endif
}

[[nodiscard]] std::vector<std::uint32_t> ownedParticleIndices(std::size_t count, int world_size, int world_rank) {
  std::vector<std::uint32_t> owned;
  for (std::size_t i = 0; i < count; ++i) {
    if (static_cast<int>(i % static_cast<std::size_t>(world_size)) == world_rank) {
      owned.push_back(static_cast<std::uint32_t>(i));
    }
  }
  return owned;
}

[[nodiscard]] VectorErrorNorm accumulateVectorError(
    std::span<const std::uint32_t> owned,
    std::span<const double> got_x,
    std::span<const double> got_y,
    std::span<const double> got_z,
    std::span<const double> ref_x,
    std::span<const double> ref_y,
    std::span<const double> ref_z) {
  VectorErrorNorm norm;
  for (std::size_t slot = 0; slot < owned.size(); ++slot) {
    const std::size_t global = owned[slot];
    const double dx = got_x[slot] - ref_x[global];
    const double dy = got_y[slot] - ref_y[global];
    const double dz = got_z[slot] - ref_z[global];
    norm.diff_l2 += dx * dx + dy * dy + dz * dz;
    const double ref_mag = std::sqrt(
        ref_x[global] * ref_x[global] + ref_y[global] * ref_y[global] + ref_z[global] * ref_z[global]);
    const double diff_mag = std::sqrt(dx * dx + dy * dy + dz * dz);
    norm.ref_l2 += ref_mag * ref_mag;
    norm.max_rel = std::max(norm.max_rel, diff_mag / std::max(ref_mag, 1.0e-24));
  }
  return norm;
}

void testDistributedPmEquivalence(int world_size, int world_rank) {
  const cosmosim::gravity::PmGridShape shape = hasFastSpectralPmBackend()
      ? cosmosim::gravity::PmGridShape{32, 24, 20}
      : cosmosim::gravity::PmGridShape{8, 6, 5};
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, world_rank);

  std::vector<double> pos_x;
  std::vector<double> pos_y;
  std::vector<double> pos_z;
  std::vector<double> mass;
  buildDeterministicParticles(hasFastSpectralPmBackend() ? 512U : 96U, 1.0, pos_x, pos_y, pos_z, mass);

  const std::vector<std::uint32_t> owned = ownedParticleIndices(pos_x.size(), world_size, world_rank);
  std::vector<double> local_x;
  std::vector<double> local_y;
  std::vector<double> local_z;
  std::vector<double> local_mass;
  local_x.reserve(owned.size());
  local_y.reserve(owned.size());
  local_z.reserve(owned.size());
  local_mass.reserve(owned.size());
  for (const std::uint32_t index : owned) {
    local_x.push_back(pos_x[index]);
    local_y.push_back(pos_y[index]);
    local_z.push_back(pos_z[index]);
    local_mass.push_back(mass[index]);
  }

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 1.0;
  options.gravitational_constant_code = 1.0;
  options.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc;
  options.enable_window_deconvolution = true;

  cosmosim::gravity::PmGridStorage local_grid(shape, layout);
  cosmosim::gravity::PmSolver local_solver(shape);
  std::vector<double> local_ax(owned.size(), 0.0);
  std::vector<double> local_ay(owned.size(), 0.0);
  std::vector<double> local_az(owned.size(), 0.0);
  local_solver.solveForParticles(local_grid, local_x, local_y, local_z, local_mass, local_ax, local_ay, local_az, options, nullptr);

  cosmosim::gravity::PmGridStorage reference_grid(shape);
  cosmosim::gravity::PmSolver reference_solver(shape);
  std::vector<double> ref_ax(pos_x.size(), 0.0);
  std::vector<double> ref_ay(pos_x.size(), 0.0);
  std::vector<double> ref_az(pos_x.size(), 0.0);
  reference_solver.solveForParticles(reference_grid, pos_x, pos_y, pos_z, mass, ref_ax, ref_ay, ref_az, options, nullptr);

  const VectorErrorNorm local_norm =
      accumulateVectorError(owned, local_ax, local_ay, local_az, ref_ax, ref_ay, ref_az);

#if COSMOSIM_ENABLE_MPI
  VectorErrorNorm global_norm{};
  MPI_Allreduce(&local_norm.diff_l2, &global_norm.diff_l2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_norm.ref_l2, &global_norm.ref_l2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_norm.max_rel, &global_norm.max_rel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  const VectorErrorNorm global_norm = local_norm;
#endif

  const double rel_l2 = std::sqrt(global_norm.diff_l2 / std::max(global_norm.ref_l2, 1.0e-30));
  std::ostringstream msg;
  msg << "distributed PM equivalence failed: world_size=" << world_size << ", rel_l2=" << rel_l2
      << " (required <= 1e-10), max_rel=" << global_norm.max_rel;
  requireOrThrow(rel_l2 <= 1.0e-10, msg.str());
}

void runTreePmEquivalenceCase(
    int world_size,
    int world_rank,
    std::size_t case_index,
    const TreePmEquivalenceCase& test_case) {
  const cosmosim::gravity::PmGridShape pm_shape = hasFastSpectralPmBackend()
      ? cosmosim::gravity::PmGridShape{32, 32, 32}
      : cosmosim::gravity::PmGridShape{8, 8, 8};
  const auto layout = cosmosim::parallel::makePmSlabLayout(pm_shape.nx, pm_shape.ny, pm_shape.nz, world_size, world_rank);

  std::vector<double> pos_x;
  std::vector<double> pos_y;
  std::vector<double> pos_z;
  std::vector<double> mass;
  const std::size_t particle_count = hasFastSpectralPmBackend()
      ? test_case.fast_particle_count
      : test_case.fallback_particle_count;
  switch (test_case.particle_pattern) {
    case ParticlePattern::kUniform:
      buildDeterministicParticles(particle_count, 1.0, pos_x, pos_y, pos_z, mass);
      break;
    case ParticlePattern::kClustered:
      buildClusteredParticles(particle_count, 1.0, pos_x, pos_y, pos_z, mass);
      break;
    case ParticlePattern::kPeriodicSeamCluster:
      buildPeriodicSeamCluster(particle_count, 1.0, pos_x, pos_y, pos_z, mass);
      break;
  }

  std::vector<std::uint32_t> owned_sources;
  if (!test_case.disjoint_uneven_ownership) {
    owned_sources = ownedParticleIndices(pos_x.size(), world_size, world_rank);
  } else if (world_size == 1 || world_rank == 0) {
    owned_sources.resize(pos_x.size());
    for (std::size_t i = 0; i < owned_sources.size(); ++i) {
      owned_sources[i] = static_cast<std::uint32_t>(i);
    }
  }

  std::vector<double> local_x;
  std::vector<double> local_y;
  std::vector<double> local_z;
  std::vector<double> local_mass;
  local_x.reserve(owned_sources.size());
  local_y.reserve(owned_sources.size());
  local_z.reserve(owned_sources.size());
  local_mass.reserve(owned_sources.size());
  for (const std::uint32_t index : owned_sources) {
    local_x.push_back(pos_x[index]);
    local_y.push_back(pos_y[index]);
    local_z.push_back(pos_z[index]);
    local_mass.push_back(mass[index]);
  }

  std::vector<std::uint32_t> owned_targets;
  if (!test_case.disjoint_uneven_ownership) {
    owned_targets = owned_sources;
  } else if (world_size == 1 || world_rank == world_size - 1) {
    owned_targets.resize(pos_x.size());
    for (std::size_t i = 0; i < owned_targets.size(); ++i) {
      owned_targets[i] = static_cast<std::uint32_t>(i);
    }
  }

  std::vector<std::uint32_t> local_active(owned_targets.size(), 0U);
  std::vector<double> local_target_x;
  std::vector<double> local_target_y;
  std::vector<double> local_target_z;
  if (test_case.disjoint_uneven_ownership) {
    std::fill(local_active.begin(), local_active.end(), std::numeric_limits<std::uint32_t>::max());
    local_target_x.reserve(owned_targets.size());
    local_target_y.reserve(owned_targets.size());
    local_target_z.reserve(owned_targets.size());
    for (const std::uint32_t index : owned_targets) {
      local_target_x.push_back(pos_x[index]);
      local_target_y.push_back(pos_y[index]);
      local_target_z.push_back(pos_z[index]);
    }
  } else {
    for (std::size_t i = 0; i < local_active.size(); ++i) {
      local_active[i] = static_cast<std::uint32_t>(i);
    }
  }

  std::vector<double> local_previous_acceleration;
  if (test_case.opening_criterion == cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError) {
    local_previous_acceleration.reserve(owned_targets.size());
    for (const std::uint32_t global_index : owned_targets) {
      local_previous_acceleration.push_back(8.0 + 0.125 * static_cast<double>(global_index % 17U));
    }
  }

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.assignment_scheme = test_case.assignment_scheme;
  options.pm_options.enable_window_deconvolution = true;
  options.tree_options.opening_theta = 0.55;
  options.tree_options.opening_criterion = test_case.opening_criterion;
  options.tree_options.multipole_order = test_case.multipole_order;
  options.tree_options.relative_force_tolerance = 0.005;
  options.tree_options.relative_force_acceleration_floor_code = 1.0e-12;
  options.tree_options.max_leaf_size = 8;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 0.008;
  // Keep the direct-DFT fallback mesh small while respecting the periodic
  // one-minimum-image cutoff; the FFTW validation retains the production value.
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25,
      hasFastSpectralPmBackend() ? 6.25 : 3.9,
      1.0 / static_cast<double>(pm_shape.nx));
  options.tree_exchange_batch_bytes = test_case.exchange_batch_bytes;

  std::vector<double> dist_ax(local_active.size(), 0.0);
  std::vector<double> dist_ay(local_active.size(), 0.0);
  std::vector<double> dist_az(local_active.size(), 0.0);
  const cosmosim::gravity::TreePmForceAccumulatorView dist_acc{
      .active_particle_index = local_active,
      .accel_x_comoving = dist_ax,
      .accel_y_comoving = dist_ay,
      .accel_z_comoving = dist_az,
      .previous_acceleration_magnitude_code = local_previous_acceleration,
      .target_pos_x_comoving = local_target_x,
      .target_pos_y_comoving = local_target_y,
      .target_pos_z_comoving = local_target_z,
  };
  cosmosim::gravity::TreePmCoordinator dist_coordinator(pm_shape, layout);
  cosmosim::gravity::TreePmDiagnostics dist_diag;
  ExchangeMetricTotals local_exchange_totals;
  int pm_refresh_count = 0;
  int pm_reuse_count = 0;
  std::uint64_t observed_pm_solve_count = 0;
  std::uint64_t observed_pm_reuse_count = 0;

  cosmosim::gravity::TreePmOptions last_options = options;
  for (int solve_index = 0; solve_index < test_case.solve_count; ++solve_index) {
    last_options = options;
    const bool refresh_long_range = test_case.solve_count == 1 ||
        solve_index == 0 || solve_index == test_case.solve_count - 1;
    const std::uint64_t field_generation =
        (refresh_long_range && solve_index > 0) ? 1U : 0U;
    last_options.decomposition_epoch =
        1000ULL + static_cast<std::uint64_t>(case_index) * 100ULL + field_generation;
    last_options.force_epoch =
        10000ULL + static_cast<std::uint64_t>(case_index) * 100ULL + field_generation;
    if (test_case.solve_count > 1) {
      dist_coordinator.solveActiveSetWithPmCadence(
          local_x,
          local_y,
          local_z,
          local_mass,
          dist_acc,
          last_options,
          refresh_long_range,
          nullptr,
          &dist_diag,
          {});
      pm_refresh_count += refresh_long_range ? 1 : 0;
      pm_reuse_count += refresh_long_range ? 0 : 1;
    } else {
      dist_coordinator.solveActiveSet(
          local_x, local_y, local_z, local_mass, dist_acc, last_options, nullptr, &dist_diag);
      ++pm_refresh_count;
    }
    observed_pm_solve_count += dist_diag.pm_solve_count;
    observed_pm_reuse_count += dist_diag.pm_reuse_count;
    accumulateExchangeMetrics(dist_diag, local_exchange_totals);
  }

  std::vector<std::uint32_t> full_active(pos_x.size(), 0U);
  if (test_case.disjoint_uneven_ownership) {
    std::fill(full_active.begin(), full_active.end(), std::numeric_limits<std::uint32_t>::max());
  } else {
    for (std::size_t i = 0; i < full_active.size(); ++i) {
      full_active[i] = static_cast<std::uint32_t>(i);
    }
  }
  std::vector<double> reference_previous_acceleration;
  if (test_case.opening_criterion == cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError) {
    reference_previous_acceleration.reserve(pos_x.size());
    for (std::size_t i = 0; i < pos_x.size(); ++i) {
      reference_previous_acceleration.push_back(8.0 + 0.125 * static_cast<double>(i % 17U));
    }
  }
  const std::span<const double> reference_target_x = test_case.disjoint_uneven_ownership
      ? std::span<const double>(pos_x)
      : std::span<const double>{};
  const std::span<const double> reference_target_y = test_case.disjoint_uneven_ownership
      ? std::span<const double>(pos_y)
      : std::span<const double>{};
  const std::span<const double> reference_target_z = test_case.disjoint_uneven_ownership
      ? std::span<const double>(pos_z)
      : std::span<const double>{};
  std::vector<double> ref_ax(pos_x.size(), 0.0);
  std::vector<double> ref_ay(pos_x.size(), 0.0);
  std::vector<double> ref_az(pos_x.size(), 0.0);
  const cosmosim::gravity::TreePmForceAccumulatorView ref_acc{
      .active_particle_index = full_active,
      .accel_x_comoving = ref_ax,
      .accel_y_comoving = ref_ay,
      .accel_z_comoving = ref_az,
      .previous_acceleration_magnitude_code = reference_previous_acceleration,
      .target_pos_x_comoving = reference_target_x,
      .target_pos_y_comoving = reference_target_y,
      .target_pos_z_comoving = reference_target_z,
  };
  cosmosim::gravity::TreePmCoordinator ref_coordinator(pm_shape);
  ref_coordinator.solveActiveSet(pos_x, pos_y, pos_z, mass, ref_acc, last_options, nullptr, nullptr);

  const VectorErrorNorm local_norm =
      accumulateVectorError(owned_targets, dist_ax, dist_ay, dist_az, ref_ax, ref_ay, ref_az);

#if COSMOSIM_ENABLE_MPI
  VectorErrorNorm global_norm{};
  MPI_Allreduce(&local_norm.diff_l2, &global_norm.diff_l2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_norm.ref_l2, &global_norm.ref_l2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_norm.max_rel, &global_norm.max_rel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  const VectorErrorNorm global_norm = local_norm;
#endif

  const double rel_l2 = std::sqrt(global_norm.diff_l2 / std::max(global_norm.ref_l2, 1.0e-30));
  const ExchangeMetricTotals global_exchange_totals{
      .request_packets = globalSumU64(local_exchange_totals.request_packets),
      .response_packets = globalSumU64(local_exchange_totals.response_packets),
      .request_bytes = globalSumU64(local_exchange_totals.request_bytes),
      .response_bytes = globalSumU64(local_exchange_totals.response_bytes),
      .batches = globalSumU64(local_exchange_totals.batches),
      .peer_participations = globalSumU64(local_exchange_totals.peer_participations),
      .pruned_nodes = globalSumU64(local_exchange_totals.pruned_nodes),
      .cutoff_pair_skips = globalSumU64(local_exchange_totals.cutoff_pair_skips),
      .remote_pairs_pruned_by_bounds = globalSumU64(local_exchange_totals.remote_pairs_pruned_by_bounds),
  };
  const std::uint64_t global_source_count = globalSumU64(static_cast<std::uint64_t>(local_x.size()));
  const std::uint64_t global_target_count = globalSumU64(static_cast<std::uint64_t>(owned_targets.size()));
  const std::uint64_t empty_source_rank_count = globalSumU64(local_x.empty() ? 1ULL : 0ULL);
  const std::uint64_t empty_target_rank_count = globalSumU64(owned_targets.empty() ? 1ULL : 0ULL);

  std::ostringstream msg;
  msg << "distributed TreePM equivalence failed: case=" << test_case.name
      << ", world_size=" << world_size
      << ", rel_l2=" << rel_l2 << " (required <= 5e-6)"
      << ", max_rel=" << global_norm.max_rel << " (required <= 5e-5)"
      << ", residual_pruned_nodes=" << global_exchange_totals.pruned_nodes
      << ", residual_pair_skips_cutoff=" << global_exchange_totals.cutoff_pair_skips
      << ", force_l2_pm_global=" << dist_diag.force_l2_pm_global
      << ", force_l2_tree_short_range_local=" << dist_diag.force_l2_tree_short_range_local
      << ", force_l2_tree_short_range_remote=" << dist_diag.force_l2_tree_short_range_remote
      << ", residual_remote_request_packets=" << global_exchange_totals.request_packets
      << ", residual_remote_response_packets=" << global_exchange_totals.response_packets
      << ", residual_remote_request_bytes=" << global_exchange_totals.request_bytes
      << ", residual_remote_response_bytes=" << global_exchange_totals.response_bytes
      << ", residual_remote_request_batches=" << global_exchange_totals.batches
      << ", residual_remote_peer_participations=" << global_exchange_totals.peer_participations
      << ", pm_refresh=" << pm_refresh_count
      << ", pm_reuse=" << pm_reuse_count;
  requireOrThrow(rel_l2 <= 5.0e-6, msg.str());
  requireOrThrow(global_norm.max_rel <= 5.0e-5, msg.str());
  requireOrThrow(global_source_count == particle_count, msg.str());
  requireOrThrow(global_target_count == particle_count, msg.str());
  requireOrThrow(global_exchange_totals.request_packets == global_exchange_totals.response_packets, msg.str());
  requireOrThrow(observed_pm_solve_count == static_cast<std::uint64_t>(pm_refresh_count), msg.str());
  requireOrThrow(observed_pm_reuse_count == static_cast<std::uint64_t>(pm_reuse_count), msg.str());
  if (test_case.solve_count > 1) {
    requireOrThrow(pm_refresh_count >= 2, msg.str());
    requireOrThrow(pm_reuse_count >= 1, msg.str());
  }
  if (test_case.disjoint_uneven_ownership && world_size > 1) {
    requireOrThrow(empty_source_rank_count == static_cast<std::uint64_t>(world_size - 1), msg.str());
    requireOrThrow(empty_target_rank_count == static_cast<std::uint64_t>(world_size - 1), msg.str());
  }
  if (test_case.require_multiple_batches && world_size > 1) {
    requireOrThrow(
        global_exchange_totals.batches >
            static_cast<std::uint64_t>(world_size * test_case.solve_count),
        msg.str());
    requireOrThrow(global_exchange_totals.request_packets > 0U, msg.str());
    requireOrThrow(global_exchange_totals.peer_participations > 0U, msg.str());
  }
  if (test_case.require_cutoff_evidence && hasFastSpectralPmBackend()) {
    // Depending on particle layout and MAC acceptance, valid cutoff work can be
    // accounted as internal-node pruning before any leaf pair is inspected. The
    // validation contract is cutoff traversal evidence, not a specific counter.
    requireOrThrow(
        (global_exchange_totals.pruned_nodes + global_exchange_totals.cutoff_pair_skips) > 0,
        msg.str());
    if (world_size > 1) {
      requireOrThrow(global_exchange_totals.remote_pairs_pruned_by_bounds > 0, msg.str());
    }
  }

  if (world_rank == 0) {
    std::ostringstream metric;
    metric << std::scientific << std::setprecision(8)
           << "TreePM equivalence: case=" << test_case.name
           << " ranks=" << world_size
           << " relL2=" << rel_l2
           << " maxRel=" << global_norm.max_rel
           << " packets=" << global_exchange_totals.request_packets
           << " bytes=" << (global_exchange_totals.request_bytes + global_exchange_totals.response_bytes)
           << " peers=" << global_exchange_totals.peer_participations
           << " batches=" << global_exchange_totals.batches
           << " solves=" << test_case.solve_count
           << " refresh=" << pm_refresh_count
           << " reuse=" << pm_reuse_count
           << " emptySourceRanks=" << empty_source_rank_count
           << " emptyTargetRanks=" << empty_target_rank_count;
    std::cout << metric.str() << '\n';
  }
}

void testRestartRoundtripContinuationContract(int world_size, int world_rank) {
#if COSMOSIM_ENABLE_HDF5
  std::stringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[mode]\n";
  stream << "mode = zoom_in\n";
  stream << "ic_file = generated\n";
  stream << "zoom_high_res_region = false\n\n";
  stream << "[numerics]\n";
  stream << "time_begin_code = 0.01\n";
  stream << "time_end_code = 0.0103\n";
  stream << "max_global_steps = 3\n";
  stream << "hierarchical_max_rung = 0\n";
  stream << "treepm_pm_grid = 16\n";
  stream << "treepm_asmth_cells = 1.25\n";
  stream << "treepm_rcut_cells = 6.25\n";
  stream << "treepm_update_cadence_steps = 1\n";
  stream << "treepm_tree_exchange_batch_bytes = 256\n\n";
  stream << "[physics]\n";
  stream << "enable_cooling = false\n";
  stream << "enable_star_formation = false\n";
  stream << "enable_feedback = false\n";
  stream << "enable_stellar_evolution = false\n\n";
  stream << "[output]\n";
  stream << "run_name = validation_phase2_restart\n";
  stream << "output_directory = integration_outputs\n";
  stream << "output_stem = snapshot\n";
  stream << "restart_stem = restart\n";
  stream << "snapshot_interval_steps = 1\n";
  stream << "write_restarts = true\n\n";
  stream << "[parallel]\n";
  stream << "mpi_ranks_expected = " << world_size << "\n";
  stream << "omp_threads = 1\n";
  stream << "gpu_devices = 0\n";

  const auto frozen = cosmosim::core::loadFrozenConfigFromString(stream.str(), "validation_phase2_mpi_gravity_restart");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);

  const std::filesystem::path out_dir =
      std::filesystem::temp_directory_path() / "cosmosim_validation_phase2_restart";
  const auto report = runner.run(out_dir, cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = true});
  requireOrThrow(report.restart_roundtrip_executed, "restart roundtrip was not executed");
  requireOrThrow(report.restart_roundtrip_ok, "restart roundtrip failed verification");
  requireOrThrow(!report.restart_path.empty(), "restart path missing from workflow report");
  requireOrThrow(report.world_size == world_size, "workflow report world_size mismatch");
  requireOrThrow(report.world_rank == world_rank, "workflow report world_rank mismatch");
  requireOrThrow(report.global_particle_count == 42ULL, "distributed workflow global particle count mismatch");
  requireOrThrow(report.global_cell_count == 6ULL, "distributed workflow global cell count mismatch");
  requireOrThrow(report.global_particle_id_sum == (42ULL * 43ULL) / 2ULL, "distributed workflow global particle-id sum mismatch");
  requireOrThrow(report.global_particle_id_xor == xorRangeOneToN(42ULL), "distributed workflow global particle-id xor mismatch");
  if (world_size > 1) {
    requireOrThrow(report.restart_path.filename().string().find("_rank") != std::string::npos, "distributed restart filename must be rank-qualified");
  }

#if COSMOSIM_ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if (world_rank == 0) {
    std::filesystem::remove_all(out_dir);
  }
#if COSMOSIM_ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#else
  (void)world_size;
  (void)world_rank;
#endif
}

void testExplicitFailureContracts(int world_size, int world_rank) {
  {
    bool threw = false;
    try {
      const cosmosim::parallel::MpiContext runtime(
          /*is_enabled=*/world_size > 1,
          world_size,
          world_rank);
      (void)cosmosim::parallel::buildDistributedExecutionTopology(
          16,
          16,
          16,
          runtime,
          world_size + 1,
          /*configured_gpu_devices=*/0,
          /*cuda_runtime_available=*/false,
          /*visible_device_count=*/0);
    } catch (const std::runtime_error&) {
      threw = true;
    }
    requireOrThrow(threw, "rank-count/config mismatch must throw");
  }

  {
    const auto frozen = cosmosim::core::loadFrozenConfigFromString(
        "schema_version = 1\n\n[mode]\nmode = zoom_in\nic_file = generated\n\n"
        "[numerics]\ntreepm_pm_decomposition_mode = pencil\n",
        "validation_phase2_pencil_decomposition");
    requireOrThrow(
        frozen.config.numerics.treepm_pm_decomposition_mode == cosmosim::core::PmDecompositionMode::kPencil,
        "pencil decomposition mode must parse to typed enum");
  }

#if COSMOSIM_ENABLE_MPI
  if (world_size > 1) {
    bool threw = false;
    bool carried_layout_diagnostic = false;
    try {
      const cosmosim::gravity::PmGridShape shape{16, 12, 10};
      const int mismatched_rank = (world_rank + 1) % world_size;
      const auto mismatched_layout = cosmosim::parallel::makePmSlabLayout(
          shape.nx, shape.ny, shape.nz, world_size, mismatched_rank);
      cosmosim::gravity::PmGridStorage grid(shape, mismatched_layout);
      cosmosim::gravity::PmSolver solver(shape);
      cosmosim::gravity::PmSolveOptions options;
      options.box_size_mpc_comoving = 1.0;
      options.scale_factor = 1.0;
      options.gravitational_constant_code = 1.0;
      std::fill(grid.density().begin(), grid.density().end(), 0.0);
      solver.solvePoissonPeriodic(grid, options, nullptr);
    } catch (const std::exception& error) {
      threw = true;
      carried_layout_diagnostic =
          std::string(error.what()).find("slab layout world metadata") !=
          std::string::npos;
    }
    requireOrThrow(threw, "communicator/layout mismatch must throw in PM solve");
    requireOrThrow(
        carried_layout_diagnostic,
        "communicator/layout mismatch must retain its coordinated PM diagnostic");
  }
#endif

  {
    bool threw = false;
    try {
      (void)cosmosim::parallel::DistributedRestartState::deserialize(
          "schema_version=2\ndecomposition_epoch=1\nworld_size=2\npm_grid_nx=8\npm_grid_ny=8\npm_grid_nz=8\n"
          "pm_decomposition_mode=slab\ngravity_kick_opportunity=0\npm_update_cadence_steps=1\n"
          "long_range_field_version=0\nlast_long_range_refresh_opportunity=0\n"
          "long_range_field_built_step_index=0\nlong_range_field_built_scale_factor=1\n"
          "long_range_restart_policy=deterministic_rebuild\nitem_count=1\nrank[0]=0\n"
          "pm_slab_rank_count=2\npm_slab_begin_x[0]=0\npm_slab_end_x[0]=4\n");
    } catch (const std::runtime_error&) {
      threw = true;
    }
    requireOrThrow(threw, "missing distributed restart metadata must throw");
  }

  {
    bool threw = false;
    try {
      (void)cosmosim::parallel::DistributedRestartState::deserialize(
          "schema_version=2\ndecomposition_epoch=1\nworld_size=2\npm_grid_nx=8\npm_grid_ny=8\npm_grid_nz=8\n"
          "pm_decomposition_mode=slab\ngravity_kick_opportunity=2\npm_update_cadence_steps=2\n"
          "long_range_field_version=0\nlast_long_range_refresh_opportunity=1\n"
          "long_range_field_built_step_index=0\nlong_range_field_built_scale_factor=1\n"
          "long_range_restart_policy=deterministic_rebuild\nitem_count=1\nrank[0]=0\n"
          "pm_slab_rank_count=2\npm_slab_begin_x[0]=0\npm_slab_end_x[0]=4\npm_slab_begin_x[1]=4\npm_slab_end_x[1]=8\n");
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    requireOrThrow(threw, "inconsistent cadence restart state must throw");
  }
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#else
  const int world_size = 1;
  const int world_rank = 0;
#endif

  testDistributedPmEquivalence(world_size, world_rank);
  const std::vector<TreePmEquivalenceCase> equivalence_cases{
      TreePmEquivalenceCase{
          .name = "uniform_tsc_quadrupole_com_balanced",
          .particle_pattern = ParticlePattern::kUniform,
          .assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc,
          .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole,
          .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance,
          .fast_particle_count = 192U,
          .fallback_particle_count = 64U,
      },
      TreePmEquivalenceCase{
          .name = "uniform_cic_monopole_geometric_batched_cadence",
          .particle_pattern = ParticlePattern::kUniform,
          .assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kCic,
          .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kMonopole,
          .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric,
          .fast_particle_count = 192U,
          .fallback_particle_count = 64U,
          .solve_count = 4,
          .exchange_batch_bytes = 384U,
          .require_multiple_batches = true,
          .require_cutoff_evidence = true,
      },
      TreePmEquivalenceCase{
          .name = "clustered_tsc_quadrupole_com_uneven_empty_ranks",
          .particle_pattern = ParticlePattern::kClustered,
          .disjoint_uneven_ownership = true,
          .assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc,
          .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole,
          .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance,
          .fast_particle_count = 96U,
          .fallback_particle_count = 40U,
      },
      TreePmEquivalenceCase{
          .name = "periodic_seam_tsc_quadrupole_relative_uneven_empty_ranks",
          .particle_pattern = ParticlePattern::kPeriodicSeamCluster,
          .disjoint_uneven_ownership = true,
          .assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc,
          .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole,
          .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError,
          .fast_particle_count = 96U,
          .fallback_particle_count = 40U,
          .solve_count = 3,
          .exchange_batch_bytes = 768U,
          .require_multiple_batches = true,
      },
  };
  for (std::size_t case_index = 0; case_index < equivalence_cases.size(); ++case_index) {
    runTreePmEquivalenceCase(world_size, world_rank, case_index, equivalence_cases[case_index]);
  }
  testRestartRoundtripContinuationContract(world_size, world_rank);
  testExplicitFailureContracts(world_size, world_rank);

#if COSMOSIM_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
