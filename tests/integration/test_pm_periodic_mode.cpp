#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <span>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void testPeriodicSinusoidalForceResponse() {
  const cosmosim::gravity::PmGridShape shape{32, 8, 8};
  cosmosim::gravity::PmGridStorage grid(shape);
  cosmosim::gravity::PmSolver solver(shape);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_x_mpc_comoving = 1.5;
  options.box_size_y_mpc_comoving = 0.75;
  options.box_size_z_mpc_comoving = 2.25;
  options.scale_factor = 0.8;
  options.gravitational_constant_code = 1.0;

  const double lx = options.box_size_x_mpc_comoving;
  const double amplitude = 0.1;
  const double kx = 2.0 * k_pi / lx;
  for (std::size_t ix = 0; ix < shape.nx; ++ix) {
    const double x = (static_cast<double>(ix) + 0.5) / static_cast<double>(shape.nx) * lx;
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        grid.density()[grid.linearIndex(ix, iy, iz)] = amplitude * std::sin(kx * x);
      }
    }
  }

  cosmosim::gravity::PmProfileEvent profile;
  solver.solvePoissonPeriodic(grid, options, &profile);

  const double expected_force_amp = 4.0 * k_pi * options.gravitational_constant_code *
      amplitude / kx;
  const double expected_phi_amp = -4.0 * k_pi * options.gravitational_constant_code *
      amplitude / (kx * kx);

  double corr_force = 0.0;
  double norm_force_expected = 0.0;
  double norm_force_got = 0.0;
  double corr_phi = 0.0;
  double norm_phi_expected = 0.0;
  double norm_phi_got = 0.0;
  double consistency_rms = 0.0;
  double consistency_ref_rms = 0.0;
  double max_transverse = 0.0;
  const double dx = lx / static_cast<double>(shape.nx);
  for (std::size_t ix = 0; ix < shape.nx; ++ix) {
    const double x = (static_cast<double>(ix) + 0.5) / static_cast<double>(shape.nx) * lx;
    const double expected_force = expected_force_amp * std::cos(kx * x);
    const double expected_phi = expected_phi_amp * std::sin(kx * x);
    const double consistency_force = expected_force;
    const std::size_t ix_prev = (ix + shape.nx - 1U) % shape.nx;
    const std::size_t ix_next = (ix + 1U) % shape.nx;
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        const std::size_t index = grid.linearIndex(ix, iy, iz);
        const std::size_t prev = grid.linearIndex(ix_prev, iy, iz);
        const std::size_t next = grid.linearIndex(ix_next, iy, iz);
        const double accel_x = grid.force_x()[index];
        const double accel_y = grid.force_y()[index];
        const double accel_z = grid.force_z()[index];
        const double phi = grid.potential()[index];
        const double accel_from_phi = -(grid.potential()[next] - grid.potential()[prev]) / (2.0 * dx);

        corr_force += expected_force * accel_x;
        norm_force_expected += expected_force * expected_force;
        norm_force_got += accel_x * accel_x;
        corr_phi += expected_phi * phi;
        norm_phi_expected += expected_phi * expected_phi;
        norm_phi_got += phi * phi;
        const double consistency_error = accel_from_phi - consistency_force;
        consistency_rms += consistency_error * consistency_error;
        consistency_ref_rms += consistency_force * consistency_force;
        max_transverse = std::max(max_transverse, std::max(std::abs(accel_y), std::abs(accel_z)));
      }
    }
  }

  const double force_cosine_similarity = corr_force / std::sqrt(std::max(norm_force_expected * norm_force_got, 1.0e-20));
  const double phi_cosine_similarity = corr_phi / std::sqrt(std::max(norm_phi_expected * norm_phi_got, 1.0e-20));
  const double force_consistency_rel = std::sqrt(consistency_rms / std::max(consistency_ref_rms, 1.0e-20));
#if COSMOSIM_ENABLE_FFTW
  const double min_cosine_similarity = 0.9;
  const double max_consistency_rel = 0.08;
#else
  const double min_cosine_similarity = 0.05;
  const double max_consistency_rel = 5.0;
#endif

  std::ostringstream diag;
  diag << "PM periodic sinusoidal response validation failed: backend='" << cosmosim::gravity::PmSolver::fftBackendName()
       << "', force_cosine_similarity=" << force_cosine_similarity
       << " (required >= " << min_cosine_similarity << ")"
       << ", potential_cosine_similarity=" << phi_cosine_similarity
       << ", force_from_potential_rel_l2=" << force_consistency_rel
       << " (required <= " << max_consistency_rel << ")"
       << ", transverse_max=" << max_transverse
       << ", force_signal_norm=" << std::sqrt(norm_force_got)
       << ", potential_signal_norm=" << std::sqrt(norm_phi_got)
       << ", build_flag.COSMOSIM_ENABLE_FFTW=" << (COSMOSIM_ENABLE_FFTW ? "ON" : "OFF");
  requireOrThrow(std::isfinite(force_cosine_similarity), diag.str());
  requireOrThrow(std::isfinite(phi_cosine_similarity), diag.str());
  requireOrThrow(std::isfinite(force_consistency_rel), diag.str());
  requireOrThrow(norm_force_got > 0.0, diag.str());
  requireOrThrow(norm_phi_got > 0.0, diag.str());
  requireOrThrow(std::abs(force_cosine_similarity) > min_cosine_similarity, diag.str());
  requireOrThrow(std::abs(phi_cosine_similarity) > min_cosine_similarity, diag.str());
  requireOrThrow(force_consistency_rel <= max_consistency_rel, diag.str());
  requireOrThrow(max_transverse < 5.0e-2, diag.str());

  requireOrThrow(profile.assign_ms >= 0.0, "PM profile.assign_ms must be non-negative");
  requireOrThrow(profile.fft_forward_ms >= 0.0, "PM profile.fft_forward_ms must be non-negative");
  requireOrThrow(profile.fft_inverse_ms >= 0.0, "PM profile.fft_inverse_ms must be non-negative");
  requireOrThrow(profile.bytes_moved > 0, "PM profile.bytes_moved must be positive");
}

#if COSMOSIM_ENABLE_MPI
void runDistributedDensityAssignmentCase(cosmosim::gravity::PmAssignmentScheme scheme) {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    return;
  }

  const cosmosim::gravity::PmGridShape shape{8, 6, 6};
  const auto local_layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, world_rank);
  cosmosim::gravity::PmGridStorage local_grid(shape, local_layout);
  cosmosim::gravity::PmSolver local_solver(shape);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 1.0;
  options.assignment_scheme = scheme;

  // Deterministic distributed-owner particle partition:
  // owner rank = particle index % world_size.
  const std::vector<double> all_x{0.499999, 0.500001, -0.0001, 1.0002, 0.125, 0.875};
  const std::vector<double> all_y{0.10, 0.20, 0.30, 0.40, 0.55, 0.65};
  const std::vector<double> all_z{0.15, 0.25, 0.35, 0.45, 0.75, 0.85};
  const std::vector<double> all_mass{1.0, 2.0, 1.5, 2.5, 3.0, 4.0};

  std::vector<double> local_x;
  std::vector<double> local_y;
  std::vector<double> local_z;
  std::vector<double> local_mass;
  for (std::size_t i = 0; i < all_x.size(); ++i) {
    if (static_cast<int>(i % static_cast<std::size_t>(world_size)) == world_rank) {
      local_x.push_back(all_x[i]);
      local_y.push_back(all_y[i]);
      local_z.push_back(all_z[i]);
      local_mass.push_back(all_mass[i]);
    }
  }

  cosmosim::gravity::PmProfileEvent profile;
  local_solver.assignDensity(local_grid, local_x, local_y, local_z, local_mass, options, &profile);
  requireOrThrow(profile.routed_density_records > 0, "Distributed PM density route did not report routed records");
  requireOrThrow(profile.routed_density_peer_count > 0, "Distributed PM density route did not report participating peers");
  requireOrThrow(profile.routed_mpi_bytes_sent > 0, "Distributed PM density route did not report sent wire bytes");
  requireOrThrow(
      profile.routed_mpi_bytes_received > 0,
      "Distributed PM density route did not report received wire bytes");
  requireOrThrow(profile.routed_mpi_wait_ms >= 0.0, "Distributed PM density route reported invalid MPI wait time");

  const double cell_volume = std::pow(options.box_size_mpc_comoving, 3.0) / static_cast<double>(shape.cellCount());
  const double local_mass_deposited =
      std::accumulate(local_grid.density().begin(), local_grid.density().end(), 0.0) * cell_volume;
  double global_mass_deposited = 0.0;
  MPI_Allreduce(&local_mass_deposited, &global_mass_deposited, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double global_particle_mass = std::accumulate(all_mass.begin(), all_mass.end(), 0.0);
  requireOrThrow(
      std::abs(global_mass_deposited - global_particle_mass) < 1.0e-12,
      "Distributed PM density assignment mass conservation failed");

  // Public PM entry is collective in an MPI world, including a full-domain
  // serial reference solve.  Every rank builds the same independent reference
  // and selects only its authoritative comparison slab.
  cosmosim::gravity::PmGridStorage reference_grid(shape);
  cosmosim::gravity::PmSolver reference_solver(shape);
  reference_solver.assignDensity(reference_grid, all_x, all_y, all_z, all_mass, options, nullptr);
  std::vector<double> reference_local_density(local_layout.localCellCount(), 0.0);
  for (std::size_t local_ix = 0; local_ix < local_layout.local_nx(); ++local_ix) {
    const std::size_t ix = local_layout.globalXFromLocal(local_ix);
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        reference_local_density[(local_ix * shape.ny + iy) * shape.nz + iz] =
            reference_grid.density()[reference_grid.linearIndex(ix, iy, iz)];
      }
    }
  }

  double local_max_abs = 0.0;
  double local_ref_l2 = 0.0;
  double local_diff_l2 = 0.0;
  for (std::size_t i = 0; i < local_grid.localCellCount(); ++i) {
    const double diff = local_grid.density()[i] - reference_local_density[i];
    local_max_abs = std::max(local_max_abs, std::abs(diff));
    local_ref_l2 += reference_local_density[i] * reference_local_density[i];
    local_diff_l2 += diff * diff;
  }
  double global_max_abs = 0.0;
  double global_ref_l2 = 0.0;
  double global_diff_l2 = 0.0;
  MPI_Allreduce(&local_max_abs, &global_max_abs, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&local_ref_l2, &global_ref_l2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_diff_l2, &global_diff_l2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double rel_l2 = std::sqrt(global_diff_l2 / std::max(global_ref_l2, 1.0e-30));
  requireOrThrow(global_max_abs <= 1.0e-12, "Distributed PM density assignment max-abs drift exceeds tolerance");
  requireOrThrow(rel_l2 <= 1.0e-12, "Distributed PM density assignment relative L2 drift exceeds tolerance");
}

void testDistributedDensityAssignmentMatchesReference() {
  runDistributedDensityAssignmentCase(cosmosim::gravity::PmAssignmentScheme::kCic);
  runDistributedDensityAssignmentCase(cosmosim::gravity::PmAssignmentScheme::kTsc);
}

void testDistributedDensityRejectsNonFiniteSourceCollectively() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    return;
  }

  const cosmosim::gravity::PmGridShape shape{8, 4, 4};
  const auto layout =
      cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, world_rank);
  cosmosim::gravity::PmGridStorage grid(shape, layout);
  cosmosim::gravity::PmSolver solver(shape);
  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 1.0;

  const std::vector<double> x = world_rank == 0
      ? std::vector<double>{std::numeric_limits<double>::quiet_NaN()}
      : std::vector<double>{};
  const std::vector<double> y(x.size(), 0.25);
  const std::vector<double> z(x.size(), 0.75);
  const std::vector<double> mass(x.size(), 1.0);

  int local_rejected = 0;
  try {
    solver.assignDensity(grid, x, y, z, mass, options, nullptr);
  } catch (const std::invalid_argument&) {
    local_rejected = 1;
  }
  int all_rejected = 0;
  MPI_Allreduce(&local_rejected, &all_rejected, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  requireOrThrow(
      all_rejected == 1,
      "Distributed PM density validation did not reject a rank-local non-finite source collectively");
}

void testDistributedApiPreflightRejectsRankLocalErrorsCollectively() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    return;
  }

  const cosmosim::gravity::PmGridShape shape{8, 4, 4};
  const auto layout =
      cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, world_rank);
  cosmosim::gravity::PmGridStorage grid(shape, layout);
  cosmosim::gravity::PmSolver solver(shape);
  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;

  const auto require_collective_rejection = [&](auto&& call, const std::string& message) {
    int local_rejected = 0;
    try {
      call();
    } catch (const std::exception&) {
      local_rejected = 1;
    }
    int all_rejected = 0;
    MPI_Allreduce(&local_rejected, &all_rejected, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    requireOrThrow(all_rejected == 1, message);
  };

  const std::vector<double> one{0.25};
  const std::vector<double> empty;
  const auto divergent_layout = world_rank == 0
      ? cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, /*world_size=*/1, /*world_rank=*/0)
      : layout;
  cosmosim::gravity::PmGridStorage divergent_grid(shape, divergent_layout);
  require_collective_rejection(
      [&]() {
        solver.assignDensity(divergent_grid, empty, empty, empty, empty, options, nullptr);
      },
      "Divergent serial/distributed PM layout metadata was not rejected collectively");
  require_collective_rejection(
      [&]() {
        if (world_rank == 0) {
          solver.assignDensity(grid, empty, empty, empty, empty, options, nullptr);
        } else {
          solver.solvePoissonPeriodic(grid, options, nullptr);
        }
      },
      "Divergent PM public API ordering was not rejected collectively");

  const cosmosim::gravity::PmGridShape rank_local_shape =
      world_rank == 0 ? cosmosim::gravity::PmGridShape{8, 4, 4}
                      : cosmosim::gravity::PmGridShape{10, 4, 4};
  const auto rank_local_shape_layout = cosmosim::parallel::makePmSlabLayout(
      rank_local_shape.nx,
      rank_local_shape.ny,
      rank_local_shape.nz,
      world_size,
      world_rank);
  cosmosim::gravity::PmGridStorage rank_local_shape_grid(
      rank_local_shape, rank_local_shape_layout);
  cosmosim::gravity::PmSolver rank_local_shape_solver(rank_local_shape);
  require_collective_rejection(
      [&]() {
        rank_local_shape_solver.solvePoissonPeriodic(
            rank_local_shape_grid, options, nullptr);
      },
      "Divergent PM mesh shapes were not rejected before FFTW-MPI plan collectives");

  require_collective_rejection(
      [&]() {
        solver.assignDensity(
            grid,
            world_rank == 0 ? std::span<const double>(one) : std::span<const double>(empty),
            std::span<const double>(empty),
            world_rank == 0 ? std::span<const double>(one) : std::span<const double>(empty),
            world_rank == 0 ? std::span<const double>(one) : std::span<const double>(empty),
            options,
            nullptr);
      },
      "Rank-local PM density span error was not rejected collectively");

  std::vector<double> one_output(1, 0.0);
  std::vector<double> empty_output;
  require_collective_rejection(
      [&]() {
        const std::span<const double> coordinates =
            world_rank == 0 ? std::span<const double>(one) : std::span<const double>(empty);
        solver.interpolateForces(
            grid,
            coordinates,
            coordinates,
            coordinates,
            std::span<double>(empty_output),
            world_rank == 0 ? std::span<double>(one_output) : std::span<double>(empty_output),
            world_rank == 0 ? std::span<double>(one_output) : std::span<double>(empty_output),
            options,
            nullptr);
      },
      "Rank-local PM force span error was not rejected collectively");

  auto invalid_options = options;
  if (world_rank == 0) {
    invalid_options.scale_factor = std::numeric_limits<double>::quiet_NaN();
  }
  require_collective_rejection(
      [&]() { solver.solvePoissonPeriodic(grid, invalid_options, nullptr); },
      "Rank-local distributed PM solve option error was not rejected collectively");
  require_collective_rejection(
      [&]() {
        solver.interpolatePotential(
            grid,
            empty,
            empty,
            empty,
            empty_output,
            invalid_options,
            nullptr);
      },
      "Rank-local PM potential option error was not rejected collectively");
}

void testDistributedPmAcceptsZeroRecordRank() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    return;
  }

  const cosmosim::gravity::PmGridShape shape{8, 4, 4};
  const auto layout =
      cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, world_rank);
  cosmosim::gravity::PmGridStorage grid(shape, layout);
  cosmosim::gravity::PmSolver solver(shape);
  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 1.0;
  options.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kCic;

  const std::vector<double> x = world_rank == 0 ? std::vector<double>{0.25} : std::vector<double>{};
  const std::vector<double> y(x.size(), 0.25);
  const std::vector<double> z(x.size(), 0.25);
  const std::vector<double> mass(x.size(), 1.0);
  std::vector<double> ax(x.size(), 0.0);
  std::vector<double> ay(x.size(), 0.0);
  std::vector<double> az(x.size(), 0.0);
  std::vector<double> potential(x.size(), 0.0);
  cosmosim::gravity::PmProfileEvent profile;

  solver.solveForParticles(grid, x, y, z, mass, ax, ay, az, options, &profile);
  solver.interpolatePotential(grid, x, y, z, potential, options, &profile);
  requireOrThrow(profile.routed_mpi_wait_ms >= 0.0, "Zero-record PM rank reported invalid MPI wait time");
  if (world_rank == 1) {
    requireOrThrow(profile.routed_mpi_bytes_sent == 0U, "Zero-record PM rank unexpectedly sent wire payload");
    requireOrThrow(profile.routed_mpi_bytes_received == 0U, "Zero-record PM rank unexpectedly received wire payload");
  }
}

void testDistributedOpenBoundaryRoutingMatchesRectangularReference() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    return;
  }

  const cosmosim::gravity::PmGridShape shape{4, 3, 2};
  const auto layout =
      cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, world_rank);
  cosmosim::gravity::PmGridStorage local_grid(shape, layout);
  cosmosim::gravity::PmSolver local_solver(shape);
  cosmosim::gravity::PmSolveOptions options;
  options.box_size_x_mpc_comoving = 2.0;
  options.box_size_y_mpc_comoving = 3.0;
  options.box_size_z_mpc_comoving = 4.0;
  options.boundary_condition = cosmosim::gravity::PmBoundaryCondition::kIsolatedOpen;

  const std::vector<double> local_source_x = world_rank == 0 ? std::vector<double>{1.95} : std::vector<double>{};
  const std::vector<double> local_source_y(local_source_x.size(), 0.0);
  const std::vector<double> local_source_z(local_source_x.size(), 0.0);
  const std::vector<double> local_mass(local_source_x.size(), 1.0);
  local_solver.assignDensity(
      local_grid,
      local_source_x,
      local_source_y,
      local_source_z,
      local_mass,
      options,
      nullptr);

  cosmosim::gravity::PmGridStorage reference_grid(shape);
  cosmosim::gravity::PmSolver reference_solver(shape);
  const std::vector<double> reference_source_x{1.95};
  const std::vector<double> reference_source_y{0.0};
  const std::vector<double> reference_source_z{0.0};
  const std::vector<double> reference_mass{1.0};
  reference_solver.assignDensity(
      reference_grid,
      reference_source_x,
      reference_source_y,
      reference_source_z,
      reference_mass,
      options,
      nullptr);
  for (std::size_t local_ix = 0; local_ix < layout.local_nx(); ++local_ix) {
    const std::size_t global_ix = layout.globalXFromLocal(local_ix);
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        const std::size_t local_index = local_grid.linearIndex(global_ix, iy, iz);
        const std::size_t reference_index = reference_grid.linearIndex(global_ix, iy, iz);
        requireOrThrow(
            std::abs(local_grid.density()[local_index] - reference_grid.density()[reference_index]) < 1.0e-12,
            "Distributed isolated/open density routing wrapped a rectangular boundary stencil");
      }
    }
  }

  for (std::size_t ix = 0; ix < shape.nx; ++ix) {
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        const double value = 1.0 + static_cast<double>(ix + 10U * iy + 100U * iz);
        const std::size_t reference_index = reference_grid.linearIndex(ix, iy, iz);
        reference_grid.force_x()[reference_index] = value;
        reference_grid.force_y()[reference_index] = -2.0 * value;
        reference_grid.force_z()[reference_index] = 0.5 * value;
        reference_grid.potential()[reference_index] = 3.0 * value;
        if (layout.ownsGlobalX(ix)) {
          const std::size_t local_index = local_grid.linearIndex(ix, iy, iz);
          local_grid.force_x()[local_index] = value;
          local_grid.force_y()[local_index] = -2.0 * value;
          local_grid.force_z()[local_index] = 0.5 * value;
          local_grid.potential()[local_index] = 3.0 * value;
        }
      }
    }
  }

  const std::vector<double> target_x = world_rank == 0
      ? std::vector<double>{0.95, -0.1}
      : std::vector<double>{1.95, 2.1};
  const std::vector<double> target_y(target_x.size(), 1.0);
  const std::vector<double> target_z(target_x.size(), 1.0);
  std::vector<double> local_ax(target_x.size(), 0.0);
  std::vector<double> local_ay(target_x.size(), 0.0);
  std::vector<double> local_az(target_x.size(), 0.0);
  std::vector<double> local_potential(target_x.size(), 0.0);
  std::vector<double> reference_ax(target_x.size(), 0.0);
  std::vector<double> reference_ay(target_x.size(), 0.0);
  std::vector<double> reference_az(target_x.size(), 0.0);
  std::vector<double> reference_potential(target_x.size(), 0.0);
  local_solver.interpolateForces(
      local_grid, target_x, target_y, target_z, local_ax, local_ay, local_az, options, nullptr);
  local_solver.interpolatePotential(
      local_grid, target_x, target_y, target_z, local_potential, options, nullptr);
  reference_solver.interpolateForces(
      reference_grid,
      target_x,
      target_y,
      target_z,
      reference_ax,
      reference_ay,
      reference_az,
      options,
      nullptr);
  reference_solver.interpolatePotential(
      reference_grid,
      target_x,
      target_y,
      target_z,
      reference_potential,
      options,
      nullptr);
  for (std::size_t p = 0; p < target_x.size(); ++p) {
    requireOrThrow(std::abs(local_ax[p] - reference_ax[p]) < 1.0e-12, "Distributed open force-x routing mismatch");
    requireOrThrow(std::abs(local_ay[p] - reference_ay[p]) < 1.0e-12, "Distributed open force-y routing mismatch");
    requireOrThrow(std::abs(local_az[p] - reference_az[p]) < 1.0e-12, "Distributed open force-z routing mismatch");
    requireOrThrow(
        std::abs(local_potential[p] - reference_potential[p]) < 1.0e-12,
        "Distributed open potential routing wrapped a rectangular boundary stencil");
  }
}

void testDistributedPmSolveWithZeroWidthSlabRanks() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 4) {
    return;
  }

  const cosmosim::gravity::PmGridShape shape{2, 4, 4};
  const auto layout =
      cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, world_rank);
  cosmosim::gravity::PmGridStorage grid(shape, layout);
  cosmosim::gravity::PmSolver solver(shape);
  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;

  const std::vector<double> x = world_rank == 0 ? std::vector<double>{0.25} : std::vector<double>{};
  const std::vector<double> y(x.size(), 0.25);
  const std::vector<double> z(x.size(), 0.25);
  const std::vector<double> mass(x.size(), 1.0);
  std::vector<double> ax(x.size(), 0.0);
  std::vector<double> ay(x.size(), 0.0);
  std::vector<double> az(x.size(), 0.0);
  std::vector<double> potential(x.size(), 0.0);
  solver.solveForParticles(grid, x, y, z, mass, ax, ay, az, options, nullptr);
  solver.interpolatePotential(grid, x, y, z, potential, options, nullptr);

  const int local_zero_width = layout.local_nx() == 0U ? 1 : 0;
  int global_zero_width = 0;
  MPI_Allreduce(&local_zero_width, &global_zero_width, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  requireOrThrow(global_zero_width == 2, "Expected two zero-width PM slab ranks for nx=2,np=4");
  requireOrThrow(
      std::all_of(grid.force_x().begin(), grid.force_x().end(), [](double value) { return std::isfinite(value); }),
      "Zero-width PM slab solve produced a non-finite force");
  requireOrThrow(
      std::all_of(potential.begin(), potential.end(), [](double value) { return std::isfinite(value); }),
      "Zero-width PM slab solve produced a non-finite interpolated potential");
}


void runDistributedInterpolationAgreementCase(cosmosim::gravity::PmAssignmentScheme scheme, bool gather_potential) {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    return;
  }

  const cosmosim::gravity::PmGridShape shape{16, 8, 8};
  const auto local_layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, world_rank);
  cosmosim::gravity::PmGridStorage local_grid(shape, local_layout);
  cosmosim::gravity::PmSolver local_solver(shape);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 1.0;
  options.gravitational_constant_code = 1.0;
  options.assignment_scheme = scheme;

  const std::vector<double> all_x{0.499999, 0.500001, -0.0001, 1.0002, 0.125, 0.875, 0.24999, 0.75001};
  const std::vector<double> all_y{0.10, 0.20, 0.30, 0.40, 0.55, 0.65, 0.35, 0.85};
  const std::vector<double> all_z{0.15, 0.25, 0.35, 0.45, 0.75, 0.85, 0.95, 0.05};
  const std::vector<double> all_mass{1.0, 2.0, 1.5, 2.5, 3.0, 4.0, 1.25, 0.75};

  // Deliberately uneven owner-local populations: this catches slab owners
  // that mistake an origin-local return token for a receiver-local row.
  const auto particle_owner_rank = [](std::size_t particle_index) {
    return particle_index < 7U ? 0 : 1;
  };

  std::vector<double> local_x;
  std::vector<double> local_y;
  std::vector<double> local_z;
  std::vector<double> local_mass;
  for (std::size_t i = 0; i < all_x.size(); ++i) {
    if (particle_owner_rank(i) == world_rank) {
      local_x.push_back(all_x[i]);
      local_y.push_back(all_y[i]);
      local_z.push_back(all_z[i]);
      local_mass.push_back(all_mass[i]);
    }
  }

  std::vector<double> local_ax(local_x.size(), 0.0);
  std::vector<double> local_ay(local_x.size(), 0.0);
  std::vector<double> local_az(local_x.size(), 0.0);
  cosmosim::gravity::PmProfileEvent profile;
  local_solver.solveForParticles(local_grid, local_x, local_y, local_z, local_mass, local_ax, local_ay, local_az, options, &profile);
  requireOrThrow(profile.routed_density_records > 0, "Distributed PM solve did not report routed density records");
  const std::uint64_t local_force_route_count = profile.routed_force_requests + profile.force_halo_cache_hits;
  std::uint64_t global_force_route_count = 0U;
  MPI_Allreduce(
      &local_force_route_count,
      &global_force_route_count,
      1,
      MPI_UINT64_T,
      MPI_SUM,
      MPI_COMM_WORLD);
  requireOrThrow(
      global_force_route_count > 0,
      "Distributed PM solve did not report any global remote force route or halo-cache hit");
  requireOrThrow(profile.routed_mpi_bytes_sent > 0, "Distributed PM solve did not report sent wire bytes");
  requireOrThrow(profile.routed_mpi_bytes_received > 0, "Distributed PM solve did not report received wire bytes");
  requireOrThrow(profile.routed_mpi_wait_ms >= 0.0, "Distributed PM solve reported invalid MPI wait time");

  std::vector<double> local_phi(local_x.size(), 0.0);
  if (gather_potential) {
    local_solver.interpolatePotential(local_grid, local_x, local_y, local_z, local_phi, options, &profile);
    requireOrThrow(
        profile.routed_potential_requests > 0,
        "Distributed PM potential gather did not report routed requests");
  }

  cosmosim::gravity::PmGridStorage reference_grid(shape);
  cosmosim::gravity::PmSolver reference_solver(shape);
  std::vector<double> ref_ax(all_x.size(), 0.0);
  std::vector<double> ref_ay(all_x.size(), 0.0);
  std::vector<double> ref_az(all_x.size(), 0.0);
  reference_solver.solveForParticles(
      reference_grid, all_x, all_y, all_z, all_mass, ref_ax, ref_ay, ref_az, options, nullptr);
  std::vector<double> ref_phi(all_x.size(), 0.0);
  if (gather_potential) {
    reference_solver.interpolatePotential(reference_grid, all_x, all_y, all_z, ref_phi, options, nullptr);
  }

  std::vector<double> expected_local_ax;
  std::vector<double> expected_local_ay;
  std::vector<double> expected_local_az;
  std::vector<double> expected_local_phi;
  expected_local_ax.reserve(local_x.size());
  expected_local_ay.reserve(local_x.size());
  expected_local_az.reserve(local_x.size());
  expected_local_phi.reserve(local_x.size());
  for (std::size_t i = 0; i < all_x.size(); ++i) {
    if (particle_owner_rank(i) != world_rank) {
      continue;
    }
    expected_local_ax.push_back(ref_ax[i]);
    expected_local_ay.push_back(ref_ay[i]);
    expected_local_az.push_back(ref_az[i]);
    if (gather_potential) {
      expected_local_phi.push_back(ref_phi[i]);
    }
  }

  double local_diff2 = 0.0;
  double local_ref2 = 0.0;
  for (std::size_t i = 0; i < local_ax.size(); ++i) {
    const double dx = local_ax[i] - expected_local_ax[i];
    const double dy = local_ay[i] - expected_local_ay[i];
    const double dz = local_az[i] - expected_local_az[i];
    local_diff2 += dx * dx + dy * dy + dz * dz;
    local_ref2 += expected_local_ax[i] * expected_local_ax[i] + expected_local_ay[i] * expected_local_ay[i] +
        expected_local_az[i] * expected_local_az[i];
  }
  double global_diff2 = 0.0;
  double global_ref2 = 0.0;
  MPI_Allreduce(&local_diff2, &global_diff2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_ref2, &global_ref2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double force_rel_l2 = std::sqrt(global_diff2 / std::max(global_ref2, 1.0e-30));
  requireOrThrow(force_rel_l2 <= 1.0e-10, "Distributed PM interpolation force agreement drift exceeds tolerance");

  if (gather_potential) {
    double local_phi_diff2 = 0.0;
    double local_phi_ref2 = 0.0;
    for (std::size_t i = 0; i < local_phi.size(); ++i) {
      const double diff = local_phi[i] - expected_local_phi[i];
      local_phi_diff2 += diff * diff;
      local_phi_ref2 += expected_local_phi[i] * expected_local_phi[i];
    }
    double global_phi_diff2 = 0.0;
    double global_phi_ref2 = 0.0;
    MPI_Allreduce(&local_phi_diff2, &global_phi_diff2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_phi_ref2, &global_phi_ref2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    const double phi_rel_l2 = std::sqrt(global_phi_diff2 / std::max(global_phi_ref2, 1.0e-30));
    requireOrThrow(phi_rel_l2 <= 1.0e-10, "Distributed PM interpolation potential agreement drift exceeds tolerance");
  }
}

void testDistributedInterpolationAgreement() {
  runDistributedInterpolationAgreementCase(cosmosim::gravity::PmAssignmentScheme::kCic, true);
  runDistributedInterpolationAgreementCase(cosmosim::gravity::PmAssignmentScheme::kTsc, true);
}

void testDistributedTwoRankMatchesSingleRankReference() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    return;
  }

  const cosmosim::gravity::PmGridShape shape{16, 8, 8};
  const auto local_layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, world_rank);
  cosmosim::gravity::PmGridStorage local_grid(shape, local_layout);
  cosmosim::gravity::PmSolver local_solver(shape);
  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 0.8;
  options.gravitational_constant_code = 1.0;

  const double amplitude = 0.1;
  const double kx = 2.0 * k_pi / options.box_size_mpc_comoving;
  cosmosim::gravity::PmGridStorage reference_grid(shape);
  cosmosim::gravity::PmSolver reference_solver(shape);
  for (std::size_t ix = 0; ix < shape.nx; ++ix) {
    const double x =
        (static_cast<double>(ix) + 0.5) / static_cast<double>(shape.nx) * options.box_size_mpc_comoving;
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        reference_grid.density()[reference_grid.linearIndex(ix, iy, iz)] = amplitude * std::sin(kx * x);
      }
    }
  }
  reference_solver.solvePoissonPeriodic(reference_grid, options, nullptr);

  for (std::size_t local_ix = 0; local_ix < local_layout.local_nx(); ++local_ix) {
    const std::size_t ix = local_layout.globalXFromLocal(local_ix);
    const double x = (static_cast<double>(ix) + 0.5) / static_cast<double>(shape.nx) * options.box_size_mpc_comoving;
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        local_grid.density()[local_grid.linearIndex(ix, iy, iz)] = amplitude * std::sin(kx * x);
      }
    }
  }
  local_solver.solvePoissonPeriodic(local_grid, options, nullptr);

  std::vector<double> reference_force_x(local_grid.localCellCount(), 0.0);
  for (std::size_t local_ix = 0; local_ix < local_layout.local_nx(); ++local_ix) {
    const std::size_t ix = local_layout.globalXFromLocal(local_ix);
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        reference_force_x[(local_ix * shape.ny + iy) * shape.nz + iz] =
            reference_grid.force_x()[reference_grid.linearIndex(ix, iy, iz)];
      }
    }
  }

  constexpr double k_rel_l2_tol = 2.0e-6;
  double diff2 = 0.0;
  double ref2 = 0.0;
  for (std::size_t i = 0; i < local_grid.localCellCount(); ++i) {
    const double diff = local_grid.force_x()[i] - reference_force_x[i];
    diff2 += diff * diff;
    ref2 += reference_force_x[i] * reference_force_x[i];
  }
  double global_diff2 = 0.0;
  double global_ref2 = 0.0;
  MPI_Allreduce(&diff2, &global_diff2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ref2, &global_ref2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double rel_l2 = std::sqrt(global_diff2 / std::max(global_ref2, 1.0e-30));
  requireOrThrow(rel_l2 <= k_rel_l2_tol, "Distributed PM force_x drift vs one-rank reference exceeds tolerance");

  requireOrThrow(local_solver.cachedPlanCount() == 1U, "Distributed PM solver should cache one plan per slab layout");
  const std::size_t plan_builds_before = local_solver.planBuildCount();
  local_solver.solvePoissonPeriodic(local_grid, options, nullptr);
  requireOrThrow(local_solver.planBuildCount() == plan_builds_before, "Distributed PM solver should reuse FFT plan cache");
}
#endif

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
#endif
  testPeriodicSinusoidalForceResponse();
#if COSMOSIM_ENABLE_MPI
  testDistributedDensityRejectsNonFiniteSourceCollectively();
  testDistributedApiPreflightRejectsRankLocalErrorsCollectively();
  testDistributedPmAcceptsZeroRecordRank();
  testDistributedOpenBoundaryRoutingMatchesRectangularReference();
  testDistributedPmSolveWithZeroWidthSlabRanks();
  testDistributedDensityAssignmentMatchesReference();
  testDistributedTwoRankMatchesSingleRankReference();
  testDistributedInterpolationAgreement();
  MPI_Finalize();
#endif
  return 0;
}
