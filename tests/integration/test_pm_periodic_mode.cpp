#include <algorithm>
#include <cmath>
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
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 0.8;
  options.gravitational_constant_code = 1.0;

  const double amplitude = 0.1;
  const double kx = 2.0 * k_pi / options.box_size_mpc_comoving;
  for (std::size_t ix = 0; ix < shape.nx; ++ix) {
    const double x = (static_cast<double>(ix) + 0.5) / static_cast<double>(shape.nx) * options.box_size_mpc_comoving;
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        grid.density()[grid.linearIndex(ix, iy, iz)] = amplitude * std::sin(kx * x);
      }
    }
  }

  cosmosim::gravity::PmProfileEvent profile;
  solver.solvePoissonPeriodic(grid, options, &profile);

  const double expected_force_amp = 4.0 * k_pi * options.gravitational_constant_code *
      options.scale_factor * options.scale_factor * amplitude / kx;
  const double expected_phi_amp = -4.0 * k_pi * options.gravitational_constant_code *
      options.scale_factor * options.scale_factor * amplitude / (kx * kx);

  double corr_force = 0.0;
  double norm_force_expected = 0.0;
  double norm_force_got = 0.0;
  double corr_phi = 0.0;
  double norm_phi_expected = 0.0;
  double norm_phi_got = 0.0;
  double consistency_rms = 0.0;
  double consistency_ref_rms = 0.0;
  double max_transverse = 0.0;
  const double dx = options.box_size_mpc_comoving / static_cast<double>(shape.nx);
  for (std::size_t ix = 0; ix < shape.nx; ++ix) {
    const double x = (static_cast<double>(ix) + 0.5) / static_cast<double>(shape.nx) * options.box_size_mpc_comoving;
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
  if (world_rank == 0) {
    cosmosim::gravity::PmGridStorage reference_grid(shape);
    cosmosim::gravity::PmSolver reference_solver(shape);
    for (std::size_t ix = 0; ix < shape.nx; ++ix) {
      const double x = (static_cast<double>(ix) + 0.5) / static_cast<double>(shape.nx) * options.box_size_mpc_comoving;
      for (std::size_t iy = 0; iy < shape.ny; ++iy) {
        for (std::size_t iz = 0; iz < shape.nz; ++iz) {
          reference_grid.density()[reference_grid.linearIndex(ix, iy, iz)] = amplitude * std::sin(kx * x);
        }
      }
    }
    reference_solver.solvePoissonPeriodic(reference_grid, options, nullptr);
    for (std::size_t local_ix = 0; local_ix < local_layout.local_nx(); ++local_ix) {
      const std::size_t ix = local_layout.globalXFromLocal(local_ix);
      for (std::size_t iy = 0; iy < shape.ny; ++iy) {
        for (std::size_t iz = 0; iz < shape.nz; ++iz) {
          reference_force_x[(local_ix * shape.ny + iy) * shape.nz + iz] =
              reference_grid.force_x()[reference_grid.linearIndex(ix, iy, iz)];
        }
      }
    }
  }

  if (world_rank == 0) {
    for (int target_rank = 1; target_rank < world_size; ++target_rank) {
      const auto target_layout =
          cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, target_rank);
      std::vector<double> target_ref(target_layout.localCellCount(), 0.0);
      cosmosim::gravity::PmGridStorage reference_grid(shape);
      cosmosim::gravity::PmSolver reference_solver(shape);
      for (std::size_t ix = 0; ix < shape.nx; ++ix) {
        const double x = (static_cast<double>(ix) + 0.5) / static_cast<double>(shape.nx) * options.box_size_mpc_comoving;
        for (std::size_t iy = 0; iy < shape.ny; ++iy) {
          for (std::size_t iz = 0; iz < shape.nz; ++iz) {
            reference_grid.density()[reference_grid.linearIndex(ix, iy, iz)] = amplitude * std::sin(kx * x);
          }
        }
      }
      reference_solver.solvePoissonPeriodic(reference_grid, options, nullptr);
      for (std::size_t local_ix = 0; local_ix < target_layout.local_nx(); ++local_ix) {
        const std::size_t ix = target_layout.globalXFromLocal(local_ix);
        for (std::size_t iy = 0; iy < shape.ny; ++iy) {
          for (std::size_t iz = 0; iz < shape.nz; ++iz) {
            target_ref[(local_ix * shape.ny + iy) * shape.nz + iz] =
                reference_grid.force_x()[reference_grid.linearIndex(ix, iy, iz)];
          }
        }
      }
      MPI_Send(
          target_ref.data(),
          static_cast<int>(target_ref.size()),
          MPI_DOUBLE,
          target_rank,
          /*tag=*/41,
          MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(
        reference_force_x.data(),
        static_cast<int>(reference_force_x.size()),
        MPI_DOUBLE,
        0,
        /*tag=*/41,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);
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
  testDistributedTwoRankMatchesSingleRankReference();
  MPI_Finalize();
#endif
  return 0;
}
