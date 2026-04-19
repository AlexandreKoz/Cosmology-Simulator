#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

void buildParticles(
    std::size_t count,
    std::vector<double>& x,
    std::vector<double>& y,
    std::vector<double>& z,
    std::vector<double>& mass) {
  x.resize(count);
  y.resize(count);
  z.resize(count);
  mass.resize(count, 1.0);
  for (std::size_t i = 0; i < count; ++i) {
    x[i] = std::fmod((31.0 * static_cast<double>(i) + 7.0) * 0.00071, 1.0);
    y[i] = std::fmod((47.0 * static_cast<double>(i) + 11.0) * 0.00053, 1.0);
    z[i] = std::fmod((59.0 * static_cast<double>(i) + 13.0) * 0.00037, 1.0);
    mass[i] = 0.9 + 0.03 * static_cast<double>(i % 5U);
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

  constexpr std::size_t particle_count = 120000;
  const cosmosim::gravity::PmGridShape shape{64, 64, 64};
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, world_rank);

  std::vector<double> all_x;
  std::vector<double> all_y;
  std::vector<double> all_z;
  std::vector<double> all_mass;
  buildParticles(particle_count, all_x, all_y, all_z, all_mass);

  std::vector<double> local_x;
  std::vector<double> local_y;
  std::vector<double> local_z;
  std::vector<double> local_mass;
  for (std::size_t i = 0; i < particle_count; ++i) {
    if (static_cast<int>(i % static_cast<std::size_t>(world_size)) == world_rank) {
      local_x.push_back(all_x[i]);
      local_y.push_back(all_y[i]);
      local_z.push_back(all_z[i]);
      local_mass.push_back(all_mass[i]);
    }
  }

  std::vector<double> local_ax(local_x.size(), 0.0);
  std::vector<double> local_ay(local_x.size(), 0.0);
  std::vector<double> local_az(local_x.size(), 0.0);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 1.0;
  options.gravitational_constant_code = 1.0;
  options.enable_window_deconvolution = true;
  options.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc;

  cosmosim::gravity::PmGridStorage grid(shape, layout);
  cosmosim::gravity::PmSolver solver(shape);

  constexpr int warmup = 1;
  constexpr int measured = 5;
  cosmosim::gravity::PmProfileEvent profile;
  for (int i = 0; i < warmup; ++i) {
    solver.solveForParticles(grid, local_x, local_y, local_z, local_mass, local_ax, local_ay, local_az, options, &profile);
  }

#if COSMOSIM_ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  const auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < measured; ++i) {
    solver.solveForParticles(grid, local_x, local_y, local_z, local_mass, local_ax, local_ay, local_az, options, &profile);
  }
  const auto stop = std::chrono::steady_clock::now();
  const double local_ms = std::chrono::duration<double, std::milli>(stop - start).count() / static_cast<double>(measured);

  double max_rank_ms = local_ms;
  double avg_rank_ms = local_ms;
#if COSMOSIM_ENABLE_MPI
  MPI_Allreduce(&local_ms, &max_rank_ms, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&local_ms, &avg_rank_ms, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  avg_rank_ms /= static_cast<double>(world_size);
#endif

  std::size_t min_local = local_x.size();
  std::size_t max_local = local_x.size();
#if COSMOSIM_ENABLE_MPI
  const unsigned long long local_count = static_cast<unsigned long long>(local_x.size());
  unsigned long long global_min = local_count;
  unsigned long long global_max = local_count;
  MPI_Allreduce(&local_count, &global_min, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_count, &global_max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
  min_local = static_cast<std::size_t>(global_min);
  max_local = static_cast<std::size_t>(global_max);
#endif

  if (world_rank == 0) {
    const std::filesystem::path out_dir = std::filesystem::path("validation") / "artifacts";
    std::filesystem::create_directories(out_dir);
    const std::filesystem::path out_path = out_dir / ("pm_only_scaling_np" + std::to_string(world_size) + ".csv");
    std::ofstream out(out_path);
    out << "world_size,particle_count,pm_grid_nx,avg_rank_step_ms,max_rank_step_ms,local_particle_balance_min,local_particle_balance_max\n";
    out << world_size << ',' << particle_count << ',' << shape.nx << ',' << avg_rank_ms << ',' << max_rank_ms << ','
        << min_local << ',' << max_local << '\n';
    out.flush();
    std::cout << "bench_pm_only_scaling_mpi artifact=" << out_path.string()
              << " max_rank_step_ms=" << max_rank_ms << '\n';
  }

#if COSMOSIM_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
