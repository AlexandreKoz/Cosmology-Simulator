#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include "cosmosim/gravity/tree_gravity.hpp"

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
    x[i] = std::fmod((17.0 * static_cast<double>(i) + 3.0) * 0.00083, 1.0);
    y[i] = std::fmod((23.0 * static_cast<double>(i) + 5.0) * 0.00061, 1.0);
    z[i] = std::fmod((41.0 * static_cast<double>(i) + 7.0) * 0.00043, 1.0);
    mass[i] = 0.7 + 0.02 * static_cast<double>(i % 9U);
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

  constexpr std::size_t particle_count = 140000;
  std::vector<double> pos_x;
  std::vector<double> pos_y;
  std::vector<double> pos_z;
  std::vector<double> mass;
  buildParticles(particle_count, pos_x, pos_y, pos_z, mass);

  std::vector<std::uint32_t> active;
  active.reserve(particle_count / static_cast<std::size_t>(world_size) + 1U);
  for (std::size_t i = 0; i < particle_count; ++i) {
    if (static_cast<int>(i % static_cast<std::size_t>(world_size)) == world_rank) {
      active.push_back(static_cast<std::uint32_t>(i));
    }
  }

  std::vector<double> ax(active.size(), 0.0);
  std::vector<double> ay(active.size(), 0.0);
  std::vector<double> az(active.size(), 0.0);

  cosmosim::gravity::TreeGravityOptions options;
  options.opening_theta = 0.55;
  options.max_leaf_size = 16;
  options.gravitational_constant_code = 1.0;
  options.softening.epsilon_comoving = 1.0e-3;

  cosmosim::gravity::TreeGravitySolver solver;
  solver.build(pos_x, pos_y, pos_z, mass, options, nullptr);

  constexpr int warmup = 1;
  constexpr int measured = 5;
  for (int i = 0; i < warmup; ++i) {
    solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, ax, ay, az, options, nullptr);
  }

#if COSMOSIM_ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  const auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < measured; ++i) {
    solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, ax, ay, az, options, nullptr);
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

  std::size_t min_local = active.size();
  std::size_t max_local = active.size();
#if COSMOSIM_ENABLE_MPI
  const unsigned long long local_count = static_cast<unsigned long long>(active.size());
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
    const std::filesystem::path out_path = out_dir / ("tree_only_scaling_np" + std::to_string(world_size) + ".csv");
    std::ofstream out(out_path);
    out << "world_size,particle_count,active_count_total,avg_rank_step_ms,max_rank_step_ms,local_active_min,local_active_max\n";
    out << world_size << ',' << particle_count << ',' << particle_count << ',' << avg_rank_ms << ',' << max_rank_ms << ','
        << min_local << ',' << max_local << '\n';
    out.flush();
    std::cout << "bench_tree_only_scaling_mpi artifact=" << out_path.string()
              << " max_rank_step_ms=" << max_rank_ms << '\n';
  }

#if COSMOSIM_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
