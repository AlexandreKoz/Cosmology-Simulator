#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
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
  const cosmosim::gravity::PmGridShape shape{32, 24, 20};
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, world_rank);

  std::vector<double> pos_x;
  std::vector<double> pos_y;
  std::vector<double> pos_z;
  std::vector<double> mass;
  buildDeterministicParticles(512, 1.0, pos_x, pos_y, pos_z, mass);

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

void runTreePmCase(int world_size, int world_rank, bool communication_stress) {
  const cosmosim::gravity::PmGridShape pm_shape{32, 32, 32};
  const auto layout = cosmosim::parallel::makePmSlabLayout(pm_shape.nx, pm_shape.ny, pm_shape.nz, world_size, world_rank);

  std::vector<double> pos_x;
  std::vector<double> pos_y;
  std::vector<double> pos_z;
  std::vector<double> mass;
  buildDeterministicParticles(384, 1.0, pos_x, pos_y, pos_z, mass);
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
  std::vector<std::uint32_t> local_active(local_x.size(), 0U);
  for (std::size_t i = 0; i < local_active.size(); ++i) {
    local_active[i] = static_cast<std::uint32_t>(i);
  }

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc;
  options.pm_options.enable_window_deconvolution = true;
  options.tree_options.opening_theta = 0.55;
  options.tree_options.max_leaf_size = 8;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 0.008;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25,
      4.5,
      1.0 / static_cast<double>(pm_shape.nx));
  if (communication_stress) {
    options.tree_exchange_batch_bytes = 64;
  }

  std::vector<double> dist_ax(local_active.size(), 0.0);
  std::vector<double> dist_ay(local_active.size(), 0.0);
  std::vector<double> dist_az(local_active.size(), 0.0);
  cosmosim::gravity::TreePmForceAccumulatorView dist_acc{local_active, dist_ax, dist_ay, dist_az};
  cosmosim::gravity::TreePmCoordinator dist_coordinator(pm_shape, layout);
  cosmosim::gravity::TreePmDiagnostics dist_diag;

  const int iterations = communication_stress ? 4 : 1;
  for (int i = 0; i < iterations; ++i) {
    if (communication_stress) {
      dist_coordinator.solveActiveSetWithPmCadence(
          local_x,
          local_y,
          local_z,
          local_mass,
          dist_acc,
          options,
          (i % 2) == 0,
          nullptr,
          &dist_diag);
    } else {
      dist_coordinator.solveActiveSet(local_x, local_y, local_z, local_mass, dist_acc, options, nullptr, &dist_diag);
    }
  }

  std::vector<std::uint32_t> full_active(pos_x.size(), 0U);
  for (std::size_t i = 0; i < full_active.size(); ++i) {
    full_active[i] = static_cast<std::uint32_t>(i);
  }
  std::vector<double> ref_ax(pos_x.size(), 0.0);
  std::vector<double> ref_ay(pos_x.size(), 0.0);
  std::vector<double> ref_az(pos_x.size(), 0.0);
  cosmosim::gravity::TreePmForceAccumulatorView ref_acc{full_active, ref_ax, ref_ay, ref_az};
  cosmosim::gravity::TreePmCoordinator ref_coordinator(pm_shape);
  ref_coordinator.solveActiveSet(pos_x, pos_y, pos_z, mass, ref_acc, options, nullptr, nullptr);

  const VectorErrorNorm local_norm =
      accumulateVectorError(owned, dist_ax, dist_ay, dist_az, ref_ax, ref_ay, ref_az);

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
  msg << "distributed TreePM equivalence failed: world_size=" << world_size
      << ", communication_stress=" << (communication_stress ? "true" : "false")
      << ", rel_l2=" << rel_l2 << " (required <= 5e-6)"
      << ", max_rel=" << global_norm.max_rel << " (required <= 5e-5)"
      << ", residual_pair_skips_cutoff=" << dist_diag.residual_pair_skips_cutoff;
  requireOrThrow(rel_l2 <= 5.0e-6, msg.str());
  requireOrThrow(global_norm.max_rel <= 5.0e-5, msg.str());
  if (communication_stress) {
    requireOrThrow(dist_diag.residual_pair_skips_cutoff > 0, msg.str());
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
  stream << "hierarchical_max_rung = 1\n";
  stream << "treepm_pm_grid = 16\n";
  stream << "treepm_asmth_cells = 1.25\n";
  stream << "treepm_rcut_cells = 4.5\n";
  stream << "treepm_update_cadence_steps = 2\n";
  stream << "treepm_tree_exchange_batch_bytes = 256\n\n";
  stream << "[output]\n";
  stream << "run_name = validation_phase2_restart\n";
  stream << "output_directory = integration_outputs\n";
  stream << "output_stem = snapshot\n";
  stream << "restart_stem = restart\n";
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
      std::vector<double> density(grid.cellCount(), 0.0);
      (void)solver.solvePoissonPeriodic(grid, density, options, nullptr);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    requireOrThrow(threw, "communicator/layout mismatch must throw in PM solve");
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
  runTreePmCase(world_size, world_rank, /*communication_stress=*/false);
  runTreePmCase(world_size, world_rank, /*communication_stress=*/true);
  testRestartRoundtripContinuationContract(world_size, world_rank);
  testExplicitFailureContracts(world_size, world_rank);

#if COSMOSIM_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
