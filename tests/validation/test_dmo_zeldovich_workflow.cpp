#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;
constexpr std::size_t k_lattice_n = 4U;
constexpr std::size_t k_particle_count = k_lattice_n * k_lattice_n * k_lattice_n;
constexpr double k_box_size_mpc_comoving = 1.0;
constexpr double k_initial_scale_factor = 0.1;
constexpr double k_linear_density_amplitude = 5.0e-3;
constexpr double k_step_dt_code = 5.0e-6;
constexpr std::uint64_t k_direct_steps = 2U;
constexpr std::size_t k_power_spectrum_mesh_n = 12U;
constexpr std::size_t k_power_spectrum_bin_count = 32U;

struct ParallelRuntime {
  int world_size = 1;
  int world_rank = 0;
};

[[nodiscard]] ParallelRuntime parallelRuntime() {
  ParallelRuntime runtime;
#if COSMOSIM_ENABLE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &runtime.world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &runtime.world_rank);
#endif
  return runtime;
}

void barrier() {
#if COSMOSIM_ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

[[nodiscard]] double globalSum(double local_value) {
#if COSMOSIM_ENABLE_MPI
  double global_value = 0.0;
  MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global_value;
#else
  return local_value;
#endif
}

[[nodiscard]] std::uint64_t globalSum(std::uint64_t local_value) {
#if COSMOSIM_ENABLE_MPI
  std::uint64_t global_value = 0U;
  MPI_Allreduce(&local_value, &global_value, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  return global_value;
#else
  return local_value;
#endif
}

[[nodiscard]] std::uint64_t globalXor(std::uint64_t local_value) {
#if COSMOSIM_ENABLE_MPI
  std::uint64_t global_value = 0U;
  MPI_Allreduce(&local_value, &global_value, 1, MPI_UINT64_T, MPI_BXOR, MPI_COMM_WORLD);
  return global_value;
#else
  return local_value;
#endif
}

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

[[maybe_unused]] void requireGlobally(
    bool local_condition,
    const ParallelRuntime& runtime,
    const std::string& message) {
  const std::uint64_t passing_ranks = globalSum(
      local_condition ? std::uint64_t{1} : std::uint64_t{0});
  requireOrThrow(
      passing_ranks == static_cast<std::uint64_t>(runtime.world_size),
      message + "; passing_ranks=" + std::to_string(passing_ranks) +
          "/" + std::to_string(runtime.world_size));
}

[[nodiscard]] double wrapPeriodic(double value, double box_size) {
  double wrapped = std::fmod(value, box_size);
  if (wrapped < 0.0) {
    wrapped += box_size;
  }
  if (wrapped >= box_size) {
    wrapped -= box_size;
  }
  return wrapped;
}

[[nodiscard]] double minimumImage(double delta, double box_size) {
  return delta - std::round(delta / box_size) * box_size;
}

[[nodiscard]] cosmosim::core::SimulationState gatherParticlesByStableId(
    const cosmosim::core::SimulationState& local_state,
    const ParallelRuntime& runtime,
    bool contribute_local_state) {
  const std::size_t local_size = contribute_local_state ? local_state.particles.size() : 0U;
  requireOrThrow(
      local_size <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
      "DMO spectrum gather local particle count exceeds MPI int range");

#if COSMOSIM_ENABLE_MPI
  const int local_count = static_cast<int>(local_size);
  std::vector<int> counts(runtime.world_rank == 0 ? static_cast<std::size_t>(runtime.world_size) : 0U, 0);
  MPI_Gather(
      &local_count,
      1,
      MPI_INT,
      runtime.world_rank == 0 ? counts.data() : nullptr,
      1,
      MPI_INT,
      0,
      MPI_COMM_WORLD);

  std::vector<int> displacements;
  std::size_t global_size = 0U;
  if (runtime.world_rank == 0) {
    displacements.resize(counts.size(), 0);
    for (std::size_t rank = 0; rank < counts.size(); ++rank) {
      requireOrThrow(counts[rank] >= 0, "DMO spectrum gather received a negative count");
      requireOrThrow(
          global_size <= static_cast<std::size_t>(std::numeric_limits<int>::max() - counts[rank]),
          "DMO spectrum gather global displacement exceeds MPI int range");
      displacements[rank] = static_cast<int>(global_size);
      global_size += static_cast<std::size_t>(counts[rank]);
    }
  }

  std::vector<double> gathered_x(global_size, 0.0);
  std::vector<double> gathered_y(global_size, 0.0);
  std::vector<double> gathered_z(global_size, 0.0);
  std::vector<double> gathered_vx(global_size, 0.0);
  std::vector<double> gathered_vy(global_size, 0.0);
  std::vector<double> gathered_vz(global_size, 0.0);
  std::vector<double> gathered_mass(global_size, 0.0);
  std::vector<std::uint8_t> gathered_time_bin(global_size, 0U);
  std::vector<std::uint64_t> gathered_id(global_size, 0U);
  const auto gather_double_lane = [&](std::span<const double> lane, std::vector<double>& gathered) {
    MPI_Gatherv(
        local_count > 0 ? lane.data() : nullptr,
        local_count,
        MPI_DOUBLE,
        runtime.world_rank == 0 ? gathered.data() : nullptr,
        runtime.world_rank == 0 ? counts.data() : nullptr,
        runtime.world_rank == 0 ? displacements.data() : nullptr,
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD);
  };
  gather_double_lane(local_state.particles.position_x_comoving, gathered_x);
  gather_double_lane(local_state.particles.position_y_comoving, gathered_y);
  gather_double_lane(local_state.particles.position_z_comoving, gathered_z);
  gather_double_lane(local_state.particles.velocity_x_peculiar, gathered_vx);
  gather_double_lane(local_state.particles.velocity_y_peculiar, gathered_vy);
  gather_double_lane(local_state.particles.velocity_z_peculiar, gathered_vz);
  gather_double_lane(local_state.particles.mass_code, gathered_mass);
  MPI_Gatherv(
      local_count > 0 ? local_state.particles.time_bin.data() : nullptr,
      local_count,
      MPI_UINT8_T,
      runtime.world_rank == 0 ? gathered_time_bin.data() : nullptr,
      runtime.world_rank == 0 ? counts.data() : nullptr,
      runtime.world_rank == 0 ? displacements.data() : nullptr,
      MPI_UINT8_T,
      0,
      MPI_COMM_WORLD);
  MPI_Gatherv(
      local_count > 0 ? local_state.particle_sidecar.particle_id.data() : nullptr,
      local_count,
      MPI_UINT64_T,
      runtime.world_rank == 0 ? gathered_id.data() : nullptr,
      runtime.world_rank == 0 ? counts.data() : nullptr,
      runtime.world_rank == 0 ? displacements.data() : nullptr,
      MPI_UINT64_T,
      0,
      MPI_COMM_WORLD);
#else
  (void)runtime;
  const std::size_t global_size = local_size;
  std::vector<double> gathered_x;
  std::vector<double> gathered_y;
  std::vector<double> gathered_z;
  std::vector<double> gathered_vx;
  std::vector<double> gathered_vy;
  std::vector<double> gathered_vz;
  std::vector<double> gathered_mass;
  std::vector<std::uint8_t> gathered_time_bin;
  std::vector<std::uint64_t> gathered_id;
  if (contribute_local_state) {
    gathered_x.assign(
        local_state.particles.position_x_comoving.begin(),
        local_state.particles.position_x_comoving.end());
    gathered_y.assign(
        local_state.particles.position_y_comoving.begin(),
        local_state.particles.position_y_comoving.end());
    gathered_z.assign(
        local_state.particles.position_z_comoving.begin(),
        local_state.particles.position_z_comoving.end());
    gathered_vx.assign(
        local_state.particles.velocity_x_peculiar.begin(),
        local_state.particles.velocity_x_peculiar.end());
    gathered_vy.assign(
        local_state.particles.velocity_y_peculiar.begin(),
        local_state.particles.velocity_y_peculiar.end());
    gathered_vz.assign(
        local_state.particles.velocity_z_peculiar.begin(),
        local_state.particles.velocity_z_peculiar.end());
    gathered_mass.assign(
        local_state.particles.mass_code.begin(),
        local_state.particles.mass_code.end());
    gathered_time_bin.assign(
        local_state.particles.time_bin.begin(),
        local_state.particles.time_bin.end());
    gathered_id.assign(
        local_state.particle_sidecar.particle_id.begin(),
        local_state.particle_sidecar.particle_id.end());
  }
#endif

  cosmosim::core::SimulationState gathered_state;
  if (runtime.world_rank != 0) {
    return gathered_state;
  }
  requireOrThrow(global_size > 0U, "DMO spectrum gather produced no particles on the root rank");
  std::vector<std::size_t> order(global_size, 0U);
  std::iota(order.begin(), order.end(), 0U);
  std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    return gathered_id[lhs] < gathered_id[rhs];
  });

  gathered_state.resizeParticles(global_size);
  gathered_state.species.count_by_species.fill(0U);
  gathered_state.species.count_by_species[
      static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] = global_size;
  for (std::size_t destination = 0; destination < global_size; ++destination) {
    const std::size_t source = order[destination];
    requireOrThrow(
        destination == 0U || gathered_id[source] != gathered_state.particle_sidecar.particle_id[destination - 1U],
        "DMO spectrum gather found duplicate stable particle IDs");
    gathered_state.particles.position_x_comoving[destination] = gathered_x[source];
    gathered_state.particles.position_y_comoving[destination] = gathered_y[source];
    gathered_state.particles.position_z_comoving[destination] = gathered_z[source];
    gathered_state.particles.velocity_x_peculiar[destination] = gathered_vx[source];
    gathered_state.particles.velocity_y_peculiar[destination] = gathered_vy[source];
    gathered_state.particles.velocity_z_peculiar[destination] = gathered_vz[source];
    gathered_state.particles.mass_code[destination] = gathered_mass[source];
    gathered_state.particles.time_bin[destination] = gathered_time_bin[source];
    gathered_state.particle_sidecar.particle_id[destination] = gathered_id[source];
    gathered_state.particle_sidecar.species_tag[destination] =
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
  }
  gathered_state.rebuildSpeciesIndex();
  requireOrThrow(
      gathered_state.validateUniqueParticleIds(),
      "DMO spectrum gather did not produce unique stable particle IDs");
  return gathered_state;
}

[[nodiscard, maybe_unused]] std::uint64_t xorRangeOneToN(std::uint64_t n) {
  switch (n & 3U) {
    case 0U: return n;
    case 1U: return 1U;
    case 2U: return n + 1U;
    default: return 0U;
  }
}

[[nodiscard]] std::string configText(int world_size, std::string_view run_name) {
  std::ostringstream stream;
  stream.precision(17);
  stream << "schema_version = 1\n\n";
  stream << "[units]\n";
  stream << "length_unit = mpc\n";
  stream << "mass_unit = msun\n";
  stream << "velocity_unit = km_s\n";
  stream << "coordinate_frame = comoving\n\n";
  stream << "[mode]\n";
  stream << "mode = cosmo_cube\n";
  stream << "ic_file = generated\n";
  stream << "hydro_boundary = periodic\n";
  stream << "gravity_boundary = periodic\n\n";
  stream << "[cosmology]\n";
  stream << "omega_matter = 0.315\n";
  stream << "omega_lambda = 0.685\n";
  stream << "omega_baryon = 0.049\n";
  stream << "hubble_param = 0.674\n";
  stream << "box_size_x = " << k_box_size_mpc_comoving << " mpc\n";
  stream << "box_size_y = " << k_box_size_mpc_comoving << " mpc\n";
  stream << "box_size_z = " << k_box_size_mpc_comoving << " mpc\n\n";
  stream << "[numerics]\n";
  stream << "a_begin = " << k_initial_scale_factor << "\n";
  stream << "a_end = 0.2\n";
  stream << "t_code_begin = 0.0\n";
  stream << "t_code_end = 0.001\n";
  stream << "integrator_time_variable = scale_factor\n";
  stream << "cosmology_max_delta_ln_a = 0.02\n";
  stream << "cosmology_max_hubble_time_fraction = 0.02\n";
  stream << "max_global_steps = 2\n";
  stream << "hierarchical_max_rung = 0\n";
  stream << "gravity_softening = 0.5 kpc\n";
  stream << "treepm_pm_grid_nx = 16\n";
  stream << "treepm_pm_grid_ny = 16\n";
  stream << "treepm_pm_grid_nz = 16\n";
  stream << "treepm_asmth_cells = 1.25\n";
  stream << "treepm_rcut_cells = 6.25\n";
  stream << "treepm_tree_opening_criterion = com_distance\n";
  stream << "treepm_tree_opening_theta = 0.7\n";
  stream << "treepm_tree_relative_force_tolerance = 0.005\n";
  stream << "treepm_tree_relative_force_acceleration_floor = 1e-30\n";
  stream << "treepm_update_cadence_steps = 1\n";
  stream << "treepm_tree_exchange_batch_bytes = 4096\n\n";
  stream << "[physics]\n";
  stream << "enable_cooling = false\n";
  stream << "enable_star_formation = false\n";
  stream << "enable_feedback = false\n";
  stream << "enable_stellar_evolution = false\n\n";
  stream << "[analysis]\n";
  stream << "power_spectrum_mesh_n = " << k_power_spectrum_mesh_n << "\n";
  stream << "power_spectrum_bin_count = " << k_power_spectrum_bin_count << "\n\n";
  stream << "[output]\n";
  stream << "run_name = " << run_name << "\n";
  stream << "output_directory = validation_outputs\n";
  stream << "output_stem = snapshot\n";
  stream << "restart_stem = restart\n";
  stream << "snapshot_interval_steps = 1\n";
  stream << "write_restarts = true\n\n";
  stream << "[parallel]\n";
  stream << "mpi_ranks_expected = " << world_size << "\n";
  stream << "omp_threads = 1\n";
  stream << "gpu_devices = 0\n";
  stream << "deterministic_reduction = true\n";
  stream << "decomposition_runtime_rebalance_enabled = false\n";
  return stream.str();
}

[[nodiscard]] double initialTotalMassCode() {
  const cosmosim::core::UnitSystem units = cosmosim::core::makeUnitSystem("mpc", "msun", "km_s");
  const cosmosim::core::LambdaCdmBackground background(cosmosim::core::CosmologyBackgroundConfig{
      .hubble_param = 0.674,
      .omega_matter = 0.315,
      .omega_lambda = 0.685,
  });
  const double hubble0_code = background.hubble0Si() * units.timeSiPerCode();
  const double newton_g_code = cosmosim::core::newtonGravitationalConstantCode(units);
  requireOrThrow(
      std::isfinite(newton_g_code) && newton_g_code > 0.0,
      "DMO initial condition derived an invalid code-unit Newton constant");
  // Friedmann mean matter density in the configured Mpc/Msun/(km/s) code units.
  const double mean_matter_density_code =
      3.0 * hubble0_code * hubble0_code * background.config().omega_matter /
      (8.0 * k_pi * newton_g_code);
  return mean_matter_density_code *
      k_box_size_mpc_comoving * k_box_size_mpc_comoving * k_box_size_mpc_comoving;
}

[[nodiscard, maybe_unused]] std::array<double, 3> initialPositionForId(std::uint64_t particle_id) {
  requireOrThrow(
      particle_id >= 1U && particle_id <= k_particle_count,
      "Zel'dovich stable particle ID is outside the deterministic lattice");
  const std::size_t linear_index = static_cast<std::size_t>(particle_id - 1U);
  const std::size_t ix = linear_index / (k_lattice_n * k_lattice_n);
  const std::size_t iy = (linear_index / k_lattice_n) % k_lattice_n;
  const std::size_t iz = linear_index % k_lattice_n;
  const double spacing = k_box_size_mpc_comoving / static_cast<double>(k_lattice_n);
  const double qx = (static_cast<double>(ix) + 0.5) * spacing;
  const double qy = (static_cast<double>(iy) + 0.5) * spacing;
  const double qz = (static_cast<double>(iz) + 0.5) * spacing;
  const double wave_number = 2.0 * k_pi / k_box_size_mpc_comoving;
  const double displacement_amplitude = k_linear_density_amplitude / wave_number;
  return {
      wrapPeriodic(qx + displacement_amplitude * std::sin(wave_number * qx), k_box_size_mpc_comoving),
      qy,
      qz,
  };
}

[[nodiscard, maybe_unused]] cosmosim::core::SimulationState makeZeldovichState() {
  namespace core = cosmosim::core;
  core::SimulationState state;
  state.resizeParticles(k_particle_count);
  state.metadata.run_name = "dmo_zeldovich_workflow";
  state.metadata.scale_factor = k_initial_scale_factor;
  state.metadata.snapshot_stem = "snapshot";
  state.metadata.restart_stem = "restart";
  state.species.count_by_species.fill(0U);
  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)] =
      k_particle_count;

  const core::UnitSystem units = core::makeUnitSystem("mpc", "msun", "km_s");
  const core::LambdaCdmBackground background(core::CosmologyBackgroundConfig{
      .hubble_param = 0.674,
      .omega_matter = 0.315,
      .omega_lambda = 0.685,
  });
  const double hubble_code = background.hubbleSi(k_initial_scale_factor) * units.timeSiPerCode();
  const double e2 = std::pow(background.eFactor(k_initial_scale_factor), 2.0);
  const double omega_matter_a =
      background.config().omega_matter /
      (std::pow(k_initial_scale_factor, 3.0) * e2);
  const double logarithmic_growth_rate = std::pow(omega_matter_a, 0.55);
  const double particle_mass_code = initialTotalMassCode() / static_cast<double>(k_particle_count);
  const double wave_number = 2.0 * k_pi / k_box_size_mpc_comoving;
  const double displacement_amplitude = k_linear_density_amplitude / wave_number;

  for (std::size_t ix = 0; ix < k_lattice_n; ++ix) {
    for (std::size_t iy = 0; iy < k_lattice_n; ++iy) {
      for (std::size_t iz = 0; iz < k_lattice_n; ++iz) {
        const std::size_t row = (ix * k_lattice_n + iy) * k_lattice_n + iz;
        const std::uint64_t particle_id = static_cast<std::uint64_t>(row + 1U);
        const std::array<double, 3> position = initialPositionForId(particle_id);
        const double spacing = k_box_size_mpc_comoving / static_cast<double>(k_lattice_n);
        const double qx = (static_cast<double>(ix) + 0.5) * spacing;
        const double displacement_x = displacement_amplitude * std::sin(wave_number * qx);
        state.particles.position_x_comoving[row] = position[0];
        state.particles.position_y_comoving[row] = position[1];
        state.particles.position_z_comoving[row] = position[2];
        // Growing-mode Zel'dovich velocity: u=a dx/dt=a H(a) f(a) Psi.
        state.particles.velocity_x_peculiar[row] =
            k_initial_scale_factor * hubble_code * logarithmic_growth_rate * displacement_x;
        state.particles.velocity_y_peculiar[row] = 0.0;
        state.particles.velocity_z_peculiar[row] = 0.0;
        state.particles.mass_code[row] = particle_mass_code;
        state.particles.time_bin[row] = 0U;
        state.particle_sidecar.particle_id[row] = particle_id;
        state.particle_sidecar.sfc_key[row] = static_cast<std::uint64_t>(row);
        state.particle_sidecar.species_tag[row] =
            static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
        state.particle_sidecar.owning_rank[row] = 0U;
      }
    }
  }
  state.rebuildSpeciesIndex();
  requireOrThrow(state.validateOwnershipInvariants(), "initial Zel'dovich state violates ownership invariants");
  return state;
}

[[nodiscard]] cosmosim::core::SimulationState makeSparseEmptyRankState() {
  namespace core = cosmosim::core;
  constexpr std::size_t sparse_count = 3U;
  core::SimulationState state;
  state.resizeParticles(sparse_count);
  state.metadata.run_name = "dmo_empty_rank_smoke";
  state.metadata.snapshot_stem = "snapshot";
  state.metadata.restart_stem = "restart";
  state.species.count_by_species.fill(0U);
  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)] =
      sparse_count;
  for (std::size_t row = 0; row < sparse_count; ++row) {
    state.particles.position_x_comoving[row] =
        (static_cast<double>(row) + 0.5) / static_cast<double>(sparse_count);
    state.particles.position_y_comoving[row] = 0.5;
    state.particles.position_z_comoving[row] = 0.5;
    state.particles.mass_code[row] = initialTotalMassCode() / static_cast<double>(sparse_count);
    state.particles.time_bin[row] = 0U;
    state.particle_sidecar.particle_id[row] = 10001U + row;
    state.particle_sidecar.sfc_key[row] = row;
    state.particle_sidecar.species_tag[row] =
        static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[row] = 0U;
  }
  state.rebuildSpeciesIndex();
  requireOrThrow(state.validateOwnershipInvariants(), "sparse DMO state violates ownership invariants");
  return state;
}

struct StateSummary {
  double mass_code = 0.0;
  std::array<double, 3> momentum_code{};
  double rms_velocity_code = 0.0;
  std::uint64_t particle_count = 0U;
  std::uint64_t cell_count = 0U;
  std::uint64_t particle_id_sum = 0U;
  std::uint64_t particle_id_xor = 0U;
  bool finite = true;
  bool dm_only = true;
};

[[nodiscard, maybe_unused]] StateSummary summarizeState(
    const cosmosim::core::SimulationState& state,
    const ParallelRuntime& runtime,
    bool contribute_local_state) {
  StateSummary local;
  if (contribute_local_state) {
    local.particle_count = state.particles.size();
    local.cell_count = state.cells.size();
    local.dm_only = state.cells.size() == 0U;
    for (std::size_t row = 0; row < state.particles.size(); ++row) {
      const double mass = state.particles.mass_code[row];
      const double vx = state.particles.velocity_x_peculiar[row];
      const double vy = state.particles.velocity_y_peculiar[row];
      const double vz = state.particles.velocity_z_peculiar[row];
      local.finite = local.finite && std::isfinite(state.particles.position_x_comoving[row]) &&
          std::isfinite(state.particles.position_y_comoving[row]) &&
          std::isfinite(state.particles.position_z_comoving[row]) && std::isfinite(mass) &&
          std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz) && mass > 0.0;
      local.dm_only = local.dm_only &&
          state.particle_sidecar.species_tag[row] ==
              static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
      local.mass_code += mass;
      local.momentum_code[0] += mass * vx;
      local.momentum_code[1] += mass * vy;
      local.momentum_code[2] += mass * vz;
      local.rms_velocity_code += mass * (vx * vx + vy * vy + vz * vz);
      local.particle_id_sum += state.particle_sidecar.particle_id[row];
      local.particle_id_xor ^= state.particle_sidecar.particle_id[row];
    }
  }

  StateSummary global;
  global.mass_code = globalSum(local.mass_code);
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    global.momentum_code[axis] = globalSum(local.momentum_code[axis]);
  }
  const double global_mass_weighted_velocity2 = globalSum(local.rms_velocity_code);
  global.rms_velocity_code = global.mass_code > 0.0
      ? std::sqrt(global_mass_weighted_velocity2 / global.mass_code)
      : 0.0;
  global.particle_count = globalSum(local.particle_count);
  global.cell_count = globalSum(local.cell_count);
  global.particle_id_sum = globalSum(local.particle_id_sum);
  global.particle_id_xor = globalXor(local.particle_id_xor);
  global.finite = globalSum(local.finite ? std::uint64_t{1} : std::uint64_t{0}) ==
      static_cast<std::uint64_t>(runtime.world_size);
  global.dm_only = globalSum(local.dm_only ? std::uint64_t{1} : std::uint64_t{0}) ==
      static_cast<std::uint64_t>(runtime.world_size);
  return global;
}

struct DirectFourierMetrics {
  double fundamental_amplitude = 0.0;
  double fundamental_real = 0.0;
  double fundamental_imag = 0.0;
  double high_k_rms_amplitude = 0.0;
  double high_k_max_amplitude = 0.0;
};

[[nodiscard, maybe_unused]] DirectFourierMetrics directFourierMetrics(
    const cosmosim::core::SimulationState& state,
    bool contribute_local_state) {
  double local_mass = 0.0;
  if (contribute_local_state) {
    for (const double mass : state.particles.mass_code) {
      local_mass += mass;
    }
  }
  const double total_mass = globalSum(local_mass);
  requireOrThrow(total_mass > 0.0 && std::isfinite(total_mass), "direct Fourier estimator has no finite mass");

  std::vector<double> mode_amplitudes;
  std::vector<double> mode_real;
  std::vector<double> mode_imag;
  mode_amplitudes.reserve(k_lattice_n / 2U);
  mode_real.reserve(k_lattice_n / 2U);
  mode_imag.reserve(k_lattice_n / 2U);
  for (std::size_t mode = 1U; mode <= k_lattice_n / 2U; ++mode) {
    double local_real = 0.0;
    double local_imag = 0.0;
    if (contribute_local_state) {
      const double wave_number =
          2.0 * k_pi * static_cast<double>(mode) / k_box_size_mpc_comoving;
      for (std::size_t row = 0; row < state.particles.size(); ++row) {
        const double phase = wave_number * wrapPeriodic(
            state.particles.position_x_comoving[row], k_box_size_mpc_comoving);
        const double mass = state.particles.mass_code[row];
        local_real += mass * std::cos(phase);
        local_imag -= mass * std::sin(phase);
      }
    }
    const double real = globalSum(local_real);
    const double imag = globalSum(local_imag);
    // This is a direct particle sum, independent of the production mesh/FFT
    // diagnostics. The factor two is the real-field one-sided amplitude.
    const double normalization = 2.0 / total_mass;
    mode_real.push_back(normalization * real);
    mode_imag.push_back(normalization * imag);
    mode_amplitudes.push_back(normalization * std::hypot(real, imag));
  }

  DirectFourierMetrics metrics;
  metrics.fundamental_amplitude = mode_amplitudes.front();
  metrics.fundamental_real = mode_real.front();
  metrics.fundamental_imag = mode_imag.front();
  double high_k_square_sum = 0.0;
  for (std::size_t mode_index = 1U; mode_index < mode_amplitudes.size(); ++mode_index) {
    high_k_square_sum += mode_amplitudes[mode_index] * mode_amplitudes[mode_index];
    metrics.high_k_max_amplitude = std::max(metrics.high_k_max_amplitude, mode_amplitudes[mode_index]);
  }
  if (mode_amplitudes.size() > 1U) {
    metrics.high_k_rms_amplitude =
        std::sqrt(high_k_square_sum / static_cast<double>(mode_amplitudes.size() - 1U));
  }
  return metrics;
}

[[nodiscard]] double fourierPhaseDrift(
    const DirectFourierMetrics& reference,
    const DirectFourierMetrics& candidate) {
  const double dot = reference.fundamental_real * candidate.fundamental_real +
      reference.fundamental_imag * candidate.fundamental_imag;
  const double cross = reference.fundamental_real * candidate.fundamental_imag -
      reference.fundamental_imag * candidate.fundamental_real;
  return std::abs(std::atan2(cross, dot));
}

[[nodiscard]] double fourierPhaseCoherence(
    const DirectFourierMetrics& reference,
    const DirectFourierMetrics& candidate) {
  const double denominator =
      reference.fundamental_amplitude * candidate.fundamental_amplitude;
  requireOrThrow(
      denominator > 0.0 && std::isfinite(denominator),
      "Fourier phase-coherence comparison has zero or non-finite mode amplitude");
  return (reference.fundamental_real * candidate.fundamental_real +
          reference.fundamental_imag * candidate.fundamental_imag) /
      denominator;
}

[[nodiscard]] std::size_t validatePowerSpectrumEstimate(
    const cosmosim::analysis::PowerSpectrumEstimate& estimate,
    std::string_view label) {
  using cosmosim::analysis::PowerSpectrumMassAssignment;
  using cosmosim::analysis::PowerSpectrumShotNoisePolicy;
  using cosmosim::analysis::PowerSpectrumWindowCorrection;
  requireOrThrow(
      estimate.options.mesh_n == k_power_spectrum_mesh_n,
      std::string(label) + ": unexpected power-spectrum mesh size");
  requireOrThrow(
      estimate.options.bin_count == k_power_spectrum_bin_count &&
          estimate.bins.size() == k_power_spectrum_bin_count,
      std::string(label) + ": requested power-spectrum bins were not retained");
  requireOrThrow(
      estimate.options.mass_assignment == PowerSpectrumMassAssignment::kCloudInCell &&
          cosmosim::analysis::powerSpectrumMassAssignmentLabel(estimate.options.mass_assignment) == "cic",
      std::string(label) + ": power-spectrum assignment policy is not explicit CIC");
  requireOrThrow(
      estimate.options.window_correction ==
              PowerSpectrumWindowCorrection::kDeconvolveAssignmentWindow &&
          cosmosim::analysis::powerSpectrumWindowCorrectionLabel(
              estimate.options.window_correction) == "deconvolve_assignment_window",
      std::string(label) + ": CIC assignment-window correction is not explicit");
  requireOrThrow(
      estimate.options.shot_noise_policy ==
              PowerSpectrumShotNoisePolicy::kReportWithoutSubtraction &&
          cosmosim::analysis::powerSpectrumShotNoisePolicyLabel(
              estimate.options.shot_noise_policy) == "poisson_reported_not_subtracted",
      std::string(label) + ": shot-noise policy is not report-without-subtraction");
  requireOrThrow(
      cosmosim::analysis::powerSpectrumKUnits() == "code_length^-1" &&
          cosmosim::analysis::powerSpectrumPowerUnits() == "code_length^3",
      std::string(label) + ": power-spectrum units are not the documented code units");
  requireOrThrow(
      cosmosim::analysis::powerSpectrumFourierNormalization() ==
          "delta_k=sum(delta_grid*exp(-i*k*x_grid))/N_mesh;P=V*abs(delta_k)^2",
      std::string(label) + ": Fourier normalization contract changed");
  requireOrThrow(
      estimate.box_size_code == k_box_size_mpc_comoving &&
          std::isfinite(estimate.k_fundamental_code) && estimate.k_fundamental_code > 0.0 &&
          std::isfinite(estimate.k_axis_nyquist_code) && estimate.k_axis_nyquist_code > 0.0,
      std::string(label) + ": invalid box or Fourier scale metadata");
  requireOrThrow(
      std::isfinite(estimate.poisson_shot_noise_code_volume) &&
          estimate.poisson_shot_noise_code_volume > 0.0,
      std::string(label) + ": invalid reported Poisson shot-noise level");

  std::uint64_t total_mode_count = 0U;
  std::size_t empty_bin_count = 0U;
  double previous_upper = 0.0;
  for (std::size_t bin_index = 0; bin_index < estimate.bins.size(); ++bin_index) {
    const auto& bin = estimate.bins[bin_index];
    requireOrThrow(bin.bin_index == bin_index, std::string(label) + ": unstable bin index");
    requireOrThrow(
        std::isfinite(bin.k_lower_code) && std::isfinite(bin.k_upper_code) &&
            std::isfinite(bin.k_center_code) && bin.k_lower_code == previous_upper &&
            bin.k_lower_code < bin.k_upper_code && bin.k_center_code >= bin.k_lower_code &&
            bin.k_center_code <= bin.k_upper_code,
        std::string(label) + ": invalid deterministic bin geometry");
    requireOrThrow(
        bin.empty == (bin.mode_count == 0U),
        std::string(label) + ": empty-bin flag disagrees with mode count");
    requireOrThrow(
        std::isfinite(bin.power_code_volume),
        std::string(label) + ": non-finite power-spectrum bin");
    if (bin.empty) {
      ++empty_bin_count;
      requireOrThrow(
          bin.power_code_volume == 0.0,
          std::string(label) + ": empty bin does not carry the explicit zero sentinel");
    } else {
      requireOrThrow(
          bin.power_code_volume >= 0.0,
          std::string(label) + ": unsubtracted spectrum contains negative power");
    }
    total_mode_count += bin.mode_count;
    previous_upper = bin.k_upper_code;
  }
  const std::uint64_t expected_mode_count =
      static_cast<std::uint64_t>(k_power_spectrum_mesh_n) *
          static_cast<std::uint64_t>(k_power_spectrum_mesh_n) *
          static_cast<std::uint64_t>(k_power_spectrum_mesh_n) -
      1U;
  requireOrThrow(
      total_mode_count == expected_mode_count,
      std::string(label) + ": discrete mode counts do not cover mesh^3-1 modes");
  requireOrThrow(empty_bin_count > 0U, std::string(label) + ": fixture did not exercise empty bins");
  return empty_bin_count;
}

[[nodiscard]] const cosmosim::analysis::PowerSpectrumEstimateBin& fundamentalPowerBin(
    const cosmosim::analysis::PowerSpectrumEstimate& estimate) {
  const double fundamental = estimate.k_fundamental_code;
  for (std::size_t bin_index = 0; bin_index < estimate.bins.size(); ++bin_index) {
    const auto& bin = estimate.bins[bin_index];
    const bool is_last = bin_index + 1U == estimate.bins.size();
    if (fundamental >= bin.k_lower_code &&
        (fundamental < bin.k_upper_code || (is_last && fundamental <= bin.k_upper_code))) {
      requireOrThrow(!bin.empty, "fundamental power-spectrum shell is unexpectedly empty");
      return bin;
    }
  }
  throw std::runtime_error("fundamental Fourier mode is outside the power-spectrum bins");
}

void requirePowerSpectraEquivalent(
    const cosmosim::analysis::PowerSpectrumEstimate& lhs,
    const cosmosim::analysis::PowerSpectrumEstimate& rhs,
    std::string_view label) {
  requireOrThrow(lhs.bins.size() == rhs.bins.size(), std::string(label) + ": bin-count mismatch");
  for (std::size_t bin_index = 0; bin_index < lhs.bins.size(); ++bin_index) {
    const auto& lhs_bin = lhs.bins[bin_index];
    const auto& rhs_bin = rhs.bins[bin_index];
    requireOrThrow(
        lhs_bin.mode_count == rhs_bin.mode_count && lhs_bin.empty == rhs_bin.empty &&
            lhs_bin.k_lower_code == rhs_bin.k_lower_code &&
            lhs_bin.k_upper_code == rhs_bin.k_upper_code &&
            lhs_bin.k_center_code == rhs_bin.k_center_code,
        std::string(label) + ": bin geometry or mode-count mismatch");
    const double power_scale = std::max({
        std::abs(lhs_bin.power_code_volume),
        std::abs(rhs_bin.power_code_volume),
        1.0e-30,
    });
    requireOrThrow(
        std::abs(lhs_bin.power_code_volume - rhs_bin.power_code_volume) <= 1.0e-12 * power_scale,
        std::string(label) + ": bin power mismatch");
  }
}

[[nodiscard]] std::string serializePowerSpectrumArtifact(
    const ParallelRuntime& runtime,
    double evolved_scale_factor,
    double expected_linear_growth,
    double measured_power_growth,
    double initial_recovered_amplitude,
    const cosmosim::analysis::PowerSpectrumEstimate& initial,
    const cosmosim::analysis::PowerSpectrumEstimate& evolved) {
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  out << "{\n";
  out << "  \"schema_version\": \"dmo_zeldovich_power_spectrum_v1\",\n";
  out << "  \"world_size\": " << runtime.world_size << ",\n";
  out << "  \"mesh\": {\"nx\": " << initial.options.mesh_n
      << ", \"ny\": " << initial.options.mesh_n
      << ", \"nz\": " << initial.options.mesh_n
      << ", \"requested_bin_count\": " << initial.options.bin_count
      << ", \"box_size_code\": " << initial.box_size_code
      << ", \"k_fundamental_code\": " << initial.k_fundamental_code
      << ", \"k_axis_nyquist_code\": " << initial.k_axis_nyquist_code << "},\n";
  out << "  \"estimator\": {\n";
  out << "    \"mass_assignment\": \""
      << cosmosim::analysis::powerSpectrumMassAssignmentLabel(initial.options.mass_assignment)
      << "\",\n";
  out << "    \"window_correction\": \""
      << cosmosim::analysis::powerSpectrumWindowCorrectionLabel(initial.options.window_correction)
      << "\",\n";
  out << "    \"fourier_normalization\": \""
      << cosmosim::analysis::powerSpectrumFourierNormalization() << "\",\n";
  out << "    \"binning\": \"linear_k_to_3d_corner_dc_excluded\",\n";
  out << "    \"mode_count_convention\": \"full_discrete_mesh_modes\",\n";
  out << "    \"empty_bin_policy\": \"retain_with_null_power\",\n";
  out << "    \"k_units\": \"" << cosmosim::analysis::powerSpectrumKUnits() << "\",\n";
  out << "    \"power_units\": \"" << cosmosim::analysis::powerSpectrumPowerUnits() << "\",\n";
  out << "    \"shot_noise_policy\": \""
      << cosmosim::analysis::powerSpectrumShotNoisePolicyLabel(initial.options.shot_noise_policy)
      << "\",\n";
  out << "    \"poisson_shot_noise_code_volume\": "
      << initial.poisson_shot_noise_code_volume << "\n";
  out << "  },\n";
  out << "  \"linear_growth\": {\n";
  out << "    \"initial_scale_factor\": " << k_initial_scale_factor << ",\n";
  out << "    \"evolved_scale_factor\": " << evolved_scale_factor << ",\n";
  out << "    \"expected_scale_factor_growth\": " << expected_linear_growth << ",\n";
  out << "    \"measured_sqrt_power_growth\": " << measured_power_growth << ",\n";
  out << "    \"initial_recovered_density_amplitude\": " << initial_recovered_amplitude << "\n";
  out << "  },\n";

  const auto write_spectrum = [&](std::string_view name,
                                  const cosmosim::analysis::PowerSpectrumEstimate& estimate,
                                  bool trailing_comma) {
    out << "  \"" << name << "\": {\n";
    out << "    \"scale_factor\": "
        << (name == "initial_spectrum" ? k_initial_scale_factor : evolved_scale_factor) << ",\n";
    out << "    \"bins\": [\n";
    for (std::size_t bin_index = 0; bin_index < estimate.bins.size(); ++bin_index) {
      const auto& bin = estimate.bins[bin_index];
      out << "      {\"bin_index\": " << bin.bin_index
          << ", \"k_lower_code\": " << bin.k_lower_code
          << ", \"k_upper_code\": " << bin.k_upper_code
          << ", \"k_center_code\": " << bin.k_center_code
          << ", \"power_code_volume\": ";
      if (bin.empty) {
        out << "null";
      } else {
        out << bin.power_code_volume;
      }
      out << ", \"mode_count\": " << bin.mode_count
          << ", \"empty\": " << (bin.empty ? "true" : "false") << "}";
      if (bin_index + 1U < estimate.bins.size()) {
        out << ',';
      }
      out << '\n';
    }
    out << "    ]\n";
    out << "  }" << (trailing_comma ? "," : "") << "\n";
  };
  write_spectrum("initial_spectrum", initial, true);
  write_spectrum("evolved_spectrum", evolved, false);
  out << "}\n";
  return out.str();
}

void writeAndVerifyPowerSpectrumArtifact(
    const std::filesystem::path& path,
    const std::string& serialized) {
  std::filesystem::create_directories(path.parent_path());
  {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    requireOrThrow(output.is_open(), "failed opening DMO power-spectrum artifact for writing");
    output.write(serialized.data(), static_cast<std::streamsize>(serialized.size()));
    requireOrThrow(output.good(), "failed writing DMO power-spectrum artifact");
  }
  std::ifstream input(path, std::ios::binary);
  requireOrThrow(input.is_open(), "failed reopening DMO power-spectrum artifact");
  const std::string roundtrip(
      (std::istreambuf_iterator<char>(input)),
      std::istreambuf_iterator<char>());
  requireOrThrow(
      roundtrip == serialized,
      "DMO power-spectrum artifact byte round-trip is not reproducible");
}

constexpr std::string_view k_physical_state_artifact_schema =
    "dmo_zeldovich_physical_state_v2";

[[nodiscard]] constexpr std::string_view currentExecutionBackend() noexcept {
#if COSMOSIM_ENABLE_MPI
  return "mpi";
#else
  return "serial";
#endif
}

struct PhysicalStateArtifact {
  std::string execution_backend;
  int world_size = 0;
  double current_time_code = 0.0;
  double current_scale_factor = 0.0;
  std::uint64_t step_index = 0U;
  cosmosim::core::SimulationState state;
};

[[nodiscard]] std::filesystem::path dmoWorkflowRoot(
    std::string_view execution_backend,
    int world_size) {
  requireOrThrow(
      execution_backend == "serial" || execution_backend == "mpi",
      "unsupported DMO execution-backend tag");
  requireOrThrow(world_size >= 1 && world_size <= 4, "invalid DMO artifact world size");
  if (execution_backend == "serial") {
    requireOrThrow(world_size == 1, "serial DMO artifact must have world_size=1");
    return std::filesystem::temp_directory_path() /
        "cosmosim_dmo_zeldovich_workflow_serial";
  }
  return std::filesystem::temp_directory_path() /
      ("cosmosim_dmo_zeldovich_workflow_mpi_np" + std::to_string(world_size));
}

[[nodiscard]] std::filesystem::path physicalStateArtifactPath(
    std::string_view execution_backend,
    int world_size) {
  return dmoWorkflowRoot(execution_backend, world_size) /
      "dmo_zeldovich_physical_state.tsv";
}

void writePhysicalStateArtifact(
    const std::filesystem::path& path,
    const PhysicalStateArtifact& artifact) {
  requireOrThrow(
      artifact.execution_backend == "serial" || artifact.execution_backend == "mpi",
      "invalid DMO artifact execution backend");
  requireOrThrow(artifact.world_size >= 1 && artifact.world_size <= 4, "invalid DMO artifact world size");
  requireOrThrow(
      artifact.execution_backend == "mpi" || artifact.world_size == 1,
      "serial DMO artifact must have world_size=1");
  requireOrThrow(
      std::isfinite(artifact.current_time_code) &&
          std::isfinite(artifact.current_scale_factor) &&
          artifact.current_scale_factor > 0.0,
      "invalid DMO artifact integration epoch");
  requireOrThrow(
      artifact.state.particles.isConsistent() && artifact.state.particle_sidecar.isConsistent() &&
          artifact.state.validateUniqueParticleIds(),
      "invalid DMO artifact particle state");

  std::filesystem::create_directories(path.parent_path());
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  requireOrThrow(output.is_open(), "failed opening DMO physical-state artifact for writing");
  output << std::setprecision(std::numeric_limits<double>::max_digits10);
  output << k_physical_state_artifact_schema << '\n';
  output << "execution_backend " << artifact.execution_backend << '\n';
  output << "world_size " << artifact.world_size << '\n';
  output << "current_time_code " << artifact.current_time_code << '\n';
  output << "current_scale_factor " << artifact.current_scale_factor << '\n';
  output << "step_index " << artifact.step_index << '\n';
  output << "particle_count " << artifact.state.particles.size() << '\n';
  output << "columns particle_id position_x_comoving position_y_comoving position_z_comoving "
            "velocity_x_peculiar velocity_y_peculiar velocity_z_peculiar mass_code time_bin\n";
  for (std::size_t row = 0; row < artifact.state.particles.size(); ++row) {
    output << artifact.state.particle_sidecar.particle_id[row] << '\t'
           << artifact.state.particles.position_x_comoving[row] << '\t'
           << artifact.state.particles.position_y_comoving[row] << '\t'
           << artifact.state.particles.position_z_comoving[row] << '\t'
           << artifact.state.particles.velocity_x_peculiar[row] << '\t'
           << artifact.state.particles.velocity_y_peculiar[row] << '\t'
           << artifact.state.particles.velocity_z_peculiar[row] << '\t'
           << artifact.state.particles.mass_code[row] << '\t'
           << static_cast<unsigned int>(artifact.state.particles.time_bin[row]) << '\n';
  }
  requireOrThrow(output.good(), "failed writing DMO physical-state artifact");
}

[[nodiscard]] PhysicalStateArtifact readPhysicalStateArtifact(
    const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  requireOrThrow(
      input.is_open(),
      "missing DMO physical-state artifact: " + path.string() +
          "; run the complete np1-np4 CTest matrix before the equivalence test");

  std::string schema;
  std::getline(input, schema);
  requireOrThrow(schema == k_physical_state_artifact_schema, "unsupported DMO physical-state artifact schema");

  PhysicalStateArtifact artifact;
  std::string key;
  std::size_t particle_count = 0U;
  input >> key >> artifact.execution_backend;
  requireOrThrow(key == "execution_backend", "DMO physical-state artifact is missing execution_backend");
  input >> key >> artifact.world_size;
  requireOrThrow(key == "world_size", "DMO physical-state artifact is missing world_size");
  input >> key >> artifact.current_time_code;
  requireOrThrow(key == "current_time_code", "DMO physical-state artifact is missing current_time_code");
  input >> key >> artifact.current_scale_factor;
  requireOrThrow(key == "current_scale_factor", "DMO physical-state artifact is missing current_scale_factor");
  input >> key >> artifact.step_index;
  requireOrThrow(key == "step_index", "DMO physical-state artifact is missing step_index");
  input >> key >> particle_count;
  requireOrThrow(key == "particle_count", "DMO physical-state artifact is missing particle_count");
  input >> key;
  requireOrThrow(key == "columns", "DMO physical-state artifact is missing its column declaration");
  std::string column_declaration;
  std::getline(input, column_declaration);
  requireOrThrow(
      column_declaration ==
          " particle_id position_x_comoving position_y_comoving position_z_comoving "
          "velocity_x_peculiar velocity_y_peculiar velocity_z_peculiar mass_code time_bin",
      "DMO physical-state artifact column declaration changed");

  artifact.state.resizeParticles(particle_count);
  artifact.state.species.count_by_species.fill(0U);
  artifact.state.species.count_by_species[
      static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] = particle_count;
  for (std::size_t row = 0; row < particle_count; ++row) {
    unsigned int time_bin = 0U;
    input >> artifact.state.particle_sidecar.particle_id[row]
          >> artifact.state.particles.position_x_comoving[row]
          >> artifact.state.particles.position_y_comoving[row]
          >> artifact.state.particles.position_z_comoving[row]
          >> artifact.state.particles.velocity_x_peculiar[row]
          >> artifact.state.particles.velocity_y_peculiar[row]
          >> artifact.state.particles.velocity_z_peculiar[row]
          >> artifact.state.particles.mass_code[row]
          >> time_bin;
    requireOrThrow(input.good(), "truncated DMO physical-state artifact particle row");
    requireOrThrow(
        time_bin <= static_cast<unsigned int>(std::numeric_limits<std::uint8_t>::max()),
        "DMO physical-state artifact time bin exceeds uint8 range");
    artifact.state.particles.time_bin[row] = static_cast<std::uint8_t>(time_bin);
    artifact.state.particle_sidecar.species_tag[row] =
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    if (row > 0U) {
      requireOrThrow(
          artifact.state.particle_sidecar.particle_id[row - 1U] <
              artifact.state.particle_sidecar.particle_id[row],
          "DMO physical-state artifact is not strictly sorted by stable particle ID");
    }
  }
  std::string trailing_token;
  requireOrThrow(!(input >> trailing_token), "DMO physical-state artifact contains trailing data");
  requireOrThrow(
      (artifact.execution_backend == "serial" || artifact.execution_backend == "mpi") &&
          artifact.world_size >= 1 && artifact.world_size <= 4 &&
          (artifact.execution_backend == "mpi" || artifact.world_size == 1) &&
          std::isfinite(artifact.current_time_code) &&
          std::isfinite(artifact.current_scale_factor) &&
          artifact.current_scale_factor > 0.0 &&
          artifact.state.validateUniqueParticleIds(),
      "DMO physical-state artifact failed validation");
  return artifact;
}

struct PhysicalStateEquivalenceMetrics {
  double max_periodic_position_error_code = 0.0;
  double max_velocity_error_code = 0.0;
  double max_mass_relative_error = 0.0;
};

[[nodiscard]] PhysicalStateEquivalenceMetrics requirePhysicalStateEquivalent(
    const PhysicalStateArtifact& reference,
    const PhysicalStateArtifact& candidate,
    std::string_view comparison_label) {
  requireOrThrow(
      reference.current_time_code == candidate.current_time_code &&
          reference.current_scale_factor == candidate.current_scale_factor &&
          reference.step_index == candidate.step_index,
      "np1/npN DMO integrator epochs differ");
  requireOrThrow(
      reference.state.particles.size() == candidate.state.particles.size(),
      "np1/npN DMO particle counts differ");

  double reference_velocity_scale = 0.0;
  for (std::size_t row = 0; row < reference.state.particles.size(); ++row) {
    reference_velocity_scale = std::max({
        reference_velocity_scale,
        std::abs(reference.state.particles.velocity_x_peculiar[row]),
        std::abs(reference.state.particles.velocity_y_peculiar[row]),
        std::abs(reference.state.particles.velocity_z_peculiar[row]),
    });
  }
  requireOrThrow(
      reference_velocity_scale > 0.0 && std::isfinite(reference_velocity_scale),
      "np1 DMO artifact has no finite nonzero velocity scale");

  PhysicalStateEquivalenceMetrics metrics;
  for (std::size_t row = 0; row < reference.state.particles.size(); ++row) {
    requireOrThrow(
        reference.state.particle_sidecar.particle_id[row] ==
            candidate.state.particle_sidecar.particle_id[row],
        "np1/npN DMO stable particle-ID ordering differs");
    requireOrThrow(
        reference.state.particles.time_bin[row] == candidate.state.particles.time_bin[row],
        "np1/npN DMO particle time bins differ");
    const std::array<double, 3> reference_position{
        reference.state.particles.position_x_comoving[row],
        reference.state.particles.position_y_comoving[row],
        reference.state.particles.position_z_comoving[row],
    };
    const std::array<double, 3> candidate_position{
        candidate.state.particles.position_x_comoving[row],
        candidate.state.particles.position_y_comoving[row],
        candidate.state.particles.position_z_comoving[row],
    };
    const std::array<double, 3> reference_velocity{
        reference.state.particles.velocity_x_peculiar[row],
        reference.state.particles.velocity_y_peculiar[row],
        reference.state.particles.velocity_z_peculiar[row],
    };
    const std::array<double, 3> candidate_velocity{
        candidate.state.particles.velocity_x_peculiar[row],
        candidate.state.particles.velocity_y_peculiar[row],
        candidate.state.particles.velocity_z_peculiar[row],
    };
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      metrics.max_periodic_position_error_code = std::max(
          metrics.max_periodic_position_error_code,
          std::abs(minimumImage(
              candidate_position[axis] - reference_position[axis],
              k_box_size_mpc_comoving)));
      metrics.max_velocity_error_code = std::max(
          metrics.max_velocity_error_code,
          std::abs(candidate_velocity[axis] - reference_velocity[axis]));
    }
    const double mass_scale = std::max(
        std::abs(reference.state.particles.mass_code[row]),
        std::numeric_limits<double>::min());
    metrics.max_mass_relative_error = std::max(
        metrics.max_mass_relative_error,
        std::abs(candidate.state.particles.mass_code[row] -
                 reference.state.particles.mass_code[row]) /
            mass_scale);
  }

  constexpr double k_position_tolerance_code = 1.0e-13;
  constexpr double k_velocity_relative_tolerance = 1.0e-12;
  constexpr double k_mass_relative_tolerance = 1.0e-14;
  std::ostringstream message;
  message << std::setprecision(17)
          << comparison_label << " DMO stable-ID physical-state mismatch: "
          << "max_periodic_position_error_code=" << metrics.max_periodic_position_error_code
          << ", max_velocity_error_code=" << metrics.max_velocity_error_code
          << ", velocity_scale_code=" << reference_velocity_scale
          << ", max_mass_relative_error=" << metrics.max_mass_relative_error;
  requireOrThrow(
      metrics.max_periodic_position_error_code <= k_position_tolerance_code &&
          metrics.max_velocity_error_code <=
              k_velocity_relative_tolerance * reference_velocity_scale &&
          metrics.max_mass_relative_error <= k_mass_relative_tolerance,
      message.str());
  return metrics;
}

void compareRankCountPhysicalStateArtifacts() {
  const PhysicalStateArtifact reference =
      readPhysicalStateArtifact(physicalStateArtifactPath("mpi", 1));
  requireOrThrow(
      reference.execution_backend == "mpi" && reference.world_size == 1,
      "DMO rank-count reference is not the MPI np1 artifact");
  PhysicalStateEquivalenceMetrics worst;
  for (int world_size = 2; world_size <= 4; ++world_size) {
    const PhysicalStateArtifact candidate =
        readPhysicalStateArtifact(physicalStateArtifactPath("mpi", world_size));
    requireOrThrow(
        candidate.execution_backend == "mpi" && candidate.world_size == world_size,
        "DMO physical-state artifact world_size does not match its np-specific path");
    const PhysicalStateEquivalenceMetrics metrics =
        requirePhysicalStateEquivalent(
            reference,
            candidate,
            "MPI np1/np" + std::to_string(world_size));
    worst.max_periodic_position_error_code = std::max(
        worst.max_periodic_position_error_code,
        metrics.max_periodic_position_error_code);
    worst.max_velocity_error_code =
        std::max(worst.max_velocity_error_code, metrics.max_velocity_error_code);
    worst.max_mass_relative_error =
        std::max(worst.max_mass_relative_error, metrics.max_mass_relative_error);
  }
  std::cout << std::setprecision(17)
            << "DMO_RANK_EQUIVALENCE_PASS compared=np1,np2,np3,np4"
            << " stable_id_particles=" << reference.state.particles.size()
            << " max_periodic_position_error_code=" << worst.max_periodic_position_error_code
            << " max_velocity_error_code=" << worst.max_velocity_error_code
            << " max_mass_relative_error=" << worst.max_mass_relative_error << '\n';
}

void compareSerialAndSingleRankMpiPhysicalStateArtifacts() {
  const PhysicalStateArtifact serial =
      readPhysicalStateArtifact(physicalStateArtifactPath("serial", 1));
  const PhysicalStateArtifact mpi =
      readPhysicalStateArtifact(physicalStateArtifactPath("mpi", 1));
  requireOrThrow(
      serial.execution_backend == "serial" && serial.world_size == 1,
      "DMO serial reference artifact has the wrong backend or world size");
  requireOrThrow(
      mpi.execution_backend == "mpi" && mpi.world_size == 1,
      "DMO MPI np1 artifact has the wrong backend or world size");
  const PhysicalStateEquivalenceMetrics metrics =
      requirePhysicalStateEquivalent(serial, mpi, "serial/MPI np1");
  std::cout << std::setprecision(17)
            << "DMO_SERIAL_MPI_EQUIVALENCE_PASS compared=serial,mpi_np1"
            << " stable_id_particles=" << serial.state.particles.size()
            << " max_periodic_position_error_code="
            << metrics.max_periodic_position_error_code
            << " max_velocity_error_code=" << metrics.max_velocity_error_code
            << " max_mass_relative_error=" << metrics.max_mass_relative_error << '\n';
}

[[nodiscard, maybe_unused]] std::array<double, 3> periodicCenterOfMassDrift(
    const cosmosim::core::SimulationState& final_state) {
  std::array<double, 3> local_mass_weighted_drift{};
  double local_mass = 0.0;
  for (std::size_t row = 0; row < final_state.particles.size(); ++row) {
    const std::array<double, 3> initial =
        initialPositionForId(final_state.particle_sidecar.particle_id[row]);
    const std::array<double, 3> final{
        final_state.particles.position_x_comoving[row],
        final_state.particles.position_y_comoving[row],
        final_state.particles.position_z_comoving[row],
    };
    const double mass = final_state.particles.mass_code[row];
    local_mass += mass;
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      local_mass_weighted_drift[axis] +=
          mass * minimumImage(final[axis] - initial[axis], k_box_size_mpc_comoving);
    }
  }
  const double mass = globalSum(local_mass);
  std::array<double, 3> drift{};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    drift[axis] = globalSum(local_mass_weighted_drift[axis]) / mass;
  }
  return drift;
}

void runEmptyRankSmoke(const ParallelRuntime& runtime, const std::filesystem::path& root) {
  const cosmosim::core::SimulationState state = makeSparseEmptyRankState();
  const cosmosim::core::FrozenConfig frozen = cosmosim::core::loadFrozenConfigFromString(
      configText(runtime.world_size, "dmo_empty_rank_smoke"),
      "test_dmo_zeldovich_workflow_empty_rank");
  const cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  const auto report = runner.run(
      root / ("rank_" + std::to_string(runtime.world_rank) + "_empty_rank"),
      cosmosim::workflows::ReferenceWorkflowOptions{
          .step_index = 0U,
          .dt_time_code = k_step_dt_code,
          .write_outputs = false,
          .initial_state_override = &state,
          .restart_state_override = nullptr,
          .initial_particle_scheduler_identity_records = {},
          .max_steps_override = 1U,
      });
  requireOrThrow(report.completed_steps == 1U, "empty-rank DMO smoke did not complete one production step");
  requireOrThrow(report.global_particle_count == 3U, "empty-rank DMO smoke changed global particle count");
  requireOrThrow(report.global_cell_count == 0U, "empty-rank DMO smoke manufactured gas cells");
  requireOrThrow(report.global_particle_partition_identity_match, "empty-rank DMO ownership identity failed");
  const std::uint64_t empty_rank_count = globalSum(
      report.local_particle_count == 0U ? std::uint64_t{1} : std::uint64_t{0});
  if (runtime.world_size == 4) {
    requireOrThrow(empty_rank_count >= 1U, "four-rank sparse DMO fixture did not exercise an empty rank");
  }
}

#if COSMOSIM_ENABLE_HDF5
void runScientificValidation(const ParallelRuntime& runtime, const std::filesystem::path& root) {
  const cosmosim::core::SimulationState initial_state = makeZeldovichState();
  const StateSummary initial_summary = summarizeState(
      initial_state, runtime, runtime.world_rank == 0);
  const DirectFourierMetrics initial_fourier = directFourierMetrics(
      initial_state, runtime.world_rank == 0);
  const cosmosim::core::SimulationState gathered_initial_state = gatherParticlesByStableId(
      initial_state, runtime, runtime.world_rank == 0);

  const cosmosim::core::FrozenConfig frozen = cosmosim::core::loadFrozenConfigFromString(
      configText(runtime.world_size, "dmo_zeldovich_workflow"),
      "test_dmo_zeldovich_workflow");
  const cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);

  const auto direct_report = runner.run(
      root / ("rank_" + std::to_string(runtime.world_rank) + "_direct"),
      cosmosim::workflows::ReferenceWorkflowOptions{
          .step_index = 0U,
          .dt_time_code = k_step_dt_code,
          .write_outputs = true,
          .initial_state_override = &initial_state,
          .restart_state_override = nullptr,
          .initial_particle_scheduler_identity_records = {},
          .max_steps_override = k_direct_steps,
      });
  requireOrThrow(direct_report.completed_steps == k_direct_steps, "direct DMO workflow segment was incomplete");
  requireOrThrow(direct_report.restart_roundtrip_ok, "direct DMO workflow restart round-trip failed");
  requireOrThrow(!direct_report.restart_path.empty(), "direct DMO workflow did not expose its checkpoint path");
  const cosmosim::io::RestartReadResult direct_restart =
      cosmosim::io::readRestartCheckpointHdf5(direct_report.restart_path);

  const auto first_report = runner.run(
      root / ("rank_" + std::to_string(runtime.world_rank) + "_first"),
      cosmosim::workflows::ReferenceWorkflowOptions{
          .step_index = 0U,
          .dt_time_code = k_step_dt_code,
          .write_outputs = true,
          .initial_state_override = &initial_state,
          .restart_state_override = nullptr,
          .initial_particle_scheduler_identity_records = {},
          .max_steps_override = 1U,
      });
  requireOrThrow(first_report.completed_steps == 1U, "first DMO restart segment was incomplete");
  requireOrThrow(first_report.restart_roundtrip_ok, "first DMO restart segment did not round-trip");
  const cosmosim::io::RestartReadResult first_restart =
      cosmosim::io::readRestartCheckpointHdf5(first_report.restart_path);

  const auto resumed_report = runner.run(
      root / ("rank_" + std::to_string(runtime.world_rank) + "_resumed"),
      cosmosim::workflows::ReferenceWorkflowOptions{
          .step_index = 0U,
          .dt_time_code = 0.0,
          .write_outputs = true,
          .initial_state_override = nullptr,
          .restart_state_override = &first_restart,
          .initial_particle_scheduler_identity_records = {},
          .max_steps_override = 1U,
      });
  requireOrThrow(resumed_report.completed_steps == 1U, "resumed DMO workflow segment was incomplete");
  requireOrThrow(resumed_report.restart_roundtrip_ok, "resumed DMO workflow restart round-trip failed");
  const cosmosim::io::RestartReadResult resumed_restart =
      cosmosim::io::readRestartCheckpointHdf5(resumed_report.restart_path);

  requireGlobally(
      direct_restart.state.validateOwnershipInvariants() &&
          resumed_restart.state.validateOwnershipInvariants(),
      runtime,
      "final DMO state violates ownership invariants");
  const StateSummary direct_summary = summarizeState(direct_restart.state, runtime, true);
  const StateSummary resumed_summary = summarizeState(resumed_restart.state, runtime, true);
  const DirectFourierMetrics direct_fourier = directFourierMetrics(direct_restart.state, true);
  const DirectFourierMetrics resumed_fourier = directFourierMetrics(resumed_restart.state, true);
  const cosmosim::core::SimulationState gathered_direct_state = gatherParticlesByStableId(
      direct_restart.state, runtime, true);
  const cosmosim::core::SimulationState gathered_resumed_state = gatherParticlesByStableId(
      resumed_restart.state, runtime, true);
  const std::array<double, 3> com_drift = periodicCenterOfMassDrift(direct_restart.state);

  // A zero-force trajectory is not static in a cosmological KDK update: the
  // stored physical peculiar velocity obeys the exactly integrated homogeneous
  // Hubble drag u(a)=u(a_initial)*a_initial/a.  Compare against that defensible
  // ballistic solution so the growing-mode checks below cannot pass merely
  // because the initial condition was seeded with a Zel'dovich velocity.
  const double ballistic_hubble_drag =
      k_initial_scale_factor / direct_restart.integrator_state.current_scale_factor;
  const double scale_factor_growth =
      direct_restart.integrator_state.current_scale_factor / k_initial_scale_factor;
  const double expected_linear_growth_increment = scale_factor_growth - 1.0;
  requireOrThrow(
      std::isfinite(expected_linear_growth_increment) && expected_linear_growth_increment > 0.0,
      "DMO validation requires a finite positive linear-growth increment");
  double local_gravity_response_velocity2 = 0.0;
  double local_ballistic_velocity2 = 0.0;
  double local_gravity_response_ballistic_dot = 0.0;
  for (std::size_t row = 0; row < direct_restart.state.particles.size(); ++row) {
    const std::uint64_t particle_id = direct_restart.state.particle_sidecar.particle_id[row];
    requireOrThrow(
        particle_id >= 1U && particle_id <= k_particle_count,
        "final DMO state contains a particle ID outside the deterministic initial condition");
    const std::size_t initial_row = static_cast<std::size_t>(particle_id - 1U);
    const std::array<double, 3> ballistic_velocity{
        ballistic_hubble_drag * initial_state.particles.velocity_x_peculiar[initial_row],
        ballistic_hubble_drag * initial_state.particles.velocity_y_peculiar[initial_row],
        ballistic_hubble_drag * initial_state.particles.velocity_z_peculiar[initial_row],
    };
    const std::array<double, 3> final_velocity{
        direct_restart.state.particles.velocity_x_peculiar[row],
        direct_restart.state.particles.velocity_y_peculiar[row],
        direct_restart.state.particles.velocity_z_peculiar[row],
    };
    const double mass = direct_restart.state.particles.mass_code[row];
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const double gravity_response = final_velocity[axis] - ballistic_velocity[axis];
      local_gravity_response_velocity2 += mass * gravity_response * gravity_response;
      local_ballistic_velocity2 += mass * ballistic_velocity[axis] * ballistic_velocity[axis];
      local_gravity_response_ballistic_dot += mass * gravity_response * ballistic_velocity[axis];
    }
  }
  const double gravity_response_velocity2 = globalSum(local_gravity_response_velocity2);
  const double ballistic_velocity2 = globalSum(local_ballistic_velocity2);
  const double gravity_response_ballistic_dot = globalSum(local_gravity_response_ballistic_dot);
  const double gravity_response_rms_velocity =
      std::sqrt(gravity_response_velocity2 / direct_summary.mass_code);
  const double ballistic_rms_velocity =
      std::sqrt(ballistic_velocity2 / direct_summary.mass_code);
  const double gravity_response_to_ballistic_ratio =
      gravity_response_rms_velocity / std::max(ballistic_rms_velocity, 1.0e-300);
  const double gravity_response_alignment = gravity_response_ballistic_dot /
      std::sqrt(std::max(gravity_response_velocity2 * ballistic_velocity2, 1.0e-300));
  const double gravity_response_projection_ratio =
      gravity_response_ballistic_dot / std::max(ballistic_velocity2, 1.0e-300);
  std::ostringstream gravity_response_message;
  gravity_response_message << std::setprecision(17)
                           << "production DMO gravity/ballistic response gate failed: response_rms="
                           << gravity_response_rms_velocity
                           << ", ballistic_rms=" << ballistic_rms_velocity
                           << ", response_to_ballistic=" << gravity_response_to_ballistic_ratio
                           << ", growing_mode_projection=" << gravity_response_projection_ratio
                           << ", alignment=" << gravity_response_alignment
                           << ", hubble_drag=" << ballistic_hubble_drag
                           << ", expected_linear_growth_increment="
                           << expected_linear_growth_increment;

  requireOrThrow(direct_summary.finite && resumed_summary.finite, "DMO workflow produced non-finite state");
  requireOrThrow(direct_summary.dm_only && resumed_summary.dm_only, "DMO workflow produced non-DM state");
  requireOrThrow(
      direct_summary.particle_count == k_particle_count && direct_summary.cell_count == 0U,
      "DMO workflow changed global particle/cell counts");
  requireOrThrow(
      direct_summary.particle_id_sum == (k_particle_count * (k_particle_count + 1U)) / 2U &&
          direct_summary.particle_id_xor == xorRangeOneToN(k_particle_count),
      "DMO workflow changed the stable particle-ID set");
  requireOrThrow(
      direct_report.global_particle_partition_identity_match &&
          resumed_report.global_particle_partition_identity_match,
      "DMO workflow distributed identity summary failed");

  const double mass_relative_error =
      std::abs(direct_summary.mass_code - initial_summary.mass_code) / initial_summary.mass_code;
  requireOrThrow(mass_relative_error <= 5.0e-13, "DMO workflow did not conserve mass");
  const double momentum_delta = std::sqrt(
      std::pow(direct_summary.momentum_code[0] - initial_summary.momentum_code[0], 2.0) +
      std::pow(direct_summary.momentum_code[1] - initial_summary.momentum_code[1], 2.0) +
      std::pow(direct_summary.momentum_code[2] - initial_summary.momentum_code[2], 2.0));
  const double momentum_scale = initial_summary.mass_code *
      std::max({initial_summary.rms_velocity_code, direct_summary.rms_velocity_code, 1.0});
  requireOrThrow(
      momentum_delta / momentum_scale <= 1.0e-9,
      "periodic DMO workflow violated global momentum conservation");
  // These response windows are regression sentinels around an independently
  // constructed zero-force/Hubble-drag trajectory.  The exact initial 4^3
  // lattice force is separately certified against the periodic Ewald
  // reference; its short-step kick predicts a response of about 0.019 times
  // the ballistic RMS.  Keep a wider integration envelope here while still
  // rejecting missing, sign-flipped, and materially mis-normalized gravity.
  requireOrThrow(
      std::isfinite(gravity_response_to_ballistic_ratio) &&
          gravity_response_to_ballistic_ratio >= 0.01 &&
          gravity_response_to_ballistic_ratio <= 0.03,
      gravity_response_message.str());
  requireOrThrow(
      std::isfinite(gravity_response_projection_ratio) &&
          gravity_response_projection_ratio >= 0.01 &&
          gravity_response_projection_ratio <= 0.03,
      gravity_response_message.str());
  requireOrThrow(
      std::isfinite(gravity_response_alignment) && gravity_response_alignment >= 0.99,
      gravity_response_message.str());
  const double com_drift_norm =
      std::sqrt(com_drift[0] * com_drift[0] + com_drift[1] * com_drift[1] + com_drift[2] * com_drift[2]);
  requireOrThrow(
      com_drift_norm / k_box_size_mpc_comoving <= 1.0e-7,
      "periodic center-of-mass drift exceeded tolerance");

  requireOrThrow(
      std::abs(initial_fourier.fundamental_amplitude - k_linear_density_amplitude) /
              k_linear_density_amplitude <=
          2.0e-2,
      "direct Fourier estimator does not recover the imposed linear mode");
  const double measured_growth =
      direct_fourier.fundamental_amplitude / initial_fourier.fundamental_amplitude;
  const double imposed_phase_coherence =
      -initial_fourier.fundamental_real / initial_fourier.fundamental_amplitude;
  const double initial_quadrature_fraction =
      std::abs(initial_fourier.fundamental_imag) / initial_fourier.fundamental_amplitude;
  const double growth_phase_drift = fourierPhaseDrift(initial_fourier, direct_fourier);
  const double growth_phase_coherence = fourierPhaseCoherence(initial_fourier, direct_fourier);
  const double restart_phase_drift = fourierPhaseDrift(direct_fourier, resumed_fourier);
  requireOrThrow(
      measured_growth > 1.001,
      "Zel'dovich fundamental mode did not grow over the production workflow");
  requireOrThrow(
      std::abs((measured_growth - 1.0) - expected_linear_growth_increment) /
              expected_linear_growth_increment <=
          0.075,
      "Zel'dovich fundamental growth increment is inconsistent with short-step linear growth");
  requireOrThrow(
      std::isfinite(imposed_phase_coherence) && imposed_phase_coherence >= 1.0 - 1.0e-12 &&
          std::isfinite(initial_quadrature_fraction) && initial_quadrature_fraction <= 1.0e-10,
      "initial Zel'dovich density mode does not have the imposed cosine phase");
  requireOrThrow(
      std::isfinite(growth_phase_drift) && growth_phase_drift <= 1.0e-8 &&
          std::isfinite(growth_phase_coherence) && growth_phase_coherence >= 1.0 - 1.0e-12,
      "Zel'dovich fundamental growth changed phase or lost coherence");
  requireOrThrow(
      std::isfinite(restart_phase_drift) && restart_phase_drift <= 1.0e-12,
      "uninterrupted/resumed Zel'dovich fundamental phases differ");
  const double high_k_contamination =
      direct_fourier.high_k_max_amplitude / direct_fourier.fundamental_amplitude;
  requireOrThrow(
      std::isfinite(high_k_contamination) && high_k_contamination <= 5.0e-2,
      "DMO workflow generated excessive high-k contamination");

  requireGlobally(
      direct_report.final_state_digest != 0U &&
          direct_report.final_state_digest == resumed_report.final_state_digest,
      runtime,
      "uninterrupted/resumed deterministic state digest mismatch");
  requireOrThrow(
      globalXor(direct_report.final_state_digest) == globalXor(resumed_report.final_state_digest),
      "global uninterrupted/resumed digest reduction mismatch");
  requireOrThrow(
      direct_restart.integrator_state.current_scale_factor ==
          resumed_restart.integrator_state.current_scale_factor,
      "uninterrupted/resumed scale-factor authority mismatch");
  requireOrThrow(
      std::abs(direct_summary.mass_code - resumed_summary.mass_code) <=
          5.0e-13 * direct_summary.mass_code,
      "uninterrupted/resumed mass mismatch");
  requireOrThrow(
      std::abs(direct_fourier.fundamental_amplitude - resumed_fourier.fundamental_amplitude) <=
          1.0e-12 * direct_fourier.fundamental_amplitude,
      "uninterrupted/resumed Fourier-mode mismatch");

  double measured_power_growth = 0.0;
  double initial_power_amplitude = 0.0;
  double initial_fundamental_power = 0.0;
  double evolved_fundamental_power = 0.0;
  double shot_noise_power = 0.0;
  std::size_t empty_power_bin_count = 0U;
  const std::uint64_t power_mode_count =
      static_cast<std::uint64_t>(k_power_spectrum_mesh_n) *
          static_cast<std::uint64_t>(k_power_spectrum_mesh_n) *
          static_cast<std::uint64_t>(k_power_spectrum_mesh_n) -
      1U;
  const std::filesystem::path power_artifact_path =
      root / "dmo_zeldovich_power_spectrum.json";
  if (runtime.world_rank == 0) {
    cosmosim::analysis::DiagnosticsEngine diagnostics(frozen.config);
    const cosmosim::analysis::PowerSpectrumEstimateOptions spectrum_options{
        .mesh_n = static_cast<std::size_t>(frozen.config.analysis.power_spectrum_mesh_n),
        .bin_count = static_cast<std::size_t>(frozen.config.analysis.power_spectrum_bin_count),
        .mass_assignment = cosmosim::analysis::PowerSpectrumMassAssignment::kCloudInCell,
        .window_correction =
            cosmosim::analysis::PowerSpectrumWindowCorrection::kDeconvolveAssignmentWindow,
        .shot_noise_policy =
            cosmosim::analysis::PowerSpectrumShotNoisePolicy::kReportWithoutSubtraction,
    };
    const cosmosim::analysis::PowerSpectrumEstimate initial_spectrum =
        diagnostics.computePowerSpectrumEstimate(gathered_initial_state, spectrum_options);
    const cosmosim::analysis::PowerSpectrumEstimate evolved_spectrum =
        diagnostics.computePowerSpectrumEstimate(gathered_direct_state, spectrum_options);
    const cosmosim::analysis::PowerSpectrumEstimate resumed_spectrum =
        diagnostics.computePowerSpectrumEstimate(gathered_resumed_state, spectrum_options);

    empty_power_bin_count = validatePowerSpectrumEstimate(initial_spectrum, "initial spectrum");
    requireOrThrow(
        validatePowerSpectrumEstimate(evolved_spectrum, "evolved spectrum") == empty_power_bin_count &&
            validatePowerSpectrumEstimate(resumed_spectrum, "resumed spectrum") == empty_power_bin_count,
        "initial/evolved/resumed spectra disagree on empty-bin handling");
    requirePowerSpectraEquivalent(
        evolved_spectrum,
        resumed_spectrum,
        "uninterrupted/resumed production power spectrum");

    const double expected_shot_noise =
        std::pow(k_box_size_mpc_comoving, 3.0) / static_cast<double>(k_particle_count);
    const double shot_noise_scale = std::max(expected_shot_noise, 1.0e-30);
    requireOrThrow(
        std::abs(initial_spectrum.poisson_shot_noise_code_volume - expected_shot_noise) <=
                1.0e-14 * shot_noise_scale &&
            std::abs(evolved_spectrum.poisson_shot_noise_code_volume - expected_shot_noise) <=
                1.0e-14 * shot_noise_scale,
        "reported equal-mass Poisson shot-noise level is inconsistent with V/N");

    const auto& initial_fundamental = fundamentalPowerBin(initial_spectrum);
    const auto& evolved_fundamental = fundamentalPowerBin(evolved_spectrum);
    requireOrThrow(
        initial_fundamental.mode_count == 6U && evolved_fundamental.mode_count == 6U,
        "fundamental shell is not isolated to the six Cartesian modes");
    requireOrThrow(
        initial_fundamental.power_code_volume > 0.0 &&
            evolved_fundamental.power_code_volume > initial_fundamental.power_code_volume,
        "production binned power spectrum did not grow in the fundamental shell");
    initial_fundamental_power = initial_fundamental.power_code_volume;
    evolved_fundamental_power = evolved_fundamental.power_code_volume;
    measured_power_growth = std::sqrt(evolved_fundamental_power / initial_fundamental_power);
    initial_power_amplitude = std::sqrt(
        12.0 * initial_fundamental_power /
        (k_box_size_mpc_comoving * k_box_size_mpc_comoving * k_box_size_mpc_comoving));
    shot_noise_power = initial_spectrum.poisson_shot_noise_code_volume;
    {
      std::ostringstream message;
      message << std::setprecision(17)
              << "CIC/deconvolved production spectrum does not recover the imposed linear amplitude: measured="
              << initial_power_amplitude << " expected=" << k_linear_density_amplitude;
      requireOrThrow(
          std::abs(initial_power_amplitude - k_linear_density_amplitude) /
                  k_linear_density_amplitude <=
              2.0e-2,
          message.str());
    }
    requireOrThrow(
        std::abs((measured_power_growth - 1.0) - expected_linear_growth_increment) /
                expected_linear_growth_increment <=
            0.075,
        "production power-spectrum growth increment is inconsistent with the expected linear solution");
    requireOrThrow(
        std::abs((measured_power_growth - 1.0) / (measured_growth - 1.0) - 1.0) <= 1.0e-3,
        "production power-spectrum growth increment disagrees with the independent particle-mode estimator");

    const std::string artifact = serializePowerSpectrumArtifact(
        runtime,
        direct_restart.integrator_state.current_scale_factor,
        scale_factor_growth,
        measured_power_growth,
        initial_power_amplitude,
        initial_spectrum,
        evolved_spectrum);
    writeAndVerifyPowerSpectrumArtifact(power_artifact_path, artifact);
  }

  if (runtime.world_rank == 0) {
    writePhysicalStateArtifact(
        physicalStateArtifactPath(currentExecutionBackend(), runtime.world_size),
        PhysicalStateArtifact{
            .execution_backend = std::string(currentExecutionBackend()),
            .world_size = runtime.world_size,
            .current_time_code = direct_restart.integrator_state.current_time_code,
            .current_scale_factor = direct_restart.integrator_state.current_scale_factor,
            .step_index = direct_restart.integrator_state.step_index,
            .state = gathered_direct_state,
        });
    std::cout << std::setprecision(17)
              << "DMO_ZELDOVICH_PASS"
              << " world_size=" << runtime.world_size
              << " particles=" << direct_summary.particle_count
              << " initial_mode=" << initial_fourier.fundamental_amplitude
              << " final_mode=" << direct_fourier.fundamental_amplitude
              << " measured_growth=" << measured_growth
              << " scale_factor_growth=" << scale_factor_growth
              << " initial_phase_coherence=" << imposed_phase_coherence
              << " initial_quadrature_fraction=" << initial_quadrature_fraction
              << " growth_phase_drift_rad=" << growth_phase_drift
              << " growth_phase_coherence=" << growth_phase_coherence
              << " restart_phase_drift_rad=" << restart_phase_drift
              << " pk_initial_amplitude=" << initial_power_amplitude
              << " pk_initial_power=" << initial_fundamental_power
              << " pk_evolved_power=" << evolved_fundamental_power
              << " pk_sqrt_growth=" << measured_power_growth
              << " pk_modes=" << power_mode_count
              << " pk_empty_bins=" << empty_power_bin_count
              << " pk_shot_noise=" << shot_noise_power
              << " high_k_ratio=" << high_k_contamination
              << " ballistic_hubble_drag=" << ballistic_hubble_drag
              << " gravity_response_rms_velocity=" << gravity_response_rms_velocity
              << " ballistic_rms_velocity=" << ballistic_rms_velocity
              << " gravity_response_to_ballistic_ratio=" << gravity_response_to_ballistic_ratio
              << " gravity_response_projection_ratio=" << gravity_response_projection_ratio
              << " gravity_response_alignment=" << gravity_response_alignment
              << " mass_relative_error=" << mass_relative_error
              << " momentum_relative_error=" << momentum_delta / momentum_scale
              << " periodic_com_drift=" << com_drift_norm
              << " power_artifact=" << power_artifact_path.string()
              << " physical_state_artifact="
              << physicalStateArtifactPath(currentExecutionBackend(), runtime.world_size).string()
              << " restart_equivalent=true\n";
  }
}
#endif

}  // namespace

int main(int argc, char** argv) {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
#endif
  const ParallelRuntime runtime = parallelRuntime();
  try {
    if (argc == 2 && std::string_view(argv[1]) == "--compare-rank-artifacts") {
      requireOrThrow(
          runtime.world_size == 1,
          "DMO rank-artifact comparison must run as a single process");
#if COSMOSIM_ENABLE_HDF5
      compareRankCountPhysicalStateArtifacts();
#else
      throw std::runtime_error(
          "DMO rank-artifact comparison requires COSMOSIM_ENABLE_HDF5=ON");
#endif
#if COSMOSIM_ENABLE_MPI
      MPI_Finalize();
#endif
      return 0;
    }
    if (argc == 2 && std::string_view(argv[1]) == "--compare-serial-mpi-artifacts") {
      requireOrThrow(
          runtime.world_size == 1,
          "DMO serial/MPI artifact comparison must run as a single process");
#if COSMOSIM_ENABLE_HDF5
      compareSerialAndSingleRankMpiPhysicalStateArtifacts();
#else
      throw std::runtime_error(
          "DMO serial/MPI artifact comparison requires COSMOSIM_ENABLE_HDF5=ON");
#endif
#if COSMOSIM_ENABLE_MPI
      MPI_Finalize();
#endif
      return 0;
    }
    if (argc == 2 && std::string_view(argv[1]) == "--empty-rank-smoke") {
      const std::filesystem::path smoke_root =
          std::filesystem::temp_directory_path() /
          ("cosmosim_dmo_empty_rank_smoke_" + std::string(currentExecutionBackend()) +
           "_np" + std::to_string(runtime.world_size));
      if (runtime.world_rank == 0) {
        std::filesystem::remove_all(smoke_root);
      }
      barrier();
      runEmptyRankSmoke(runtime, smoke_root);
      barrier();
#if COSMOSIM_ENABLE_MPI
      MPI_Finalize();
#endif
      return 0;
    }
    requireOrThrow(argc == 1, "unsupported DMO validation command-line argument");
    if (runtime.world_size < 1 || runtime.world_size > 4) {
      if (runtime.world_rank == 0) {
        std::cout << "DMO_ZELDOVICH_SKIP unsupported_world_size=" << runtime.world_size
                  << " supported=1,2,3,4\n";
      }
#if COSMOSIM_ENABLE_MPI
      MPI_Finalize();
#endif
      return 0;
    }

    const std::filesystem::path root =
        dmoWorkflowRoot(currentExecutionBackend(), runtime.world_size);
    if (runtime.world_rank == 0) {
      std::filesystem::remove_all(root);
    }
    barrier();

#if COSMOSIM_ENABLE_HDF5
    runScientificValidation(runtime, root);
#else
    throw std::runtime_error(
        "scientific DMO final-state and restart validation requires COSMOSIM_ENABLE_HDF5=ON");
#endif

    barrier();
#if COSMOSIM_ENABLE_MPI
    MPI_Finalize();
#endif
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "DMO_ZELDOVICH_FAIL rank=" << runtime.world_rank
              << " world_size=" << runtime.world_size
              << " reason=" << error.what() << '\n';
#if COSMOSIM_ENABLE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    return 1;
  }
}
