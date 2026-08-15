#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "cosmosim/cosmosim.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/workflows/gravity_source_ownership.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

std::string configText(int world_size, std::string_view run_name) {
  std::ostringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[units]\nlength_unit = kpc\nmass_unit = msun\nvelocity_unit = km_s\n"
            "coordinate_frame = physical\n\n";
  stream << "[mode]\nmode = cosmo_cube\nic_file = generated\nhydro_boundary = periodic\n"
            "gravity_boundary = periodic\n\n";
  stream << "[cosmology]\nbox_size_x = 1.0\nbox_size_y = 1.0\nbox_size_z = 1.0\n\n";
  stream << "[numerics]\ntime_begin_code = 0.01\ntime_end_code = 0.010000001\n"
            "max_global_steps = 1\nhierarchical_max_rung = 0\ntreepm_pm_grid = 8\n"
            "treepm_asmth_cells = 1.0\ntreepm_rcut_cells = 3.0\n"
            "treepm_enable_window_deconvolution = false\ntreepm_update_cadence_steps = 1\n\n";
  stream << "[physics]\nenable_star_formation = true\n"
            "star_formation_model = effective_multiphase_tng_like\n"
            "sf_effective_t_star_at_threshold = 1.0e-9 code\n"
            "sf_effective_n_h_threshold = 0.13 cm^-3\n"
            "sf_star_particle_mass_policy = fixed\n"
            "sf_target_star_particle_mass_code = 100.0\n"
            "sf_min_star_particle_mass_code = 1.0\n"
            "sf_max_star_particle_mass_code = 1000000.0\n"
            "sf_max_spawn_particles_per_cell_step = 1\n"
            "sf_max_fractional_mass_conversion = 0.25\n"
            "sf_min_remaining_gas_fraction = 0.01\n"
            "sf_min_remaining_gas_mass_code = 1.0\n"
            "sf_stochastic_spawning = false\n"
            "sf_random_seed = 424242\n\n";
  stream << "[output]\nrun_name = " << run_name
         << "\noutput_directory = integration_outputs\noutput_stem = snapshot\n"
            "restart_stem = restart\nsnapshot_interval_steps = 1\nwrite_restarts = true\n\n";
  stream << "[parallel]\nmpi_ranks_expected = " << world_size
         << "\ndecomposition_runtime_rebalance_enabled = false\n"
            "deterministic_reduction = true\n";
  return stream.str();
}

cosmosim::core::SimulationState makeGlobalState() {
  namespace core = cosmosim::core;
  core::SimulationState state;
  state.resizeParticles(2);
  state.resizeCells(2);
  state.resizePatches(2);
  state.metadata.run_name = "star_formation_source_runtime_mpi";
  state.metadata.snapshot_stem = "snapshot";
  state.metadata.restart_stem = "restart";

  std::vector<core::GasCellIdentityRecord> records;
  for (std::uint32_t row = 0; row < 2U; ++row) {
    const std::uint32_t owner = row;
    const std::uint64_t particle_id = 1001U + row;
    const std::uint64_t gas_cell_id = 2001U + row;
    state.patches.patch_id[row] = 3001U + row;
    state.patches.level[row] = 0U;
    state.patches.first_cell[row] = row;
    state.patches.cell_count[row] = 1U;
    state.patches.owning_rank[row] = owner;
    state.patches.origin_x_comoving[row] = 0.5 * static_cast<double>(row);
    state.patches.origin_y_comoving[row] = 0.0;
    state.patches.origin_z_comoving[row] = 0.0;
    state.patches.extent_x_comoving[row] = 0.5;
    state.patches.extent_y_comoving[row] = 1.0;
    state.patches.extent_z_comoving[row] = 1.0;
    state.patches.cell_dim_x[row] = 1U;
    state.patches.cell_dim_y[row] = 1U;
    state.patches.cell_dim_z[row] = 1U;

    state.particles.position_x_comoving[row] = 0.25 + 0.5 * static_cast<double>(row);
    state.particles.position_y_comoving[row] = 0.5;
    state.particles.position_z_comoving[row] = 0.5;
    state.particles.mass_code[row] = 1.0e6;
    state.particle_sidecar.particle_id[row] = particle_id;
    state.particle_sidecar.sfc_key[row] = particle_id;
    state.particle_sidecar.species_tag[row] =
        static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    state.particle_sidecar.owning_rank[row] = owner;

    state.cells.center_x_comoving[row] = state.particles.position_x_comoving[row];
    state.cells.center_y_comoving[row] = 0.5;
    state.cells.center_z_comoving[row] = 0.5;
    state.cells.mass_code[row] = 1.0e6;
    state.cells.patch_index[row] = row;
    state.cells.time_bin[row] = 0U;
    state.gas_cells.gas_cell_id[row] = gas_cell_id;
    state.gas_cells.parent_particle_id[row] = particle_id;
    state.gas_cells.density_code[row] = row == 0U ? 1.0e15 : 1.0e-20;
    state.gas_cells.pressure_code[row] = 1.0e3;
    state.gas_cells.internal_energy_code[row] = 1.0e-3;
    state.gas_cells.temperature_code[row] = 100.0;
    state.gas_cells.sound_speed_code[row] = 1.0e-3;
    state.gas_cells.metal_mass_code[row] = 2.0e4;
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = gas_cell_id,
        .parent_particle_id = particle_id,
        .owning_patch_id = 3001U + row,
        .local_cell_row = row,
    });
  }
  state.species.count_by_species.fill(0U);
  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)] = 2U;
  state.rebuildSpeciesIndex();
  state.replaceGasCellIdentityRecords(std::move(records));
  assert(state.validateOwnershipInvariants());
  return state;
}

void runGravityGasOwnershipContract(int world_size, int world_rank) {
  auto state = makeGlobalState();
  // Exercise parentless and shared-lineage legality through the same selector
  // and accessor used by GravityRuntime. Generic gas particles remain mirrors.
  state.gas_cells.parent_particle_id[0] = 0U;
  state.gas_cells.parent_particle_id[1] = 0U;
  auto rows = cosmosim::workflows::internal::selectAuthoritativeGravitySourceRows(
      state, static_cast<std::uint32_t>(world_rank),
      static_cast<std::uint32_t>(world_size), "mpi authoritative gas gravity contract");
  assert(rows.particle_rows.empty());
  double local_mass = 0.0;
  for (const std::uint32_t row : rows.gas_cell_rows) {
    local_mass += cosmosim::workflows::internal::authoritativeGasGravitySource(
        state, row, "mpi authoritative gas gravity contract").mass_code;
  }
  std::uint64_t local_count = static_cast<std::uint64_t>(rows.gas_cell_rows.size());
  const cosmosim::parallel::MpiContext mpi_context{};
  const std::uint64_t global_count = mpi_context.allreduceSumUint64(local_count);
  const double global_mass = mpi_context.allreduceSumDouble(local_mass);
  assert(global_count == 2U);
  assert(std::abs(global_mass - 2.0e6) <= 1.0e-9);
  if (world_size == 3 && world_rank == 2) {
    assert(rows.gas_cell_rows.empty());
  }
}


struct DistributedGasForceResult {
  std::vector<double> ax;
  std::vector<double> ay;
  std::vector<double> az;
  std::vector<std::size_t> target_global_cell_row;
  bool has_independent_target = false;
};

void setGasLineage(
    cosmosim::core::SimulationState& state,
    bool shared_lineage) {
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(state.cells.size());
  const std::uint64_t shared_parent = 1001U;
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const std::optional<std::uint64_t> parent = shared_lineage
        ? std::optional<std::uint64_t>{shared_parent}
        : std::nullopt;
    state.gas_cells.parent_particle_id[row] = parent.value_or(0U);
    records.push_back(cosmosim::core::GasCellIdentityRecord{
        .gas_cell_id = state.gas_cells.gas_cell_id[row],
        .parent_particle_id = parent,
        .owning_patch_id = state.patches.patch_id[state.cells.patch_index[row]],
        .local_cell_row = row,
    });
  }
  state.replaceGasCellIdentityRecords(std::move(records));
}

cosmosim::gravity::TreePmOptions gasForceOptions(
    std::span<const cosmosim::parallel::TopDomainLeaf> local_domain_leaves) {
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.opening_theta = 0.5;
  options.tree_options.max_leaf_size = 1U;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.0, 3.0, 1.0 / 8.0);
  options.source_generation = cosmosim::gravity::GravitySourceGeneration{1U};
  options.decomposition_epoch = cosmosim::gravity::DecompositionEpoch{1U};
  options.force_epoch.sequence = 1U;
  options.authoritative_domain_leaves = local_domain_leaves;
  return options;
}

DistributedGasForceResult solveDistributedGasTreePm(
    const cosmosim::core::SimulationState& state,
    int world_size,
    int world_rank) {
  const auto selected = cosmosim::workflows::internal::selectAuthoritativeGravitySourceRows(
      state, static_cast<std::uint32_t>(world_rank),
      static_cast<std::uint32_t>(world_size), "distributed gas numerical force acceptance");
  assert(selected.particle_rows.empty());

  std::vector<double> local_x;
  std::vector<double> local_y;
  std::vector<double> local_z;
  std::vector<double> local_mass;
  std::vector<std::size_t> local_to_global;
  std::vector<cosmosim::parallel::DecompositionItem> local_domain_items;
  for (const std::uint32_t row : selected.gas_cell_rows) {
    const auto source = cosmosim::workflows::internal::authoritativeGasGravitySource(
        state, row, "distributed gas numerical force acceptance");
    local_x.push_back(source.x_comoving);
    local_y.push_back(source.y_comoving);
    local_z.push_back(source.z_comoving);
    local_mass.push_back(source.mass_code);
    local_to_global.push_back(row);
    cosmosim::parallel::DecompositionItem item;
    item.entity_id = state.gas_cells.gas_cell_id[row];
    item.kind = cosmosim::parallel::DecompositionEntityKind::kAmrPatch;
    item.current_owner_rank = world_rank;
    item.x_comov = source.x_comoving;
    item.y_comov = source.y_comoving;
    item.z_comov = source.z_comoving;
    item.has_spatial_bounds = true;
    item.min_x_comov = item.max_x_comov = item.x_comov;
    item.min_y_comov = item.max_y_comov = item.y_comov;
    item.min_z_comov = item.max_z_comov = item.z_comov;
    item.work_units = 1.0;
    local_domain_items.push_back(item);
  }
  cosmosim::parallel::DecompositionConfig domain_config;
  domain_config.world_size = world_size;
  domain_config.work_weight = 1.0;
  const auto local_domain_leaves = cosmosim::parallel::buildAuthoritativeTopDomainLeaves(
      local_domain_items, domain_config, world_rank, 1U, 4U);

  std::vector<std::uint32_t> active;
  std::vector<double> independent_x;
  std::vector<double> independent_y;
  std::vector<double> independent_z;
  DistributedGasForceResult result;
  if (world_size == 3 && world_rank == 1) {
    // Source-only rank: it contributes authoritative gas but owns no active target.
  } else if (world_size == 3 && world_rank == 2) {
    // Target-only rank: no local authoritative source, but remote gas must still
    // produce the same force as the serial reference at an independent target.
    active.push_back(std::numeric_limits<std::uint32_t>::max());
    independent_x.push_back(0.50);
    independent_y.push_back(0.25);
    independent_z.push_back(0.75);
    result.has_independent_target = true;
  } else {
    active.resize(local_x.size());
    for (std::size_t i = 0; i < active.size(); ++i) {
      active[i] = static_cast<std::uint32_t>(i);
      result.target_global_cell_row.push_back(local_to_global[i]);
    }
  }
  result.ax.assign(active.size(), 0.0);
  result.ay.assign(active.size(), 0.0);
  result.az.assign(active.size(), 0.0);
  cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      .active_particle_index = active,
      .accel_x_comoving = result.ax,
      .accel_y_comoving = result.ay,
      .accel_z_comoving = result.az,
      .target_pos_x_comoving = independent_x,
      .target_pos_y_comoving = independent_y,
      .target_pos_z_comoving = independent_z,
  };
  const cosmosim::gravity::PmGridShape shape{8U, 8U, 8U};
  const auto layout = cosmosim::parallel::makePmSlabLayout(
      shape.nx, shape.ny, shape.nz, world_size, world_rank);
  cosmosim::gravity::TreePmCoordinator distributed(shape, layout);
  auto options = gasForceOptions(local_domain_leaves);
  distributed.solveActiveSet(
      local_x, local_y, local_z, local_mass, accumulator, options, nullptr, nullptr);
  return result;
}

std::array<double, 3> serialGasTreePmReference(
    const cosmosim::core::SimulationState& state,
    std::optional<std::size_t> target_cell_row,
    std::optional<std::array<double, 3>> independent_target) {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> mass;
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const auto source = cosmosim::workflows::internal::authoritativeGasGravitySource(
        state, row, "serial gas numerical force reference");
    x.push_back(source.x_comoving);
    y.push_back(source.y_comoving);
    z.push_back(source.z_comoving);
    mass.push_back(source.mass_code);
  }
  std::vector<std::uint32_t> active;
  std::vector<double> tx;
  std::vector<double> ty;
  std::vector<double> tz;
  if (target_cell_row.has_value()) {
    active.push_back(static_cast<std::uint32_t>(*target_cell_row));
  } else {
    assert(independent_target.has_value());
    active.push_back(std::numeric_limits<std::uint32_t>::max());
    tx.push_back((*independent_target)[0]);
    ty.push_back((*independent_target)[1]);
    tz.push_back((*independent_target)[2]);
  }
  std::vector<double> ax(1U, 0.0), ay(1U, 0.0), az(1U, 0.0);
  cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      .active_particle_index = active,
      .accel_x_comoving = ax,
      .accel_y_comoving = ay,
      .accel_z_comoving = az,
      .target_pos_x_comoving = tx,
      .target_pos_y_comoving = ty,
      .target_pos_z_comoving = tz,
  };
  const cosmosim::gravity::PmGridShape shape{8U, 8U, 8U};
  cosmosim::gravity::TreePmCoordinator reference(shape);
  auto options = gasForceOptions({});
  reference.solveActiveSet(x, y, z, mass, accumulator, options, nullptr, nullptr);
  return {ax[0], ay[0], az[0]};
}

void runGravityGasNumericalForceContract(int world_size, int world_rank) {
  auto parentless = makeGlobalState();
  setGasLineage(parentless, false);
  const auto parentless_force = solveDistributedGasTreePm(parentless, world_size, world_rank);
  for (std::size_t i = 0; i < parentless_force.ax.size(); ++i) {
    const auto reference = parentless_force.has_independent_target
        ? serialGasTreePmReference(parentless, std::nullopt, std::array<double, 3>{0.50, 0.25, 0.75})
        : serialGasTreePmReference(parentless, parentless_force.target_global_cell_row[i], std::nullopt);
    assert(std::abs(parentless_force.ax[i] - reference[0]) < 1.0e-9);
    assert(std::abs(parentless_force.ay[i] - reference[1]) < 1.0e-9);
    assert(std::abs(parentless_force.az[i] - reference[2]) < 1.0e-9);
  }

  auto shared_lineage = makeGlobalState();
  setGasLineage(shared_lineage, true);
  const auto shared_force = solveDistributedGasTreePm(shared_lineage, world_size, world_rank);
  assert(shared_force.ax.size() == parentless_force.ax.size());
  for (std::size_t i = 0; i < shared_force.ax.size(); ++i) {
    assert(std::abs(shared_force.ax[i] - parentless_force.ax[i]) < 1.0e-12);
    assert(std::abs(shared_force.ay[i] - parentless_force.ay[i]) < 1.0e-12);
    assert(std::abs(shared_force.az[i] - parentless_force.az[i]) < 1.0e-12);
  }

  auto no_mirrors = parentless;
  no_mirrors.resizeParticles(0U);
  no_mirrors.species.count_by_species.fill(0U);
  no_mirrors.rebuildSpeciesIndex();
  const auto no_mirror_force = solveDistributedGasTreePm(no_mirrors, world_size, world_rank);
  assert(no_mirror_force.ax.size() == parentless_force.ax.size());
  for (std::size_t i = 0; i < no_mirror_force.ax.size(); ++i) {
    assert(std::abs(no_mirror_force.ax[i] - parentless_force.ax[i]) < 1.0e-12);
    assert(std::abs(no_mirror_force.ay[i] - parentless_force.ay[i]) < 1.0e-12);
    assert(std::abs(no_mirror_force.az[i] - parentless_force.az[i]) < 1.0e-12);
  }
}

std::uint64_t localStarCount(const cosmosim::core::SimulationState& state) {
  return static_cast<std::uint64_t>(state.star_particles.size());
}

std::string sidecarText(const cosmosim::core::SimulationState& state, const std::string& name) {
  const auto* block = state.sidecars.find(name);
  if (block == nullptr) return {};
  return std::string(reinterpret_cast<const char*>(block->payload.data()), block->payload.size());
}

}  // namespace

int main(int argc, char** argv) {
#if COSMOSIM_ENABLE_MPI && COSMOSIM_ENABLE_HDF5
  MPI_Init(&argc, &argv);
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  const std::string mode = argc > 1 ? argv[1] : "source_runtime";
  if (mode == "gravity_gas_ownership" || mode == "gravity_gas_force") {
    if (world_size == 2 || world_size == 3) {
      runGravityGasOwnershipContract(world_size, world_rank);
      if (mode == "gravity_gas_force") {
        runGravityGasNumericalForceContract(world_size, world_rank);
      }
    }
    MPI_Finalize();
    return 0;
  }
  if (world_size != 2) {
    MPI_Finalize();
    return 0;
  }
  const auto initial_state = makeGlobalState();
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      configText(world_size, "star_formation_source_runtime_mpi_" + mode),
      "test_star_formation_source_runtime_mpi");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  const auto root = std::filesystem::temp_directory_path() /
      ("cosmosim_star_formation_source_runtime_mpi_" + mode);
  const auto report = runner.run(
      root / ("rank_" + std::to_string(world_rank)),
      cosmosim::workflows::ReferenceWorkflowOptions{
          .dt_time_code = 1.0e-9,
          .write_outputs = true,
          .initial_state_override = &initial_state,
          .max_steps_override = 1});
  assert(report.restart_roundtrip_ok);
  const auto restart = cosmosim::io::readRestartCheckpointHdf5(report.restart_path);
  const std::uint64_t global_stars = cosmosim::parallel::MpiContext{}.allreduceSumUint64(
      localStarCount(restart.state));
  assert(global_stars == 1U);
  assert(restart.state.validateUniqueParticleIds());

  const std::string eos_sidecar = sidecarText(restart.state, "effective_multiphase_ism");
  assert(eos_sidecar.find("table_hash=") != std::string::npos);
  const std::uint64_t local_hash = cosmosim::core::stableConfigHash(eos_sidecar);
  const std::uint64_t hash_xor = cosmosim::parallel::MpiContext{}.allreduceXorUint64(local_hash);
  assert(hash_xor == 0U);

  if (mode == "restart") {
    const auto resumed = runner.run(
        root / ("rank_" + std::to_string(world_rank) + "_resumed"),
        cosmosim::workflows::ReferenceWorkflowOptions{
            .dt_time_code = 1.0e-9,
            .write_outputs = true,
            .restart_state_override = &restart,
            .max_steps_override = 1});
    assert(resumed.restart_roundtrip_ok);
  }
  MPI_Finalize();
#else
  (void)argc;
  (void)argv;
#endif
  return 0;
}
