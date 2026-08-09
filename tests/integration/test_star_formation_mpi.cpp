#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/physics/star_formation.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

constexpr std::array<std::uint64_t, 6> kGasCellIds{
    81001U, 81002U, 81003U, 81004U, 81005U, 81006U};
constexpr double kInitialGasMassCode = 1000.0;
constexpr double kInitialMetalFraction = 0.02;
constexpr double kDtCode = 1.0e-5;

struct BirthRecord {
  std::uint64_t birth_key = 0;
  std::uint64_t particle_id = 0;
  std::uint64_t parent_gas_cell_id = 0;
  std::uint64_t birth_tick = 0;
  std::uint32_t birth_ordinal = 0;
  double birth_mass_code = 0.0;
};

[[nodiscard]] bool operator<(const BirthRecord& lhs, const BirthRecord& rhs) {
  if (lhs.birth_key != rhs.birth_key) {
    return lhs.birth_key < rhs.birth_key;
  }
  return lhs.particle_id < rhs.particle_id;
}

[[nodiscard]] bool sameRecord(const BirthRecord& lhs, const BirthRecord& rhs) {
  return lhs.birth_key == rhs.birth_key &&
      lhs.particle_id == rhs.particle_id &&
      lhs.parent_gas_cell_id == rhs.parent_gas_cell_id &&
      lhs.birth_tick == rhs.birth_tick &&
      lhs.birth_ordinal == rhs.birth_ordinal &&
      lhs.birth_mass_code == rhs.birth_mass_code;
}

[[nodiscard]] cosmosim::physics::StarFormationConfig modelConfig() {
  cosmosim::physics::StarFormationConfig config;
  config.enabled = true;
  config.model = cosmosim::core::StarFormationModelKind::kAdaptiveBoundJeans;
  config.newton_g_code = 1.0;
  config.density_is_comoving = false;
  config.geometry_is_comoving = false;
  config.epsilon_ff = 1.0;
  config.bound_alpha_vir_max = 1.0;
  config.require_converging_flow = true;
  config.collapse_timescale =
      cosmosim::core::StarFormationCollapseTimescale::kMinimumFreeFallOrCompression;
  config.jeans_mass_floor_code = 0.0;
  config.star_particle_mass_policy = cosmosim::core::StarParticleMassPolicy::kFixed;
  config.target_star_particle_mass_code = 10.0;
  config.min_star_particle_mass_code = 1.0;
  config.max_star_particle_mass_code = 100.0;
  config.max_spawn_particles_per_cell_step = 8U;
  config.max_fractional_mass_conversion = 0.25;
  config.min_remaining_gas_fraction = 0.01;
  config.min_remaining_gas_mass_code = 1.0;
  config.stochastic_spawning = true;
  config.random_seed = 0x1234abcd5678ef90ULL;
  config.metadata_schema_version = 3U;
  return config;
}

[[nodiscard]] int ownerForCell(
    std::uint64_t gas_cell_id,
    int world_size,
    int policy) {
  if (world_size <= 1) {
    return 0;
  }
  if (policy == 0) {
    return 0;  // explicit zero-local-eligible-cells exercise on non-root ranks
  }
  const int parity_owner = static_cast<int>(gas_cell_id % static_cast<std::uint64_t>(world_size));
  return policy == 1 ? parity_owner : (world_size - 1 - parity_owner);
}

[[nodiscard]] cosmosim::core::SimulationState makeLocalState(
    int world_rank,
    int world_size,
    int ownership_policy) {
  namespace core = cosmosim::core;
  std::vector<std::uint64_t> local_ids;
  for (const std::uint64_t gas_cell_id : kGasCellIds) {
    if (ownerForCell(gas_cell_id, world_size, ownership_policy) == world_rank) {
      local_ids.push_back(gas_cell_id);
    }
  }

  core::SimulationState state;
  state.resizeParticles(local_ids.size());
  state.resizeCells(local_ids.size());
  state.resizePatches(local_ids.empty() ? 0U : 1U);
  state.metadata.run_name = "star_formation_mpi";
  state.metadata.snapshot_stem = "snapshot";
  state.metadata.restart_stem = "restart";

  if (!local_ids.empty()) {
    state.patches.patch_id[0] = 90000U + static_cast<std::uint64_t>(world_rank);
    state.patches.level[0] = 0U;
    state.patches.first_cell[0] = 0U;
    state.patches.cell_count[0] = static_cast<std::uint32_t>(local_ids.size());
    state.patches.owning_rank[0] = static_cast<std::uint32_t>(world_rank);
    state.patches.origin_x_comoving[0] = 0.0;
    state.patches.origin_y_comoving[0] = 0.0;
    state.patches.origin_z_comoving[0] = 0.0;
    state.patches.extent_x_comoving[0] = 0.1 * static_cast<double>(local_ids.size());
    state.patches.extent_y_comoving[0] = 0.1;
    state.patches.extent_z_comoving[0] = 0.1;
    state.patches.cell_dim_x[0] = static_cast<std::uint32_t>(local_ids.size());
    state.patches.cell_dim_y[0] = 1U;
    state.patches.cell_dim_z[0] = 1U;
  }

  std::vector<core::GasCellIdentityRecord> identity_records;
  identity_records.reserve(local_ids.size());
  state.species.count_by_species.fill(0U);
  state.species.count_by_species[
      static_cast<std::size_t>(core::ParticleSpecies::kGas)] = local_ids.size();
  for (std::size_t row = 0; row < local_ids.size(); ++row) {
    const std::uint64_t gas_cell_id = local_ids[row];
    const std::uint64_t particle_id = 100000U + gas_cell_id;
    const double velocity_x = 1.0e-3 * static_cast<double>((gas_cell_id % 5U) + 1U);
    const double velocity_y = -0.5 * velocity_x;
    const double velocity_z = 0.25 * velocity_x;

    state.particles.position_x_comoving[row] = 0.05 + 0.1 * static_cast<double>(row);
    state.particles.position_y_comoving[row] = 0.05;
    state.particles.position_z_comoving[row] = 0.05;
    state.particles.velocity_x_peculiar[row] = velocity_x;
    state.particles.velocity_y_peculiar[row] = velocity_y;
    state.particles.velocity_z_peculiar[row] = velocity_z;
    state.particles.mass_code[row] = kInitialGasMassCode;
    state.particles.time_bin[row] = 0U;
    state.particle_sidecar.particle_id[row] = particle_id;
    state.particle_sidecar.sfc_key[row] = gas_cell_id;
    state.particle_sidecar.species_tag[row] =
        static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    state.particle_sidecar.owning_rank[row] = static_cast<std::uint32_t>(world_rank);

    state.cells.center_x_comoving[row] = state.particles.position_x_comoving[row];
    state.cells.center_y_comoving[row] = state.particles.position_y_comoving[row];
    state.cells.center_z_comoving[row] = state.particles.position_z_comoving[row];
    state.cells.mass_code[row] = kInitialGasMassCode;
    state.cells.patch_index[row] = 0U;
    state.cells.time_bin[row] = 0U;
    state.gas_cells.gas_cell_id[row] = gas_cell_id;
    state.gas_cells.parent_particle_id[row] = particle_id;
    state.gas_cells.velocity_x_peculiar[row] = velocity_x;
    state.gas_cells.velocity_y_peculiar[row] = velocity_y;
    state.gas_cells.velocity_z_peculiar[row] = velocity_z;
    state.gas_cells.density_code[row] = 1.0e6;
    state.gas_cells.pressure_code[row] = 1.0;
    state.gas_cells.internal_energy_code[row] = 1.0e-6;
    state.gas_cells.temperature_code[row] = 100.0;
    state.gas_cells.sound_speed_code[row] = 1.0e-3;
    state.gas_cells.metal_mass_code[row] =
        kInitialMetalFraction * kInitialGasMassCode;
    identity_records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = gas_cell_id,
        .parent_particle_id = particle_id,
        .owning_patch_id = state.patches.patch_id[0],
        .local_cell_row = static_cast<std::uint32_t>(row),
    });
  }
  state.rebuildSpeciesIndex();
  state.replaceGasCellIdentityRecords(std::move(identity_records));
  assert(state.validateOwnershipInvariants());
  return state;
}

[[nodiscard]] std::vector<cosmosim::physics::StarFormationCellInput> makeInputs(
    const cosmosim::core::SimulationState& state,
    int world_rank) {
  std::vector<cosmosim::physics::StarFormationCellInput> inputs;
  inputs.reserve(state.cells.size());
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    cosmosim::physics::StarFormationCellInput input;
    input.cell_index = row;
    input.gas_cell_id = state.gas_cells.gas_cell_id[row];
    input.owning_rank = static_cast<std::uint32_t>(world_rank);
    input.is_active = true;
    input.is_owned = true;
    input.is_leaf = true;
    input.is_ghost = false;
    input.gas_mass_code = state.cells.mass_code[row];
    input.gas_density_code = state.gas_cells.density_code[row];
    input.cell_volume_code = input.gas_mass_code / input.gas_density_code;
    input.gas_temperature_k = state.gas_cells.temperature_code[row];
    input.gas_sound_speed_code = state.gas_cells.sound_speed_code[row];
    input.velocity_x_peculiar = state.gas_cells.velocity_x_peculiar[row];
    input.velocity_y_peculiar = state.gas_cells.velocity_y_peculiar[row];
    input.velocity_z_peculiar = state.gas_cells.velocity_z_peculiar[row];
    input.velocity_divergence_code = -100.0;
    input.velocity_gradient_frobenius_sq_code = 1.0e-6;
    input.gas_metal_mass_code = state.gas_cells.metal_mass_code[row];
    input.metallicity_mass_fraction = input.gas_metal_mass_code / input.gas_mass_code;
    input.center_x_comoving = state.cells.center_x_comoving[row];
    input.center_y_comoving = state.cells.center_y_comoving[row];
    input.center_z_comoving = state.cells.center_z_comoving[row];
    inputs.push_back(input);
  }
  return inputs;
}

void applyStep(
    cosmosim::core::SimulationState& state,
    int world_rank,
    std::uint64_t tick) {
  const cosmosim::physics::StarFormationModel model(modelConfig());
  const auto inputs = makeInputs(state, world_rank);
  (void)model.applyFromInputs(state, inputs, kDtCode, 1.0, tick);
  assert(state.validateOwnershipInvariants());
}

[[nodiscard]] std::vector<BirthRecord> localBirthRecords(
    const cosmosim::core::SimulationState& state) {
  std::vector<BirthRecord> records;
  records.reserve(state.star_particles.size());
  for (std::size_t star_row = 0; star_row < state.star_particles.size(); ++star_row) {
    const std::uint32_t particle_row = state.star_particles.particle_index[star_row];
    records.push_back(BirthRecord{
        .birth_key = state.star_particles.birth_key[star_row],
        .particle_id = state.particle_sidecar.particle_id[particle_row],
        .parent_gas_cell_id = state.star_particles.parent_gas_cell_id[star_row],
        .birth_tick = state.star_particles.birth_tick[star_row],
        .birth_ordinal = state.star_particles.birth_ordinal[star_row],
        .birth_mass_code = state.star_particles.birth_mass_code[star_row],
    });
  }
  std::sort(records.begin(), records.end());
  return records;
}

#if COSMOSIM_ENABLE_MPI
[[nodiscard]] std::vector<BirthRecord> gatherBirthRecords(
    const cosmosim::core::SimulationState& state,
    int world_size) {
  const std::vector<BirthRecord> local = localBirthRecords(state);
  if (local.size() > static_cast<std::size_t>(std::numeric_limits<int>::max() / sizeof(BirthRecord))) {
    throw std::overflow_error("local star-formation record count exceeds MPI int byte count");
  }
  const int local_bytes = static_cast<int>(local.size() * sizeof(BirthRecord));
  std::vector<int> byte_counts(static_cast<std::size_t>(world_size), 0);
  MPI_Allgather(
      &local_bytes, 1, MPI_INT,
      byte_counts.data(), 1, MPI_INT,
      MPI_COMM_WORLD);
  std::vector<int> displacements(static_cast<std::size_t>(world_size), 0);
  int total_bytes = 0;
  for (int rank = 0; rank < world_size; ++rank) {
    displacements[static_cast<std::size_t>(rank)] = total_bytes;
    total_bytes += byte_counts[static_cast<std::size_t>(rank)];
  }
  if (total_bytes % static_cast<int>(sizeof(BirthRecord)) != 0) {
    throw std::runtime_error("MPI star-formation record payload is misaligned");
  }
  std::vector<BirthRecord> gathered(
      static_cast<std::size_t>(total_bytes) / sizeof(BirthRecord));
  MPI_Allgatherv(
      local.empty() ? nullptr : static_cast<const void*>(local.data()),
      local_bytes,
      MPI_BYTE,
      gathered.empty() ? nullptr : static_cast<void*>(gathered.data()),
      byte_counts.data(),
      displacements.data(),
      MPI_BYTE,
      MPI_COMM_WORLD);
  std::sort(gathered.begin(), gathered.end());
  return gathered;
}
#endif

void assertSameRecords(
    std::span<const BirthRecord> lhs,
    std::span<const BirthRecord> rhs) {
  assert(lhs.size() == rhs.size());
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    assert(sameRecord(lhs[i], rhs[i]));
  }
}

void assertUniqueParticleIds(std::span<const BirthRecord> records) {
  std::vector<std::uint64_t> ids;
  ids.reserve(records.size());
  for (const BirthRecord& record : records) {
    assert(record.birth_key != 0U);
    assert(record.particle_id != 0U);
    ids.push_back(record.particle_id);
  }
  std::sort(ids.begin(), ids.end());
  assert(std::adjacent_find(ids.begin(), ids.end()) == ids.end());
}

#if COSMOSIM_ENABLE_MPI && COSMOSIM_ENABLE_HDF5
[[nodiscard]] cosmosim::io::RestartReadResult restartRoundTrip(
    const cosmosim::core::SimulationState& state,
    int world_rank,
    int world_size) {
  cosmosim::core::IntegratorState integrator_state;
  integrator_state.step_index = 41U;
  integrator_state.current_time_code = kDtCode;
  integrator_state.current_scale_factor = 1.0;
  integrator_state.last_completed_restart_safe = true;
  integrator_state.last_completed_boundary_kind =
      cosmosim::core::StepBoundaryKind::kGlobalSynchronizationPoint;

  cosmosim::core::HierarchicalTimeBinScheduler particle_scheduler(0U);
  particle_scheduler.reset(
      static_cast<std::uint32_t>(state.particles.size()), 0U, 0U);
  cosmosim::core::HierarchicalTimeBinScheduler gas_scheduler(0U);
  gas_scheduler.reset(
      static_cast<std::uint32_t>(state.cells.size()), 0U, 0U);

  cosmosim::io::RestartWritePayload payload;
  payload.persistent_state.simulation_state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &particle_scheduler;
  payload.gas_cell_scheduler = &gas_scheduler;
  payload.normalized_config_text =
      "schema_version = 1\n[physics]\nstar_formation_model = adaptive_bound_jeans\n";
  payload.normalized_config_hash_hex =
      cosmosim::core::stableConfigHashHex(payload.normalized_config_text);
  payload.provenance = cosmosim::core::makeProvenanceRecord(
      payload.normalized_config_hash_hex, "star_formation_mpi_test", 0, payload.normalized_config_text);
  payload.distributed_gravity_state.schema_version = 2U;
  payload.distributed_gravity_state.world_size = world_size;
  payload.distributed_gravity_state.pm_grid_nx = 4U;
  payload.distributed_gravity_state.pm_grid_ny = 4U;
  payload.distributed_gravity_state.pm_grid_nz = 4U;
  payload.distributed_gravity_state.owning_rank_by_item.assign(
      state.particle_sidecar.owning_rank.begin(),
      state.particle_sidecar.owning_rank.end());
  payload.distributed_gravity_state.pm_slab_begin_x_by_rank.resize(
      static_cast<std::size_t>(world_size), 0U);
  payload.distributed_gravity_state.pm_slab_end_x_by_rank.resize(
      static_cast<std::size_t>(world_size), 0U);
  for (int rank = 0; rank < world_size; ++rank) {
    payload.distributed_gravity_state.pm_slab_begin_x_by_rank[
        static_cast<std::size_t>(rank)] = static_cast<std::uint64_t>(rank * 2);
    payload.distributed_gravity_state.pm_slab_end_x_by_rank[
        static_cast<std::size_t>(rank)] = static_cast<std::uint64_t>((rank + 1) * 2);
  }
  payload.stochastic_state.modules.push_back(
      cosmosim::io::StochasticModulePersistentState{
          .module_name = "star_formation",
          .schema_version = 2U,
          .rng_policy = "splitmix64_counter_keyed(seed,gas_cell_id,tick,ordinal,schema)",
          .random_seed = modelConfig().random_seed,
          .rank_local_seed_offset = 0U,
          .last_committed_step_index = 41U,
          .deterministic_from_serialized_inputs = true,
      });

  const std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      ("cosmosim_star_formation_mpi_rank_" + std::to_string(world_rank) + ".hdf5");
  std::filesystem::remove(path);
  cosmosim::io::writeRestartCheckpointHdf5(path, payload);
  cosmosim::io::RestartReadResult restored =
      cosmosim::io::readRestartCheckpointHdf5(path);
  std::filesystem::remove(path);
  return restored;
}
#endif

}  // namespace

int main(int argc, char** argv) {
#if COSMOSIM_ENABLE_MPI && COSMOSIM_ENABLE_HDF5
  MPI_Init(&argc, &argv);
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    MPI_Finalize();
    return 0;
  }
  const std::string_view mode = argc > 1 ? std::string_view(argv[1]) : "base";

  if (mode == "base") {
    cosmosim::core::SimulationState distributed =
        makeLocalState(world_rank, world_size, 0);
    applyStep(distributed, world_rank, 17U);
    const std::vector<BirthRecord> gathered =
        gatherBirthRecords(distributed, world_size);

    cosmosim::core::SimulationState serial = makeLocalState(0, 1, 0);
    applyStep(serial, 0, 17U);
    const std::vector<BirthRecord> expected = localBirthRecords(serial);
    assertSameRecords(gathered, expected);
    assertUniqueParticleIds(gathered);
    if (world_rank == 1) {
      assert(distributed.star_particles.size() == 0U);
    }
  } else if (mode == "repartition") {
    cosmosim::core::SimulationState parity =
        makeLocalState(world_rank, world_size, 1);
    cosmosim::core::SimulationState swapped =
        makeLocalState(world_rank, world_size, 2);
    applyStep(parity, world_rank, 23U);
    applyStep(swapped, world_rank, 23U);
    const std::vector<BirthRecord> parity_records =
        gatherBirthRecords(parity, world_size);
    const std::vector<BirthRecord> swapped_records =
        gatherBirthRecords(swapped, world_size);
    assertSameRecords(parity_records, swapped_records);
    assertUniqueParticleIds(parity_records);
  } else if (mode == "restart") {
    cosmosim::core::SimulationState direct =
        makeLocalState(world_rank, world_size, 1);
    applyStep(direct, world_rank, 41U);
    cosmosim::core::SimulationState checkpoint_source = direct;
    applyStep(direct, world_rank, 42U);

    cosmosim::io::RestartReadResult restored =
        restartRoundTrip(checkpoint_source, world_rank, world_size);
    assert(restored.state.validateOwnershipInvariants());
    applyStep(restored.state, world_rank, 42U);

    const std::vector<BirthRecord> direct_records =
        gatherBirthRecords(direct, world_size);
    const std::vector<BirthRecord> resumed_records =
        gatherBirthRecords(restored.state, world_size);
    assertSameRecords(direct_records, resumed_records);
    assertUniqueParticleIds(resumed_records);

    const double direct_local_gas_mass = std::accumulate(
        direct.cells.mass_code.begin(), direct.cells.mass_code.end(), 0.0);
    const double resumed_local_gas_mass = std::accumulate(
        restored.state.cells.mass_code.begin(),
        restored.state.cells.mass_code.end(),
        0.0);
    double direct_global_gas_mass = 0.0;
    double resumed_global_gas_mass = 0.0;
    MPI_Allreduce(
        &direct_local_gas_mass, &direct_global_gas_mass,
        1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(
        &resumed_local_gas_mass, &resumed_global_gas_mass,
        1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    assert(direct_global_gas_mass == resumed_global_gas_mass);
  } else {
    MPI_Finalize();
    throw std::invalid_argument("unknown star-formation MPI test mode");
  }

  MPI_Finalize();
#else
  (void)argc;
  (void)argv;
#endif
  return 0;
}
