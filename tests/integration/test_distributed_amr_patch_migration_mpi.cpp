#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_scheduler.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "workflows/internal/migration_wire.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#include "../support/mpi_test_workspace.hpp"
#endif

namespace {

using cosmosim::core::AmrPatchMigrationCommit;
using cosmosim::core::AmrPatchMigrationRecord;
using cosmosim::core::SimulationState;

void seedPatch(
    SimulationState& state,
    std::uint64_t patch_id,
    std::uint32_t owner_rank,
    std::span<const std::uint64_t> gas_cell_ids,
    double origin_x_comoving) {
  state.resizePatches(1U);
  state.patches.patch_id[0] = patch_id;
  state.patches.level[0] = 1;
  state.patches.first_cell[0] = 0U;
  state.patches.cell_count[0] = static_cast<std::uint32_t>(gas_cell_ids.size());
  state.patches.parent_patch_id[0] = 0U;
  state.patches.morton_key[0] = patch_id;
  state.patches.owning_rank[0] = owner_rank;
  state.patches.origin_x_comoving[0] = origin_x_comoving;
  state.patches.origin_y_comoving[0] = 0.0;
  state.patches.origin_z_comoving[0] = 0.0;
  state.patches.extent_x_comoving[0] = static_cast<double>(gas_cell_ids.size());
  state.patches.extent_y_comoving[0] = 1.0;
  state.patches.extent_z_comoving[0] = 1.0;
  state.patches.cell_dim_x[0] = static_cast<std::uint16_t>(gas_cell_ids.size());
  state.patches.cell_dim_y[0] = 1U;
  state.patches.cell_dim_z[0] = 1U;

  state.resizeCells(gas_cell_ids.size());
  std::vector<cosmosim::core::GasCellIdentityRecord> identity_records;
  identity_records.reserve(gas_cell_ids.size());
  for (std::uint32_t row = 0U; row < gas_cell_ids.size(); ++row) {
    state.cells.center_x_comoving[row] = origin_x_comoving + 0.5 + static_cast<double>(row);
    state.cells.center_y_comoving[row] = 0.5;
    state.cells.center_z_comoving[row] = 0.5;
    state.cells.mass_code[row] = 1.0 + 0.25 * static_cast<double>(row);
    state.cells.patch_index[row] = 0U;
    state.cells.time_bin[row] = 0U;
    state.gas_cells.gas_cell_id[row] = gas_cell_ids[row];
    state.gas_cells.parent_particle_id[row] = 0U;
    state.gas_cells.velocity_x_peculiar[row] = 0.1 * static_cast<double>(row + 1U);
    state.gas_cells.velocity_y_peculiar[row] = 0.0;
    state.gas_cells.velocity_z_peculiar[row] = 0.0;
    state.gas_cells.density_code[row] = 1.0 + 0.1 * static_cast<double>(row);
    state.gas_cells.pressure_code[row] = 1.0;
    state.gas_cells.internal_energy_code[row] = 2.5;
    state.gas_cells.temperature_code[row] = 1.0;
    state.gas_cells.sound_speed_code[row] = 1.0;
    state.gas_cells.metal_mass_code[row] = 0.01 * state.cells.mass_code[row];
    identity_records.push_back(cosmosim::core::GasCellIdentityRecord{
        .gas_cell_id = gas_cell_ids[row],
        .parent_particle_id = std::nullopt,
        .owning_patch_id = patch_id,
        .local_cell_row = row,
    });
  }
  state.replaceGasCellIdentityRecords(std::move(identity_records));
  state.rebuildSpeciesIndex();
  assert(state.validateOwnershipInvariants());
}

SimulationState makeSourceState() {
  SimulationState state;
  constexpr std::array<std::uint64_t, 2> k_gas_ids{7001U, 7002U};
  seedPatch(state, 101U, 0U, k_gas_ids, 0.0);

  cosmosim::core::PendingFluxRegisterRecord pending;
  pending.register_key = 0xA501U;
  pending.coarse_patch_id = 101U;
  pending.coarse_gas_cell_id = 7001U;
  pending.coarse_cell_index = 0U;
  pending.level = 1U;
  pending.axis = 0U;
  pending.orientation = 1U;
  pending.expected_area_comov = 0.5;
  pending.interval_start_code = 1.0;
  pending.interval_end_code = 1.25;
  pending.coarse_dt_code = 0.25;
  pending.expected_fine_substeps = 2U;
  pending.completed_fine_substeps = 1U;
  pending.fine_substep_coverage_mask = 1U;
  pending.coarse_face_count = 1U;
  pending.fine_face_count = 1U;
  pending.gas_cell_identity_generation = state.gasCellIdentityGeneration();
  pending.patch_geometry_generation = state.cellIndexGeneration();
  pending.coarse_mass_flux_integral_code = 0.125;
  pending.fine_mass_flux_integral_code = 0.25;
  state.pending_flux_registers.assign({pending});

  cosmosim::core::AmrTemporalBoundaryHistoryRecord history;
  history.patch_id = 101U;
  history.patch_level = 1U;
  history.patch_geometry_fingerprint = 0xCAFEU;
  history.gas_cell_identity_generation = state.gasCellIdentityGeneration();
  history.interval_start_code = 1.0;
  history.interval_end_code = 1.25;
  history.end_state_valid = true;
  for (std::size_t offset = 0U; offset < k_gas_ids.size(); ++offset) {
    cosmosim::core::AmrTemporalBoundaryHistoryCellRecord cell;
    cell.gas_cell_id = k_gas_ids[offset];
    cell.patch_local_cell = offset;
    cell.start_mass_density_comoving = 10.0 + static_cast<double>(offset);
    cell.end_mass_density_comoving = 11.0 + static_cast<double>(offset);
    history.cells.push_back(cell);
  }
  state.amr_temporal_boundary_history.assign({history});
  return state;
}

SimulationState makeDestinationState() {
  SimulationState state;
  constexpr std::array<std::uint64_t, 1> k_gas_ids{8001U};
  seedPatch(state, 202U, 1U, k_gas_ids, 4.0);
  return state;
}

std::uint32_t patchRow(const SimulationState& state, std::uint64_t patch_id) {
  const auto it = std::find(state.patches.patch_id.begin(), state.patches.patch_id.end(), patch_id);
  assert(it != state.patches.patch_id.end());
  return static_cast<std::uint32_t>(std::distance(state.patches.patch_id.begin(), it));
}

void assertMigratedSynchronizationState(const SimulationState& state) {
  const std::uint32_t migrated_patch_row = patchRow(state, 101U);
  assert(state.patches.owning_rank[migrated_patch_row] == 1U);
  assert(state.patches.cell_count[migrated_patch_row] == 2U);

  const auto gas_row = state.rowForGasCellId(7001U);
  assert(gas_row.has_value());
  assert(state.owningPatchIdForGasCellId(7001U).value() == 101U);
  assert(state.owningPatchIdForGasCellId(7002U).value() == 101U);

  assert(state.pending_flux_registers.size() == 1U);
  const auto& pending = state.pending_flux_registers.records().front();
  assert(pending.register_key == 0xA501U);
  assert(pending.coarse_patch_id == 101U);
  assert(pending.coarse_gas_cell_id == 7001U);
  assert(pending.coarse_cell_index == *gas_row);
  assert(pending.gas_cell_identity_generation == state.gasCellIdentityGeneration());
  assert(pending.patch_geometry_generation == state.cellIndexGeneration());

  assert(state.amr_temporal_boundary_history.size() == 1U);
  const auto& history = state.amr_temporal_boundary_history.records().front();
  assert(history.patch_id == 101U);
  assert(history.gas_cell_identity_generation == state.gasCellIdentityGeneration());
  assert(history.cells.size() == 2U);
  assert(history.cells[0].gas_cell_id == 7001U);
  assert(history.cells[1].gas_cell_id == 7002U);
  assert(state.validateOwnershipInvariants());
  assert(cosmosim::amr::hasProductionAmrHydroCoverage(state));
}

cosmosim::io::RestartReadResult checkpointAndReload(
    SimulationState& state,
    const std::filesystem::path& path) {
  cosmosim::core::IntegratorState integrator_state;
  integrator_state.current_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
  integrator_state.last_completed_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
  integrator_state.last_completed_restart_safe = true;
  integrator_state.step_index = 7U;

  cosmosim::core::HierarchicalTimeBinScheduler particle_scheduler(0U);
  particle_scheduler.reset(0U, 0U, 0U);
  cosmosim::core::HierarchicalTimeBinScheduler gas_cell_scheduler(0U);
  gas_cell_scheduler.reset(static_cast<std::uint32_t>(state.cells.size()), 0U, 0U);
  cosmosim::core::syncGasCellTimeBinMirrorsFromGasCellScheduler(gas_cell_scheduler, state);

  cosmosim::io::RestartWritePayload payload;
  payload.persistent_state.simulation_state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &particle_scheduler;
  payload.gas_cell_scheduler = &gas_cell_scheduler;
  payload.normalized_config_text = "schema_version = 1\nmode = cosmo_cube\n";
  payload.normalized_config_hash_hex =
      cosmosim::core::stableConfigHashHex(payload.normalized_config_text);
  payload.provenance = cosmosim::core::makeProvenanceRecord(
      payload.normalized_config_hash_hex, "m2b2-test", 1, payload.normalized_config_text);
  payload.distributed_gravity_state.schema_version = 2U;
  payload.distributed_gravity_state.decomposition_epoch = 7U;
  payload.distributed_gravity_state.world_size = 2;
  payload.distributed_gravity_state.pm_grid_nx = 4U;
  payload.distributed_gravity_state.pm_grid_ny = 4U;
  payload.distributed_gravity_state.pm_grid_nz = 4U;
  payload.distributed_gravity_state.pm_decomposition_mode = "slab";
  payload.distributed_gravity_state.pm_slab_begin_x_by_rank = {0U, 2U};
  payload.distributed_gravity_state.pm_slab_end_x_by_rank = {2U, 4U};

  cosmosim::io::writeRestartCheckpointHdf5(path, payload);
  auto restored = cosmosim::io::readRestartCheckpointHdf5(path);
  std::filesystem::remove(path);
  return restored;
}

void continueMigrationContractAfterRestart(SimulationState& state) {
  const std::uint32_t migrated_patch_row = patchRow(state, 101U);
  auto records = state.packAmrPatchMigrationRecords(
      std::span<const std::uint32_t>(&migrated_patch_row, 1U));
  assert(records.size() == 1U);
  assert(records[0].pending_flux_register_records.size() == 1U);
  assert(records[0].temporal_boundary_history_records.size() == 1U);
  records[0].patch.owning_rank = 1U;
  for (auto& gas_record : records[0].gas_cell_records) {
    gas_record.owning_rank = 1U;
  }

  const std::uint64_t generation_before = state.gasCellIdentityGeneration();
  AmrPatchMigrationCommit commit;
  commit.world_rank = 1;
  commit.expected_gas_cell_identity_generation = generation_before;
  commit.outbound_local_patch_indices = {migrated_patch_row};
  commit.inbound_records = std::move(records);
  state.commitAmrPatchMigration(commit);
  assert(state.gasCellIdentityGeneration() == generation_before + 1U);
  assertMigratedSynchronizationState(state);
}

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

  assert(::setenv("COSMOSIM_MPI_TEST_TRANSPORT_LIMIT_BYTES", "8", 1) == 0);
  cosmosim::parallel::MpiContext mpi_context;
  SimulationState state = world_rank == 0 ? makeSourceState() : makeDestinationState();

  std::vector<std::vector<std::uint8_t>> send_payloads(static_cast<std::size_t>(world_size));
  if (world_rank == 0) {
    const std::uint32_t source_patch_row = patchRow(state, 101U);
    auto records = state.packAmrPatchMigrationRecords(
        std::span<const std::uint32_t>(&source_patch_row, 1U));
    assert(records.size() == 1U);
    records[0].patch.owning_rank = 1U;
    for (auto& gas_record : records[0].gas_cell_records) {
      gas_record.owning_rank = 1U;
    }
    send_payloads[1] = cosmosim::workflows::internal::migration_wire::encodeAmrPatchMigrationRecord(records[0]);
  }

  const auto received_payloads =
      cosmosim::parallel::exchangeBoundedAlltoallBytes(mpi_context, send_payloads);
  std::vector<AmrPatchMigrationRecord> inbound_records;
  if (world_rank == 1) {
    assert(!received_payloads[0].empty());
    inbound_records.push_back(
        cosmosim::workflows::internal::migration_wire::decodeAmrPatchMigrationRecord(received_payloads[0]));
    assert(inbound_records[0].patch.owning_rank == 1U);
    assert(inbound_records[0].pending_flux_register_records.size() == 1U);
    assert(inbound_records[0].temporal_boundary_history_records.size() == 1U);
  } else {
    assert(received_payloads[1].empty());
  }

  AmrPatchMigrationCommit commit;
  commit.world_rank = world_rank;
  commit.expected_gas_cell_identity_generation = state.gasCellIdentityGeneration();
  if (world_rank == 0) {
    commit.outbound_local_patch_indices = {patchRow(state, 101U)};
  } else {
    commit.inbound_records = std::move(inbound_records);
  }
  state.commitAmrPatchMigration(commit);

  if (world_rank == 0) {
    assert(state.patches.size() == 0U);
    assert(state.cells.size() == 0U);
    assert(state.pending_flux_registers.empty());
    assert(state.amr_temporal_boundary_history.empty());
  } else {
    assert(state.patches.size() == 2U);
    assert(state.cells.size() == 3U);
    assertMigratedSynchronizationState(state);
  }

  const std::uint64_t local_migrated_patch_count =
      std::count(state.patches.patch_id.begin(), state.patches.patch_id.end(), 101U);
  const std::uint64_t local_pending_count = state.pending_flux_registers.size();
  const std::uint64_t local_history_count = state.amr_temporal_boundary_history.size();
  assert(mpi_context.allreduceSumUint64(local_migrated_patch_count) == 1U);
  assert(mpi_context.allreduceSumUint64(local_pending_count) == 1U);
  assert(mpi_context.allreduceSumUint64(local_history_count) == 1U);

  auto workspace = cosmosim::test_support::createMpiSharedWorkspace(
      "cosmosim_distributed_amr_patch_migration_mpi");
  if (world_rank == 1) {
    auto restored = checkpointAndReload(state, workspace.root() / "rank_1_migrated_restart.hdf5");
    assertMigratedSynchronizationState(restored.state);
    continueMigrationContractAfterRestart(restored.state);
  }

  assert(::unsetenv("COSMOSIM_MPI_TEST_TRANSPORT_LIMIT_BYTES") == 0);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#else
  (void)argc;
  (void)argv;
#endif
  return 0;
}
