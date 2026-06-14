#include <cassert>
#include <array>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"

namespace {

using cosmosim::core::GasCellIdentityRecord;
using cosmosim::core::GasCellMigrationCommit;
using cosmosim::core::GasCellMigrationRecord;
using cosmosim::core::GasCellStaleGhostRecord;
using cosmosim::core::ParticleSpecies;
using cosmosim::core::SimulationState;

constexpr std::uint32_t speciesTag(ParticleSpecies species) {
  return static_cast<std::uint32_t>(species);
}

void mirrorIdentityToSidecars(SimulationState& state) {
  for (const GasCellIdentityRecord& record : state.gas_cell_identity.records()) {
    const std::uint32_t row = record.local_cell_row;
    state.gas_cells.gas_cell_id[row] = record.gas_cell_id;
    state.gas_cells.parent_particle_id[row] = record.parent_particle_id.value_or(0U);
  }
}

void seedGasFields(SimulationState& state) {
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    state.cells.center_x_comoving[row] = 10.0 + static_cast<double>(row);
    state.cells.center_y_comoving[row] = 20.0 + static_cast<double>(row);
    state.cells.center_z_comoving[row] = 30.0 + static_cast<double>(row);
    state.cells.mass_code[row] = 40.0 + static_cast<double>(row);
    state.cells.time_bin[row] = static_cast<std::uint8_t>(row + 1U);
    state.gas_cells.velocity_x_peculiar[row] = 50.0 + static_cast<double>(row);
    state.gas_cells.velocity_y_peculiar[row] = 60.0 + static_cast<double>(row);
    state.gas_cells.velocity_z_peculiar[row] = 70.0 + static_cast<double>(row);
    state.gas_cells.density_code[row] = 80.0 + static_cast<double>(row);
    state.gas_cells.pressure_code[row] = 90.0 + static_cast<double>(row);
    state.gas_cells.internal_energy_code[row] = 100.0 + static_cast<double>(row);
    state.gas_cells.temperature_code[row] = 110.0 + static_cast<double>(row);
    state.gas_cells.sound_speed_code[row] = 120.0 + static_cast<double>(row);
  }
}

void test_parentless_gas_cell_migrates_without_particle() {
  SimulationState source;
  source.resizeCells(2);
  source.gas_cell_identity.assign({
      {.gas_cell_id = 8101, .parent_particle_id = 42, .owning_patch_id = 0, .local_cell_row = 0},
      {.gas_cell_id = 8102, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 1},
  });
  mirrorIdentityToSidecars(source);
  seedGasFields(source);
  assert(source.validateOwnershipInvariants());

  const std::array<std::uint32_t, 1> rows{1};
  auto records = source.packGasCellMigrationRecords(rows, 77);
  assert(records.size() == 1);
  assert(records[0].fields.gas_cell_id == 8102U);
  assert(records[0].fields.has_parent_particle == 0U);
  assert(records[0].fields.parent_particle_id == 0U);
  assert(records[0].fields.ghost_hydro_epoch == 77U);

  SimulationState destination;
  destination.resizeCells(0);
  records[0].owning_rank = 0;
  GasCellMigrationCommit commit;
  commit.world_rank = 0;
  commit.inbound_records = records;
  destination.commitGasCellMigration(commit);

  assert(destination.cells.size() == 1);
  assert(destination.gas_cells.gas_cell_id[0] == 8102U);
  assert(destination.gas_cells.parent_particle_id[0] == 0U);
  assert(destination.gas_cells.density_code[0] == 81.0);
  assert(!destination.parentParticleIdForGasCellRow(0).has_value());
  assert(destination.validateOwnershipInvariants());
}

void test_multi_parent_cells_migrate_by_gas_cell_id_not_parent_id() {
  SimulationState source;
  source.resizeCells(2);
  source.gas_cell_identity.assign({
      {.gas_cell_id = 8201, .parent_particle_id = 9001, .owning_patch_id = 0, .local_cell_row = 0},
      {.gas_cell_id = 8202, .parent_particle_id = 9001, .owning_patch_id = 0, .local_cell_row = 1},
  });
  mirrorIdentityToSidecars(source);
  seedGasFields(source);
  assert(source.gas_cell_identity.rowsForParentParticleId(9001) == std::vector<std::uint32_t>({0, 1}));

  const std::array<std::uint32_t, 2> rows{0, 1};
  auto records = source.packGasCellMigrationRecords(rows, 3);
  std::swap(records[0], records[1]);
  for (auto& record : records) {
    record.owning_rank = 0;
  }

  SimulationState destination;
  destination.resizeCells(0);
  GasCellMigrationCommit commit;
  commit.world_rank = 0;
  commit.inbound_records = records;
  destination.commitGasCellMigration(commit);

  assert(destination.cells.size() == 2);
  assert(destination.gas_cells.gas_cell_id[0] == 8202U);
  assert(destination.gas_cells.gas_cell_id[1] == 8201U);
  assert(destination.gas_cell_identity.rowsForParentParticleId(9001) == std::vector<std::uint32_t>({0, 1}));
  assert(destination.gas_cells.density_code[destination.rowForGasCellId(8201).value()] == 80.0);
  assert(destination.gas_cells.density_code[destination.rowForGasCellId(8202).value()] == 81.0);
  assert(destination.validateOwnershipInvariants());
}

void test_row_remap_preserves_host_cell_sidecars() {
  SimulationState state;
  state.resizeParticles(1);
  state.particle_sidecar.particle_id[0] = 9301;
  state.particle_sidecar.species_tag[0] = speciesTag(ParticleSpecies::kTracer);
  state.particle_sidecar.owning_rank[0] = 0;
  state.species.count_by_species = {0, 0, 0, 0, 1};
  state.rebuildSpeciesIndex();
  state.tracers.resize(1);
  state.tracers.particle_index[0] = 0;
  state.tracers.host_cell_index[0] = 2;

  state.resizeCells(3);
  state.gas_cell_identity.assign({
      {.gas_cell_id = 8301, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 0},
      {.gas_cell_id = 8302, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 1},
      {.gas_cell_id = 8303, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 2},
  });
  mirrorIdentityToSidecars(state);
  seedGasFields(state);
  assert(state.validateOwnershipInvariants());

  GasCellMigrationCommit commit;
  commit.world_rank = 0;
  commit.outbound_local_cell_rows = {0};
  state.commitGasCellMigration(commit);

  assert(state.cells.size() == 2);
  assert(state.gas_cells.gas_cell_id[0] == 8302U);
  assert(state.gas_cells.gas_cell_id[1] == 8303U);
  assert(state.tracers.host_cell_index[0] == 1U);
  assert(state.validateOwnershipInvariants());
}

void test_stale_gas_ghost_rejection_checks_id_generation_and_epoch() {
  SimulationState state;
  state.resizeCells(1);
  state.gas_cell_identity.assign({
      {.gas_cell_id = 8401, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 0},
  });
  mirrorIdentityToSidecars(state);
  seedGasFields(state);
  const std::uint64_t generation = state.gasCellIdentityGeneration();

  GasCellMigrationCommit commit;
  commit.world_rank = 0;
  commit.expected_gas_cell_identity_generation = generation;
  commit.expected_ghost_hydro_epoch = 9;
  commit.stale_local_ghost_records = {GasCellStaleGhostRecord{
      .gas_cell_id = 8401,
      .local_cell_row = 0,
      .gas_cell_identity_generation = generation + 1U,
      .ghost_hydro_epoch = 9,
  }};
  bool threw = false;
  try {
    state.commitGasCellMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);

  commit.stale_local_ghost_records[0].gas_cell_identity_generation = generation;
  commit.stale_local_ghost_records[0].ghost_hydro_epoch = 8;
  threw = false;
  try {
    state.commitGasCellMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);

  commit.stale_local_ghost_records[0].ghost_hydro_epoch = 9;
  commit.stale_local_ghost_records[0].gas_cell_id = 9999;
  threw = false;
  try {
    state.commitGasCellMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
  assert(state.cells.size() == 1);
  assert(state.validateOwnershipInvariants());
}

}  // namespace

int main() {
  test_parentless_gas_cell_migrates_without_particle();
  test_multi_parent_cells_migrate_by_gas_cell_id_not_parent_id();
  test_row_remap_preserves_host_cell_sidecars();
  test_stale_gas_ghost_rejection_checks_id_generation_and_epoch();
  return 0;
}
