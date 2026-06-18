#include <cassert>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace {

using cosmosim::core::AmrPatchMigrationCommit;
using cosmosim::core::ParticleSpecies;
using cosmosim::core::SimulationState;

constexpr std::uint32_t speciesTag(ParticleSpecies species) {
  return static_cast<std::uint32_t>(species);
}

void seedPatchState(SimulationState& state) {
  state.resizeParticles(1);
  state.particle_sidecar.particle_id[0] = 9001;
  state.particle_sidecar.species_tag[0] = speciesTag(ParticleSpecies::kGas);
  state.particle_sidecar.owning_rank[0] = 0;
  state.particle_sidecar.sfc_key[0] = 11;
  state.particles.mass_code[0] = 4.0;
  state.particles.time_bin[0] = 2;
  state.species.count_by_species = {0, 1, 0, 0, 0};
  state.rebuildSpeciesIndex();

  state.resizePatches(2);
  state.patches.patch_id[0] = 101;
  state.patches.level[0] = 1;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 2;
  state.patches.parent_patch_id[0] = 10;
  state.patches.morton_key[0] = 0xAA;
  state.patches.origin_x_comoving[0] = 0.0;
  state.patches.origin_y_comoving[0] = 0.0;
  state.patches.origin_z_comoving[0] = 0.0;
  state.patches.extent_x_comoving[0] = 2.0;
  state.patches.extent_y_comoving[0] = 1.0;
  state.patches.extent_z_comoving[0] = 1.0;
  state.patches.cell_dim_x[0] = 2;
  state.patches.cell_dim_y[0] = 1;
  state.patches.cell_dim_z[0] = 1;
  state.patches.owning_rank[0] = 0;

  state.patches.patch_id[1] = 202;
  state.patches.level[1] = 0;
  state.patches.first_cell[1] = 2;
  state.patches.cell_count[1] = 1;
  state.patches.parent_patch_id[1] = 0;
  state.patches.morton_key[1] = 0xBB;
  state.patches.origin_x_comoving[1] = 2.0;
  state.patches.origin_y_comoving[1] = 0.0;
  state.patches.origin_z_comoving[1] = 0.0;
  state.patches.extent_x_comoving[1] = 1.0;
  state.patches.extent_y_comoving[1] = 1.0;
  state.patches.extent_z_comoving[1] = 1.0;
  state.patches.cell_dim_x[1] = 1;
  state.patches.cell_dim_y[1] = 1;
  state.patches.cell_dim_z[1] = 1;
  state.patches.owning_rank[1] = 0;

  state.resizeCells(3);
  const std::array<std::uint64_t, 3> gas_ids{7001, 7002, 8001};
  const std::array<std::uint64_t, 3> parent_ids{9001, 0, 0};
  const std::array<std::uint32_t, 3> patch_indices{0, 0, 1};
  for (std::uint32_t row = 0; row < 3; ++row) {
    state.cells.center_x_comoving[row] = 0.25 + static_cast<double>(row);
    state.cells.center_y_comoving[row] = 0.5;
    state.cells.center_z_comoving[row] = 0.75;
    state.cells.mass_code[row] = 10.0 + static_cast<double>(row);
    state.cells.time_bin[row] = static_cast<std::uint8_t>(row + 1);
    state.cells.patch_index[row] = patch_indices[row];
    state.gas_cells.gas_cell_id[row] = gas_ids[row];
    state.gas_cells.parent_particle_id[row] = parent_ids[row];
    state.gas_cells.velocity_x_peculiar[row] = 1.0 + static_cast<double>(row);
    state.gas_cells.velocity_y_peculiar[row] = 2.0 + static_cast<double>(row);
    state.gas_cells.velocity_z_peculiar[row] = 3.0 + static_cast<double>(row);
    state.gas_cells.density_code[row] = 100.0 + static_cast<double>(row);
    state.gas_cells.pressure_code[row] = 200.0 + static_cast<double>(row);
    state.gas_cells.internal_energy_code[row] = 300.0 + static_cast<double>(row);
    state.gas_cells.temperature_code[row] = 400.0 + static_cast<double>(row);
    state.gas_cells.sound_speed_code[row] = 500.0 + static_cast<double>(row);
  }
  state.refreshGasCellIdentityMapFromSidecarLanes();
}

std::uint32_t rowForGasCell(const SimulationState& state, std::uint64_t gas_cell_id) {
  const auto row = state.rowForGasCellId(gas_cell_id);
  assert(row.has_value());
  return *row;
}

void test_patch_payload_carries_descriptor_and_gas_sidecars_atomically() {
  SimulationState state;
  seedPatchState(state);
  const std::uint64_t generation_before = state.gasCellIdentityGeneration();

  const std::array<std::uint32_t, 1> outbound_patch{0};
  auto records = state.packAmrPatchMigrationRecords(outbound_patch, 42);
  assert(records.size() == 1);
  assert(records[0].patch.patch_id == 101);
  assert(records[0].patch.parent_patch_id == 10);
  assert(records[0].patch.morton_key == 0xAA);
  assert(records[0].patch.origin_x_comoving == 0.0);
  assert(records[0].patch.origin_y_comoving == 0.0);
  assert(records[0].patch.origin_z_comoving == 0.0);
  assert(records[0].patch.extent_x_comoving == 2.0);
  assert(records[0].patch.extent_y_comoving == 1.0);
  assert(records[0].patch.extent_z_comoving == 1.0);
  assert(records[0].patch.cell_dim_x == 2U);
  assert(records[0].patch.cell_dim_y == 1U);
  assert(records[0].patch.cell_dim_z == 1U);
  assert(records[0].patch.cell_count == 2);
  assert(records[0].gas_cell_records.size() == 2);
  assert(records[0].gas_cell_records[0].fields.gas_cell_id == 7001);
  assert(records[0].gas_cell_records[0].fields.parent_particle_id == 9001);
  assert(records[0].gas_cell_records[1].fields.gas_cell_id == 7002);
  assert(records[0].gas_cell_records[1].fields.has_parent_particle == 0);
  assert(records[0].gas_cell_records[1].fields.ghost_hydro_epoch == 42);

  AmrPatchMigrationCommit commit;
  commit.world_rank = 0;
  commit.expected_gas_cell_identity_generation = state.gasCellIdentityGeneration();
  commit.expected_ghost_hydro_epoch = 42;
  commit.outbound_local_patch_indices = {0};
  commit.inbound_records = records;
  state.commitAmrPatchMigration(commit);

  assert(state.gasCellIdentityGeneration() == generation_before + 1);
  assert(state.patches.size() == 2);
  assert(state.patches.patch_id[0] == 202);
  assert(state.patches.patch_id[1] == 101);
  assert(state.patches.first_cell[0] == 0);
  assert(state.patches.first_cell[1] == 1);
  assert(state.patches.cell_count[1] == 2);
  assert(state.patches.parent_patch_id[1] == 10);
  assert(state.patches.morton_key[1] == 0xAA);
  assert(state.patches.origin_x_comoving[1] == 0.0);
  assert(state.patches.extent_x_comoving[1] == 2.0);
  assert(state.patches.cell_dim_x[1] == 2U);

  const std::uint32_t parented_row = rowForGasCell(state, 7001);
  const std::uint32_t parentless_row = rowForGasCell(state, 7002);
  assert(parented_row == 1U);
  assert(parentless_row == 2U);
  assert(state.cells.patch_index[parented_row] == 1);
  assert(state.cells.patch_index[parentless_row] == 1);
  assert(state.parentParticleIdForGasCellId(7001).value() == 9001);
  assert(!state.parentParticleIdForGasCellId(7002).has_value());
  assert(state.owningPatchIdForGasCellId(7001).value() == 101);
  assert(state.gas_cells.density_code[parented_row] == 100.0);
  assert(state.gas_cells.internal_energy_code[parentless_row] == 301.0);
  assert(state.validateOwnershipInvariants());
  assert(cosmosim::amr::hasProductionAmrHydroCoverage(state));
}

void test_patch_commit_rejects_stale_gas_ghost_epoch() {
  SimulationState state;
  seedPatchState(state);

  AmrPatchMigrationCommit commit;
  commit.world_rank = 0;
  commit.expected_gas_cell_identity_generation = state.gasCellIdentityGeneration();
  commit.expected_ghost_hydro_epoch = 7;
  commit.stale_local_ghost_records.push_back(cosmosim::core::GasCellStaleGhostRecord{
      .gas_cell_id = 7002,
      .local_cell_row = 1,
      .gas_cell_identity_generation = state.gasCellIdentityGeneration(),
      .ghost_hydro_epoch = 6,
  });

  bool threw = false;
  try {
    state.commitAmrPatchMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
  assert(state.rowForGasCellId(7002).has_value());
}

}  // namespace

int main() {
  test_patch_payload_carries_descriptor_and_gas_sidecars_atomically();
  test_patch_commit_rejects_stale_gas_ghost_epoch();
  return 0;
}
