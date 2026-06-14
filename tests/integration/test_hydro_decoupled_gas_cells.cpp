#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace {

constexpr double k_gamma = 1.4;

cosmosim::hydro::HydroPatchGeometry makeGeometry() {
  return cosmosim::hydro::makeCartesianPatchGeometry(cosmosim::hydro::HydroCartesianPatchSpec{
      .nx = 3,
      .ny = 1,
      .nz = 1,
      .origin_x_comoving = 0.0,
      .origin_y_comoving = 0.0,
      .origin_z_comoving = 0.0,
      .cell_width_x_comoving = 1.0,
      .cell_width_y_comoving = 1.0,
      .cell_width_z_comoving = 1.0});
}

std::unordered_map<std::uint64_t, std::uint32_t> particleRowById(
    const cosmosim::core::SimulationState& state) {
  std::unordered_map<std::uint64_t, std::uint32_t> rows;
  for (std::uint32_t row = 0; row < state.particles.size(); ++row) {
    rows.emplace(state.particle_sidecar.particle_id[row], row);
  }
  return rows;
}

std::optional<std::uint32_t> parentRowForCell(
    const cosmosim::core::SimulationState& state,
    std::uint32_t cell_row,
    const std::unordered_map<std::uint64_t, std::uint32_t>& particle_row_by_id) {
  const auto* record = state.gas_cell_identity.findByLocalRow(cell_row);
  if (record == nullptr) {
    throw std::runtime_error("parentRowForCell: missing gas-cell identity row");
  }
  if (!record->parent_particle_id.has_value()) {
    return std::nullopt;
  }
  const auto found = particle_row_by_id.find(*record->parent_particle_id);
  if (found == particle_row_by_id.end()) {
    throw std::runtime_error("parentRowForCell: parent particle is not local");
  }
  return found->second;
}

cosmosim::core::SimulationState makeDecoupledGasState() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(2);
  state.resizeCells(3);
  state.resizePatches(1);
  state.patches.patch_id[0] = 9001;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 3;
  state.patches.owning_rank[0] = 0;

  state.particle_sidecar.particle_id = {5001, 7001};
  state.particle_sidecar.species_tag = {
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter)};
  state.particle_sidecar.sfc_key = {20, 10};
  state.particle_sidecar.owning_rank = {0, 0};
  state.species.count_by_species = {1, 1, 0, 0, 0};
  state.rebuildSpeciesIndex();

  const auto geometry = makeGeometry();
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    state.cells.center_x_comoving[row] = static_cast<double>(row) + 0.5;
    state.cells.center_y_comoving[row] = 0.5;
    state.cells.center_z_comoving[row] = 0.5;
    state.cells.mass_code[row] = 1.0 + 0.25 * static_cast<double>(row);
    state.cells.patch_index[row] = 0;
    state.gas_cells.density_code[row] = state.cells.mass_code[row] / geometry.cell_volume_comoving;
    state.gas_cells.pressure_code[row] = 1.0 + 0.1 * static_cast<double>(row);
    state.gas_cells.velocity_x_peculiar[row] = 0.15 + 0.05 * static_cast<double>(row);
    state.gas_cells.velocity_y_peculiar[row] = -0.02 * static_cast<double>(row + 1);
    state.gas_cells.velocity_z_peculiar[row] = 0.01 * static_cast<double>(row + 1);
    state.gas_cells.internal_energy_code[row] =
        state.gas_cells.pressure_code[row] / ((k_gamma - 1.0) * state.gas_cells.density_code[row]);
    state.gas_cells.temperature_code[row] = state.gas_cells.pressure_code[row] / state.gas_cells.density_code[row];
    state.gas_cells.sound_speed_code[row] =
        std::sqrt(k_gamma * state.gas_cells.pressure_code[row] / state.gas_cells.density_code[row]);
  }
  state.particles.mass_code[0] = state.cells.mass_code[1];
  state.particles.velocity_x_peculiar[0] = state.gas_cells.velocity_x_peculiar[1];
  state.particles.velocity_y_peculiar[0] = state.gas_cells.velocity_y_peculiar[1];
  state.particles.velocity_z_peculiar[0] = state.gas_cells.velocity_z_peculiar[1];

  state.gas_cells.gas_cell_id = {8001, 8002, 8003};
  state.gas_cells.parent_particle_id = {0, 5001, 5001};
  state.gas_cell_identity.assign({
      {.gas_cell_id = 8001, .parent_particle_id = std::nullopt, .owning_patch_id = 9001, .local_cell_row = 0},
      {.gas_cell_id = 8002, .parent_particle_id = 5001, .owning_patch_id = 9001, .local_cell_row = 1},
      {.gas_cell_id = 8003, .parent_particle_id = 5001, .owning_patch_id = 9001, .local_cell_row = 2},
  });
  assert(state.validateOwnershipInvariants());
  return state;
}

std::array<std::uint32_t, 3> runHydroStep(cosmosim::core::SimulationState& state) {
  const auto geometry = makeGeometry();
  cosmosim::hydro::HydroConservedStateSoa conserved(state.cells.size());
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const cosmosim::hydro::HydroPrimitiveState primitive{
        .rho_comoving = state.gas_cells.density_code[row],
        .vel_x_peculiar = state.gas_cells.velocity_x_peculiar[row],
        .vel_y_peculiar = state.gas_cells.velocity_y_peculiar[row],
        .vel_z_peculiar = state.gas_cells.velocity_z_peculiar[row],
        .pressure_comoving = state.gas_cells.pressure_code[row],
    };
    conserved.storeCell(row, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
  }

  const cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-3, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code / geometry.cell_width_x_comoving,
      .dt_over_cell_width_code = {
          update.dt_code / geometry.cell_width_x_comoving,
          update.dt_code / geometry.cell_width_y_comoving,
          update.dt_code / geometry.cell_width_z_comoving},
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = true});
  cosmosim::hydro::HllcRiemannSolver riemann_solver;
  solver.advancePatch(conserved, geometry, update, reconstruction, riemann_solver, {}, source_context, nullptr);

  const auto parent_rows = particleRowById(state);
  std::array<std::uint32_t, 3> parent_mirror_write_count{0, 0, 0};
  std::vector<std::uint8_t> parent_mirror_updated(state.particles.size(), 0U);
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const auto primitive =
        cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(row), k_gamma);
    state.gas_cells.density_code[row] = primitive.rho_comoving;
    state.gas_cells.pressure_code[row] = primitive.pressure_comoving;
    state.gas_cells.velocity_x_peculiar[row] = primitive.vel_x_peculiar;
    state.gas_cells.velocity_y_peculiar[row] = primitive.vel_y_peculiar;
    state.gas_cells.velocity_z_peculiar[row] = primitive.vel_z_peculiar;
    state.gas_cells.internal_energy_code[row] =
        primitive.pressure_comoving / ((k_gamma - 1.0) * std::max(primitive.rho_comoving, 1.0e-10));
    state.gas_cells.temperature_code[row] = primitive.pressure_comoving / std::max(primitive.rho_comoving, 1.0e-10);
    state.gas_cells.sound_speed_code[row] =
        std::sqrt(k_gamma * primitive.pressure_comoving / std::max(primitive.rho_comoving, 1.0e-10));
    state.cells.mass_code[row] = primitive.rho_comoving * geometry.cell_volume_comoving;

    const auto parent_row = parentRowForCell(state, row, parent_rows);
    if (parent_row.has_value() && parent_mirror_updated[*parent_row] == 0U) {
      state.particles.mass_code[*parent_row] = state.cells.mass_code[row];
      state.particles.velocity_x_peculiar[*parent_row] = primitive.vel_x_peculiar;
      state.particles.velocity_y_peculiar[*parent_row] = primitive.vel_y_peculiar;
      state.particles.velocity_z_peculiar[*parent_row] = primitive.vel_z_peculiar;
      parent_mirror_updated[*parent_row] = 1U;
      ++parent_mirror_write_count[*parent_row];
    }
  }
  return parent_mirror_write_count;
}

void reorderGasRows(cosmosim::core::SimulationState& state, const std::array<std::uint32_t, 3>& new_to_old) {
  auto old_cells = state.cells;
  auto old_gas = state.gas_cells;
  const auto old_map = state.gas_cell_identity;
  for (std::uint32_t new_row = 0; new_row < new_to_old.size(); ++new_row) {
    const std::uint32_t old_row = new_to_old[new_row];
    state.cells.center_x_comoving[new_row] = old_cells.center_x_comoving[old_row];
    state.cells.center_y_comoving[new_row] = old_cells.center_y_comoving[old_row];
    state.cells.center_z_comoving[new_row] = old_cells.center_z_comoving[old_row];
    state.cells.mass_code[new_row] = old_cells.mass_code[old_row];
    state.cells.time_bin[new_row] = old_cells.time_bin[old_row];
    state.cells.patch_index[new_row] = old_cells.patch_index[old_row];
    state.gas_cells.gas_cell_id[new_row] = old_gas.gas_cell_id[old_row];
    state.gas_cells.parent_particle_id[new_row] = old_gas.parent_particle_id[old_row];
    state.gas_cells.velocity_x_peculiar[new_row] = old_gas.velocity_x_peculiar[old_row];
    state.gas_cells.velocity_y_peculiar[new_row] = old_gas.velocity_y_peculiar[old_row];
    state.gas_cells.velocity_z_peculiar[new_row] = old_gas.velocity_z_peculiar[old_row];
    state.gas_cells.density_code[new_row] = old_gas.density_code[old_row];
    state.gas_cells.pressure_code[new_row] = old_gas.pressure_code[old_row];
    state.gas_cells.internal_energy_code[new_row] = old_gas.internal_energy_code[old_row];
    state.gas_cells.temperature_code[new_row] = old_gas.temperature_code[old_row];
    state.gas_cells.sound_speed_code[new_row] = old_gas.sound_speed_code[old_row];
  }

  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  for (std::uint32_t new_row = 0; new_row < new_to_old.size(); ++new_row) {
    const auto* old_record = old_map.findByLocalRow(new_to_old[new_row]);
    assert(old_record != nullptr);
    records.push_back(cosmosim::core::GasCellIdentityRecord{
        .gas_cell_id = old_record->gas_cell_id,
        .parent_particle_id = old_record->parent_particle_id,
        .owning_patch_id = old_record->owning_patch_id,
        .local_cell_row = new_row,
    });
  }
  state.gas_cell_identity.assign(std::move(records));
  assert(state.gasCellIdentityMapMatchesSidecarLanes());
}

void testProductionHydroStateSupportsDecoupledGasCellIdentity() {
  cosmosim::core::SimulationState state = makeDecoupledGasState();
  bool legacy_contract_rejected = false;
  try {
    cosmosim::core::legacyRequireParticleBoundGasCellContract(state, "decoupled hydro regression");
  } catch (const std::runtime_error&) {
    legacy_contract_rejected = true;
  }
  assert(legacy_contract_rejected);

  const auto by_sfc = cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kBySfcKey);
  cosmosim::core::reorderParticles(state, by_sfc);
  assert(state.gas_cell_identity.coversDenseLocalRows(state.cells.size()));

  const auto first_mirror_writes = runHydroStep(state);
  std::uint32_t gas_parent_row = particleRowById(state).at(5001);
  assert(first_mirror_writes[gas_parent_row] == 1);
  assert(std::abs(state.gas_cells.velocity_x_peculiar[0]) > 0.0);
  assert(state.gas_cell_identity.rowForGasCellId(8001).value() == 0);

  reorderGasRows(state, {2, 0, 1});
  assert(state.gas_cell_identity.rowForGasCellId(8001).value() == 1);
  assert(state.gas_cell_identity.rowForGasCellId(8002).value() == 2);
  assert(state.gas_cell_identity.rowForGasCellId(8003).value() == 0);

  const auto second_mirror_writes = runHydroStep(state);
  gas_parent_row = particleRowById(state).at(5001);
  assert(second_mirror_writes[gas_parent_row] == 1);
  assert(state.gas_cell_identity.findByGasCellId(8001)->parent_particle_id == std::nullopt);
  assert(state.gas_cell_identity.rowsForParentParticleId(5001).size() == 2);
  assert(state.validateOwnershipInvariants());
}

}  // namespace

int main() {
  testProductionHydroStateSupportsDecoupledGasCellIdentity();
  return 0;
}
