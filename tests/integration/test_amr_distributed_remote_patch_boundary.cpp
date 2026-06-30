#include <cassert>
#include <cmath>
#include <cstdint>
#include <optional>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace {

constexpr double k_gamma = 1.4;

void setPatch(
    cosmosim::core::SimulationState& state,
    std::size_t patch_index,
    const cosmosim::amr::PatchDescriptor& descriptor,
    std::uint32_t first_cell,
    std::uint32_t cell_count,
    std::uint32_t owner_rank) {
  cosmosim::amr::writePatchDescriptorToStateRow(state, patch_index, descriptor);
  state.patches.first_cell[patch_index] = first_cell;
  state.patches.cell_count[patch_index] = cell_count;
  state.patches.owning_rank[patch_index] = owner_rank;
}

void setCell(
    cosmosim::core::SimulationState& state,
    std::uint32_t row,
    double x,
    double rho,
    double vx,
    double pressure,
    std::uint32_t patch_index,
    std::uint64_t patch_id,
    std::uint64_t gas_cell_id,
    std::vector<cosmosim::core::GasCellIdentityRecord>& records) {
  state.cells.center_x_comoving[row] = x;
  state.cells.center_y_comoving[row] = 0.5;
  state.cells.center_z_comoving[row] = 0.5;
  state.cells.patch_index[row] = patch_index;
  state.cells.mass_code[row] = rho;
  state.gas_cells.gas_cell_id[row] = gas_cell_id;
  state.gas_cells.parent_particle_id[row] = 0U;
  state.gas_cells.density_code[row] = rho;
  state.gas_cells.velocity_x_peculiar[row] = vx;
  state.gas_cells.velocity_y_peculiar[row] = 0.0;
  state.gas_cells.velocity_z_peculiar[row] = 0.0;
  state.gas_cells.pressure_code[row] = pressure;
  state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
  state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
  state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * pressure / rho);
  records.push_back(cosmosim::core::GasCellIdentityRecord{
      .gas_cell_id = gas_cell_id,
      .parent_particle_id = std::nullopt,
      .owning_patch_id = patch_id,
      .local_cell_row = row});
}

[[nodiscard]] cosmosim::amr::PatchDescriptor coarsePatch() {
  return cosmosim::amr::PatchDescriptor{
      .patch_id = 101,
      .level = 0,
      .morton_key = 101,
      .origin_comov = {0.0, 0.0, 0.0},
      .extent_comov = {1.0, 1.0, 1.0},
      .cell_dims = {2, 1, 1}};
}

[[nodiscard]] cosmosim::amr::PatchDescriptor finePatch() {
  return cosmosim::amr::PatchDescriptor{
      .patch_id = 201,
      .parent_patch_id = 101,
      .level = 1,
      .morton_key = 201,
      .origin_comov = {1.0, 0.0, 0.0},
      .extent_comov = {0.5, 1.0, 1.0},
      .cell_dims = {2, 1, 1}};
}

[[nodiscard]] cosmosim::core::SimulationState makeFineRankState() {
  cosmosim::core::SimulationState state;
  state.resizeCells(2);
  state.resizePatches(1);
  setPatch(state, 0, finePatch(), 0, 2, 1);
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  setCell(state, 0, 1.125, 0.8, -0.25, 0.8, 0, 201, 9101, records);
  setCell(state, 1, 1.375, 0.7, -0.10, 0.7, 0, 201, 9102, records);
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

[[nodiscard]] cosmosim::amr::DistributedAmrRemotePatch makeRemoteCoarsePatch() {
  cosmosim::amr::DistributedAmrRemotePatch remote;
  remote.patch = coarsePatch();
  remote.owner_rank = 0;
  remote.ghost_hydro_epoch = 7;
  remote.expected_ghost_hydro_epoch = 7;
  remote.gas_cell_ids = {9001, 9002};
  remote.conserved_cells.push_back(cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(
      {.rho_comoving = 1.0, .vel_x_peculiar = 0.20, .pressure_comoving = 1.0}, k_gamma));
  remote.conserved_cells.push_back(cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(
      {.rho_comoving = 0.9, .vel_x_peculiar = 0.35, .pressure_comoving = 0.9}, k_gamma));
  return remote;
}

void testRemoteCoarseGhostProducesRemoteRefluxRegister() {
  cosmosim::core::SimulationState state = makeFineRankState();
  const cosmosim::amr::DistributedAmrRemotePatch remote = makeRemoteCoarsePatch();
  std::vector<cosmosim::amr::FluxRegisterEntry> outbound;
  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::HlleRiemannSolver riemann;
  const cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-5, .scale_factor = 1.0};
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  const std::vector<std::uint32_t> active_rows{0, 1};
  const auto diagnostics = cosmosim::amr::advanceDistributedProductionAmrHydro(
      state,
      active_rows,
      update,
      source_context,
      solver,
      riemann,
      {},
      cosmosim::amr::ProductionAmrHydroOptions{
          .physical_boundary_kind = cosmosim::hydro::HydroBoundaryKind::kOpen,
          .adiabatic_index = k_gamma,
          .density_floor = 1.0e-12,
          .pressure_floor = 1.0e-12,
          .state_time_code = 0.0,
          .ghost_fill_time_code = 0.0},
      cosmosim::amr::DistributedAmrHydroExchange{
          .local_rank = 1,
          .ghost_hydro_epoch = 7,
          .expected_ghost_hydro_epoch = 7,
          .remote_patches = std::span<const cosmosim::amr::DistributedAmrRemotePatch>(&remote, 1),
          .outbound_remote_flux_registers = &outbound});

  assert(diagnostics.advanced_patch_count == 1U);
  assert(diagnostics.ghost_fill.coarse_to_fine_ghosts_filled > 0U);
  assert(!outbound.empty());
  bool saw_remote_coarse_target = false;
  for (const auto& entry : outbound) {
    if (entry.coarse_patch_id == 101 && entry.coarse_gas_cell_id == 9002 &&
        entry.fine_face_count > 0U && entry.coarse_face_count == 0U) {
      saw_remote_coarse_target = true;
    }
  }
  assert(saw_remote_coarse_target);
}

void testDistributedSubcyclingIsGated() {
  cosmosim::core::SimulationState state = makeFineRankState();
  const cosmosim::amr::DistributedAmrRemotePatch remote = makeRemoteCoarsePatch();
  std::vector<cosmosim::amr::FluxRegisterEntry> outbound;
  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::HlleRiemannSolver riemann;
  const cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-5, .scale_factor = 1.0};
  bool threw = false;
  try {
    (void)cosmosim::amr::advanceDistributedProductionAmrHydro(
        state,
        {},
        update,
        {.update = update},
        solver,
        riemann,
        {},
        cosmosim::amr::ProductionAmrHydroOptions{
            .adiabatic_index = k_gamma,
            .sweep_mode = cosmosim::amr::ProductionAmrHydroSweepMode::kLocalSubcycling},
        cosmosim::amr::DistributedAmrHydroExchange{
            .local_rank = 1,
            .remote_patches = std::span<const cosmosim::amr::DistributedAmrRemotePatch>(&remote, 1),
            .outbound_remote_flux_registers = &outbound});
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}

}  // namespace

int main() {
  testRemoteCoarseGhostProducesRemoteRefluxRegister();
  testDistributedSubcyclingIsGated();
  return 0;
}
