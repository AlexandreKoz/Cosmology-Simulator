#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "restart_equivalence_harness.hpp"
#include "restart_equivalence_scenarios.hpp"

namespace {

constexpr double k_gamma = 1.4;

void setPatch(
    cosmosim::core::SimulationState& state,
    std::size_t patch_index,
    const cosmosim::amr::PatchDescriptor& descriptor,
    std::uint32_t first_cell,
    std::uint32_t cell_count) {
  cosmosim::amr::writePatchDescriptorToStateRow(state, patch_index, descriptor);
  state.patches.first_cell[patch_index] = first_cell;
  state.patches.cell_count[patch_index] = cell_count;
  state.patches.owning_rank[patch_index] = 0;
}

void setGasParticle(
    cosmosim::core::SimulationState& state,
    std::uint32_t row,
    std::uint64_t particle_id,
    double x,
    double y,
    double z,
    double mass,
    double vx,
    double vy,
    double vz) {
  state.particle_sidecar.particle_id[row] = particle_id;
  state.particle_sidecar.sfc_key[row] = 0xB000U + row;
  state.particle_sidecar.species_tag[row] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
  state.particle_sidecar.owning_rank[row] = 0;
  state.particles.position_x_comoving[row] = x;
  state.particles.position_y_comoving[row] = y;
  state.particles.position_z_comoving[row] = z;
  state.particles.velocity_x_peculiar[row] = vx;
  state.particles.velocity_y_peculiar[row] = vy;
  state.particles.velocity_z_peculiar[row] = vz;
  state.particles.mass_code[row] = mass;
}

void setCell(
    cosmosim::core::SimulationState& state,
    std::uint32_t row,
    double x,
    double y,
    double z,
    double cell_volume,
    double rho,
    double velocity_x,
    double velocity_y,
    double velocity_z,
    double pressure,
    std::uint32_t patch_index,
    std::uint64_t patch_id,
    std::uint64_t gas_cell_id,
    std::uint64_t parent_particle_id,
    std::vector<cosmosim::core::GasCellIdentityRecord>& records) {
  state.cells.center_x_comoving[row] = x;
  state.cells.center_y_comoving[row] = y;
  state.cells.center_z_comoving[row] = z;
  state.cells.patch_index[row] = patch_index;
  state.cells.time_bin[row] = 0;
  state.cells.mass_code[row] = rho * cell_volume;
  state.gas_cells.gas_cell_id[row] = gas_cell_id;
  state.gas_cells.parent_particle_id[row] = parent_particle_id;
  state.gas_cells.density_code[row] = rho;
  state.gas_cells.pressure_code[row] = pressure;
  state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
  state.gas_cells.velocity_x_peculiar[row] = velocity_x;
  state.gas_cells.velocity_y_peculiar[row] = velocity_y;
  state.gas_cells.velocity_z_peculiar[row] = velocity_z;
  state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
  state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * pressure / rho);
  records.push_back(cosmosim::core::GasCellIdentityRecord{
      .gas_cell_id = gas_cell_id,
      .parent_particle_id = parent_particle_id,
      .owning_patch_id = patch_id,
      .local_cell_row = row});
}

[[nodiscard]] cosmosim::core::SimulationState makeProductionAmrHydroRestartState(const std::string& run_name) {
  cosmosim::core::SimulationState state;
  state.metadata.run_name = run_name;
  state.metadata.normalized_config_hash_hex = run_name;
  state.metadata.snapshot_stem = "snapshot";
  state.metadata.restart_stem = "restart";

  state.resizeParticles(4);
  state.resizeCells(4);
  state.resizePatches(2);
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)] = 4;

  setPatch(state, 0, cosmosim::amr::PatchDescriptor{
      .patch_id = 101,
      .parent_patch_id = 0,
      .level = 0,
      .morton_key = 0x101,
      .origin_comov = {0.0, 0.0, 0.0},
      .extent_comov = {1.0, 1.0, 1.0},
      .cell_dims = {2, 1, 1}}, 0, 2);
  setPatch(state, 1, cosmosim::amr::PatchDescriptor{
      .patch_id = 201,
      .parent_patch_id = 101,
      .level = 1,
      .morton_key = 0x201,
      .origin_comov = {1.0, 0.0, 0.0},
      .extent_comov = {0.5, 1.0, 1.0},
      .cell_dims = {2, 1, 1}}, 2, 2);

  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(4);
  const double coarse_volume = 0.5;
  const double fine_volume = 0.25;
  setCell(state, 0, 0.25, 0.5, 0.5, coarse_volume, 1.00, 0.20, 0.03, -0.02, 1.00, 0, 101, 9001, 7001, records);
  setCell(state, 1, 0.75, 0.5, 0.5, coarse_volume, 0.92, 0.33, -0.04, 0.01, 0.93, 0, 101, 9002, 7002, records);
  setCell(state, 2, 1.125, 0.5, 0.5, fine_volume, 0.84, -0.24, 0.02, 0.04, 0.82, 1, 201, 9101, 7003, records);
  setCell(state, 3, 1.375, 0.5, 0.5, fine_volume, 0.73, -0.11, -0.03, -0.01, 0.74, 1, 201, 9102, 7004, records);

  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    setGasParticle(
        state,
        row,
        7001U + row,
        state.cells.center_x_comoving[row],
        state.cells.center_y_comoving[row],
        state.cells.center_z_comoving[row],
        state.cells.mass_code[row],
        state.gas_cells.velocity_x_peculiar[row],
        state.gas_cells.velocity_y_peculiar[row],
        state.gas_cells.velocity_z_peculiar[row]);
  }

  state.gas_cell_identity.assign(std::move(records));
  state.rebuildSpeciesIndex();
  assert(state.validateOwnershipInvariants());
  assert(cosmosim::amr::hasProductionAmrHydroCoverage(state));
  return state;
}

void applyProductionAmrHydroStep(cosmosim::tests::RestartEquivalenceStepContext& context) {
  assert(cosmosim::amr::hasProductionAmrHydroCoverage(context.state));

  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::HlleRiemannSolver riemann;
  const cosmosim::hydro::HydroUpdateContext update{
      .dt_code = 1.0e-4,
      .scale_factor = context.integrator_state.current_scale_factor,
      .hubble_rate_code = 0.0};
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  const std::vector<std::uint32_t> active_rows{0, 1, 2, 3};
  const cosmosim::amr::ProductionAmrHydroOptions options{
      .physical_boundary_kind = cosmosim::hydro::HydroBoundaryKind::kOpen,
      .adiabatic_index = k_gamma,
      .density_floor = 1.0e-12,
      .pressure_floor = 1.0e-12};

  const auto diagnostics = cosmosim::amr::advanceProductionAmrHydro(
      context.state,
      active_rows,
      update,
      source_context,
      solver,
      riemann,
      {},
      options);

  assert(diagnostics.patch_count == 2U);
  assert(diagnostics.advanced_patch_count == 2U);
  assert(diagnostics.flux_register_entry_count > 0U);
  assert(
      diagnostics.ghost_fill.same_level_ghosts_filled + diagnostics.ghost_fill.coarse_to_fine_ghosts_filled +
          diagnostics.ghost_fill.fine_to_coarse_ghosts_filled >
      0U);
  assert(context.state.gas_cell_identity.coversDenseLocalRows(context.state.cells.size()));
  assert(context.state.gasCellIdentityMapMatchesSidecarLanes());

  for (std::uint32_t row = 0; row < context.state.cells.size(); ++row) {
    assert(context.state.gas_cells.gas_cell_id[row] != 0U);
    assert(context.state.gas_cells.density_code[row] > 0.0);
    assert(context.state.gas_cells.pressure_code[row] > 0.0);
    assert(context.state.gas_cells.internal_energy_code[row] > 0.0);
    assert(std::isfinite(context.state.gas_cells.velocity_x_peculiar[row]));
    assert(std::isfinite(context.state.gas_cells.velocity_y_peculiar[row]));
    assert(std::isfinite(context.state.gas_cells.velocity_z_peculiar[row]));
  }
}

void assertExplicitPatchGeometrySurvived(const cosmosim::core::SimulationState& state) {
  assert(cosmosim::amr::hasProductionAmrHydroCoverage(state));
  assert(state.patches.size() == 2U);
  assert(state.patches.patch_id[0] == 101U);
  assert(state.patches.patch_id[1] == 201U);
  assert(state.patches.parent_patch_id[1] == 101U);
  assert(state.patches.morton_key[0] == 0x101U);
  assert(state.patches.morton_key[1] == 0x201U);
  assert(state.patches.origin_x_comoving[1] == 1.0);
  assert(state.patches.extent_x_comoving[1] == 0.5);
  assert(state.patches.cell_dim_x[0] == 2U);
  assert(state.patches.cell_dim_y[0] == 1U);
  assert(state.patches.cell_dim_z[0] == 1U);
  assert(state.patches.cell_dim_x[1] == 2U);
  assert(state.patches.cell_dim_y[1] == 1U);
  assert(state.patches.cell_dim_z[1] == 1U);
  assert(state.gas_cell_identity.coversDenseLocalRows(state.cells.size()));
  assert(state.gasCellIdentityMapMatchesSidecarLanes());
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_HDF5
  const auto restart_path = cosmosim::tests::stage8RestartPath("restart_equivalence_amr_hydro");
  auto state = makeProductionAmrHydroRestartState("restart_equivalence_amr_hydro");
  auto scheduler = cosmosim::tests::makeStage8Scheduler(static_cast<std::uint32_t>(state.particles.size()), 2);
  scheduler.setElementBin(2, 1, scheduler.currentTick());
  auto integrator_state = cosmosim::tests::makeStage8IntegratorState(2, 2);
  auto output_state = cosmosim::tests::makeStage8OutputCadenceState(false);
  auto scenario = cosmosim::tests::makeStage8Scenario(
      std::move(state), integrator_state, std::move(scheduler), std::move(output_state), restart_path, 8, 3);
  scenario.step_kernel = applyProductionAmrHydroStep;
  scenario.tolerances.position_abs = 1.0e-12;
  scenario.tolerances.velocity_abs = 1.0e-12;
  scenario.tolerances.scalar_abs = 1.0e-12;

  const auto result = cosmosim::tests::runRestartEquivalenceScenario(std::move(scenario));
  assertExplicitPatchGeometrySurvived(result.direct_state);
  assertExplicitPatchGeometrySurvived(result.restarted_state);
  assert(result.direct_scheduler_state.bin_index == result.restarted_scheduler_state.bin_index);
  assert(result.direct_integrator_state.step_index == result.restarted_integrator_state.step_index);
  assert(result.direct_state.gas_cell_identity.records().size() == 4U);
  assert(result.restarted_state.gas_cell_identity.records().size() == 4U);
  std::filesystem::remove(restart_path);
#endif
  return 0;
}
