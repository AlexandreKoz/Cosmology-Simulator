#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/build_config.hpp"
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
    double mass) {
  state.particle_sidecar.particle_id[row] = particle_id;
  state.particle_sidecar.sfc_key[row] = 0xCAFE000U + row;
  state.particle_sidecar.species_tag[row] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
  state.particle_sidecar.owning_rank[row] = 0;
  state.particles.position_x_comoving[row] = x;
  state.particles.position_y_comoving[row] = y;
  state.particles.position_z_comoving[row] = z;
  state.particles.velocity_x_peculiar[row] = 0.0;
  state.particles.velocity_y_peculiar[row] = 0.0;
  state.particles.velocity_z_peculiar[row] = 0.0;
  state.particles.mass_code[row] = mass;
}

void setCell(
    cosmosim::core::SimulationState& state,
    std::uint32_t row,
    double x,
    double rho,
    double pressure,
    std::uint32_t patch_index,
    std::uint64_t patch_id,
    std::uint64_t gas_cell_id,
    std::uint64_t parent_particle_id,
    std::vector<cosmosim::core::GasCellIdentityRecord>& records) {
  state.cells.center_x_comoving[row] = x;
  state.cells.center_y_comoving[row] = 0.5;
  state.cells.center_z_comoving[row] = 0.5;
  state.cells.patch_index[row] = patch_index;
  state.cells.mass_code[row] = rho;
  state.cells.time_bin[row] = 0;
  state.gas_cells.gas_cell_id[row] = gas_cell_id;
  state.gas_cells.parent_particle_id[row] = parent_particle_id;
  state.gas_cells.density_code[row] = rho;
  state.gas_cells.pressure_code[row] = pressure;
  state.gas_cells.internal_energy_code[row] = pressure / ((k_gamma - 1.0) * rho);
  state.gas_cells.velocity_x_peculiar[row] = 0.0;
  state.gas_cells.velocity_y_peculiar[row] = 0.0;
  state.gas_cells.velocity_z_peculiar[row] = 0.0;
  state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
  state.gas_cells.sound_speed_code[row] = std::sqrt(k_gamma * pressure / rho);
  records.push_back(cosmosim::core::GasCellIdentityRecord{
      .gas_cell_id = gas_cell_id,
      .parent_particle_id = parent_particle_id,
      .owning_patch_id = patch_id,
      .local_cell_row = row});
}

[[nodiscard]] cosmosim::core::SimulationState makeRestartPendingFluxState(const std::string& run_name) {
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
      .morton_key = 101,
      .origin_comov = {0.0, 0.0, 0.0},
      .extent_comov = {1.0, 1.0, 1.0},
      .cell_dims = {2, 1, 1}}, 0, 2);
  setPatch(state, 1, cosmosim::amr::PatchDescriptor{
      .patch_id = 201,
      .parent_patch_id = 101,
      .level = 1,
      .morton_key = 201,
      .origin_comov = {1.0, 0.0, 0.0},
      .extent_comov = {0.5, 1.0, 1.0},
      .cell_dims = {2, 1, 1}}, 2, 2);

  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(4);
  setCell(state, 0, 0.25, 1.0, 1.0, 0, 101, 9001, 7001, records);
  setCell(state, 1, 0.75, 0.9, 0.9, 0, 101, 9002, 7002, records);
  setCell(state, 2, 1.125, 0.8, 0.8, 1, 201, 9101, 7003, records);
  setCell(state, 3, 1.375, 0.7, 0.7, 1, 201, 9102, 7004, records);
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    setGasParticle(
        state,
        row,
        7001U + row,
        state.cells.center_x_comoving[row],
        state.cells.center_y_comoving[row],
        state.cells.center_z_comoving[row],
        state.cells.mass_code[row]);
  }
  state.gas_cell_identity.assign(std::move(records));
  state.rebuildSpeciesIndex();
  assert(state.validateOwnershipInvariants());
  assert(cosmosim::amr::hasProductionAmrHydroCoverage(state));
  return state;
}

[[nodiscard]] cosmosim::amr::FluxRegisterEntry makeFluxEntry(bool coarse, bool fine) {
  cosmosim::amr::FluxRegisterEntry entry;
  entry.register_key = 7701;
  entry.coarse_patch_id = 101;
  entry.coarse_gas_cell_id = 9002;
  entry.coarse_cell_index = 1;
  entry.level = 0;
  entry.axis = cosmosim::hydro::HydroFaceAxis::kX;
  entry.orientation = cosmosim::hydro::HydroFaceSide::kUpper;
  entry.face_area_comov = 0.5;
  entry.dt_code = coarse ? 1.0e-4 : 5.0e-5;
  if (coarse) {
    entry.coarse_area_comov = 0.5;
    entry.coarse_face_count = 1;
    entry.coarse_face_flux_code.mass_code = 0.10;
    entry.coarse_face_flux_code.momentum_x_code = 0.20;
    entry.coarse_face_flux_code.momentum_y_code = 0.01;
    entry.coarse_face_flux_code.momentum_z_code = -0.02;
    entry.coarse_face_flux_code.total_energy_code = 0.30;
  }
  if (fine) {
    entry.fine_area_comov = 0.5;
    entry.fine_face_count = 1;
    entry.fine_face_flux_code.mass_code = 0.40;
    entry.fine_face_flux_code.momentum_x_code = 0.50;
    entry.fine_face_flux_code.momentum_y_code = 0.03;
    entry.fine_face_flux_code.momentum_z_code = -0.04;
    entry.fine_face_flux_code.total_energy_code = 0.60;
  }
  return entry;
}

[[nodiscard]] cosmosim::amr::ProductionAmrHydroOptions pendingOptions(std::uint32_t fine_substep_index) {
  cosmosim::amr::ProductionAmrHydroOptions options;
  options.adiabatic_index = k_gamma;
  options.persist_incomplete_flux_registers = true;
  options.reflux_coarse_dt_code = 1.0e-4;
  options.reflux_interval_start_code = 0.0;
  options.reflux_interval_end_code = 1.0e-4;
  options.expected_fine_substeps = 2;
  options.fine_substep_index = fine_substep_index;
  return options;
}

[[nodiscard]] cosmosim::amr::AmrHydroGeometryOptions coarseTemporalOptions() {
  cosmosim::amr::AmrHydroGeometryOptions options;
  options.boundary_classes[1] = cosmosim::amr::AmrHydroBoundaryClass::kCoarseFine;
  return options;
}

[[nodiscard]] cosmosim::amr::AmrHydroGeometryOptions fineTemporalOptions() {
  cosmosim::amr::AmrHydroGeometryOptions options;
  options.boundary_classes[0] = cosmosim::amr::AmrHydroBoundaryClass::kCoarseFine;
  return options;
}

void consumeTemporalGhostAt(
    cosmosim::core::SimulationState& state,
    double fill_time_code,
    bool require_endpoint) {
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  auto coarse = cosmosim::amr::buildAmrHydroPatchGeometry(state, descriptors[0], coarseTemporalOptions());
  auto fine = cosmosim::amr::buildAmrHydroPatchGeometry(state, descriptors[1], fineTemporalOptions());
  auto coarse_conserved = cosmosim::amr::loadAmrHydroConservedState(state, coarse, k_gamma);
  auto fine_conserved = cosmosim::amr::loadAmrHydroConservedState(state, fine, k_gamma);
  std::vector<cosmosim::amr::AmrHydroGhostFillPatch> patches{
      {.geometry = &coarse,
       .conserved = &coarse_conserved,
       .target_state_time_code = 0.0,
       .ghost_fill_time_code = fill_time_code,
       .source_current_state_time_code = 1.0,
       .temporal_boundary_history = &state.amr_temporal_boundary_history,
       .enable_temporal_coarse_to_fine = true},
      {.geometry = &fine,
       .conserved = &fine_conserved,
       .target_state_time_code = fill_time_code,
       .ghost_fill_time_code = fill_time_code,
       .source_current_state_time_code = fill_time_code,
       .temporal_boundary_history = &state.amr_temporal_boundary_history,
       .enable_temporal_coarse_to_fine = true}};
  const auto diagnostics = cosmosim::amr::fillAmrHydroGhostCells(patches, k_gamma);
  assert(diagnostics.temporal_coarse_to_fine_ghosts_filled > 0U);
  if (require_endpoint) {
    assert(diagnostics.temporal_endpoint_ghosts_filled > 0U);
  }
}

void applyTemporalGhostRestartStep(cosmosim::tests::RestartEquivalenceStepContext& context) {
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(context.state);
  if (context.integrator_state.step_index == 0U) {
    cosmosim::amr::captureAmrTemporalBoundaryHistoryStart(context.state, descriptors, 0.0, k_gamma);
    // This is the coarse end state of the open interval.  The first fine substep
    // still consumes the recorded start state, not this future state.
    context.state.gas_cells.density_code[1] = 3.0;
    context.state.gas_cells.pressure_code[1] = 1.0;
    context.state.gas_cells.internal_energy_code[1] = 1.0 / ((k_gamma - 1.0) * 3.0);
    context.state.gas_cells.temperature_code[1] = context.state.gas_cells.internal_energy_code[1];
    context.state.gas_cells.sound_speed_code[1] = std::sqrt(k_gamma / 3.0);
    cosmosim::amr::captureAmrTemporalBoundaryHistoryEnd(context.state, descriptors, 1.0, k_gamma);
    consumeTemporalGhostAt(context.state, 0.0, true);

    std::vector<cosmosim::amr::FluxRegisterEntry> entries{makeFluxEntry(true, false), makeFluxEntry(false, true)};
    assert(cosmosim::amr::mergeFluxRegistersIntoPendingStore(context.state, entries, pendingOptions(0)) == 1U);
    assert(!context.state.amr_temporal_boundary_history.empty());
    assert(context.state.pending_flux_registers.size() == 1U);
    return;
  }
  if (context.integrator_state.step_index == 1U) {
    // This executes after HDF5 reload on the restart branch.  It proves the
    // deserialized temporal history is the actual ghost source at mid-interval.
    consumeTemporalGhostAt(context.state, 0.5, false);
    std::vector<cosmosim::amr::FluxRegisterEntry> entries{makeFluxEntry(false, true)};
    assert(cosmosim::amr::mergeFluxRegistersIntoPendingStore(context.state, entries, pendingOptions(1)) == 0U);
    const auto diagnostics = cosmosim::amr::applyCompletePendingFluxRegistersToSimulationState(
        context.state,
        descriptors,
        k_gamma);
    assert(diagnostics.complete_register_count == 1U);
    cosmosim::amr::retireAmrTemporalBoundaryHistory(context.state);
    assert(context.state.pending_flux_registers.empty());
    assert(context.state.amr_temporal_boundary_history.empty());
    return;
  }
  assert(false);
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_HDF5
  const auto restart_path = cosmosim::tests::stage8RestartPath("restart_equivalence_amr_temporal_ghosts");
  auto state = makeRestartPendingFluxState("restart_equivalence_amr_temporal_ghosts");
  auto scheduler = cosmosim::tests::makeStage8Scheduler(static_cast<std::uint32_t>(state.particles.size()), 2);
  auto integrator_state = cosmosim::tests::makeStage8IntegratorState(2, 2);
  auto output_state = cosmosim::tests::makeStage8OutputCadenceState(false);
  auto scenario = cosmosim::tests::makeStage8Scenario(
      std::move(state), integrator_state, std::move(scheduler), std::move(output_state), restart_path, 2, 1);
  scenario.step_kernel = applyTemporalGhostRestartStep;
  scenario.tolerances.position_abs = 1.0e-14;
  scenario.tolerances.velocity_abs = 1.0e-14;
  scenario.tolerances.scalar_abs = 1.0e-14;
  const auto result = cosmosim::tests::runRestartEquivalenceScenario(std::move(scenario));
  assert(result.direct_state.amr_temporal_boundary_history.empty());
  assert(result.restarted_state.amr_temporal_boundary_history.empty());
  assert(result.direct_state.pending_flux_registers.empty());
  assert(result.restarted_state.pending_flux_registers.empty());
  std::filesystem::remove(restart_path);
#endif
  return 0;
}
