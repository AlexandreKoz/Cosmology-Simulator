#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace {

constexpr double k_gamma = 1.4;
constexpr double k_conservation_tol = 5.0e-4;

struct Totals {
  double mass = 0.0;
  double momentum_x = 0.0;
  double momentum_y = 0.0;
  double momentum_z = 0.0;
  double total_energy = 0.0;
};

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

[[nodiscard]] cosmosim::hydro::HydroPrimitiveState primitiveForRow(
    const cosmosim::core::SimulationState& state,
    std::uint32_t row) {
  return cosmosim::hydro::HydroPrimitiveState{
      .rho_comoving = state.gas_cells.density_code[row],
      .vel_x_peculiar = state.gas_cells.velocity_x_peculiar[row],
      .vel_y_peculiar = state.gas_cells.velocity_y_peculiar[row],
      .vel_z_peculiar = state.gas_cells.velocity_z_peculiar[row],
      .pressure_comoving = state.gas_cells.pressure_code[row]};
}

[[nodiscard]] Totals totalState(
    const cosmosim::core::SimulationState& state,
    const std::vector<cosmosim::amr::PatchDescriptor>& descriptors) {
  Totals totals;
  for (const cosmosim::amr::PatchDescriptor& patch : descriptors) {
    const double volume = patch.extent_comov[0] * patch.extent_comov[1] * patch.extent_comov[2] /
        static_cast<double>(static_cast<std::size_t>(patch.cell_dims[0]) * patch.cell_dims[1] * patch.cell_dims[2]);
    for (const std::uint32_t row : state.gas_cell_identity.rowsForPatch(patch.patch_id)) {
      const auto conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(
          primitiveForRow(state, row), k_gamma);
      totals.mass += conserved.mass_density_comoving * volume;
      totals.momentum_x += conserved.momentum_density_x_comoving * volume;
      totals.momentum_y += conserved.momentum_density_y_comoving * volume;
      totals.momentum_z += conserved.momentum_density_z_comoving * volume;
      totals.total_energy += conserved.total_energy_density_comoving * volume;
    }
  }
  return totals;
}

void setCell(
    cosmosim::core::SimulationState& state,
    std::uint32_t row,
    double x,
    double rho,
    double velocity_x,
    double velocity_y,
    double velocity_z,
    double pressure,
    std::uint32_t patch_index,
    std::uint64_t patch_id,
    std::uint64_t gas_cell_id,
    std::vector<cosmosim::core::GasCellIdentityRecord>& records) {
  state.cells.center_x_comoving[row] = x;
  state.cells.center_y_comoving[row] = 0.5;
  state.cells.center_z_comoving[row] = 0.5;
  state.cells.patch_index[row] = patch_index;
  state.cells.time_bin[row] = 0;
  state.cells.mass_code[row] = rho;
  state.gas_cells.gas_cell_id[row] = gas_cell_id;
  state.gas_cells.parent_particle_id[row] = 0;
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
      .parent_particle_id = std::nullopt,
      .owning_patch_id = patch_id,
      .local_cell_row = row});
}

[[nodiscard]] cosmosim::core::SimulationState makeCoarseFineState() {
  cosmosim::core::SimulationState state;
  state.resizeCells(4);
  state.resizePatches(2);
  setPatch(state, 0, cosmosim::amr::PatchDescriptor{
      .patch_id = 101,
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
  setCell(state, 0, 0.25, 1.0, 0.20, 0.04, -0.03, 1.0, 0, 101, 9001, records);
  setCell(state, 1, 0.75, 0.9, 0.35, 0.08, -0.02, 0.9, 0, 101, 9002, records);
  setCell(state, 2, 1.125, 0.8, -0.25, -0.05, 0.06, 0.8, 1, 201, 9101, records);
  setCell(state, 3, 1.375, 0.7, -0.10, -0.02, 0.03, 0.7, 1, 201, 9102, records);
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

void assertNear(double lhs, double rhs) {
  assert(std::abs(lhs - rhs) < k_conservation_tol);
}


cosmosim::amr::FluxRegisterEntry makeCompleteRegister(
    std::uint64_t coarse_gas_cell_id = 9002,
    std::size_t coarse_cell_index = 1) {
  return cosmosim::amr::FluxRegisterEntry{
      .register_key = 199,
      .coarse_patch_id = 101,
      .coarse_gas_cell_id = coarse_gas_cell_id,
      .coarse_cell_index = coarse_cell_index,
      .level = 0,
      .coarse_face_flux_code = cosmosim::amr::ConservedState{
          .mass_code = 0.10,
          .momentum_x_code = 0.20,
          .momentum_y_code = -0.10,
          .momentum_z_code = 0.05,
          .total_energy_code = 0.30},
      .fine_face_flux_code = cosmosim::amr::ConservedState{
          .mass_code = 0.40,
          .momentum_x_code = 0.70,
          .momentum_y_code = 0.20,
          .momentum_z_code = -0.15,
          .total_energy_code = 0.90},
      .face_area_comov = 0.5,
      .coarse_area_comov = 0.5,
      .fine_area_comov = 0.5,
      .dt_code = 1.0e-5,
      .coarse_face_count = 1,
      .fine_face_count = 1};
}

void rebuildIdentityFromCurrentRows(cosmosim::core::SimulationState& state) {
  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(state.cells.size());
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    records.push_back(cosmosim::core::GasCellIdentityRecord{
        .gas_cell_id = state.gas_cells.gas_cell_id[row],
        .parent_particle_id = std::nullopt,
        .owning_patch_id = state.patches.patch_id[state.cells.patch_index[row]],
        .local_cell_row = row});
  }
  state.gas_cell_identity.assign(std::move(records));
}

void swapGasRows(cosmosim::core::SimulationState& state, std::uint32_t a, std::uint32_t b) {
  using std::swap;
  swap(state.cells.center_x_comoving[a], state.cells.center_x_comoving[b]);
  swap(state.cells.center_y_comoving[a], state.cells.center_y_comoving[b]);
  swap(state.cells.center_z_comoving[a], state.cells.center_z_comoving[b]);
  swap(state.cells.mass_code[a], state.cells.mass_code[b]);
  swap(state.cells.time_bin[a], state.cells.time_bin[b]);
  swap(state.cells.patch_index[a], state.cells.patch_index[b]);
  swap(state.gas_cells.gas_cell_id[a], state.gas_cells.gas_cell_id[b]);
  swap(state.gas_cells.parent_particle_id[a], state.gas_cells.parent_particle_id[b]);
  swap(state.gas_cells.velocity_x_peculiar[a], state.gas_cells.velocity_x_peculiar[b]);
  swap(state.gas_cells.velocity_y_peculiar[a], state.gas_cells.velocity_y_peculiar[b]);
  swap(state.gas_cells.velocity_z_peculiar[a], state.gas_cells.velocity_z_peculiar[b]);
  swap(state.gas_cells.density_code[a], state.gas_cells.density_code[b]);
  swap(state.gas_cells.pressure_code[a], state.gas_cells.pressure_code[b]);
  swap(state.gas_cells.internal_energy_code[a], state.gas_cells.internal_energy_code[b]);
  swap(state.gas_cells.temperature_code[a], state.gas_cells.temperature_code[b]);
  swap(state.gas_cells.sound_speed_code[a], state.gas_cells.sound_speed_code[b]);
  rebuildIdentityFromCurrentRows(state);
}

void testAutomaticRefluxChangesCoarseOwnerAndPreservesTotals() {
  cosmosim::core::SimulationState state = makeCoarseFineState();
  const auto before_descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  const Totals before = totalState(state, before_descriptors);
  const double coarse_density_before = state.gas_cells.density_code[1];

  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::HlleRiemannSolver riemann;
  const cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-5, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
  const std::vector<std::uint32_t> active_rows{0, 1, 2, 3};
  const cosmosim::amr::ProductionAmrHydroOptions options{
      .physical_boundary_kind = cosmosim::hydro::HydroBoundaryKind::kOpen,
      .adiabatic_index = k_gamma,
      .density_floor = 1.0e-12,
      .pressure_floor = 1.0e-12};

  const auto diagnostics = cosmosim::amr::advanceProductionAmrHydro(
      state, active_rows, update, source_context, solver, riemann, {}, options);

  assert(diagnostics.flux_register_entry_count > 0U);
  assert(diagnostics.reflux.complete_register_count > 0U);
  assert(diagnostics.reflux.corrected_cells > 0U);
  assert(diagnostics.reflux.corrected_mass_code > 0.0);
  assert(diagnostics.reflux.corrected_momentum_x_code > 0.0);
  assert(diagnostics.reflux.corrected_momentum_y_code > 0.0);
  assert(diagnostics.reflux.corrected_momentum_z_code > 0.0);
  assert(diagnostics.reflux.corrected_total_energy_code > 0.0);
  assert(diagnostics.reflux.corrected_internal_energy_code >= 0.0);
  assert(std::abs(state.gas_cells.density_code[1] - coarse_density_before) > 0.0);

  const Totals after = totalState(state, cosmosim::amr::buildProductionAmrPatchDescriptors(state));
  assertNear(after.mass, before.mass);
  assertNear(after.momentum_x, before.momentum_x);
  assertNear(after.momentum_y, before.momentum_y);
  assertNear(after.momentum_z, before.momentum_z);
  assertNear(after.total_energy, before.total_energy);
}


void testManualRefluxAcceptanceSkipsUnsafeRegisters() {
  cosmosim::core::SimulationState state = makeCoarseFineState();
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);

  auto incomplete_coarse_only = makeCompleteRegister();
  incomplete_coarse_only.fine_face_count = 0;
  auto incomplete_fine_only = makeCompleteRegister();
  incomplete_fine_only.coarse_face_count = 0;
  auto area_mismatch = makeCompleteRegister();
  area_mismatch.fine_area_comov = 0.25;
  auto missing_id = makeCompleteRegister(999999, 1);
  auto wrong_patch_owner = makeCompleteRegister(9101, 0);
  std::vector<cosmosim::amr::FluxRegisterEntry> entries{
      incomplete_coarse_only,
      incomplete_fine_only,
      area_mismatch,
      missing_id,
      wrong_patch_owner};

  const double density_before = state.gas_cells.density_code[1];
  const auto diagnostics = cosmosim::amr::applyFluxRegistersToSimulationState(
      state, entries, descriptors, k_gamma);
  assert(diagnostics.corrected_cells == 0U);
  assert(diagnostics.skipped_incomplete_register_count == 2U);
  assert(diagnostics.skipped_area_mismatch_count == 1U);
  assert(diagnostics.skipped_missing_target_count == 2U);
  assert(state.gas_cells.density_code[1] == density_before);
}

void testCompleteManualRefluxAppliesOnceByStableGasCellId() {
  cosmosim::core::SimulationState state = makeCoarseFineState();
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  const double density_before = state.gas_cells.density_code[1];
  std::vector<cosmosim::amr::FluxRegisterEntry> entries{makeCompleteRegister()};
  const auto diagnostics = cosmosim::amr::applyFluxRegistersToSimulationState(
      state, entries, descriptors, k_gamma);
  assert(diagnostics.complete_register_count == 1U);
  assert(diagnostics.corrected_cells == 1U);
  assert(diagnostics.corrected_mass_code > 0.0);
  assert(diagnostics.corrected_momentum_x_code > 0.0);
  assert(diagnostics.corrected_momentum_y_code > 0.0);
  assert(diagnostics.corrected_momentum_z_code > 0.0);
  assert(diagnostics.corrected_total_energy_code > 0.0);
  assert(state.gas_cells.density_code[1] != density_before);
}

void testRefluxRowReorderTargetsStableGasCellId() {
  cosmosim::core::SimulationState state = makeCoarseFineState();
  swapGasRows(state, 0, 1);
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  const double target_density_before = state.gas_cells.density_code[0];
  const double other_density_before = state.gas_cells.density_code[1];
  std::vector<cosmosim::amr::FluxRegisterEntry> entries{makeCompleteRegister(9002, 1)};

  const auto diagnostics = cosmosim::amr::applyFluxRegistersToSimulationState(
      state, entries, descriptors, k_gamma);
  assert(diagnostics.corrected_cells == 1U);
  assert(state.gas_cells.gas_cell_id[0] == 9002U);
  assert(state.gas_cells.density_code[0] != target_density_before);
  assert(state.gas_cells.density_code[1] == other_density_before);
}

void testStaleRefluxTargetMappingFails() {
  cosmosim::core::SimulationState state = makeCoarseFineState();
  const auto descriptors = cosmosim::amr::buildProductionAmrPatchDescriptors(state);
  std::vector<cosmosim::amr::FluxRegisterEntry> entries{cosmosim::amr::FluxRegisterEntry{
      .register_key = 99,
      .coarse_patch_id = 101,
      .coarse_gas_cell_id = 9002,
      .coarse_cell_index = 0,
      .level = 0,
      .coarse_face_flux_code = cosmosim::amr::ConservedState{.mass_code = 1.0},
      .fine_face_flux_code = cosmosim::amr::ConservedState{.mass_code = 2.0},
      .face_area_comov = 0.5,
      .coarse_area_comov = 0.5,
      .fine_area_comov = 0.5,
      .dt_code = 1.0e-5,
      .coarse_face_count = 1,
      .fine_face_count = 1}};

  bool rejected = false;
  try {
    (void)cosmosim::amr::applyFluxRegistersToSimulationState(state, entries, descriptors, k_gamma);
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  assert(rejected);
}

}  // namespace

int main() {
  testAutomaticRefluxChangesCoarseOwnerAndPreservesTotals();
  testManualRefluxAcceptanceSkipsUnsafeRegisters();
  testCompleteManualRefluxAppliesOnceByStableGasCellId();
  testRefluxRowReorderTargetsStableGasCellId();
  testStaleRefluxTargetMappingFails();
  return 0;
}
