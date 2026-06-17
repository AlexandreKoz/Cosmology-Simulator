#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace cosmosim::tests::amr_hydro_validation {

constexpr double k_gamma = 1.4;

struct CellCenter {
  double x_comoving = 0.0;
  double y_comoving = 0.0;
  double z_comoving = 0.0;
};

struct Totals {
  double mass = 0.0;
  double momentum_x = 0.0;
  double momentum_y = 0.0;
  double momentum_z = 0.0;
  double total_energy = 0.0;
};

using CellInitializer = std::function<hydro::HydroPrimitiveState(
    const amr::PatchDescriptor& patch,
    std::size_t patch_cell,
    const CellCenter& center)>;

inline void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

[[nodiscard]] inline std::size_t patchCellCount(const amr::PatchDescriptor& patch) {
  return static_cast<std::size_t>(patch.cell_dims[0]) *
      static_cast<std::size_t>(patch.cell_dims[1]) *
      static_cast<std::size_t>(patch.cell_dims[2]);
}

[[nodiscard]] inline double patchCellVolume(const amr::PatchDescriptor& patch) {
  return patch.extent_comov[0] * patch.extent_comov[1] * patch.extent_comov[2] /
      static_cast<double>(patchCellCount(patch));
}

[[nodiscard]] inline CellCenter cellCenter(const amr::PatchDescriptor& patch, std::size_t patch_cell) {
  const std::size_t nx = patch.cell_dims[0];
  const std::size_t ny = patch.cell_dims[1];
  const std::size_t i = patch_cell % nx;
  const std::size_t j = (patch_cell / nx) % ny;
  const std::size_t k = patch_cell / (nx * ny);
  return CellCenter{
      .x_comoving = patch.origin_comov[0] +
          (static_cast<double>(i) + 0.5) * patch.extent_comov[0] / static_cast<double>(patch.cell_dims[0]),
      .y_comoving = patch.origin_comov[1] +
          (static_cast<double>(j) + 0.5) * patch.extent_comov[1] / static_cast<double>(patch.cell_dims[1]),
      .z_comoving = patch.origin_comov[2] +
          (static_cast<double>(k) + 0.5) * patch.extent_comov[2] / static_cast<double>(patch.cell_dims[2])};
}

inline void writePatch(
    core::SimulationState& state,
    std::size_t patch_index,
    const amr::PatchDescriptor& patch,
    std::uint32_t first_cell,
    std::uint32_t cell_count) {
  amr::writePatchDescriptorToStateRow(state, patch_index, patch);
  state.patches.first_cell[patch_index] = first_cell;
  state.patches.cell_count[patch_index] = cell_count;
  state.patches.owning_rank[patch_index] = 0;
}

inline void writeCell(
    core::SimulationState& state,
    std::uint32_t row,
    std::uint32_t patch_index,
    const amr::PatchDescriptor& patch,
    std::size_t patch_cell,
    std::uint64_t gas_cell_id,
    const hydro::HydroPrimitiveState& primitive,
    std::vector<core::GasCellIdentityRecord>& records) {
  const CellCenter center = cellCenter(patch, patch_cell);
  const double volume = patchCellVolume(patch);
  state.cells.center_x_comoving[row] = center.x_comoving;
  state.cells.center_y_comoving[row] = center.y_comoving;
  state.cells.center_z_comoving[row] = center.z_comoving;
  state.cells.patch_index[row] = patch_index;
  state.cells.time_bin[row] = 0;
  state.cells.mass_code[row] = primitive.rho_comoving * volume;
  state.gas_cells.gas_cell_id[row] = gas_cell_id;
  state.gas_cells.parent_particle_id[row] = 0;
  state.gas_cells.density_code[row] = primitive.rho_comoving;
  state.gas_cells.pressure_code[row] = primitive.pressure_comoving;
  state.gas_cells.internal_energy_code[row] =
      primitive.pressure_comoving / ((k_gamma - 1.0) * primitive.rho_comoving);
  state.gas_cells.velocity_x_peculiar[row] = primitive.vel_x_peculiar;
  state.gas_cells.velocity_y_peculiar[row] = primitive.vel_y_peculiar;
  state.gas_cells.velocity_z_peculiar[row] = primitive.vel_z_peculiar;
  state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
  state.gas_cells.sound_speed_code[row] =
      std::sqrt(k_gamma * primitive.pressure_comoving / primitive.rho_comoving);
  records.push_back(core::GasCellIdentityRecord{
      .gas_cell_id = gas_cell_id,
      .parent_particle_id = std::nullopt,
      .owning_patch_id = patch.patch_id,
      .local_cell_row = row});
}

[[nodiscard]] inline core::SimulationState makePatchState(
    const std::vector<amr::PatchDescriptor>& patches,
    const CellInitializer& initializer,
    std::uint64_t first_gas_cell_id = 10000U) {
  std::size_t cell_count = 0;
  for (const auto& patch : patches) {
    cell_count += patchCellCount(patch);
  }

  core::SimulationState state;
  state.resizeCells(cell_count);
  state.resizePatches(patches.size());

  std::vector<core::GasCellIdentityRecord> records;
  records.reserve(cell_count);
  std::uint32_t row = 0;
  std::uint64_t gas_cell_id = first_gas_cell_id;
  for (std::size_t patch_index = 0; patch_index < patches.size(); ++patch_index) {
    const amr::PatchDescriptor& patch = patches[patch_index];
    const std::uint32_t first_cell = row;
    const std::uint32_t count = static_cast<std::uint32_t>(patchCellCount(patch));
    writePatch(state, patch_index, patch, first_cell, count);
    for (std::size_t patch_cell = 0; patch_cell < count; ++patch_cell) {
      writeCell(
          state,
          row,
          static_cast<std::uint32_t>(patch_index),
          patch,
          patch_cell,
          gas_cell_id++,
          initializer(patch, patch_cell, cellCenter(patch, patch_cell)),
          records);
      ++row;
    }
  }
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

[[nodiscard]] inline hydro::HydroPrimitiveState primitiveForRow(
    const core::SimulationState& state,
    std::uint32_t row) {
  return hydro::HydroPrimitiveState{
      .rho_comoving = state.gas_cells.density_code[row],
      .vel_x_peculiar = state.gas_cells.velocity_x_peculiar[row],
      .vel_y_peculiar = state.gas_cells.velocity_y_peculiar[row],
      .vel_z_peculiar = state.gas_cells.velocity_z_peculiar[row],
      .pressure_comoving = state.gas_cells.pressure_code[row]};
}

[[nodiscard]] inline Totals totalState(
    const core::SimulationState& state,
    const std::vector<amr::PatchDescriptor>& descriptors) {
  Totals totals;
  for (const amr::PatchDescriptor& patch : descriptors) {
    const double volume = patchCellVolume(patch);
    for (const std::uint32_t row : state.gas_cell_identity.rowsForPatch(patch.patch_id)) {
      const auto conserved = hydro::HydroCoreSolver::conservedFromPrimitive(primitiveForRow(state, row), k_gamma);
      totals.mass += conserved.mass_density_comoving * volume;
      totals.momentum_x += conserved.momentum_density_x_comoving * volume;
      totals.momentum_y += conserved.momentum_density_y_comoving * volume;
      totals.momentum_z += conserved.momentum_density_z_comoving * volume;
      totals.total_energy += conserved.total_energy_density_comoving * volume;
    }
  }
  return totals;
}

inline void requireFinitePositiveState(const core::SimulationState& state, const std::string& case_name) {
  requireOrThrow(state.gas_cell_identity.coversDenseLocalRows(state.cells.size()), case_name + ": sparse identity map");
  requireOrThrow(state.gasCellIdentityMapMatchesSidecarLanes(), case_name + ": stale gas-cell identity sidecars");
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const hydro::HydroPrimitiveState primitive = primitiveForRow(state, row);
    requireOrThrow(state.gas_cells.gas_cell_id[row] != 0U, case_name + ": zero gas_cell_id");
    requireOrThrow(std::isfinite(primitive.rho_comoving), case_name + ": non-finite density");
    requireOrThrow(std::isfinite(primitive.pressure_comoving), case_name + ": non-finite pressure");
    requireOrThrow(std::isfinite(primitive.vel_x_peculiar), case_name + ": non-finite velocity x");
    requireOrThrow(std::isfinite(primitive.vel_y_peculiar), case_name + ": non-finite velocity y");
    requireOrThrow(std::isfinite(primitive.vel_z_peculiar), case_name + ": non-finite velocity z");
    requireOrThrow(std::isfinite(state.gas_cells.internal_energy_code[row]), case_name + ": non-finite internal energy");
    requireOrThrow(primitive.rho_comoving > 0.0, case_name + ": non-positive density");
    requireOrThrow(primitive.pressure_comoving > 0.0, case_name + ": non-positive pressure");
    requireOrThrow(state.gas_cells.internal_energy_code[row] > 0.0, case_name + ": non-positive internal energy");
  }
}

[[nodiscard]] inline std::vector<std::uint32_t> allCellRows(const core::SimulationState& state) {
  std::vector<std::uint32_t> rows(state.cells.size());
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    rows[row] = row;
  }
  return rows;
}

[[nodiscard]] inline amr::ProductionAmrHydroDiagnostics advanceProductionHydroSteps(
    core::SimulationState& state,
    std::size_t step_count,
    double dt_code,
    hydro::HydroBoundaryKind boundary_kind = hydro::HydroBoundaryKind::kReflective) {
  hydro::HydroCoreSolver solver(k_gamma);
  hydro::HllcRiemannSolver riemann;
  amr::ProductionAmrHydroDiagnostics last;
  const amr::ProductionAmrHydroOptions options{
      .physical_boundary_kind = boundary_kind,
      .adiabatic_index = k_gamma,
      .density_floor = 1.0e-10,
      .pressure_floor = 1.0e-10};
  for (std::size_t step = 0; step < step_count; ++step) {
    const hydro::HydroUpdateContext update{
        .dt_code = dt_code,
        .scale_factor = 1.0,
        .hubble_rate_code = 0.0};
    const hydro::HydroSourceContext source_context{.update = update};
    const std::vector<std::uint32_t> active_rows = allCellRows(state);
    last = amr::advanceProductionAmrHydro(
        state,
        active_rows,
        update,
        source_context,
        solver,
        riemann,
        {},
        options);
  }
  return last;
}

[[nodiscard]] inline double relativeDifference(double lhs, double rhs) {
  return std::abs(lhs - rhs) / std::max({1.0e-14, std::abs(lhs), std::abs(rhs)});
}

inline void requireConservedClose(
    const Totals& actual,
    const Totals& expected,
    double tolerance,
    const std::string& case_name) {
  requireOrThrow(relativeDifference(actual.mass, expected.mass) < tolerance, case_name + ": mass drift");
  requireOrThrow(
      relativeDifference(actual.momentum_x, expected.momentum_x) < tolerance,
      case_name + ": momentum-x drift");
  requireOrThrow(
      relativeDifference(actual.momentum_y, expected.momentum_y) < tolerance,
      case_name + ": momentum-y drift");
  requireOrThrow(
      relativeDifference(actual.momentum_z, expected.momentum_z) < tolerance,
      case_name + ": momentum-z drift");
  requireOrThrow(
      relativeDifference(actual.total_energy, expected.total_energy) < tolerance,
      case_name + ": total-energy drift");
}

}  // namespace cosmosim::tests::amr_hydro_validation
