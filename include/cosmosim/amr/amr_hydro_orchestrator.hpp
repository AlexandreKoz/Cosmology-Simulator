#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/amr/amr_framework.hpp"
#include "cosmosim/amr/amr_ghost_fill.hpp"
#include "cosmosim/amr/amr_hydro_geometry.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace cosmosim::amr {

struct ProductionAmrHydroOptions {
  hydro::HydroBoundaryKind physical_boundary_kind = hydro::HydroBoundaryKind::kOpen;
  double adiabatic_index = 5.0 / 3.0;
  double density_floor = 1.0e-10;
  double pressure_floor = 1.0e-10;
};

struct ProductionAmrHydroDiagnostics {
  std::size_t patch_count = 0;
  std::size_t advanced_patch_count = 0;
  std::size_t active_cell_count = 0;
  std::size_t active_face_count = 0;
  std::size_t flux_register_entry_count = 0;
  AmrHydroGhostFillDiagnostics ghost_fill;
  RefluxDiagnostics reflux;
};

struct ProductionAmrRegridDiagnostics {
  std::size_t refined_patch_count = 0;
  std::size_t derefined_patch_count = 0;
  std::size_t created_gas_cell_count = 0;
  std::size_t retired_gas_cell_count = 0;
  double conserved_mass_before = 0.0;
  double conserved_mass_after = 0.0;
  double conserved_momentum_x_before = 0.0;
  double conserved_momentum_x_after = 0.0;
  double conserved_momentum_y_before = 0.0;
  double conserved_momentum_y_after = 0.0;
  double conserved_momentum_z_before = 0.0;
  double conserved_momentum_z_after = 0.0;
  double conserved_total_energy_before = 0.0;
  double conserved_total_energy_after = 0.0;
};

[[nodiscard]] bool hasProductionAmrHydroCoverage(const core::SimulationState& state);

[[nodiscard]] bool patchStateRowHasExplicitGeometry(
    const core::SimulationState& state,
    std::size_t patch_index);

[[nodiscard]] PatchDescriptor patchDescriptorFromStateRow(
    const core::SimulationState& state,
    std::size_t patch_index);

void writePatchDescriptorToStateRow(
    core::SimulationState& state,
    std::size_t patch_index,
    const PatchDescriptor& patch);

[[nodiscard]] std::vector<PatchDescriptor> buildProductionAmrPatchDescriptors(
    const core::SimulationState& state);

void populateAmrHydroFluxRegisterFaces(
    AmrHydroPatchGeometry& patch_geometry,
    std::span<const PatchDescriptor> all_patches);

void scatterAmrHydroConservedState(
    core::SimulationState& state,
    const AmrHydroPatchGeometry& patch_geometry,
    const hydro::HydroConservedStateSoa& conserved,
    double adiabatic_index);

[[nodiscard]] RefluxDiagnostics applyFluxRegistersToSimulationState(
    core::SimulationState& state,
    std::span<const FluxRegisterEntry> entries,
    std::span<const PatchDescriptor> all_patches,
    double adiabatic_index);

[[nodiscard]] ProductionAmrHydroDiagnostics advanceProductionAmrHydro(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options = {});

[[nodiscard]] ProductionAmrRegridDiagnostics refineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    std::uint64_t first_child_patch_id,
    std::uint64_t first_child_gas_cell_id,
    const ProductionAmrHydroOptions& options = {});

[[nodiscard]] ProductionAmrRegridDiagnostics refineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    const ProductionAmrHydroOptions& options = {});

[[nodiscard]] ProductionAmrRegridDiagnostics derefineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    std::uint64_t replacement_gas_cell_id,
    const ProductionAmrHydroOptions& options = {});

[[nodiscard]] ProductionAmrRegridDiagnostics derefineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    const ProductionAmrHydroOptions& options = {});

}  // namespace cosmosim::amr
