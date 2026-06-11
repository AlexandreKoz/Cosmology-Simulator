#pragma once

#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::hydro {

// Add one ghost row per physical patch boundary face. Ghost rows are transient
// hydro scratch rows addressed after the real cell range and are not
// SimulationState cells.
void appendCartesianBoundaryGhostFaces(
    HydroPatchGeometry& geometry,
    HydroBoundaryKind boundary_kind);

// Fill physical-boundary ghost rows from authoritative real cells. Imported MPI
// ghosts are read-only and intentionally skipped by this helper.
void fillHydroBoundaryGhostCells(
    HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    double adiabatic_index);

}  // namespace cosmosim::hydro
