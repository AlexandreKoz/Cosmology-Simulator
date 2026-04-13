"""Notebook-friendly wrappers over the stable C++ binding surface."""

from __future__ import annotations

from ._cosmosim import DiagnosticsEngine, SnapshotReadResult


def summarize_run_health(snapshot: SnapshotReadResult, config) -> dict:
    """Return conservative run-health counters for a loaded snapshot."""
    diagnostics = DiagnosticsEngine(config)
    counters = diagnostics.compute_run_health(snapshot.state)
    return {
        "particle_count": int(counters.particle_count),
        "cell_count": int(counters.cell_count),
        "star_count": int(counters.star_count),
        "ownership_invariants_ok": bool(counters.ownership_invariants_ok),
        "unique_particle_ids_ok": bool(counters.unique_particle_ids_ok),
        "non_finite_particles": int(counters.non_finite_particles),
        "non_finite_cells": int(counters.non_finite_cells),
    }


def snapshot_mass_sum(snapshot: SnapshotReadResult) -> float:
    """Return total particle mass in code units."""
    return float(snapshot.state.mass_code().sum())
