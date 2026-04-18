# CosmoSim Validation Ladder

## Purpose

This document defines the routine validation ladder for gravity, hydro, AMR, and selected galaxy-formation modules.
The ladder is designed for desktop-first CI and developer loops, while preserving explicit tolerance policy and reference-data versioning.

## Validation classes and CTest mapping

- `validation_unit`: localized invariant checks (e.g., two-body tree gravity symmetry, tolerance-table integrity).
- `validation_integration`: cross-module checks through real solver paths (PM mode response, hydro mass conservation, cooling monotonicity).
- `validation_convergence`: resolution/algorithm-parameter trends (tree opening-angle convergence, hydro smooth-wave self-convergence).
- `validation_regression`: frozen behavior checks with explicit thresholds (AMR reflux correction, star-formation mass budget).
- `bench_validation_ladder`: profiling-only hook; not a correctness gate.

## Tolerance and regression storage policy

- Authoritative tolerance table: `validation/reference/validation_tolerances_v1.txt`.
- Keys are namespaced by validation case and metric.
- Changes to tolerance semantics require a new versioned file and changelog note.
- Existing versioned files are immutable in normal workflow.

## Numerical conventions and assumptions

- Gravity and hydro tests operate in **code units** and in the currently implemented frame conventions.
- Hydro tests use conservative variables and enforce mass conservation checks separately from accuracy checks.
- PM mode validation uses mode-shape correlation (cosine similarity) instead of absolute-amplitude equality to remain backend-agnostic.
- Convergence checks are deliberately lightweight and CI-sized; they are trend checks, not publication-grade convergence studies.

## Coverage ladder (current)

Implemented now:
- Gravity: two-body symmetry check, PM single-mode response, opening-angle convergence.
- TreePM split-kernel correctness: Gaussian SR/LR complementarity, mesh-cell (`asmth_cells`, `rcut_cells`) derivation audit,
  residual-cutoff pruning counters, and PM-only vs tree-only vs split consistency on small periodic setups.
- Hydro: Sod-like conservation check, smooth-wave self-convergence.
- AMR: coarse-fine reflux synchronization regression.
- Galaxy modules: cooling monotonicity check, star-formation mass-budget regression.

Planned as modules mature:
- Sedov blast convergence ladder.
- Evrard collapse regression and convergence set.
- MPI single-rank vs multi-rank consistency tests once MPI path is active in CI.

## Benchmark reporting expectations

`bench_validation_ladder` reports:
- build type,
- effective feature set,
- setup and steady-state timing separation for hydro,
- tree traversal counters,
- throughput proxy (`hydro_face_rate_mface_s`).

These measurements support profiling and regression triage only; they are not correctness evidence.
