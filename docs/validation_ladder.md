# CosmoSim Validation Ladder

## Purpose

This document defines the routine validation ladder for gravity, hydro, AMR, and selected galaxy-formation modules.
The ladder is designed for desktop-first CI and developer loops, while preserving explicit tolerance policy and reference-data versioning.

## Validation classes and CTest mapping

- `validation_unit`: localized invariant checks (e.g., two-body tree gravity symmetry, small-`N` direct-sum agreement in open boundaries, split-kernel composition identity).
- `validation_integration`: cross-module checks through real solver paths (PM analytic mode response, PM uniform-density cancellation, TreePM periodic coupling consistency, hydro mass conservation, cooling monotonicity).
- `validation_convergence`: resolution/algorithm-parameter trends and controlled dynamics checks (tree opening-angle convergence, hydro smooth-wave self-convergence, two-body orbit energy-drift monitor, static-halo radial-force profile trend).
- `validation_regression`: frozen behavior checks with explicit thresholds (AMR reflux correction, star-formation mass budget). Regression gates are stability sentinels and are not standalone scientific proof.
- `bench_validation_ladder`: profiling-only hook; not a correctness gate.

## Tolerance and regression storage policy

- Authoritative tolerance table: `validation/reference/validation_tolerances_v1.txt`.
- Keys are namespaced by validation case and metric.
- Changes to tolerance semantics require a new versioned file and changelog note.
- Existing versioned files are immutable in normal workflow.

## Numerical conventions and assumptions

- Gravity and hydro tests operate in **code units** and in the currently implemented frame conventions.
- Hydro tests use conservative variables and enforce mass conservation checks separately from accuracy checks.
- PM analytic mode validation compares force and potential to the periodic spectral solution for the selected plane-wave mode.
- Controlled two-body orbit checks use the same softened potential law as the tree solver (`epsilon_comoving`-regularized force and potential).
- Static-halo profile checks compare tree force against softened direct-summation reference on fixed shells.
- Periodic TreePM **force-error mapping** uses a **fine spectral PM + exact pairwise short-range residual periodic proxy reference**; this is explicitly not Ewald-exact.
- Periodic TreePM **integration consistency** (`test_tree_pm_coupling_periodic`) uses a **minimum-image periodic direct reference**; this is also not Ewald-exact.

## Coverage ladder (Phase 1 gravity evidence now wired)

Implemented now:
- Gravity:
  - two-body symmetry check;
  - small-`N` open-boundary tree-vs-direct summation check;
  - PM single-mode (analytic plane-wave) response;
  - PM uniform-density cancellation on both exact lattice and deterministic parity-balanced glass-like proxy placements;
  - PM-only vs tree-only vs split TreePM consistency check;
  - split-kernel composition identity check;
  - opening-angle convergence check;
  - two-body orbit energy-drift check;
  - static-halo radial-force profile check vs softened direct-sum reference.
- Hydro: Sod-like conservation check, smooth-wave self-convergence.
- AMR: coarse-fine reflux synchronization regression.
- Galaxy modules: cooling monotonicity check, star-formation mass-budget regression.

Phase 2 distributed TreePM gravity gate (implemented):
- `validation_phase2_mpi_gravity_single_rank`: one-rank reference baseline through the same distributed-code path harness.
- `validation_phase2_mpi_gravity_two_rank`: distributed PM and distributed full TreePM vs one-rank reference, with explicit thresholds:
  - PM force `rel_L2 <= 1e-10`.
  - TreePM force `rel_L2 <= 5e-6`, `max_rel <= 5e-5`.
- `validation_phase2_mpi_gravity_three_rank`: rank-count reproducibility sweep using the same invariants at `np=3`.
- Communication stress is included in the same executable using `tree_exchange_batch_bytes=64` and PM-cadence toggling to force repeated export/import rounds.
- Restart continuation gate is included through reference-workflow restart write/read verification in MPI mode (`restart_roundtrip_ok`), exact distributed ownership/slab metadata replay, and rank-qualified restart file naming.
- Distributed-workflow honesty is checked separately in `test_reference_workflow_distributed_treepm_mpi` by reducing local owned counts and particle-ID checksums; identical per-rank full-state digests are not treated as evidence of a real distributed-memory workflow.

Planned as modules mature:
- Sedov blast convergence ladder.
- Evrard collapse regression and convergence set.

## Force-error map artifact (PMGRID/ASMTH/RCUT)

- Executable: `bench_tree_pm_force_error_map`.
- Output artifact (deterministic): `validation/artifacts/tree_pm_force_error_map.csv`.
- Swept coordinates:
  - `pm_grid ∈ {16,24,32}`
  - `asmth_cells ∈ {0.8,1.25,2.0}`
  - `rcut_cells ∈ {3.0,4.5,6.0}`
- CSV columns:
  - `reference_method` (currently `periodic_spectral_direct_proxy`)
  - `pm_grid`, `asmth_cells`, `rcut_cells`
  - `relative_l2_error`
- Interpretation:
  - lower `relative_l2_error` indicates closer agreement to the fine spectral-plus-direct periodic proxy;
  - this map is a parameter-sensitivity evidence class for Phase 1 and must be read with reference-method limits.

## Benchmark reporting expectations

`bench_validation_ladder` reports:
- build type,
- effective feature set,
- setup and steady-state timing separation for hydro,
- tree traversal counters,
- throughput proxy (`hydro_face_rate_mface_s`).

These measurements support profiling and regression triage only; they are not correctness evidence.


## Phase 3 status (contract frozen, implementation pending)

- Phase 3 closure criteria are now frozen in `docs/treepm_phase3_contract.md`.
- Current ladder evidence remains Phase 1/2-focused and does **not** by itself close Phase 3.
- Existing `pm_only_scaling_*` and `tree_only_scaling_*` artifacts are performance artifacts and must not be treated as strong/weak scaling certification.

## Phase 2 scaling artifacts (separate PM-only and tree-only)

- Target: `generate_mpi_gravity_scaling_artifacts` (MPI-enabled builds).
- Produces deterministic CSV artifacts in `validation/artifacts/`:
  - `pm_only_scaling_np1.csv`, `pm_only_scaling_np2.csv`
  - `tree_only_scaling_np1.csv`, `tree_only_scaling_np2.csv`
- These are performance evidence artifacts, not correctness gates.
