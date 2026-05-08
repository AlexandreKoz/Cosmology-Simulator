# CI, regression, and reproducibility pipeline

This document defines the CI contract for CosmoSim and how it maps to reproducibility and regression discipline.

## CI workflow topology

Workflow file: `.github/workflows/ci.yml`.

### Presubmit / push matrix

`build_test_matrix` runs three representative presets:

1. `cpu-only-debug` (full test preset) for baseline architecture health.
2. `hdf5-debug` (targeted snapshot/restart/provenance tests) for schema and I/O path checks.
3. `pm-hdf5-fftw-debug` (targeted PM + validation tests) for feature-rich gravity path checks.

Each matrix leg executes `scripts/ci/run_preset_pipeline.sh` to perform:
- configure,
- build,
- test,
- feature summary / metadata archive,
- optional lightweight benchmark sentinel (`bench_layout_smoke`).

### Reproducibility gate

`reproducibility_gate` runs `scripts/ci/run_reproducibility_gate.sh`, which enforces:
- deterministic config/provenance-oriented checks (`unit_config_parser`, `integration_provenance_roundtrip`, `integration_feature_summary`),
- feature summary contract check,
- validation regression check,
- SHA-256 recording for versioned tolerance and benchmark size baselines.

### Optional MPI path

`optional_mpi_smoke` is restricted to `schedule` and manual dispatch. It exercises
`mpi-release` build + selected distributed-memory tests without inflating presubmit cost.

## Regression artifact contract

Baseline expectations are versioned in:

- `validation/reference/ci_build_metadata_expectations_v1.json`

`scripts/ci/check_build_metadata.py` validates generated `cosmosim_build_metadata.json` against
that contract using conservative checks:

- preset identity,
- enabled feature booleans,
- dependency state (when declared in baseline),
- required core targets.

Any intentional baseline drift must be committed by updating this contract with rationale in the same PR.

## Local usage

```bash
./scripts/ci/run_preset_pipeline.sh cpu-only-debug build-cpu-debug test-cpu-debug - ci_artifacts/local_cpu 1
./scripts/ci/check_build_metadata.py \
  ci_artifacts/local_cpu/cosmosim_build_metadata-cpu-only-debug.json \
  validation/reference/ci_build_metadata_expectations_v1.json \
  cpu-only-debug
./scripts/ci/run_reproducibility_gate.sh ci_artifacts/local_repro
```

## Assumptions and limitations

- MPI and GPU paths are optional in default hosted CI to keep turnaround affordable.
- Regression checks rely on curated scientific tests plus metadata contract checks; this is not a substitute for full validation campaigns.
- Benchmark outputs are archived as profiling sentinels only and must not be interpreted as correctness proof.

## Stage 1 runtime-truth CPU P0 gate

The dependency-free Stage 1 runtime-truth gate is the local and CI entry point for runtime-truth repair closure. It is intentionally CPU-only: it must configure, build, and run without HDF5, FFTW, MPI, CUDA, or Python being enabled.

One-command local gate:

```bash
./scripts/ci/run_stage1_runtime_truth_gate.sh ci_artifacts/stage1_runtime_truth
```

Equivalent expanded commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4
ctest --preset test-stage1-runtime-truth-cpu-debug --output-on-failure
```

Compatibility preset retained for earlier repair notes:

```bash
ctest --preset test-stage0-runtime-truth-cpu-debug --output-on-failure
```

The Stage 1 P0 suite is registered explicitly rather than selected by broad substrings. The suite covers:

- active views and scheduler mirrors: `unit_simulation_state`, `integration_hierarchical_time_bins`, `integration_hierarchical_timestep_regression`;
- gas identity and sidecar layout: `unit_gas_cell_identity_invariants`, `unit_hot_cold_sidecar_layout`, `integration_reorder_compaction_sidecars`;
- migration/transform invariants: `integration_species_migration_invariants`, `integration_transform_fuzz_invariants`;
- softening ownership: `integration_softening_ownership_invariants`;
- snapshot/restart/provenance contracts runnable in CPU-only mode: `unit_snapshot_hdf5_schema`, `unit_restart_checkpoint_schema`, `integration_snapshot_hdf5_roundtrip`, `integration_restart_checkpoint_roundtrip`, `integration_provenance_roundtrip`;
- CTest registration audit and optional-feature negative check: `integration_runtime_truth_ctest_labels`.

Expected gate behavior:

- `test-stage1-runtime-truth-cpu-debug` runs only the P0 runtime-truth list above.
- Every P0 runtime-truth test has the `runtime_truth` and `p0` CTest labels plus a subsystem label such as `softening`, `sidecar`, `gas`, `migration`, `restart`, `provenance`, `scheduler`, or `active_views`.
- CPU-only builds do not register HDF5 app-smoke, CUDA, Python, or MPI multi-rank tests; `integration_runtime_truth_ctest_labels` fails if those feature-gated tests appear while their feature option is disabled.
- Feature-enabled presets still fail during configure when their requested dependency is missing because the CMake feature options use required dependency discovery (`MPI`, `HDF5`, `FFTW`, `CUDAToolkit`, `Python3`, and `pybind11`).

Feature-specific gates remain explicit and dependency-correct:

```bash
./scripts/ci/run_preset_pipeline.sh hdf5-debug build-hdf5-debug test-hdf5-debug \
  unit_snapshot_hdf5_schema\|unit_restart_checkpoint_schema\|integration_snapshot_hdf5_roundtrip\|integration_restart_checkpoint_roundtrip \
  ci_artifacts/hdf5_runtime_truth 0
./scripts/ci/run_preset_pipeline.sh pm-hdf5-fftw-debug build-pm-hdf5-fftw-debug test-pm-hdf5-fftw-debug \
  unit_pm_solver\|integration_pm_periodic_mode\|integration_tree_pm_coupling_periodic \
  ci_artifacts/pm_hdf5_fftw 0
```

Do not commit `build/` directories, `ci_artifacts/`, CTest XML, or generated feature/metadata outputs. Archive them only as CI artifacts.
