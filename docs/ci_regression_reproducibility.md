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
- deterministic config/provenance-oriented checks (`unit_config_parser`, `integration_feature_summary`),
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
