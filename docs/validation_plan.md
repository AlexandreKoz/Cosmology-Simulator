# Validation plan

CosmoSim uses a four-level validation ladder: unit, integration, regression, convergence.

## Validation structure in-repo

- Unit: `tests/unit/`
- Integration: `tests/integration/`
- Validation ladder: `tests/validation/`
- Frozen tolerances: `validation/reference/validation_tolerances_v1.txt`

`validation/input_decks/` is currently documentation/reference material only. It is **not** the authoritative executable runtime path in this repository build. The real runnable path is the config-driven application entry point (`cosmosim_harness <config.param.txt>`) plus the tested release/runtime smoke configs.

## Required command sequence

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

For HDF5 snapshot/restart flows and the real runtime smoke path:

```bash
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug
ctest --preset test-hdf5-debug --output-on-failure
```

For PM/TreePM FFTW-specific validation when FFTW is available:

```bash
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
```

## Gravity Phase 1 evidence classes (explicit references)

- **Open-boundary small-`N` tree check**
  - Method: direct summation, softened with the same `epsilon_comoving` as tree solver.
  - Claim: local tree operator correctness for non-periodic pair accumulation.
- **Periodic PM plane-wave mode check**
  - Method: analytic periodic spectral solution for selected Fourier mode(s).
  - Claim: PM operator/sign/phase correctness in periodic domain.
- **Periodic PM uniform-density cancellation**
  - Method: exact lattice plus deterministic parity-balanced micro-jitter (glass-like proxy) particle placement, periodic Poisson solve with mean subtraction.
  - Claim: near-zero force under uniform density.
- **Periodic TreePM consistency check (integration test)**
  - Method: minimum-image periodic direct reference.
  - Explicit limitation: not Ewald exact.
- **Periodic TreePM force-error mapping (benchmark artifact)**
  - Method: fine spectral PM plus exact pairwise short-range residual periodic proxy reference.
  - Explicit limitation: not Ewald exact.
- **Split-kernel matching checks**
  - Method: Gaussian split SR/LR factor complementarity and TreePM PM-only/tree-only/split comparisons.
  - Claim: split coupling semantics are wired and monotone in expected direction.
- **Two-body orbit controlled run + energy drift monitor**
  - Method: leapfrog integration with tree accelerations and softened two-body energy accounting.
  - Claim: bounded long-run drift in a controlled setup.
- **Static halo radial-force profile**
  - Method: tree force on shell probes versus softened direct-summation reference.
  - Claim: radial-force profile agreement for a static mass distribution.
- **Regression-only checks (AMR reflux / SF mass budget)**
  - Claim: implementation stability against historical behavior, not standalone validation proof.

## Force-error mapping artifact requirement

Run:

```bash
./build/pm-hdf5-fftw-debug/bench_tree_pm_force_error_map
```

Expected deterministic artifact:

- `validation/artifacts/tree_pm_force_error_map.csv` (generated artifact; do not rely on a source-controlled copy)

The CSV is the required Phase 1 force-error map over PMGRID/ASMTH/RCUT and must be attached/referenced in review evidence.

## Validation invariants

- No silent comoving/proper or code/SI unit mixing.
- Mode and boundary policy must remain explicit and validated.
- The runtime application must load a real config and emit a real run directory outcome.
- Restart roundtrip must preserve state consistency and provenance.
- Snapshot schema compatibility and alias-read behavior remain tested.
- Deterministic/reproducible config normalization and hashes remain stable for unchanged inputs.

## Expected tests for feature work

- **Unit:** parser/math/layout or local invariants.
- **Integration:** real pipeline path through owning modules.
- **Regression:** preserve prior output metadata/behavior with tolerances.
- **Convergence:** when numerical order or resolution behavior is affected.

## Validation reporting in PRs

Every PR should include:

- build preset(s),
- exact test commands,
- pass/fail status,
- any intentionally deferred validation and rationale.

## Runnable smoke/config examples

Use one of these config-driven paths for honest runtime smoke checks:

- `configs/release/release_smoke_zoom_in.param.txt`
- `configs/release/release_smoke_isolated_galaxy.param.txt`
- `configs/release/release_smoke_cosmo_cube.param.txt`

Or run the built-in integration smoke gate:

- `integration_runtime_app_smoke` (HDF5-enabled builds)


## Phase 2 distributed TreePM validation gate (MPI gravity gate)

- Distributed workflow snapshot/restart paths are rank-qualified (`..._rank###.hdf5`) so MPI restart validation does not rely on unsafe multi-rank writes to a single file.

The Phase 2 gate is now a dedicated MPI validation suite:

- Executable: `test_validation_phase2_mpi_gravity`
- CTest entries:
  - `validation_phase2_mpi_gravity_single_rank`
  - `validation_phase2_mpi_gravity_two_rank`
  - `validation_phase2_mpi_gravity_three_rank`

### Numerical contracts enforced

- Distributed PM equivalence vs one-rank reference: `rel_L2 <= 1e-10`.
- Distributed full TreePM equivalence vs one-rank reference: `rel_L2 <= 5e-6` and `max_rel <= 5e-5`.
- Communication stress path: tiny tree exchange batches (`tree_exchange_batch_bytes=64`) plus PM cadence refresh/reuse toggles, checked against the same TreePM thresholds.
- Restart continuation contract in MPI mode: reference workflow restart write/read roundtrip must report `restart_roundtrip_ok=true`.

### Phase 2 scaling artifacts

Run in MPI-enabled builds:

```bash
cmake --build --preset build-mpi-hdf5-fftw-debug --target generate_mpi_gravity_scaling_artifacts
```

Outputs (artifact files):

- `validation/artifacts/pm_only_scaling_np1.csv`
- `validation/artifacts/pm_only_scaling_np2.csv`
- `validation/artifacts/tree_only_scaling_np1.csv`
- `validation/artifacts/tree_only_scaling_np2.csv`

For final integration hard-gates and exact command bundle, see:

- `docs/treepm_phase2_closeout.md`
