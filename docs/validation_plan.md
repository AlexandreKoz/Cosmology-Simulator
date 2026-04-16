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
