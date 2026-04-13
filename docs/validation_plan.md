# Validation plan

CosmoSim uses a four-level validation ladder: unit, integration, regression, convergence.

## Validation structure in-repo

- Unit: `tests/unit/`
- Integration: `tests/integration/`
- Validation ladder: `tests/validation/`
- Input decks: `validation/input_decks/`
- Frozen tolerances: `validation/reference/validation_tolerances_v1.txt`

## Required command sequence

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

For PM/TreePM and HDF5 flows:

```bash
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
```

## Validation invariants

- No silent comoving/proper or code/SI unit mixing.
- Mode and boundary policy must remain explicit and validated.
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

## Reference decks

Use `validation/input_decks/*.param.txt` for deterministic validation runs:

- `gravity_two_body.param.txt`
- `hydro_sod_1d.param.txt`
- `cooling_box.param.txt`
- `amr_reflux_sync.param.txt`
