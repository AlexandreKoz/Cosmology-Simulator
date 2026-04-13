# Repair state recap (post-repair audit snapshot)

_Date captured: 2026-04-07 (UTC)_

This recap records **current command-backed audit evidence** for the emergency repair closeout pass.

## 1) Preset and dependency gates

Commands:

```bash
cat CMakePresets.json
cmake --preset hdf5-debug
cmake --preset pm-hdf5-fftw-debug
```

Observed:

- Presets now include `hdf5-debug` and `pm-hdf5-fftw-debug`.
- Dependency failure behavior is explicit and actionable.
- In this environment, both dependency-enabled presets fail at configure with:
  - `COSMOSIM_ENABLE_HDF5=ON but HDF5 was not found.`

Interpretation:

- Preset-level closure for P05/P06 exists in-tree.
- Runtime closeout for HDF5/PM feature paths is blocked by missing dependency in this environment.

## 2) CPU-only baseline evidence

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4
ctest --preset test-cpu-debug --output-on-failure
```

Outcome:

- Configure: **PASS**
- Build: **PASS**
- Tests: **PASS** (`36/36`)

Interpretation:

- Baseline path remains healthy and no unrelated regressions were observed.

## 3) Guardrail scripts

Commands:

```bash
./scripts/ci/check_repo_hygiene.sh
./scripts/ci/guard_feature_paths.sh
```

Outcome:

- `check_repo_hygiene.sh`: **PASS**.
- `guard_feature_paths.sh`: CPU segment passes; feature segment fails at HDF5 configure gate, as expected when dependency is absent.

Interpretation:

- Hygiene and process checks are active and correctly prevent false success claims from CPU-only evidence.

## 4) IC reader gas import state evidence

Command:

```bash
nl -ba src/io/ic_reader.cpp | sed -n '562,592p'
```

Observed:

- For `PartType0`, IC import now writes:
  - `gas_cells.density_code[cell_i]`
  - `gas_cells.internal_energy_code[cell_i]`
- Unsupported reports are for unmapped metallicity/smoothing length, not thermodynamic bypass.

Related test artifact:

- `tests/unit/test_ic_reader.cpp` includes HDF5-gated checks for gas thermodynamic mapping and optional-density defaulting behavior.

Interpretation:

- P04 implementation appears closed in code/tests.
- Full runtime verification remains dependent on HDF5 availability.

## 5) Docs and naming/process hygiene

Commands:

```bash
test -f docs/build_instructions.md
test -f CONTRIBUTING.md
find . -maxdepth 1 -type f | sed 's#^./##' | sort
```

Outcome:

- `docs/build_instructions.md`: present.
- `CONTRIBUTING.md`: present.
- Root file naming is now policy-safe (no legacy violating spreadsheet at repo root).

Interpretation:

- P08/P09/P10 are closed at repository state level.

## 6) Audit conclusion pointer

See `docs/repair_closeout_report.md` for the stop/go decision.
Current state is **STOP** until HDF5 (and then PM/HDF5/FFTW) preset build+test commands pass on the intended feature paths.

## 7) Core boundary repair (reference workflow assembly)

Commands:

```bash
ctest --test-dir build/cpu-only-debug --output-on-failure -R "integration_core_dependency_direction|integration_reference_workflow"
```

Observed:

- Core dependency-direction guard passes and fails fast on forbidden upward includes in `include/cosmosim/core/**` and `src/core/**`.
- Reference workflow integration smoke test still passes after moving assembly to `workflows/`.

Interpretation:

- The architecture boundary leak from `core` into analysis/I/O/physics workflow assembly is repaired while keeping smoke behavior intact.

## 8) Configuration-contract hardening (typed freeze path)

Commands:

```bash
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_config_parser|unit_simulation_mode|integration_simulation_mode_toy_runs"
```

Observed:

- Policy-like config values are parsed from param strings, then frozen into typed enums for solver selection, coordinate frame, mode boundaries, and feedback mode/variant.
- Unknown key handling and deprecated alias mapping use a centralized key/alias registry in `src/core/config.cpp`.
- Unit/integration coverage validates unknown key rejection, alias behavior, invalid enum rejection, and deterministic normalized config/hash behavior.

Interpretation:

- String-driven runtime policy drift is reduced after freeze while preserving param-style UX and deterministic normalized-config provenance semantics.
