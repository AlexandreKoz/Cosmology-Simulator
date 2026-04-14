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

## 9) SimulationState ownership/invariants decomposition repair

_Date captured: 2026-04-13 (UTC)_

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_simulation_state
./build/cpu-only-debug/test_unit_simulation_state
```

Observed:

- `src/core/simulation_state.cpp` responsibility concentration was split by invariant boundary:
  - storage-lane consistency (`simulation_state_structures.cpp`)
  - ownership checks (`simulation_state_ownership.cpp`)
  - species indexing and transfer packing (`simulation_state_species.cpp`)
  - metadata/module-sidecar serialization (`simulation_state_metadata.cpp`)
  - active-view assembly/scatter (`simulation_state_active_views.cpp`)
  - reorder/scratch logic remains in `simulation_state.cpp`
- Header-level hot-field contract notes now explicitly name allowed gravity/hydro active-view lanes.
- `tests/unit/test_simulation_state.cpp` now covers:
  - ownership invariant failure/recovery
  - unique-ID invariant failure/recovery
  - species-index rebuild correctness
  - transfer packet pack/unpack-equivalence checks
  - metadata serialize/deserialize round-trip
  - active-view hot-field writeback behavior and cold-lane protection
  - static assertion guardrails for gravity/hydro active-view compactness.

Interpretation:

- State ownership and active-view invariants are now reviewable in focused files with explicit hot/cold contracts and command-backed checks.

## 10) Snapshot/restart I/O contract boundary hardening

_Date captured: 2026-04-13 (UTC)_

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_snapshot_hdf5_schema test_unit_restart_checkpoint_schema test_integration_restart_checkpoint_roundtrip
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_snapshot_hdf5_schema|unit_restart_checkpoint_schema|integration_restart_checkpoint_roundtrip"
```

Observed:

- Shared continuation-metadata contract names/validation are centralized via `include/cosmosim/io/io_contract.hpp` + `src/io/internal/io_contract.cpp`.
- Restart docs and API now publish an exact-restart checklist (`exactRestartCompletenessChecklist()`), separating restart-completeness obligations from snapshot interoperability.
- New negative checks cover restart schema mismatch rejection, missing required scheduler dataset rejection, and finalize-failure behavior that leaves target path untouched while preserving the temporary artifact for diagnosis.

Interpretation:

- Snapshot and restart responsibilities are now explicitly separated and reviewable, with stronger error behavior on continuation-critical contract violations.

## 11) Operational observability repair (structured run events)

_Date captured: 2026-04-13 (UTC)_

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_profiling test_integration_reference_workflow
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_profiling|integration_reference_workflow"
```

Observed:

## 12) Diagnostics maturity-tier repair (analysis honesty/scalability guard)

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_config_parser test_unit_analysis_diagnostics test_integration_analysis_bundle
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_config_parser|unit_analysis_diagnostics|integration_analysis_bundle"
```

Observed:

- Analysis config now carries typed `diagnostics_execution_policy` with normalized-config persistence.
- Diagnostics bundles now emit per-diagnostic maturity metadata (`tier`, `maturity`, `scalability`, execution policy).
- Provisional heavy reference diagnostics (power spectrum direct summation) no longer run under default policy.
- Unit/integration tests now assert that heavy provisional diagnostics only run when explicitly opted in.

Interpretation:

- Infrastructure run-health counters remain first-class and cheap.
- Validated light science diagnostics remain available in default runs.
- Reference/provisional heavy diagnostics are quarantined behind explicit non-default policy.

- `core::ProfilerSession` now carries a minimal structured runtime event model (kind, severity, subsystem, optional step/time/scale context, message, key/value payload).
- Reference workflow emits a machine-readable operational report (`reference_operational_events.json`) linked to deterministic config provenance via `provenance_config_hash_hex`.
- Key infrastructure lifecycle/failure surfaces are explicit in event records (config freeze validation, restart/snapshot write/read begin/complete/failure).

Interpretation:

- Operational troubleshooting and CI artifact review no longer depend on ad hoc text/exception surfaces alone.
- Reproducibility posture remains unchanged: operational events are additive and provenance-linked; no solver behavior or restart/snapshot schema semantics were changed.
