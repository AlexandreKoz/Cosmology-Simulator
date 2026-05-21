# Stage 5 config contract audit map

Date: 2026-05-21  
Scope: audit-only infrastructure repair note for config-contract hardening. No solver behavior changes.

## 1) Raw input → normalized config → derived/runtime state inventory

- **Raw input boundary (must remain string/alias tolerant):** `src/core/config.cpp` parser + key registry + alias translation (`parseKeyValueText`, normalization/freeze path). Keep legacy aliases input-boundary-only here.
- **Typed normalized contract (authoritative):** `core::SimulationConfig` / `core::FrozenConfig` in `include/cosmosim/core/config.hpp` and freeze/validation in `src/core/config.cpp`.
- **Derived runtime projection (must not become new truth owner):** `core::DerivedRuntimeConfig` and `deriveRuntimeConfig(...)` in `src/core/config.cpp`.
- **Runtime state metadata mirrors:** `core::SimulationState::metadata.normalized_config_hash{,_hex}` in `include/cosmosim/core/simulation_state.hpp`; populated in workflow runtime init.
- **Provenance serialization authority:** `src/core/provenance.cpp` + `include/cosmosim/core/provenance.hpp`.

## 2) Highest-risk mixing points (raw/normalized/derived/runtime)

1. `src/core/config.cpp`
   - Single file currently owns parser, alias policy, validation, and normalized text emission; safe for input-boundary aliases, but high drift risk if runtime-only keys leak back into parser conditionals.
2. `src/workflows/reference_workflow.cpp`
   - Runtime init threads `FrozenConfig`, `SimulationConfig`, normalized text/hash, and live `SimulationState` metadata through multiple setup helpers (`prepareRunDirectory`, `finalizeStateMetadata`, restart/snapshot writer payload construction).
3. `src/io/restart_checkpoint.cpp` and `src/io/snapshot_hdf5.cpp`
   - Continuation payloads carry both typed config-dependent fields and serialized normalized/provenance lanes; compatibility fallbacks exist for legacy provenance attributes and can hide boundary drift if not isolated.
4. `src/python_bindings.cpp`
   - Exposes `FrozenConfig`/config hash outward; bindings are a second entrypoint that must not bypass freeze/validation.

## 3) Runtime entrypoints that must stay frozen-config-only

- `src/core/harness_main.cpp` (`loadFrozenConfigFromFile` → `ReferenceWorkflowRunner`): canonical CLI path already frozen-config-only.
- `src/workflows/reference_workflow.cpp` (`ReferenceWorkflowRunner::runImpl` and helpers): must continue accepting typed/frozen config only; no raw key/value maps.
- `src/io/ic_reader.cpp` public readers: take `const core::SimulationConfig&`; callers must pass frozen typed config, never ad-hoc raw params.
- `src/io/snapshot_hdf5.cpp` and `src/io/restart_checkpoint.cpp`: payloads include typed config pointer/reference + normalized text/hash + provenance; construction must remain centralized from frozen config state.
- `src/python_bindings.cpp` config loaders: must continue routing through `loadFrozenConfigFromFile/String`.

## 4) Stage-5 field focus map (where to patch in follow-up prompts)

- **Cosmology/time integration fields:**
  - `include/cosmosim/core/config.hpp` (`CosmologyConfig`, `TimeConfig`) + validation/normalization in `src/core/config.cpp`.
  - Runtime consumption surfaces: `src/workflows/reference_workflow.cpp`, `src/core/time_integration.cpp`.
- **Solver/module enable flags:**
  - `PhysicsConfig`/`AnalysisConfig` booleans in `include/cosmosim/core/config.hpp`.
  - Consumption in workflow callback registration and module setup: `src/workflows/reference_workflow.cpp`, `src/analysis/diagnostics.cpp`, `src/physics/*` module callbacks.
- **Provenance hash/metadata serialization:**
  - `src/core/provenance.cpp`, plus I/O propagation in snapshot/restart writers.
- **Snapshot/restart metadata contracts:**
  - `include/cosmosim/io/io_contract.hpp`, `include/cosmosim/io/restart_checkpoint.hpp`, `include/cosmosim/io/snapshot_hdf5.hpp` and implementations.

## 5) Legacy aliases that must remain input-boundary-only

- Scalar-to-axis cosmology and PM grid compatibility aliases documented in config docs and handled in parser normalization (`box_size*`, `treepm_pm_grid*` compatibility lanes).
- Snapshot/restart provenance compatibility aliases/fallback attributes in readers must remain read-boundary compatibility only; writers should keep canonical schema attributes.

## 6) Minimal independent prompt map (implementation order)

1. **Prompt S5-01 (registry extraction prep):** Add a tiny internal checklist comment block in `src/core/config.cpp` at alias mapping section to mark “input-boundary-only aliases” and “canonical key emission only” invariants.
2. **Prompt S5-02 (runtime init seam audit):** Centralize/label one helper in `src/workflows/reference_workflow.cpp` for constructing restart/snapshot config/provenance payloads from `FrozenConfig` only.
3. **Prompt S5-03 (entrypoint guardrails):** Add targeted negative tests in `tests/unit/test_config_parser.cpp` / `tests/unit/test_restart_checkpoint_schema.cpp` proving normalized hash/provenance mismatch rejection remains strict.
4. **Prompt S5-04 (docs/contract sync):** Update `docs/configuration.md` + `docs/output_schema.md` only if any canonical key/schema behavior changes; otherwise keep no-op.

## 7) Reproducibility impact statement for this audit note

This patch is documentation-only and introduces no config-key, schema, solver, or runtime-behavior changes. Determinism, normalized config dumps, provenance hashing, output naming, and restart/snapshot compatibility behavior remain unchanged.
