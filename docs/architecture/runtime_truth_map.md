# Runtime truth map (Stage 0 freeze audit)

_Date: 2026-04-25 (UTC)_

## Scope and method

This document is a code-first runtime-state ownership audit. It was compiled from direct repository inspection and command-backed build/test checks, before runtime-state implementation edits.

### Commands used for discovery

- `rg --files -g 'AGENTS.md'`
- `rg --files`
- `rg -n "time_bin|timestep|active_set|scheduler|reorder|species|sidecar|softening|restart|checkpoint|snapshot|normalize|derived|provenance|migration|compact|view|cache|mirror|cell_id|gas" src include tests docs`
- targeted `sed -n` reads of core/workflow/io/config/provenance files listed below

### Build/test commands executed

- `cmake --preset cpu-only-debug`
- `cmake --build --preset build-cpu-debug`
- `ctest --preset test-cpu-debug --output-on-failure`

## Exact files inspected

### Core runtime state and scheduling

- `include/cosmosim/core/simulation_state.hpp`
- `src/core/simulation_state.cpp`
- `src/core/simulation_state_structures.cpp`
- `src/core/simulation_state_ownership.cpp`
- `src/core/simulation_state_species.cpp`
- `src/core/simulation_state_active_views.cpp`
- `include/cosmosim/core/time_integration.hpp`
- `src/core/time_integration.cpp`
- `src/core/simulation_state_metadata.cpp`

### Runtime orchestration / mutation authority

- `src/workflows/reference_workflow.cpp`
- `include/cosmosim/workflows/reference_workflow.hpp`

### Config / normalization / derived values

- `include/cosmosim/core/config.hpp`
- `src/core/config.cpp`
- `include/cosmosim/core/simulation_mode.hpp`
- `src/core/simulation_mode.cpp`

### Restart/checkpoint/snapshot/provenance

- `include/cosmosim/io/restart_checkpoint.hpp`
- `src/io/restart_checkpoint.cpp`
- `include/cosmosim/io/snapshot_hdf5.hpp`
- `src/io/snapshot_hdf5.cpp`
- `include/cosmosim/core/provenance.hpp`
- `src/core/provenance.cpp`
- `src/io/internal/io_contract.cpp`

### Tests inspected for coverage mapping

- `tests/unit/test_simulation_state.cpp`
- `tests/unit/test_time_integration.cpp`
- `tests/unit/test_species_state_organization.cpp`
- `tests/unit/test_hot_cold_sidecar_layout.cpp`
- `tests/unit/test_config_parser.cpp`
- `tests/unit/test_restart_checkpoint_schema.cpp`
- `tests/unit/test_snapshot_hdf5_schema.cpp`
- `tests/unit/test_units_cosmology_provenance.cpp`
- `tests/integration/test_hierarchical_time_bins.cpp`
- `tests/integration/test_hierarchical_timestep_regression.cpp`
- `tests/integration/test_time_integration_loop.cpp`
- `tests/integration/test_reorder_compaction_sidecars.cpp`
- `tests/integration/test_species_migration_invariants.cpp`
- `tests/integration/test_soa_species_pipeline.cpp`
- `tests/integration/test_species_mixed_state.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/integration/test_snapshot_hdf5_roundtrip.cpp`
- `tests/integration/test_provenance_roundtrip.cpp`
- `tests/integration/test_parallel_two_rank_restart.cpp`

---

## Domain map

## 1) Timestep and timestep-bin state

### Storage locations (current truth)

- Per-particle timestep bin: `SimulationState::particles.time_bin`.
- Per-cell timestep bin: `SimulationState::cells.time_bin`.
- Integrator coarse state: `IntegratorState::{current_time_code,current_scale_factor,dt_time_code,step_index,scheme,time_bins}`.
- Scheduler authoritative bin timeline: `HierarchicalTimeBinScheduler` internals (`m_hot`, `m_elements_by_bin`, `m_active_elements`, `m_current_tick`).

### Mutation paths

- Scheduler mutates bin/tick truth via `reset`, `setElementBin`, `requestBinTransition`, `beginSubstep` + `endSubstep`.
- Runtime bin adaptation mutates pending transitions in `updateAdaptiveTimeBins(...)`.
- State mirrors refreshed from scheduler in `syncTimeBinMirrorsFromScheduler(...)`.
- Mirror drift detection is now explicit via `timeBinMirrorsMatchScheduler(...)` / `debugAssertTimeBinMirrorAuthorityInvariant(...)`.
- Restart load mutates integrator/scheduler via `readRestartCheckpointHdf5(...)` and `importPersistentState(...)`.

### Read paths

- Active element extraction for stepping: `scheduler.beginSubstep()`.
- Stage execution uses `ActiveSetDescriptor` built from scheduler active elements.
- Diagnostics/provenance read bin distribution via scheduler diagnostics and metadata payload fields.

### Serialization paths

- Restart serializes:
  - `IntegratorState` scalar fields and `time_bins` fields.
  - Full `TimeBinPersistentState` arrays in `/scheduler/*`.
  - Particle/cell `time_bin` arrays under `/state`.
- Snapshot currently writes/reads particle data but imported particles are assigned `time_bin = 0` during read; no persistent timestep-bin continuity in snapshot I/O.

### Cache/mirror/view classification

- `SimulationState::{particles.time_bin,cells.time_bin}` is a mirror of scheduler truth during workflow execution (rebuilt from scheduler each substep).
- `IntegratorState::time_bins` is metadata-level mirror and not sufficient alone to reconstruct scheduler occupancy/membership.

### Reorder/resize/migration behavior

- Particle reorder/migration physically carries particle `time_bin` with moved rows.
- Gas cell rebuild from gas-particle IDs carries cell `time_bin` through `GasCellMigrationRecord`.

### Restart/reload behavior

- Restart restores exact scheduler persistent vectors and verifies payload integrity.
- Post-read scheduler can be reconstructed exactly (subject to schema compatibility and integrity checks).

### Ambiguity/duplication

- **Ambiguous authority by design**: runtime has three lanes for bin-like info (`scheduler hot metadata`, `state particle/cell time_bin`, `integrator_state.time_bins`); workflow currently treats scheduler as authority and syncs state arrays after each substep.

### Tests covering this domain

- `tests/unit/test_time_integration.cpp`
- `tests/integration/test_hierarchical_time_bins.cpp`
- `tests/integration/test_hierarchical_timestep_regression.cpp`
- `tests/integration/test_time_integration_loop.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/unit/test_restart_checkpoint_schema.cpp`

### Remaining tests (audit-identified)

- Snapshot read/write test asserting intentional loss (or future preservation) of timestep bins, to avoid accidental contract drift.

---

## 2) Active-set construction

### Storage locations

- Authoritative active indices: `HierarchicalTimeBinScheduler::m_active_elements` (transient per substep).
- Step-local descriptor: `ActiveSetDescriptor` in orchestrator calls.
- Workspace compact buffers: `TransientStepWorkspace` and active views (`ParticleActiveView`, `CellActiveView`, `GravityParticleKernelView`, `HydroCellKernelView`).

### Mutation/read paths

- Mutated in scheduler `rebuildActiveSet()` and workflow active-set split loop.
- Read by stage callbacks through `StepContext.active_set`.
- Gravity/hydro compact views gather/scatter from state arrays per active indices.

### Serialization

- Active-set vectors are not serialized directly.
- Restart serializes scheduler persistent state, from which active sets can be recomputed.

### Cache/mirror/view paths

- Active views are compact transient materializations (gather buffers) with explicit scatter back to state.

### Reorder/resize/migration behavior

- Active sets are index-based and therefore sensitive to reorder/migration operations; downstream paths rely on updated indices from reorder/migration routines.
- Active eligibility is species-agnostic in the scheduler contract; migration changes tags/sidecars but does not itself change active membership unless scheduler bins are reassigned through scheduler APIs.

### Restart/reload behavior

- Active sets are recomputed after restart via scheduler state rather than loaded as persisted lists.

### Ambiguity/duplication

- No second persistent active-set truth detected; active views are clearly transient.

### Tests covering this domain

- `tests/unit/test_time_integration.cpp`
- `tests/integration/test_time_integration_loop.cpp`
- `tests/integration/test_hierarchical_time_bins.cpp`

### Missing tests

- Direct tests for `build*/scatter*ActiveView` stale-index fault handling under reordered/migrated state.

---

## 3) Particle and species ordering

### Storage locations

- Particle SoA arrays in `SimulationState::particles`.
- Ordering keys/metadata in `SimulationState::particle_sidecar` (`particle_id`, `sfc_key`, `species_tag`, `owning_rank`, flags).
- Species index maps in `SimulationState::particle_species_index`.
- Species counts in `SimulationState::species.count_by_species`.

### Mutation paths

- Full reorder through `reorderParticles(...)` with `ParticleReorderMap`.
- Species map rebuild via `SimulationState::rebuildSpeciesIndex()`.
- Migration commit (`commitParticleMigration`) compacts/remaps and then rebuilds species counts/index.

### Read paths

- Species-local loops use `particle_species_index.globalIndices(...)`.
- Workflows and I/O read sidecar tags/counts for classification and packing.

### Serialization

- Restart serializes all particle arrays, particle sidecar lanes, and species counts.
- Snapshot serializes particle content by type groups; species re-derived on read from type mapping.

### Cache/mirror/view

- `ParticleSpeciesIndex` is derived/cache-like index map rebuilt from sidecar tags.

### Reorder/resize/migration behavior

- Reorder supports sidecar policies (`kMoveWithParent` or remapped parent indirection).
- Migration commit rewrites particle arrays and sidecar rows, then recomputes species ledger/index.

### Restart/reload behavior

- Restart reads species counts and sidecars, then calls `rebuildSpeciesIndex()` + invariant validation.

### Ambiguity/duplication

- Duplicated species truth exists intentionally:
  - `species.count_by_species` (ledger),
  - `particle_sidecar.species_tag` (row truth),
  - `particle_species_index` (derived index).
- Ownership invariant checks enforce consistency but mutation authority is spread across multiple functions.

### Tests covering

- `tests/unit/test_simulation_state.cpp`
- `tests/unit/test_species_state_organization.cpp`
- `tests/integration/test_reorder_compaction_sidecars.cpp`
- `tests/integration/test_soa_species_pipeline.cpp`
- `tests/integration/test_species_mixed_state.cpp`

### Missing tests

- Fault-injection coverage for partial updates where species ledger changed without index rebuild.

---

## 4) Species sidecars

### Storage locations

- `SimulationState::{star_particles, black_holes, tracers}` as index-addressed sidecars keyed by `particle_index`.

### Mutation paths

- Reorder remaps sidecar parent indices according to policy.
- Migration commit rebuilds each sidecar using kept rows plus inbound records.
- Physics callbacks mutate sidecar payload values (stellar evolution, BH, tracer updates).

### Read paths

- Physics kernels and diagnostics read sidecars by species.
- Migration pack path reads required sidecar rows by parent particle index.

### Serialization

- Restart serializes/deserializes all sidecar fields (including stellar channel arrays).
- Snapshot serializes selected sidecar-derived datasets where schema supports them (e.g., tracer host metadata, softening sidecar).

### Cache/mirror/view

- No separate persistent mirror found; sidecars are canonical for species-specific cold data.

### Reorder/resize/migration behavior

- Explicit sidecar sync policies and stale-index guards exist.
- Migration can drop/rebuild sidecar rows as ownership changes.

### Restart/reload behavior

- Restart reads sidecars then validates one-to-one invariants vs species-tagged particles.

### Ambiguity/duplication

- Parent link duplication: species tag lane and sidecar membership both encode species identity.

### Tests covering

- `tests/unit/test_simulation_state.cpp`
- `tests/integration/test_reorder_compaction_sidecars.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`

### Missing tests

- Snapshot roundtrip completeness assertions for all sidecar families (currently strongest guarantees are restart-focused).

---

## 5) Gas cell identity

### Storage locations

- Gas cell structural lanes: `SimulationState::cells`.
- Hydro/cold lanes: `SimulationState::gas_cells`.
- Identity bridge to particles: implicit mapping through `particle_species_index.globalIndices(kGas)` and particle IDs.

### Mutation paths

- Drift callback syncs cell centers from gas particle positions.
- Compaction/migration paths rebuild cell arrays from gas-particle ID keyed `GasCellMigrationRecord` map.
- Host cell indices for BH/tracers remapped after rebuild.

### Read paths

- Hydro kernels read from cell/gas arrays (often through active compact views).
- Migration and diagnostics read gas-particle↔cell correspondence.

### Serialization

- Restart serializes `cells` and `gas_cells` fully.
- Snapshot focuses on particle-centric schema; gas cell arrays are not a first-class persisted group there.

### Cache/mirror/view

- Cell centers can behave like mirror-of-gas-particles during drift-coupled paths.

### Reorder/resize/migration behavior

- Gas identity continuity across migration is preserved via particle IDs, not raw cell indices.

### Restart/reload behavior

- Restart restores explicit cell arrays and patch indices.

### Ambiguity/duplication

- **Ambiguous identity authority**: gas identity is partly index-based (cell index), partly particle-ID-based (migration rebuild), and partly species-order-based (gas global index list).

### Tests covering

- `tests/integration/test_reorder_compaction_sidecars.cpp`
- `tests/integration/test_hydro_sod_like.cpp`
- `tests/integration/test_amr_static_refinement_sync.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`

### Missing tests

- Explicit invariant tests for gas cell identity continuity across reorder + migration + restart combined path.

---

## 6) Softening state

### Storage locations

- Config-level defaults/species lanes in `SimulationConfig::numerics`.
- Optional per-particle softening value sidecar: `particle_sidecar.gravity_softening_comoving`; authoritative override membership is `particle_sidecar.has_gravity_softening_override`. A value without a mask bit is a materialized default/diagnostic mirror, not override authority.
- Provenance lanes (`ProvenanceRecord` softening policy/kernel/epsilon).

### Mutation paths

- Initial population from species policy in `maybeInitializeParticleSofteningFromSpeciesPolicy(...)`.
- Reordered/migrated alongside particle rows.
- Snapshot import may populate sidecar if `GravitySofteningComoving` dataset exists.

### Read paths

- Time-step criteria gravity hook uses sidecar if present, else species/global fallback.
- Gravity modules consume policy and sidecar values for force/timestep behavior.

### Serialization

- Restart serializes full sidecar when present.
- Snapshot writes per-particle softening additively and reads it if present.

### Cache/mirror/view

- Sidecar is canonical runtime lane for per-particle softening overrides; config values are policy defaults.

### Reorder/resize/migration behavior

- Reorder/migration preserve softening lane with particle rows.

### Restart/reload behavior

- Restart restores sidecar or empty state (policy-only mode).

### Ambiguity/duplication

- Duplicate representational lanes exist: global/species config policy plus optional per-particle sidecar overrides.

### Tests covering

- `tests/integration/test_snapshot_hdf5_roundtrip.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/unit/test_restart_checkpoint_schema.cpp`

### Missing tests

- Coverage for conflict resolution when config species policy and snapshot softening sidecar disagree.

---

## 7) Config-derived runtime values

### Storage locations

- Canonical typed config: `FrozenConfig::config`.
- Canonical normalized text/hash: `FrozenConfig::{normalized_text, provenance.config_hash_hex}`.
- Runtime derived policy objects built downstream (`ModePolicy`, per-module derived values).

### Mutation paths

- Parsing + normalization in `normalizeValidateFreeze(...)`.
- Compatibility alias/scalar lanes collapsed into axis-aware typed fields.
- Validation enforces bounds/required relationships.

### Read paths

- Reference workflow and modules consume typed config.
- Provenance generation and output schemas consume normalized hash + selected derived fields.

### Serialization

- Normalized config text written to run directory and embedded in snapshot/restart attributes.
- Hash continuity checked by shared IO contract validation.

### Cache/mirror/view

- Normalized text acts as canonical serialized mirror of typed config.

### Reorder/resize/migration behavior

- Not applicable directly; config drives policies that influence these paths.

### Restart/reload behavior

- Restart requires normalized text/hash/provenance agreement.

### Ambiguity/duplication

- Intentional backward-compatibility duplication exists (e.g., scalar and axis-aware PM/box keys); normalized output collapses to canonical axis-aware form.

### Tests covering

- `tests/unit/test_config_parser.cpp`
- `tests/unit/test_simulation_mode.cpp`
- `tests/integration/test_config_examples.cpp`
- `tests/integration/test_simulation_mode_toy_runs.cpp`

### Missing tests

- Cross-check tests asserting derived runtime constants used by workflow callbacks remain in lockstep with normalized config text for every mode.

---

## 8) Restart and checkpoint truth

### Storage locations

- Authoritative restart payload: `RestartWritePayload` (state + integrator + scheduler + provenance + normalized config + distributed gravity state).
- On-disk groups: `/state`, `/integrator`, `/scheduler`, `/distributed_gravity`, plus root hash/schema attributes.

### Mutation paths

- Restart writer assembles payload, computes integrity hash, writes temp file, finalizes atomically.
- Reader reconstructs structures, re-computes integrity hash, validates schema + invariants.

### Read paths

- Workflow restart path loads restart payload then resumes scheduler/integrator/state.

### Serialization overlap with snapshot

- Overlap domains: particles, IDs, softening sidecar (optional), normalized config, provenance.
- Non-overlap: restart carries scheduler persistent state, module sidecars, strict continuation integrity hash.

### Cache/mirror/view

- No separate restart cache detected; restart payload integrity hash is the main consistency guard.

### Reorder/resize/migration behavior

- Restart stores post-mutation state as-is; no special reorder metadata beyond arrays/index fields and sidecars.

### Restart/reload behavior

- Exact continuation target with strict schema compatibility and integrity checks.

### Ambiguity/duplication

- Snapshot and restart both carry normalized config and provenance but with different completeness guarantees; this is an intentional dual-contract surface.

### Tests covering

- `tests/unit/test_restart_checkpoint_schema.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/integration/test_parallel_two_rank_restart.cpp`
- `tests/integration/test_provenance_roundtrip.cpp`

### Missing tests

- Contract test asserting snapshot->restart->snapshot roundtrip equivalence boundaries (what is intentionally lossy vs exact).

---

## Cross-domain duplication/mirror inventory (classified)

- `scheduler hot metadata` vs `state.*.time_bin`: **mirror** (scheduler-authoritative during loop).
- `species_tag` vs `species.count_by_species` vs `ParticleSpeciesIndex`: **canonical + derived + ledger duplication** with explicit invariant checks.
- Gas identity by `cell index` vs `gas species global index` vs `gas particle_id`: **ambiguous multi-lane identity** (currently reconciled during compaction via particle IDs).
- Softening policy config vs optional per-particle softening sidecar: **policy + override dual truth**.
- Normalized config hash in metadata/provenance/restart attributes: **redundant consistency mirrors** validated by contract checks.

## Repair tickets / follow-up prompts (no broad code change in this stage)

1. **P31-STATE-TEST-DRIFT-026 (Open)**
   - Symptom: `unit_snapshot_hdf5_schema` currently asserts `gadget_arepo_v2` while runtime schema reports `gadget_arepo_v4` during `ctest --preset test-cpu-debug`.
   - Scope: test/schema contract alignment only.

2. **P32-STATE-OWNER-AMBIGUITY-027 (Open)**
   - Scope: make mutation authority explicit for scheduler bins vs state bin mirrors and gas cell identity lanes; add invariant checks/tests before behavior changes.

3. **P33-STATE-SNAPSHOT-RESTART-BOUNDARY-028 (Open)**
   - Scope: codify and test intentional overlap/loss boundaries between snapshot and restart (especially timestep bins, scheduler state, module sidecars).

4. **P34-STATE-SOFTENING-OVERRIDE-029 (Open)**
   - Scope: explicit policy for precedence/compatibility between config species softening and imported per-particle softening sidecar values.

## Build/test status for this audit

- Configure/build with `cpu-only-debug` preset succeeded.
- Full `ctest --preset test-cpu-debug --output-on-failure` did **not** complete cleanly:
  - `unit_snapshot_hdf5_schema` failed on schema name assertion (`gadget_arepo_v2` expected).
  - Test run subsequently stalled in this environment after progressing into later integration tests, so this audit records a partial run with explicit blocker instead of claiming green closure.

## Reproducibility impact statement

This stage is documentation/audit only. No runtime behavior, schema payload, deterministic scheduling mode, config normalization flow, or provenance completeness path was modified by this change.

### Stage 0 follow-up note: migration softening value vs override authority

Migration payloads now distinguish `has_gravity_softening_value` from `has_gravity_softening_override`. The former preserves a materialized numeric value when present; the latter is the only flag that carries per-particle override authority across migration.
