# Stage 2 adversarial timestep ownership review

Date: 2026-05-11

Scope: infrastructure review only; no solver math, new physics model, or timestep algorithm changes.

## Review posture

This review intentionally assumes the Stage 2 patch is hostile: comments are not accepted as ownership, mirrors are presumed fake until a stale-detection or rebuild path is found, and restart/migration paths are treated as possible split-brain timestep importers.

Evidence sweep commands used during the review:

- `rg --files -g 'AGENTS.md' -g '!build*' -g '!_deps'`
- `rg -n "(time_step|timestep|dt|delta_t|timeStep|step_size|next_step|last_step|active_set|active set|active_mask|active_indices|active_species|mirror|cache|stale|scheduler|Schedule|cadence|restart|migration)" include src tests docs -g '!build*' -g '!*.csv' -g '!*.json'`
- `rg -n "(ActiveSetDescriptor\\{|activeElements\\(|active_flag|particles\\.time_bin|cells\\.time_bin|pending_bin_index|\\.bin_index\\s*=|TimeBinPersistentState\\{|importPersistentState|exportPersistentState|makeSchedulerActiveSetDescriptor|syncTimeBinMirrorsFromScheduler|debugAssertTimeBinMirrorAuthorityInvariant|current_tick|next_activation_tick)" src include tests -g '!*.json'`

## Blocker findings

### B1 — Restart cell mirror validation had an exact-size bypass

**Status:** patched in this change.

**Evidence:** `src/io/restart_checkpoint.cpp` previously validated `CellSoa::time_bin` only when the scheduler `bin_index` array was exactly cell-sized. The v6 scheduler lane is particle-sized, and valid particle-bound gas fixtures commonly have `cells.size() != particles.size()`. That left stale gas-cell timestep mirrors able to pass `restartPayloadIntegrityHash(...)`, `writeRestartCheckpointHdf5(...)`, and `readRestartCheckpointHdf5(...)` validation whenever particle mirrors matched scheduler truth but cell mirrors did not.

**Failure mode:** a restart payload with correct particle scheduler state but stale `CellSoa::time_bin` could survive validation in mixed particle/gas states, preserving a fake mirror that later diagnostics or hydro-facing code could mistake for timestep truth.

**Required closure evidence:**

- code must validate each non-empty cell mirror through `requireParticleBoundGasCellContract(...)` and `gasParticleIndexForCellRow(...)` against the parent particle's scheduler `bin_index`;
- restart import mirror rebuild must use the same parent-particle mapping instead of exact-size fallback;
- CPU-only regression must create a particle count larger than gas-cell count, intentionally stale a cell mirror, and prove restart payload validation rejects it;
- HDF5 restart round-trip should be rerun when HDF5 is available because this is restart-schema-adjacent compatibility behavior, even though no payload field or schema version changed.

## Serious findings

### S1 — Direct public `ParticleSoa::time_bin` and `CellSoa::time_bin` writes remain easy in tests and import scaffolds

**Status:** acceptable only as a documented test-fixture risk; do not treat mirror comments as ownership.

**Evidence:** multiple tests still assign `state.particles.time_bin[...]` or `state.cells.time_bin[...]` directly to build fixtures. Production write/import paths that intentionally materialize external data also initialize mirrors before scheduler reconciliation.

**Risk:** fixture code can hide production defects by making mirrors agree manually instead of proving `syncTimeBinMirrorsFromScheduler(...)` or restart import rebuild produced agreement.

**Required follow-up prompt:** convert Stage 2 timestep-authority tests that assert scheduler behavior to build mirrors through scheduler sync helpers except where the test is explicitly negative and mutates mirrors to prove stale detection.

### S2 — `TimeBinPersistentState` is intentionally writable as a DTO

**Status:** acceptable risk for restart/migration boundaries only; every importer must remain hostile.

**Evidence:** tests and restart decode construct `TimeBinPersistentState` directly, and scheduler authority is restored through `HierarchicalTimeBinScheduler::importPersistentState(...)`.

**Risk:** any future importer that uses `TimeBinPersistentState::active_flag` or `pending_bin_index` directly as live scheduler state without import validation would reintroduce split-brain timestep truth.

**Required follow-up prompt:** keep all new restart/migration transforms behind `importPersistentState(...)`, `remapSchedulerPersistentStateByParticleId(...)`, or `rebuildSchedulerPersistentStateFromIdentityRecords(...)`; add negative tests for any new transform.

### S3 — Active-set descriptors can still be manually constructed for non-scheduler unit callbacks

**Status:** acceptable only outside hierarchical scheduler execution; solver-runtime paths must carry scheduler provenance.

**Evidence:** `StepOrchestrator::executeSingleStep(...)` enforces scheduler provenance only when an expected scheduler tick is supplied. Baseline/unit tests can still pass manually constructed descriptors for non-hierarchical paths.

**Risk:** a future hierarchical workflow could omit `expected_scheduler_tick` and downgrade active-set verification to generation/index freshness only.

**Required follow-up prompt:** any workflow enabling hierarchical bins must call `makeSchedulerActiveSetDescriptor(...)` and pass `expected_scheduler_tick`; add workflow-level negative coverage if a new hierarchical execution path is introduced.

## Minor findings

### M1 — Documentation needs to stay explicit that cell mirrors map through gas parent identity

**Status:** patched in this change.

The restart/output schema guide now records the v6 compatibility behavior: cell timestep mirrors remain payload mirrors but are validated/rebuilt against the parent gas particle's scheduler entry when cells are present.

### M2 — HDF5 gate remains feature-path evidence, not optional closure

**Status:** documented validation requirement.

The CPU-only regression is sufficient to prove the validation function rejects the loophole, but full restart closure still requires `hdf5-debug` restart round-trip coverage in an environment with HDF5 available.

## Acceptable risks

- `HierarchicalTimeBinScheduler::hotMetadata()` exposes a const reference only. The writable scheduler lanes remain private and are mutated through scheduler methods or validated import.
- `active_flag` is treated as a scheduler cache and rejected as stale outside an open substep by scheduler validation/import paths.
- `ParticleMigrationRecord::time_bin` remains a mirror carried for diagnostics/compatibility; exact migration continuation requires scheduler identity records and rebuild helpers.
- PM cadence state is separate from particle timestep bins but is scheduler-owned through `PmSynchronizationState`; this review did not find a direct mutable cadence mirror bypass in the Stage 2 scope.

## Reproducibility impact

The blocker patch does not change solver numerics, timestep scheduling, restart schema fields, or normalized config/provenance content. It tightens restart validation and import mirror rebuild behavior so existing v6 payloads with stale particle-bound gas-cell timestep mirrors fail fast instead of being accepted. Deterministic continuation improves because particle scheduler truth is now the only source used to rebuild both particle and gas-cell timestep mirrors after restart import.
