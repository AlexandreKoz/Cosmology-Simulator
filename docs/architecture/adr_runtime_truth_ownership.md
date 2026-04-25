# ADR-INFRA-OWNERSHIP-013: Single-source-of-truth runtime ownership policy

## Status

Accepted (infrastructure repair policy baseline)

Date: 2026-04-25 (UTC)

## Context

Stage 0 runtime-truth mapping identified real runtime domains where authoritative state, mirrors, and mutation paths existed but were not yet formalized as a single ownership policy. The baseline audit is captured in `docs/architecture/runtime_truth_map.md`.

This ADR defines one authoritative owner for each runtime state domain and codifies allowed mirrors/caches/views, refresh/invalidation semantics, restart/reload behavior, reorder/resize/migration rules, and forbidden duplicate-authority patterns.

This ADR is infrastructure policy only; it does not authorize solver behavior changes.

## Decision

### A. Ownership table (single authority per runtime domain)

| Runtime domain | Authoritative owner | Authoritative storage | Primary mutators | Notes |
|---|---|---|---|---|
| Current simulation time and scale factor | Integrator subsystem (`core::IntegratorState`) | `IntegratorState::{current_time_code,current_scale_factor,dt_time_code,step_index}` | `core::StepOrchestrator::executeSingleStep` | Scheduler never writes integrator time scalars directly. |
| Timestep-bin assignment (all scheduled elements) | Scheduler subsystem (`core::HierarchicalTimeBinScheduler`) | `HierarchicalTimeBinScheduler::m_hot.bin_index` + `m_elements_by_bin` + `m_current_tick` | `reset`, `setElementBin`, `requestBinTransition` + `endSubstep` (`applyPendingTransitions`) | `SimulationState::{particles.time_bin,cells.time_bin}` are mirrors only. |
| Active/inactive status | Scheduler subsystem | `HierarchicalTimeBinScheduler::m_hot.active_flag` and `m_active_elements` | `beginSubstep` (`rebuildActiveSet`), `endSubstep` | `ActiveSetDescriptor` is per-step derived data. |
| Canonical particle order | `core::SimulationState` particle SoA+sidecar lanes | `SimulationState::{particles,particle_sidecar}` row order | `reorderParticles`, `commitParticleMigration` | All sidecars must be synchronized through these APIs. |
| Species partition/grouping | `ParticleSidecar::species_tag` with ledger/index maintained by `SimulationState` | `particle_sidecar.species_tag`; `species.count_by_species`; `particle_species_index` | `reorderParticles` (then `rebuildSpeciesIndex`), `commitParticleMigration`, species-changing physics callbacks | `species.count_by_species` and `particle_species_index` are derived/validated mirrors. |
| Gas cell identity and gas-cell row order | Local gas particle identity (`particle_sidecar.particle_id` for gas species) mediated by workflow gas rebuild path | `cells` + `gas_cells` rows keyed by gas particle IDs during migration/rebuild | `rebuildLocalGasStateFromParticleIds` and related migration helpers in `workflows/reference_workflow.cpp` | Gas cell index is not globally stable; gas particle ID is stable identity anchor. |
| Softening (global/species/per-particle) | Frozen typed config for policy defaults; particle sidecar for overrides | Defaults: `SimulationConfig.numerics.gravity_softening_*`; overrides: `particle_sidecar.gravity_softening_comoving` | Config load/normalization; `materializePerParticleSoftening` | Diagnostics/provenance are observers, never authorities. |
| Raw config, normalized config, derived runtime constants, provenance | Config/provenance subsystem (`core::FrozenConfig`, `core::ProvenanceRecord`) | `FrozenConfig::{config,normalized_text,provenance}` + persisted normalized/provenance artifacts | `loadFrozenConfigFromFile/String`, `writeNormalizedConfigSnapshot`, provenance constructors and IO writers | Runtime mutation of typed config is forbidden. |
| Restart continuation truth | Restart schema payload + typed restart read result | `io::RestartReadResult::{state,integrator_state,scheduler_state,provenance,normalized_config_*}` | `writeRestartCheckpointHdf5`, `readRestartCheckpointHdf5` | Restart preserves full scheduler persistent state; snapshot does not. |
| Active-set construction | Scheduler active-elements + workflow split into particle/cell subsets | `scheduler.beginSubstep()` output then `ActiveSetDescriptor` | Reference workflow step loop | Active views are transient gather/scatter workspaces only. |

### B. Mirror/cache/view rules

| Mirror/cache/view | Owner of mirror object | Source of truth | Refresh policy | Invalidation policy | Allowed lifetime |
|---|---|---|---|---|---|
| `SimulationState::particles.time_bin` | `SimulationState` storage lane (mirror role) | `HierarchicalTimeBinScheduler` persistent/hot bin state | Must refresh via `syncTimeBinsFromScheduler(...)` after scheduler mutation cycles (`initializeSchedulerBins`, post-`endSubstep`, post-restart import path) | Invalid immediately after any scheduler bin mutation until sync | Persistent lane, but non-authoritative during stepping |
| `SimulationState::cells.time_bin` | `SimulationState` storage lane (mirror role) | Scheduler bin state (shared index-space convention in reference workflow) | Same as particle bin mirror | Same as particle bin mirror | Persistent lane, non-authoritative during stepping |
| `IntegratorState::time_bins` metadata | Integrator subsystem | Scheduler configuration/runtime context | Refresh at scheduler initialization/reload (`hierarchical_enabled`, `max_bin`, `active_bin`) | Invalid for occupancy/active membership inference at all times | Persistent metadata only |
| `ActiveSetDescriptor` | Orchestrator/workflow | Scheduler `beginSubstep()` result | Rebuild each substep from current active scheduler elements | Invalid after reorder/migration/bin mutation outside its construction point | Single substep |
| `ParticleActiveView`, `CellActiveView`, `GravityParticleKernelView`, `HydroCellKernelView` | `TransientStepWorkspace` | `SimulationState` + current `ActiveSetDescriptor` | Gather on build; push back with scatter functions where mutable views are used | Invalid on any state reorder/resize/migration before scatter; invalid after scatter for reuse | Stage-local/transient only |
| `ParticleSpeciesIndex` | `SimulationState` species index subsystem | `particle_sidecar.species_tag` | Must rebuild via `SimulationState::rebuildSpeciesIndex()` after any change affecting particle order/count/species tags | Invalid immediately after reorder/migration/species-tag edits until rebuild | Persistent derived index |
| `SpeciesContainer::count_by_species` ledger | `SimulationState::species` | `particle_sidecar.species_tag` | Updated by authoritative mutation APIs and validated by invariants | Invalid after unsynchronized direct species-tag edits (forbidden) | Persistent derived ledger |
| Gas migration lookup maps (e.g., `GasCellMigrationRecord` maps) | Reference workflow migration helpers | Gas particle ID + gas-cell row data | Recomputed at each migration/compaction event | Invalid after migration commit/reorder | Event-scoped transient |

### C. Mutation authority table

| State lane | Allowed mutators | Forbidden mutators |
|---|---|---|
| `IntegratorState::{current_time_code,current_scale_factor,step_index}` | `StepOrchestrator::executeSingleStep`; restart import/read assignment | Solver callbacks directly mutating integrator time/step counters |
| Scheduler bin/tick/active state | `HierarchicalTimeBinScheduler` public APIs only (`reset`,`setElementBin`,`requestBinTransition`,`beginSubstep`,`endSubstep`,`importPersistentState`) | Direct writes to scheduler internals; direct writes to `state.*.time_bin` as authority |
| Particle/cell time-bin mirrors | `syncTimeBinsFromScheduler`; restart read into `SimulationState` as load-time artifact | Any module treating mirrors as independent bin authority |
| Particle order and sidecar alignment | `reorderParticles`; `commitParticleMigration` | Ad-hoc partial reorders of subset arrays/sidecars |
| Species index and counts | `rebuildSpeciesIndex` + synchronized species-ledger updates in canonical mutation paths | Manual edits to `particle_species_index`/`count_by_species` detached from `species_tag` truth |
| Gas cell rows/identity mapping | Workflow gas rebuild helpers keyed by gas particle IDs | Assuming pre-migration cell index stability or mutating cells without ID remap |
| Per-particle softening overrides | Controlled writes to `particle_sidecar.gravity_softening_comoving` (materialization/init and explicitly scoped update paths) | Diagnostics/provenance or tree kernels writing override truth |
| Normalized config/provenance text/hash | Config/provenance loading and snapshot/restart IO contracts | Recomputing/rewriting normalized config or hashes in unrelated subsystems |
| Active-set buffers/views | Active-view builders and scatter functions | Caching and reusing stale active views across reorder/bin changes |

Debug/test enforcement helpers:
- `timeBinMirrorsMatchScheduler(...)` reports whether particle/cell mirror lanes still match scheduler authority.
- `debugAssertTimeBinMirrorAuthorityInvariant(...)` fails loudly when a non-owner mutation (or stale mirror) is observed.

### D. Restart/reload rules

1. **Restart checkpoint is continuation authority** for scheduler + integrator + state + normalized config + provenance (`readRestartCheckpointHdf5`).
2. `HierarchicalTimeBinScheduler::importPersistentState` is the only valid path to restore scheduler bin/tick membership exactly.
3. After restart state load, consumers must rebuild/validate derived indices (`rebuildSpeciesIndex`, ownership invariants) before stepping.
4. **Snapshot is not timestep-bin continuation authority**: snapshot read assigns particle `time_bin = 0` on import path, then rebuilds species/index; scheduler continuity must come from restart, not snapshot.
5. Restart compatibility and integrity requirements remain schema/version/hash gated by restart IO contracts.

### E. Reorder / resize / migration rules

1. Canonical particle reorder operations must flow through `buildParticleReorderMap` + `reorderParticles`.
2. Sidecar synchronization policy (`SidecarSyncPolicy`) must be explicit per sidecar lane; default is parent-indirection remap.
3. Any reorder/migration changing particle indices requires `rebuildSpeciesIndex` before species-index consumers execute.
4. `debugAssertNoStaleParticleIndices` should be used in debug/repair paths when reorder/migration risk stale sidecar indices.
5. `SidecarSyncMode::kMoveWithParent` requires full-row movement for each species sidecar payload, not only `particle_index` remapping.
6. Gas cell rows must be reconstructed by gas particle ID mapping in migration/compaction paths (`collectLocalGasCellRecords` / `rebuildLocalGasStateFromParticleIds`), then host-cell references remapped.
7. Resize operations (`resizeParticles`, `resizeCells`) are structural and must be followed by required derived-index/ledger synchronization before runtime stepping.
8. Allowed species migration path is `packParticleMigrationRecords` + `commitParticleMigration`; species-tag edits outside this path are forbidden because sidecars/count ledgers/indexes and per-particle softening overrides must stay synchronized.
9. Species migration that changes sidecar family (e.g., star->gas, BH->gas, gas->star/tracer) must drop obsolete sidecar rows and initialize required destination sidecar rows in the same commit boundary.

### F. Softening override priority and preservation

Priority order for gravity softening values:
1. Per-particle override (`particle_sidecar.gravity_softening_comoving[i]`) when sidecar is populated.
2. Species-specific config default (`SimulationConfig.numerics.gravity_softening_<species>_kpc_comoving`) when positive.
3. Global config default (`SimulationConfig.numerics.gravity_softening_kpc_comoving`).

Preservation rules:
- Reorder/migration must move/remap per-particle overrides with parent particles (`reorderParticles`, migration commit paths).
- Restart must persist and reload override sidecar values.
- Snapshot contracts must not be treated as authoritative for restart continuation-grade override semantics unless explicitly documented and tested.

### G. Config-derived runtime value ownership

1. **Raw source text owner**: input file/string passed to config loader.
2. **Normalized typed owner**: `core::FrozenConfig` (`config` + `normalized_text` + hash/provenance metadata).
3. **Derived runtime constants owner**: consuming subsystem at controlled construction points from `SimulationConfig` (e.g., softening policy materialization and mode policy validation), not ad-hoc recalculation across modules.
4. **Runtime-mutable values**: only dedicated runtime state (`SimulationState`, `IntegratorState`, scheduler state). `SimulationConfig` itself is immutable at runtime.
5. **Provenance owner**: `core::ProvenanceRecord` produced from frozen config/runtime context and serialized by snapshot/restart IO layers.

### H. Active-set construction policy

1. Active-set authority is scheduler output from `beginSubstep()`.
2. `ActiveSetDescriptor` is derived per substep by splitting scheduler element indices into particle/cell subsets in workflow orchestration.
3. Active views are transient gather/scatter workspaces and must not be persisted as authority.
4. Active-set caches may exist only within a substep/stage and must be invalidated after:
   - scheduler bin mutation (`requestBinTransition` + `endSubstep`),
   - any reorder/resize/migration affecting index spaces,
   - restart/reload of state/scheduler.
5. Mutable compact kernel views must carry captured index-space generation and fail scatter when generations mismatch.
6. Consumers: stage callbacks via `StepContext.active_set` and builders in `simulation_state_active_views.cpp`.
7. Active eligibility is scheduler/bin-driven and species-agnostic by contract; species migration alone does not authorize ad-hoc active-set mutations outside scheduler ownership APIs.
8. Forbidden: competing active-set builders that bypass scheduler authority for the same step.

## Forbidden duplicate authority patterns

The following patterns are explicitly forbidden:

1. Solver modules mutating scheduler bin truth indirectly by writing `state.particles.time_bin` or `state.cells.time_bin` and treating it as authoritative.
2. Any path where both scheduler internals and particle/cell arrays are treated as concurrent authoritative timestep owners.
3. Sidecar index mirrors surviving reorder/migration without remap/rebuild and stale-index validation.
4. Diagnostics/provenance payload fields becoming runtime softening truth.
5. Reusing active-set caches/views after bin mutation, reorder, resize, migration, or restart without rebuild.
6. Recomputing normalized/derived config values inconsistently in multiple subsystems outside typed config/provenance contracts.
7. Persisting/reloading active set lists as continuation truth independent of scheduler persistent state.

## Test obligations for follow-up prompts

Future repair prompts that touch these domains must include targeted tests (or cite exact existing coverage) for:

1. Scheduler truth vs state mirror synchronization invariants across full substep loops.
2. Reorder/migration stale-index prevention and species-index rebuild obligations.
3. Gas cell ID-based reconstruction correctness under migration/compaction.
4. Softening priority-order correctness and preservation across reorder/restart.
5. Restart-vs-snapshot timestep-bin contract boundaries.
6. Active-set cache invalidation after scheduler bin changes and reorder events.

## Migration notes for current ambiguous areas

1. **Scheduler vs state bin ambiguity** (runtime truth map): policy now formalizes scheduler as sole authority and state bins as mirrors only.
2. **Species multi-lane duplication** (`species_tag`, `count_by_species`, `particle_species_index`): policy now formalizes `species_tag` as root truth and other lanes as derived mirrors requiring explicit rebuild/sync.
3. **Gas cell identity under migration**: policy now formalizes gas particle ID as stable identity anchor and cell index as rebuildable local index.
4. **Snapshot timestep-bin continuity ambiguity**: policy now formalizes that snapshot import is non-authoritative for timestep-bin continuation; restart is required for exact continuation.

## Reproducibility impact

This ADR is documentation/policy only and does not change runtime behavior, schema payloads, or solver numerics. It strengthens deterministic continuation expectations by making ownership and invalidation contracts explicit.
