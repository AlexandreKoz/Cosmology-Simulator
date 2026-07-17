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
| Gas cell identity and gas-cell row order | `core::SimulationState::gas_cell_identity` | `GasCellIdentityMap` records keyed by stable nonzero `gas_cell_id`, with transient dense `local_cell_row` and generation | `replaceGasCellIdentityRecords`, `commitGasCellMigration`, `commitAmrPatchMigration`, restart import, legacy particle-bound refresh helpers | **H2 distributed contract:** production workflow compaction/rebalance migrates gas cells and AMR patches by stable `gas_cell_id`/`patch_id`. `parent_particle_id` is optional compatibility/provenance only and is not hydro truth. |
| Softening (global/species/per-particle) | Frozen typed config for policy defaults; particle sidecar for overrides | Defaults: `SimulationConfig.numerics.gravity_softening_*`; overrides: `particle_sidecar.gravity_softening_comoving` gated by `particle_sidecar.has_gravity_softening_override` | Config load/normalization; `materializePerParticleSoftening` | Diagnostics/provenance are observers, never authorities. |
| Raw config, normalized config, derived runtime constants, provenance | Config/provenance subsystem (`core::FrozenConfig`, `core::ProvenanceRecord`) | `FrozenConfig::{config,normalized_text,provenance}` + persisted normalized/provenance artifacts | `loadFrozenConfigFromFile/String`, `writeNormalizedConfigSnapshot`, provenance constructors and IO writers | Runtime mutation of typed config is forbidden. |
| Gravity Newton constant in code units | Frozen `core::UnitSystem` conversion | Result of `core::newtonGravitationalConstantCode(UnitSystem)` supplied identically to PM/tree options | Workflow gravity construction/restart reconstruction from frozen units | Physical-unit workflows must not hard-code `G_code=1`; option copies are derived values, not independent authorities. |
| Cosmological gravity response | Integrator cosmological timeline | In-flight `StepContext.timeline_step`; committed `IntegratorState::{current_scale_factor,current_hubble_rate_code}` | `CosmologicalTimeline::prepareStep` / stage dispatch / step commit | PM and TreePM return scale-free `A`; collisionless KDK and the gas conservative source apply the time response to their respective state. Post-drift hydro uses the in-flight step-end epoch, not the still-step-begin committed state. PM field-build scale is validity metadata, not a kernel multiplier. |
| Particle decomposition generation | Gravity workflow | `GravityStageCallback::m_decomposition_epoch` | Advance after a globally committed ownership transition; restore from distributed restart state | Tree wire protocols consume this epoch. Rank-local particle-index generation is not a distributed epoch. |
| Long-range PM field validity | `TreePmCoordinator` | `LongRangeFieldValidity` plus transient PM mesh/halo state | Refresh at an integrator force surface; explicit reuse only after full compatibility validation | Particle decomposition does not alone invalidate fixed FFT-slab field ownership; force epoch and operator inputs do. |
| In-flight PM refresh decision | Integrator-issued `PmRefreshDirective`; committed cadence remains `PmSynchronizationState` | Pending opportunity/version/build-step/build-scale values until commit | Stage dispatch creates directive; gravity solve observes it; cadence commit updates `PmSynchronizationState` | `gravity.pm_long_range_field` emitted before commit must copy the pending directive/decision, not stale committed fields; the event is diagnostic only. |
| Committed gravity force history | Gravity workflow force cache | Dense acceleration arrays, per-row validity, stable-ID restart payload | Refresh at force evaluation; invalidate on dense-row ownership/index changes | Relative MAC may observe only compatible committed history. |
| Gravity timestep acceleration | Adaptive scheduler criteria after step commit | Transient public `ComovingGravityTimeStepInput` containing scale-free `A`, comoving softening, and committed `IntegratorState.current_scale_factor` | `updateAdaptiveTimeBins` via `computeComovingGravityTimeStep` | With comoving length, the compatible coordinate acceleration is `A/a^3`; the helper validates finite positive `a` and evaluates `eta sqrt(a^3 epsilon_com/|A|)`. The generic `computeGravityTimeStep` remains coordinate-neutral. |
| Restart continuation truth | Restart schema payload + typed restart read result | `io::RestartReadResult::{state,integrator_state,scheduler_state,provenance,normalized_config_*}` | `writeRestartCheckpointHdf5`, `readRestartCheckpointHdf5` | Restart preserves full scheduler persistent state; snapshot does not. |
| Active-set construction | Scheduler active-elements + workflow split into particle/cell subsets | `scheduler.beginSubstep()` output then `ActiveSetDescriptor` | Reference workflow step loop | Active views are transient gather/scatter workspaces only. |

### B. Mirror/cache/view rules

| Mirror/cache/view | Owner of mirror object | Source of truth | Refresh policy | Invalidation policy | Allowed lifetime |
|---|---|---|---|---|---|
| `SimulationState::particles.time_bin` | `SimulationState` storage lane (mirror role) | `HierarchicalTimeBinScheduler` persistent/hot bin state | Must refresh via `syncTimeBinMirrorsFromScheduler(...)` after scheduler mutation cycles (`initializeSchedulerBins`, post-`endSubstep`, post-restart import path) | Invalid immediately after any scheduler bin mutation until sync | Persistent lane, but non-authoritative during stepping |
| `SimulationState::cells.time_bin` | `SimulationState` storage lane (mirror role) | Scheduler bin state (shared index-space convention in reference workflow) | Same as particle bin mirror | Same as particle bin mirror | Persistent lane, non-authoritative during stepping |
| `IntegratorState::time_bins` metadata | Integrator subsystem | Scheduler configuration/runtime context | Refresh at scheduler initialization/reload (`hierarchical_enabled`, `max_bin`, `active_bin`) | Invalid for occupancy/active membership inference at all times | Persistent metadata only |
| `ActiveSetDescriptor` | Orchestrator/workflow | Scheduler `beginSubstep()` result | Rebuild each substep from current active scheduler elements | Invalid after reorder/migration/bin mutation outside its construction point | Single substep |
| `ParticleActiveView`, `CellActiveView`, `GravityParticleKernelView`, `HydroCellKernelView` | `TransientStepWorkspace` | `SimulationState` + current `ActiveSetDescriptor` | Gather on build; push back with scatter functions where mutable views are used | Invalid on any state reorder/resize/migration before scatter; invalid after scatter for reuse | Stage-local/transient only |
| `ParticleSpeciesIndex` | `SimulationState` species index subsystem | `particle_sidecar.species_tag` | Must rebuild via `SimulationState::rebuildSpeciesIndex()` after any change affecting particle order/count/species tags | Invalid immediately after reorder/migration/species-tag edits until rebuild | Persistent derived index |
| `SpeciesContainer::count_by_species` ledger | `SimulationState::species` | `particle_sidecar.species_tag` | Updated by authoritative mutation APIs and validated by invariants | Invalid after unsynchronized direct species-tag edits (forbidden) | Persistent derived ledger |
| Gas migration payloads (`GasCellMigrationRecord`, `AmrPatchMigrationRecord`) | State/workflow migration transaction | `SimulationState::gas_cell_identity`, `PatchSoa`, gas hydro lanes, and gas scheduler identity records | Packed at each migration/compaction event by stable `gas_cell_id`/`patch_id` | Invalid after migration commit/reorder or gas identity-map generation change | Event-scoped transient |

Gas identity field classes:
- **Stable identity field:** `GasCellIdentityMap::gas_cell_id`, nonzero and unique. It is independent of parent particle identity; equality with a gas particle ID is allowed only as legacy compatibility data.
- **Compatibility mirror fields:** `GasCellSidecar::{gas_cell_id,parent_particle_id}` mirror the map for legacy restart/import and particle-bound callers. `parent_particle_id == 0` mirrors a parentless identity record; parent lineage is not a uniqueness authority.
- **Persistent owner-mutable hydro sidecar fields:** `velocity_[xyz]_peculiar`, `density_code`, `pressure_code`, `internal_energy_code`, `temperature_code`, `sound_speed_code`.
- **Transient reconstruction scratch (not identity/restart truth):** reconstruction gradients live in hydro/transient workspaces, not `GasCellSidecar`.
- **Scratch/transient hydro extraction fields:** `HydroCellKernelView` gathered active-cell buffers in `TransientStepWorkspace`; these are not persistence or identity authority.

### C. Mutation authority table

| State lane | Allowed mutators | Forbidden mutators |
|---|---|---|
| `IntegratorState::{current_time_code,current_scale_factor,step_index}` | `StepOrchestrator::executeSingleStep`; restart import/read assignment | Solver callbacks directly mutating integrator time/step counters |
| Scheduler bin/tick/active state | `HierarchicalTimeBinScheduler` public APIs only (`reset`,`setElementBin`,`requestBinTransition`,`beginSubstep`,`endSubstep`,`importPersistentState`) | Direct writes to scheduler internals; direct writes to `state.*.time_bin` as authority |
| Particle/cell time-bin mirrors | `syncTimeBinMirrorsFromScheduler`; restart read into `SimulationState` as load-time artifact | Any module treating mirrors as independent bin authority |
| Particle order and sidecar alignment | `reorderParticles`; `commitParticleMigration` | Ad-hoc partial reorders of subset arrays/sidecars |
| Species index and counts | `rebuildSpeciesIndex` + synchronized species-ledger updates in canonical mutation paths | Manual edits to `particle_species_index`/`count_by_species` detached from `species_tag` truth |
| Gas cell rows/identity mapping | `SimulationState::gas_cell_identity` materialization/remap helpers; legacy workflow gas rebuild helpers must keep the map and compatibility mirrors synchronized | Assuming pre-migration cell index stability, mutating cells without ID remap, or letting compatibility mirror lanes drift from the map |
| Per-particle softening overrides | Controlled writes to `particle_sidecar.setGravitySofteningOverride` / `clearGravitySofteningOverride`; materialized default values may populate `gravity_softening_comoving` only when the override mask is false | Diagnostics/provenance or tree kernels writing override truth |
| Normalized config/provenance text/hash | Config/provenance loading and snapshot/restart IO contracts | Recomputing/rewriting normalized config or hashes in unrelated subsystems |
| Active-set buffers/views | Active-view builders and scatter functions | Caching and reusing stale active views across reorder/bin changes |

Debug/test enforcement helpers:
- `timeBinMirrorsMatchScheduler(...)` reports whether particle/cell mirror lanes still match scheduler authority.
- `debugAssertTimeBinMirrorAuthorityInvariant(...)` fails loudly when a non-owner mutation (or stale mirror) is observed.

### D. Restart/reload rules

1. **Restart checkpoint is continuation authority** for scheduler + integrator + state + normalized config + provenance (`readRestartCheckpointHdf5`).
2. `HierarchicalTimeBinScheduler::importPersistentState` is the only valid path to restore scheduler bin/tick membership exactly.
3. After restart state load, consumers must rebuild/validate derived indices (`rebuildSpeciesIndex`, ownership invariants) before stepping.
4. Resume-time active sets are reconstructed from scheduler persistent state (`importPersistentState` + scheduler active-set rebuild), not from serialized active-view caches.
5. Missing required restart continuation fields must fail loudly; optional legacy compatibility is allowed only when explicitly documented and tested (current example: absent per-particle softening override dataset is treated as no overrides present).
6. **Snapshot is not timestep-bin continuation authority**: snapshot read assigns particle `time_bin = 0` on import path, then rebuilds species/index; scheduler continuity must come from restart, not snapshot.
7. Restart compatibility and integrity requirements remain schema/version/hash gated by restart IO contracts.

### E. Reorder / resize / migration rules

1. Canonical particle reorder operations must flow through `buildParticleReorderMap` + `reorderParticles`.
2. Sidecar synchronization policy (`SidecarSyncPolicy`) must be explicit per sidecar lane; default is parent-indirection remap.
3. Any reorder/migration changing particle indices requires `rebuildSpeciesIndex` before species-index consumers execute.
4. `debugAssertNoStaleParticleIndices` should be used in debug/repair paths when reorder/migration risk stale sidecar indices.
5. `SidecarSyncMode::kMoveWithParent` requires full-row movement for each species sidecar payload, not only `particle_index` remapping; after row movement, `particle_index` is remapped exactly once through `old_to_new_index`.
6. The default `SidecarSyncMode::kUseParentIndirection` keeps sidecar row order stable and remaps only `particle_index`; payload identity remains keyed by the parent particle ID, not by pre-reorder row number.
7. Reorder implementations must use typed sidecar lane visitors so scalar and array-channel lanes move together. Future module sidecars may register only an explicit row-move/remap contract; ad-hoc edits to one metadata lane are forbidden.
8. `debugAssertSpeciesSidecarOwnershipInvariants(...)` is the defensive post-reorder guard for exactly one eligible sidecar row and zero ineligible rows.
9. Gas cell rows must be reconstructed through the authoritative `SimulationState::gas_cell_identity` mapping. Production workflow compaction/rebalance uses explicit gas-cell/AMR patch migration payloads keyed by `gas_cell_id` and `patch_id`; particle migration commits on those paths set `preserve_gas_cell_state` so optional parent-particle removal cannot delete hydro truth. Code that needs a legacy row/ID association must use the named helpers (`parentParticleIdForGasCellRow`, `gasCellRowForParticleId`, or `gasParticleIndexForCellRow`) after `requireParticleBoundGasCellContract(...)` has passed.
10. Temporary safety guard: `reorderParticles` must fail loudly if it would change the relative gas-particle order while gas cells exist (unless a gas-cell ID-based rebuild follows in the same repair path). Hydro active-view construction and scatter also require this contract before solver kernels can consume cell rows.
11. Resize operations (`resizeParticles`, `resizeCells`) are structural and must be followed by required derived-index/ledger synchronization before runtime stepping.
12. Allowed species migration path is `packParticleMigrationRecords` + `commitParticleMigration`; species-tag edits outside this path are forbidden because sidecars/count ledgers/indexes and per-particle softening overrides must stay synchronized.
13. Species migration that changes sidecar family (e.g., star->gas, BH->gas, gas->star/tracer) must drop obsolete sidecar rows and initialize required destination sidecar rows in the same commit boundary.

### F. Softening override priority and preservation

Priority order for gravity softening values:
1. Per-particle override only when `particle_sidecar.has_gravity_softening_override[i] != 0`; `gravity_softening_comoving[i]` without the mask is a materialized default/diagnostic value, not override authority.
2. Species-specific config default (`SimulationConfig.numerics.gravity_softening_<species>_kpc_comoving`) when positive.
3. Global config default (`SimulationConfig.numerics.gravity_softening_kpc_comoving`).

Preservation rules:
- Reorder/migration must move/remap per-particle overrides with parent particles (`reorderParticles`, migration commit paths).
- Resize grow/shrink must preserve retained-row override values when the optional sidecar lane is populated.
- Active-set extraction lanes are mirrors only: gathered active softening slices must reflect sidecar truth and may never write back as an independent authority.
- Restart must persist and reload override sidecar values.
- Snapshot contracts must not be treated as authoritative for restart continuation-grade override semantics unless explicitly documented and tested.
- Diagnostics/provenance softening fields (`gravity_softening_policy`, `gravity_softening_kernel`, `gravity_softening_epsilon_kpc_comoving`) are descriptive mirrors and are never runtime authority.

### G. Config-derived runtime value ownership

1. **Raw source text owner**: input file/string passed to config loader.
2. **Normalized typed owner**: `core::FrozenConfig` (`config` + `normalized_text` + hash/provenance metadata).
3. **Derived runtime constants owner**: consuming subsystem at controlled construction points from `SimulationConfig` (e.g., softening policy materialization and mode policy validation), not ad-hoc recalculation across modules.
4. **Runtime-mutable values**: only dedicated runtime state (`SimulationState`, `IntegratorState`, scheduler state). `SimulationConfig` itself is immutable at runtime.
5. **Provenance owner**: `core::ProvenanceRecord` produced from frozen config/runtime context and serialized by snapshot/restart IO layers.

Ownership map for config-derived values used in runtime pathways:

| Value lane | Owner | Examples | Mutation policy |
|---|---|---|---|
| Raw config | Config parser input | `param.txt` key/value lines, legacy aliases (e.g., `mode`, `box_size`, `treepm_pm_grid`) | Parse-time only; never consumed directly by solver modules. |
| Normalized config | `core::FrozenConfig` | `FrozenConfig::config` typed enums/numerics, `FrozenConfig::normalized_text`, `FrozenConfig::provenance.config_hash_hex` | Immutable after freeze; canonical source for runtime construction and reproducible hashing. |
| Derived runtime constants | Runtime constructors + typed helper seams | TreePM mesh spacing/split/cutoff from `treepm_*` + box/grid axes, species softening policy arrays, `LambdaCdmBackground` and `UnitSystem` | Must be computed from normalized typed config through a single owning construction path per subsystem; no parallel ad-hoc recomputation with divergent formulas. |
| Runtime-mutable state | Runtime state owners | `IntegratorState::{current_time_code,current_scale_factor,step_index}`, scheduler bins/ticks, particle/cell arrays and sidecars | Mutated only by runtime evolution/restart import APIs; not fed back into config authority. |
| Diagnostic/provenance mirrors | `core::ProvenanceRecord` + IO payload wrappers | `gravity_treepm_mesh_spacing_*`, `gravity_treepm_split_scale_mpc_comoving`, normalized config hash mirror fields in snapshot/restart payloads | Observer-only for audit/continuation validation; never promoted to solver authority. |

Ambiguous legacy-name policy (must remain explicit in code/docs/tests):

- `numerics.t_code_begin` and `numerics.t_code_end` are **code-time domain scalars**, not redshift, scale factor, or physical seconds; legacy `time_begin_code`/`time_end_code` names are accepted only at the UserConfig input boundary.
- `IntegratorState.current_scale_factor` is the cosmological scale factor `a`; continuation/restart must carry this runtime lane explicitly.
- Redshift is diagnostic-only derived metadata (`z = 1/a - 1` when `a>0`), not a persisted runtime authority lane.
- In non-cosmological modes, scale factor fallback is treated as a runtime compatibility lane and must not be interpreted as physical-time truth.

### H. Active-set construction policy

1. Active-set authority is scheduler output from `beginSubstep()`.
2. `ActiveSetDescriptor` is derived per substep by splitting scheduler element indices into particle/cell subsets in workflow orchestration.
3. Solver-local active views (`GravityParticleKernelView`, `HydroCellKernelView`, and read-only active views) must be derived only from the current authoritative `ActiveSetDescriptor` indices for that same substep.
4. Active views are transient gather/scatter workspaces and must not be persisted as authority.
5. No module-local fallback builder is allowed to infer active membership directly from `state.*.time_bin` mirrors during stepping; mirror lanes are debug/IO mirrors only.
6. Scheduler-driven active extraction and solver-local derived views are a two-layer contract:
   - Layer A (authority): `HierarchicalTimeBinScheduler::beginSubstep()`.
   - Layer B (derived): workflow split -> `ActiveSetDescriptor` -> active view builders.
   Any alternate authority lane is forbidden.
7. Active-set caches/views must carry explicit source generation fields for mutable scatter paths; pointer-keyed global registries are forbidden because view lifetime must be local and auditable.
   - Particle compact mutable views capture `SimulationState::particleIndexGeneration()`.
   - Cell compact mutable views capture `SimulationState::cellIndexGeneration()` and `SimulationState::gasCellIdentityGeneration()`.
   Scatter must fail loudly on either generation mismatch.
8. Allowed invalidation events for active views/caches:
   - scheduler bin mutation (`requestBinTransition` + `endSubstep`),
   - particle/cell reorder,
   - particle/cell resize,
   - migration commit that rewires index spaces,
   - restart/reload that imports scheduler/state truth.
9. Rebuild policy:
   - rebuild authoritative active indices at each `beginSubstep()`,
   - rebuild solver-local compact views when entering the stage that uses them,
   - never reuse mutable compact views across structural mutations.
10. Lifetime:
   - authoritative active indices: one scheduler substep,
   - `ActiveSetDescriptor`: one orchestrator step call,
   - compact solver views: one callback stage region from build to scatter.
11. Consumers: stage callbacks via `StepContext.active_set` and builders in `simulation_state_active_views.cpp`.
12. Active eligibility is scheduler/bin-driven and species-agnostic by contract; species migration alone does not authorize ad-hoc active-set mutations outside scheduler ownership APIs.
13. Forbidden: competing active-set builders that bypass scheduler authority for the same step.

### I. Production gravity time, cache, and decomposition policy

1. Production `ReferenceWorkflow` requires `hierarchical_max_rung=0` until
   per-element kick and drift epochs are represented and restartable. Both
   particle and gas-cell scheduler restart payloads must remain bin-zero.
2. Production cadence one refreshes PM at every integrator-issued,
   rank-coordinated force-refresh surface. Lower-level local-boundary and PM
   reuse APIs are validation/future-integration seams, not production multirate
   authority.
3. An explicit PM reuse request must fail if any rank lacks a compatible field;
   it must not silently mutate the integrator's requested operation into a
   refresh.
4. The PM compatibility key includes force epoch, field-build scale,
   `G_code`, split/box geometry, assignment, boundary, decomposition mode, and
   deconvolution. It excludes the particle decomposition epoch because the PM
   mesh remains owned by fixed FFT slabs and interpolation routes current
   targets to those owners.
5. Tree hierarchy/request/response packets carry the actual workflow
   decomposition epoch. That epoch advances only on an actual committed
   particle-ownership transition and is restart truth.
6. Dense acceleration rows and their relative-MAC history are invalidated
   immediately on ownership/index change even when the PM mesh field remains
   compatible. Restart written at that boundary must persist an invalid cache,
   never stale row data.
7. Physical Newton `G` is converted from frozen units once. Neither PM nor the
   complementary tree inserts `a^2`. Collisionless KDK and gas
   `ComovingGravityExpansionSource` apply the `A/a^2` response to their own
   state.
8. During the post-drift hydro stage,
   `StepContext.timeline_step.{scale_factor_end,hubble_end_code}` is the
   source-epoch authority for both fixed-grid and AMR callbacks;
   `IntegratorState` is not committed until the step ends. Recomputing an
   SI-valued background rate inside a solver callback is forbidden.
9. Post-step gravity timestep criteria use the now-committed scale factor and
   call `computeComovingGravityTimeStep` with scale-free `A`, comoving
   softening, and finite positive `a`. That public helper converts to
   comoving-coordinate acceleration `A/a^3` before using the generic
   coordinate-neutral criterion. `A/a^2` is the peculiar-velocity response and
   is not dimensionally compatible with that length criterion.
10. PM and TreePM MPI wire formats are explicit versioned byte contracts.
   Native struct representation, padding, and local endianness are not wire
   authority.

## Forbidden duplicate authority patterns

The following patterns are explicitly forbidden:

1. Solver modules mutating scheduler bin truth indirectly by writing `state.particles.time_bin` or `state.cells.time_bin` and treating it as authoritative.
2. Any path where both scheduler internals and particle/cell arrays are treated as concurrent authoritative timestep owners.
3. Sidecar index mirrors surviving reorder/migration without remap/rebuild and stale-index validation.
4. Diagnostics/provenance payload fields becoming runtime softening truth.
5. Reusing active-set caches/views after bin mutation, reorder, resize, migration, or restart without rebuild.
6. Recomputing normalized/derived config values inconsistently in multiple subsystems outside typed config/provenance contracts.
7. Persisting/reloading active set lists as continuation truth independent of scheduler persistent state.
8. Treating hydro reconstruction scratch buffers or active-kernel compact views as persistent gas-cell identity.
9. Assuming gas-cell identity is `cell_index` alone across reorder/resize/migration/restart boundaries.
10. Using naked `gas particle count == cell count` checks in workflow or solver code instead of `requireParticleBoundGasCellContract(...)` with a caller-specific error prefix.

### Campaign A workflow owner clarification (2026-07-16)

- `workflows::RungZeroTimeState` is the live owner of both scheduler objects,
  `IntegratorState`, and pending output cadence. `ReferenceWorkflowRunner`
  constructs this owner but does not mutate its scheduler internals.
- `workflows::RuntimeModuleRegistry` is the live composition authority for
  production stage tasks. The frozen `RuntimeExecutionPlan`, not a callback
  list or the core layer catalog, determines which typed tasks execute.
- `RuntimeResourceLease` is transient validation metadata. Its captured
  generations/tick/step are never restart truth and must be rebuilt after
  migration, compaction, scheduler advance, or step commit.
- `GravityRuntime`, `HydroAmrRuntime`, `SourceRuntime`, `AnalysisRuntime`, and
  `OutputRestartRuntime` are stage owners reached only through their typed view
  tasks. The numerical modules remain dependency-downstream of their owners
  and do not include workflow composition.

## Test obligations for follow-up prompts

Future repair prompts that touch these domains must include targeted tests (or cite exact existing coverage) for:

1. Scheduler truth vs state mirror synchronization invariants across full substep loops.
2. Reorder/migration stale-index prevention and species-index rebuild obligations.
3. Gas cell ID-based reconstruction correctness under migration/compaction, including identity-map generation invalidation and compatibility-mirror drift rejection.
4. Softening priority-order correctness and preservation across reorder/restart.
5. Restart-vs-snapshot timestep-bin contract boundaries.
6. Active-set cache invalidation after scheduler bin changes and reorder events.

## Migration notes for current ambiguous areas

1. **Scheduler vs state bin ambiguity** (runtime truth map): policy now formalizes scheduler as sole authority and state bins as mirrors only.
2. **Species multi-lane duplication** (`species_tag`, `count_by_species`, `particle_species_index`): policy now formalizes `species_tag` as root truth and other lanes as derived mirrors requiring explicit rebuild/sync.
3. **Gas cell identity under migration**: policy now formalizes `SimulationState::gas_cell_identity` as the in-memory stable gas-cell identity authority. Distributed workflow compaction/rebalance no longer materializes gas state from a one-cell-per-gas-particle assumption; legacy particle-bound helpers remain compatibility-only seams for old tests/import adapters.
4. **Snapshot timestep-bin continuity ambiguity**: policy now formalizes that snapshot import is non-authoritative for timestep-bin continuation; restart is required for exact continuation.

## Reproducibility impact

This ADR is documentation/policy only and does not change runtime behavior, schema payloads, or solver numerics. It strengthens deterministic continuation expectations by making ownership and invalidation contracts explicit.

### Stage 0 P0-04 migration softening payload rule

`ParticleMigrationRecord::has_gravity_softening_value` records whether the migration payload carries a numeric softening value at all. `ParticleMigrationRecord::has_gravity_softening_override` records whether that numeric value is authoritative per-particle override truth. This distinction is required so cross-rank or local species migration can preserve materialized species/default values without promoting them into true overrides.

## Stage 0 P0-05..P0-08 repair addendum

The P0-05..P0-08 repair pass strengthened the ownership contract in four places.

- Gas-cell identity is now represented explicitly by `GasCellSidecar::gas_cell_id` and `GasCellSidecar::parent_particle_id` for particle-bound gas-cell layouts. These lanes are persistent runtime truth for restartable particle-bound gas cells and are refreshed from canonical gas-particle ordering only by the state ownership layer. The named contract check is `requireParticleBoundGasCellContract(...)`; the allowed row/ID seams are `parentParticleIdForGasCellRow`, `gasCellRowForParticleId`, and `gasParticleIndexForCellRow` so future AMR/moving-mesh decoupling can replace the mapping without changing general mutation paths.
- Per-particle softening override authority is no longer inferred from the existence of a numeric softening value. `has_gravity_softening_override` is the explicit override mask; unmasked numeric values are materialized cache/default values and fall back to species/global softening in gravity softening resolution.
- Config-derived runtime values are owned by `DerivedRuntimeConfig`, produced from `FrozenConfig` by `deriveRuntimeConfig()`. Runtime modules must not independently reinterpret raw config keys when a derived runtime value exists.
- Active-set descriptors are derived scheduler products with source generation metadata. A descriptor is invalid after particle/cell generation changes and must be rebuilt from scheduler-owned activity rather than repaired locally by a solver.

### Stage 2 structural-transform scheduler remap contract

- Scheduler element rows are local index-space rows, not physical identities. Any particle reorder, compaction, migration commit, gas-cell rebuild, or restart import that changes row spaces must either remap scheduler persistent lanes by stable identity in the same commit boundary or leave the scheduler unusable until an explicit identity-record rebuild/import occurs.
- Stable scheduler remap keys are `particle_sidecar.particle_id` for particle scheduler rows and `gas_cells.gas_cell_id` / `gas_cells.parent_particle_id` for gas-cell mirror validation. `state.particles.time_bin`, `state.cells.time_bin`, `ParticleMigrationRecord::time_bin`, and `GasCellMigrationRecord::time_bin` are diagnostic/transfer mirrors only and are never sufficient authority for exact continuation.
- Public scheduler remap helpers under `include/cosmosim/core/time_integration.hpp` provide the supported migration path: export full identity records, remap persistent state by particle ID, rebuild scheduler persistent state from identity records, then refresh particle and gas-cell mirrors from scheduler authority.
- Gas-cell time-bin mirror validation in production resolves cells through dense `GasCellIdentityMap` local rows and
  does not require local gas-particle count to equal local gas-cell count. Legacy particle-bound restart/import
  adapters may still resolve `parent_particle_id -> gas particle row` before handing state to production paths.
- MPI H2 contract: workflow AMR patch migration sends gas-cell scheduler identity records alongside the `gas_cell_id` payload and rebuilds the gas scheduler after the validated commit. Rank-count-changing restart and broader spawned-particle scheduler registration remain separate follow-ups; `ParticleMigrationRecord::time_bin` alone remains non-authoritative.
