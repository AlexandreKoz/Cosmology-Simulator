# Time Integration Framework

This document defines the baseline time-integration contract for CosmoSim and records conservative assumptions made in the first implementation.

## Stage contract

The authoritative per-step stage order is:

1. `gravity_kick_pre`
2. `drift`
3. `force_refresh`
4. `hydro_update`
5. `source_terms`
6. `gravity_kick_post`
7. `analysis_hooks`
8. `output_check`

The scheduler (`StageScheduler`) owns this ordering and solver modules interact via stage-bound `IntegrationCallback` implementations. Each callback must declare its exact `integrationStages()` set during registration; `StepOrchestrator` stores handlers in per-stage buckets and dispatches only the bucket for the active stage in stable registration order. Production callbacks must treat any off-stage direct invocation as a contract violation rather than silently self-filtering. This keeps gravity/hydro/source/output sequencing explicit and auditable, prevents solvers from receiving irrelevant stages, and avoids per-stage broadcast scans in the hot path.

### Callback migration note

Any new `IntegrationCallback` implementation must provide:

- `callbackName()` for profiling and diagnostics;
- `integrationStages()` returning the small typed stage set the handler is allowed to receive;
- `onStage(...)` logic that assumes the declared stage contract and throws/asserts on impossible off-stage direct calls instead of returning silently.

Existing code that previously registered a broad callback and checked `context.stage` internally must split the behavior into explicit stage declarations or return the exact small stage set needed by the handler. No config keys, snapshot/restart payloads, or solver numerics migrate for this interface repair.

### TreePM long-range PM kick operator contract (distributed TreePM runtime)

The reference workflow now applies an explicit, config-driven long-range PM kick operator inside KDK:

- operator symbol: `K_PM^LR(surface, opportunity, version)`
- synchronization surfaces: `gravity_kick_pre` and `gravity_kick_post` only
- operator action per kick surface:
  1) decide refresh/reuse for PM long-range field (`refresh_long_range_field`);
  2) evaluate PM long-range acceleration on the **active particle set only** using either refreshed or cached field;
  3) add PM+tree residual acceleration to active particles and apply half-kick update.

This operator is explicit in workflow runtime events and cadence records (`pm_sync_surface`, `gravity_kick_opportunity`, `field_version`, `last_refresh_opportunity`, `active_particles_kicked`, `inactive_particles_skipped`).

#### Refresh and reuse rule

- cadence control: `numerics.treepm_update_cadence_steps` (integer `>= 1`)
- cadence unit: **gravity kick opportunities** (`gravity_kick_pre` and `gravity_kick_post`)
- policy:
  - refresh PM mesh solve when no cached field exists, or when
    `current_kick_opportunity - last_refresh_opportunity >= cadence_steps`
  - otherwise reuse the cached long-range PM field and only evaluate interpolation + tree residual for active targets
- rank-consensus requirements in distributed runs:
  - each rank increments the kick opportunity on every gravity kick stage, even when the local active target set is empty;
  - refresh/reuse is a rank-consensus decision (all ranks must either refresh or reuse on the same opportunity);
  - cadence metadata (`gravity_kick_opportunity`, `field_version`, `last_refresh_opportunity`, refresh flag) is expected to remain coherent across ranks.

The temporal contract for reuse is explicit: a reused PM field corresponds to the particle state and scale factor from the last refresh opportunity, and that metadata is recorded by the workflow report/events.
With `treepm_update_cadence_steps = 1` (default), behavior reduces to immediate PM refresh at every kick opportunity.

#### KDK operator order with PM sync

For each global step:

1. `K_PM+Tree^pre` on active set at sync surface `kick_pre`
2. `D` drift on active set
3. `force_refresh` PM cadence/field refresh surface
4. hydro/source stages
5. `K_PM+Tree^post` on active set at sync surface `kick_post`
6. analysis/output stages

Both kick surfaces use the same PM sync contract above.

#### Inactive-particle treatment

- Inactive particles do **not** receive PM or tree kick updates on a given kick opportunity.
- Their state evolves only when their bin activates on a later synchronized kick surface.
- This is recorded per event as `inactive_particles_skipped`.

#### Restart continuation rule

- Restart payload persists cadence metadata (`gravity_kick_opportunity`, cadence steps, field version, last refresh opportunity, field-built step/scale factor).
- Restart policy is `deterministic_rebuild`: PM mesh field values are not checkpointed; after resume, the first kick opportunity refreshes/rebuilds the long-range PM field before reuse can occur.
- Cadence metadata remains auditable across restart boundaries via distributed restart state and provenance fields.


## Stage 2 timestep authority contract

Stage 2 uses a single-owner timestep model: `HierarchicalTimeBinScheduler` is the only live authority for per-element bin assignment, next activation, active flags, pending transitions, active-set construction, and PM kick cadence metadata. Solver callbacks may propose timestep candidates and consume scheduler-built active sets, but they must not treat `ParticleSoa::time_bin`, `CellSoa::time_bin`, migration records, or restart mirrors as authority.

### Scheduler authority and mirror policy

- Authoritative live lanes: `HierarchicalTimeBinScheduler` hot metadata (`bin_index`, `next_activation_tick`, `active_flag`, `pending_bin_index`) and `PmSynchronizationState` cadence fields (`gravity_kick_opportunity`, cadence steps, field version, last refresh opportunity, field-built step/scale factor).
- Derived mirrors: `ParticleSoa::time_bin`, `CellSoa::time_bin`, migration/transfer `time_bin` payloads, and serialized restart state mirrors are diagnostic or compatibility mirrors only. They are refreshed from scheduler state with `syncTimeBinMirrorsFromScheduler(...)` and may be validated for corruption, but they are never fallback scheduling truth.
- Cell mirrors in particle-bound gas states map through the parent gas particle identity (`gasParticleIndexForCellRow(...)`), so a cell row count mismatch is not a reason to trust `CellSoa::time_bin` independently.
- Active-set descriptors for hierarchical execution must be constructed from scheduler output (`makeSchedulerActiveSetDescriptor(...)`) and must carry scheduler tick/generation provenance before solver callbacks consume them.

### Candidate criteria flow

1. Physics modules compute local candidate timesteps through typed criteria hooks such as CFL, gravity acceleration, source-term, or user-clamp hooks.
2. Candidates are submitted to the scheduler with an auditable `TimeStepCandidateSource` label.
3. `mapDtToTimeBin(...)` normalizes physical `dt` proposals into the integer power-of-two bin hierarchy using the configured `TimeStepLimits`.
4. `reconcileCandidateTransitions()` conservatively keeps the finest submitted bin per element, validates synchronization legality, and records clipping/candidate counters.
5. Pending transitions are committed only at scheduler-controlled substep boundaries; invalid destination synchronization fails fast rather than being silently clipped or delegated to mirrors.

### Active-set construction flow

1. `beginSubstep()` validates scheduler internals and rebuilds the compact active list from scheduler `bin_index` and `next_activation_tick`.
2. The workflow splits the scheduler active element list into particle/cell subsets as needed for callbacks.
3. `makeSchedulerActiveSetDescriptor(...)` stamps scheduler tick and state-generation provenance and validates that descriptor indices still match scheduler activity.
4. `StepOrchestrator::executeSingleStep(...)` runs the canonical KDK stage order and checks descriptor freshness when hierarchical callers pass an expected scheduler tick.
5. `endSubstep()` reconciles candidates, applies legal pending transitions, clears active flags, advances the integer tick, and only then are public `time_bin` mirrors refreshed from scheduler truth.

### Invariant framework

The Stage 2 invariant framework intentionally traps split-brain timestep ownership before solver callbacks can consume it:

- scheduler metadata arrays must remain same-sized, in-range, and internally consistent;
- active flags are scheduler caches for an open substep, not serialized live authority;
- active-set descriptors must match scheduler active elements and source generations;
- restart import must validate scheduler lanes and reject stale particle or cell `time_bin` mirrors before rebuilding mirrors from scheduler state;
- PM refresh events must be committed before the next kick opportunity can be registered;
- distributed PM cadence decisions remain rank-consensus metadata, but current evidence does not make Phase 3 multirate TreePM synchronization production-proven.

This contract changes documentation and test guardrails only. It does not change restart schema fields, solver numerics, normalized config output, provenance format, or deterministic scheduling semantics.

## Hierarchical integer timeline bins

`HierarchicalTimeBinScheduler` implements power-of-two integer bins (`dt_bin = dt_min * 2^bin`) and maintains:

- hot sidecars (`TimeBinHotMetadata`) for `bin_index`, `next_activation_tick`, and `active_flag`
- bin-local compact ownership lists (`m_elements_by_bin`) with O(1) move/erase sidecars
- per-substep compact active lists sorted for deterministic memory access

This avoids repeated full-state scans: each substep only touches bins synchronized to the current integer tick.

## State and active sets

`IntegratorState` tracks:

- `current_time_code`
- `current_scale_factor`
- `dt_time_code`
- `step_index`
- `TimeBinContext` (hierarchical-bin scaffold)

`ActiveSetDescriptor` carries compact particle/cell subsets and subset flags, avoiding global full-state sweeps in scheduler logic.

## Timestep criteria hooks

Typed helper criteria are provided and can be combined via `TimeStepCriteriaRegistry`:

- CFL limiter: `dt_CFL <= C_cfl * (Delta x / (|u| + c_s))`
- Gravity limiter: `dt_grav <= eta * sqrt(eps / |a|)`
- source-term limiter and user clamp hooks

`combineTimeStepCriteria` takes the conservative minimum across registered hooks and a fallback dt.

## Synchronization and legal transitions

Elements may request bin promotion/demotion through `requestBinTransition`. Transitions are only applied when the element is active and the current tick is synchronized to the destination bin period. Illegal attempts are recorded in diagnostics (`illegal_transition_attempts`) instead of being silently clipped.

## Diagnostics and pathological collapse visibility

`TimeBinDiagnostics` reports:

- occupancy and active counts per bin
- active fraction
- promoted/demoted counts
- clipped-to-min/max counters from dt mapping
- illegal transition attempts
- collapse candidate count when finest-bin occupancy dominates

## Cosmology helper equations

The helper API is aligned with comoving-variable usage:

- `dx/dt = v_pec / a`
- `dv_pec/dt = -H(a) v_pec - (1/a) grad_x phi + source_terms`
- `da/dt = a H(a)`

The implementation provides midpoint-integrated drift/kick prefactors over scale factor and a closed-form Hubble drag factor `a_begin / a_end`.

## Assumptions (explicit)

- Baseline stage scheme is kick-drift-kick only.
- PM long-range cadence is now explicit and rank-consensus audited; it remains a cadence-gated approximation rather than a full asynchronous multirate PM integrator.
- Time bins use strict power-of-two integer periods in tick space.
- Scale factor update in the orchestrator uses forward Euler (`advanceScaleFactorEuler`) for conservative simplicity.
- If no cosmology background is passed to `StepOrchestrator`, the orchestrator advances `time` only and leaves `current_scale_factor` unchanged.
- Output cadence policy remains external to the core integrator; output triggering is represented by the explicit `output_check` stage callback.

## Provenance and schema implications

No snapshot schema field is changed by this patch. Hierarchical scheduler metadata remains runtime-only and can be serialized in a follow-up restart/state I/O patch while preserving current output naming and external dataset compatibility.
