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
- `stageContracts()` returning a typed `StageContract` entry for each declared stage (required inputs, mutable domains, outputs, side effects, sync requirement, active-set family, restart/output safety, and owning subsystem);
- `onStage(...)` logic that assumes the declared stage contract and throws/asserts on impossible off-stage direct calls instead of returning silently.

Existing code that previously registered a broad callback and checked `context.stage` internally must split the behavior into explicit stage declarations or return the exact small stage set needed by the handler. No config keys, snapshot/restart payloads, or solver numerics migrate for this interface repair.

### TreePM long-range force-refresh contract (distributed TreePM runtime)

The reference workflow applies an integrator-owned, config-driven TreePM
operator inside KDK. An invalid initial force cache is bootstrapped at the
globally authorized `gravity_kick_pre` boundary. After drift, the
`force_refresh` stage evaluates the force at the end-of-step source epoch; the
`gravity_kick_post` stage consumes that cache for the closing half-kick.

The synchronization decision is explicit in workflow runtime events and
cadence records (`pm_sync_surface`, `gravity_kick_opportunity`,
`field_version`, `last_refresh_opportunity`, `active_particles_kicked`,
`inactive_particles_skipped`).

For the in-flight refresh, the integrator-issued `PmRefreshDirective` is the
decision authority until `PmSynchronizationState` commits it. The
`gravity.pm_long_range_field` event is emitted in that interval, so its
opportunity, version, build step, and build scale-factor payload must come from
the directive/solver decision. Reading the still-previous committed cadence
state there would produce a self-contradictory refresh event. This event fix
does not move live cadence ownership away from `PmSynchronizationState`.

#### Production refresh and lower-level reuse rule

- `numerics.treepm_update_cadence_steps` is normalized and validated as exactly
  `1` for the production `ReferenceWorkflow`.
- Every authorized, rank-coordinated production force-refresh evaluates a new
  PM field. All ranks enter the decision and collectives, including ranks with
  zero local targets.
- Lower-level scheduler and `TreePmCoordinator` interfaces retain explicit
  refresh/reuse support for focused tests and future multirate work. Reuse is
  legal only when a compatible cached field is present on every rank. A
  missing or incompatible cache throws coherently; it is never converted into
  an implicit refresh.
- Cache compatibility covers force epoch, force-evaluation scale-factor
  metadata, `G_code`, split and box geometry, assignment, boundary,
  decomposition-layout mode, and window-deconvolution policy. Particle
  ownership/decomposition epoch is intentionally excluded because the cached
  PM field is owned by the fixed FFT slab layout, not particle rows.
- Refresh/reuse votes and cadence metadata are rank-consensus facts.

The scale factor stored with a field is source-epoch/validity metadata. The PM
and tree kernels are scale-free; cosmological time dependence is applied only
by the KDK factors described below.

#### KDK operator order with PM sync

For each global step:

1. `K_PM+Tree^pre` on the active set, bootstrapping force truth if necessary
2. `D` drift on active set
3. rank-coordinated `force_refresh` at the post-drift source epoch
4. hydro/source stages
5. `K_PM+Tree^post` on the active set using the refreshed force cache
6. analysis/output stages

The integrator owns refresh authorization and commits the PM synchronization
event only after the TreePM callback consumes it successfully.

#### Inactive-particle treatment

- The production workflow currently requires `hierarchical_max_rung = 0`, so
  mixed-rung inactive-particle kicks are not a supported production path.
- The lower-level scheduler still records compact active sets and
  `inactive_particles_skipped` for infrastructure validation. Those lanes do
  not constitute a production-certified multirate TreePM integrator.

#### Restart continuation rule

- Restart payload persists cadence metadata (`gravity_kick_opportunity`,
  cadence steps, field version, last refresh opportunity, and field-build
  step/scale factor).
- Restart policy is `deterministic_rebuild`: PM mesh/tree scratch is not
  checkpointed and cannot be treated as reusable field truth. Restart schema
  v20 does persist the committed gravity-force history keyed by stable particle
  ID for exact KDK continuation; import remaps compatible rows and invalidates
  stale or migrated rows. The next integrator-authorized force surface
  reconstructs transient PM/tree state.
- Cadence metadata remains auditable across restart boundaries via distributed restart state and provenance fields.


## Stage 2 timestep authority contract

Stage 2 uses a single-owner timestep model: `HierarchicalTimeBinScheduler` is the only live authority for per-element bin assignment, next activation, active flags, pending transitions, active-set construction, and PM kick cadence metadata. Solver callbacks may propose timestep candidates and consume scheduler-built active sets, but they must not treat `ParticleSoa::time_bin`, `CellSoa::time_bin`, migration records, or restart mirrors as authority.

Production currently fails closed at `numerics.hierarchical_max_rung = 0`.
Nonzero rungs are rejected by typed config validation and by the production
workflow because per-element kick/drift epochs required for mixed-rung KDK are
not yet authoritative. The scheduler interfaces below remain infrastructure
for future multirate work, not an enabled production capability.

### Scheduler authority and mirror policy

- Authoritative live lanes: `HierarchicalTimeBinScheduler` hot metadata (`bin_index`, `next_activation_tick`, `active_flag`, `pending_bin_index`) and `PmSynchronizationState` cadence fields (`gravity_kick_opportunity`, cadence steps, field version, last refresh opportunity, field-built step/scale factor).
- Derived mirrors: `ParticleSoa::time_bin`, `CellSoa::time_bin`, migration/transfer `time_bin` payloads, and serialized restart state mirrors are diagnostic or compatibility mirrors only. They are refreshed from scheduler state with `syncTimeBinMirrorsFromScheduler(...)` and may be validated for corruption, but they are never fallback scheduling truth.
- Cell mirrors in production gas-cell states map through dense local gas-cell rows validated by `GasCellIdentityMap`.
  Legacy particle-bound import paths may still map through a parent gas particle, but production scheduler invariants
  do not require local gas-particle count to equal local gas-cell count.
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
- distributed PM decisions remain rank-consensus metadata;
- production restart import requires particle and gas schedulers to use
  `max_bin = 0`, with all current/pending bins at zero or unset as appropriate;
- current evidence does not make mixed-rung TreePM synchronization
  production-proven.

The fail-closed rung-zero policy prevents the scheduler scaffold from being
mistaken for production multirate integration.

## Hierarchical integer timeline bins (infrastructure)

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
- directional hydro CFL limiter:
  `dt_hydro_CFL = min_axes(C_cfl * cell_width_axis_code / (abs(v_axis_code) + c_s))`
- Gravity limiter for comoving softening and scale-free TreePM `A`:
  `dt_grav = eta * sqrt(a^3 eps_com / |A|)`
- source-term limiter and user clamp hooks

`computeGravityTimeStep(...)` itself remains coordinate-system neutral and
evaluates `eta sqrt(length/acceleration)`. Cosmological callers should instead
use the public `ComovingGravityTimeStepInput` and
`computeComovingGravityTimeStep(...)` interface. It accepts `eps_com`, the
scale-free TreePM magnitude `|A|`, and `a`, validates a finite positive scale
factor, converts to comoving-coordinate acceleration `|A|/a^3`, and delegates
the final length/acceleration calculation to the generic helper. The production
workflow uses the committed `IntegratorState.current_scale_factor` because
adaptive-bin proposals run after step commit. Using the peculiar-velocity
response `A/a^2` directly with a comoving length would mix coordinate systems.

### Public API migration: cosmological gravity timestep

Code that previously constructed `GravityTimeStepInput` from a comoving
softening and a TreePM force lane must migrate to
`ComovingGravityTimeStepInput` from
`include/cosmosim/core/time_integration.hpp`. Pass the unscaled TreePM
`|A|` and the force-epoch scale factor; do not pre-divide by `a^2` or `a^3`.
The existing `GravityTimeStepInput`/`computeGravityTimeStep(...)` API remains
available for callers whose length and acceleration are already expressed in
one coordinate system.

`combineTimeStepCriteria` takes the conservative minimum across registered hooks and a fallback dt.

Hydro CFL candidates submitted by the reference workflow use local Cartesian patch widths where gas-cell
center metadata identifies a row-ordered patch, otherwise the same near-cubic fixed-grid fallback used by the
hydro patch builder. The scheduler remains the only authority for the accepted bin; hydro submits
`TimeStepCandidateSource::kHydroCfl` candidates and the callback verifies the accepted `dt_time_code` against the
current directional CFL bound before invoking the finite-volume update. A violation is a hard runtime error, not a
silent clamp inside the hydro solver.

`HydroCflDiagnostics` records the worst checked gas cell with local row, gas-cell ID when available, patch ID/row
when available, proposed and accepted `dt`, CFL number, safety factor, velocity components, sound speed, cell widths,
and limiting axis. The reference workflow emits this through the profiler event `hydro.cfl_guard` for reporting.

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

The helper API uses physical peculiar velocity `u` and comoving position `x`:

- `dx/dt = u / a`
- `du/dt + H(a) u = A / a^2 + source_terms`
- `da/dt = a H(a)`

Here `A = -grad_x Psi` is the scale-free comoving Newtonian kernel returned by
both PM and the complementary tree residual. Neither gravity branch inserts an
additional factor of `a^2`. Collisionless particles receive the response
through KDK Hubble-drag/kick factors. Gas receives the same `A/a^2` response in
`ComovingGravityExpansionSource`, which also applies the conservative
expansion sources using the integrator-owned `a` and `H` in code units.

For hydro conserved variables `m=rho_com u`,
`E=rho_com(e+|u|^2/2)`, and kinetic energy density `K`, the explicit source is

```text
S_m = rho_com A/a^2 - H m
S_E = rho_com (u dot A)/a^2 - H (2 K + 3 P_com).
```

Because hydro runs after drift/force refresh but before the step commit,
fixed-grid and AMR callbacks consume
`StepContext.timeline_step.{scale_factor_end,hubble_end_code}`. The committed
`IntegratorState` still represents the step-begin epoch at that point; solver
callbacks must not use it or independently recompute an SI-valued Hubble rate.

The implementation provides midpoint-integrated drift/kick prefactors over scale factor and a closed-form Hubble drag factor `a_begin / a_end`.

## Assumptions (explicit)

- Baseline stage scheme is kick-drift-kick only.
- Production PM refresh cadence is exactly one at every rank-coordinated force
  evaluation. Lower-level reuse is explicit and fail-closed, but is not a
  production asynchronous multirate PM integrator.
- Time bins use strict power-of-two integer periods in tick space.
- Production config currently permits only rung zero despite the lower-level
  multi-bin scheduler implementation.
- Scale factor update in the orchestrator uses forward Euler (`advanceScaleFactorEuler`) for conservative simplicity.
- If no cosmology background is passed to `StepOrchestrator`, the orchestrator advances `time` only and leaves `current_scale_factor` unchanged.
- Output cadence policy remains external to the core integrator; output triggering is represented by the explicit `output_check` stage callback.

## Provenance and schema implications

No snapshot schema field is changed by this gravity repair. Restart schema
version 20 already persists scheduler/cadence truth needed for deterministic
continuation plus committed stable-ID-keyed force history. PM mesh/tree scratch
and coordinator-local cache structures remain transient and are
deterministically rebuilt.
