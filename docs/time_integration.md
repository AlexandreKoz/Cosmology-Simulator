# Time Integration Framework

This document defines the baseline time-integration contract for CosmoSim and records conservative assumptions made in the first implementation.

## Stage contract

The authoritative per-step stage order is:

1. `gravity_kick_pre`
2. `drift`
3. `hydro_update`
4. `source_terms`
5. `gravity_kick_post`
6. `analysis_hooks`
7. `output_check`

The scheduler (`StageScheduler`) owns this ordering and solver modules interact via `IntegrationCallback` implementations. This keeps gravity/hydro/source/output sequencing explicit and auditable.

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
3. hydro/source stages
4. `K_PM+Tree^post` on active set at sync surface `kick_post`
5. analysis/output stages

Both kick surfaces use the same PM sync contract above.

#### Inactive-particle treatment

- Inactive particles do **not** receive PM or tree kick updates on a given kick opportunity.
- Their state evolves only when their bin activates on a later synchronized kick surface.
- This is recorded per event as `inactive_particles_skipped`.

#### Restart continuation rule

- Restart payload persists cadence metadata (`gravity_kick_opportunity`, cadence steps, field version, last refresh opportunity, field-built step/scale factor).
- Restart policy is `deterministic_rebuild`: PM mesh field values are not checkpointed; after resume, the first kick opportunity refreshes/rebuilds the long-range PM field before reuse can occur.
- Cadence metadata remains auditable across restart boundaries via distributed restart state and provenance fields.

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
