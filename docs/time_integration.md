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

### TreePM long-range cadence contract (Phase 1)

The reference workflow now applies an explicit, config-driven long-range PM refresh cadence:

- cadence control: `numerics.treepm_update_cadence_steps` (integer `>= 1`)
- cadence unit in Phase 1: **gravity kick opportunities** (`gravity_kick_pre` and `gravity_kick_post`)
- policy:
  - refresh PM mesh solve when no cached field exists, or when
    `current_kick_opportunity - last_refresh_opportunity >= cadence_steps`
  - otherwise reuse the cached long-range field and only rebuild/evaluate the short-range tree residual

The temporal contract for reuse is explicit: a reused PM field corresponds to the particle state and scale factor from the last refresh opportunity, and that metadata is recorded by the workflow report/events.
With `treepm_update_cadence_steps = 1` (default), behavior reduces to immediate PM refresh at every kick opportunity.

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
- PM long-range cadence in Phase 1 is single-rank and intentionally conservative; it is not a full multirate integrator redesign.
- Time bins use strict power-of-two integer periods in tick space.
- Scale factor update in the orchestrator uses forward Euler (`advanceScaleFactorEuler`) for conservative simplicity.
- If no cosmology background is passed to `StepOrchestrator`, the orchestrator advances `time` only and leaves `current_scale_factor` unchanged.
- Output cadence policy remains external to the core integrator; output triggering is represented by the explicit `output_check` stage callback.

## Provenance and schema implications

No snapshot schema field is changed by this patch. Hierarchical scheduler metadata remains runtime-only and can be serialized in a follow-up restart/state I/O patch while preserving current output naming and external dataset compatibility.
