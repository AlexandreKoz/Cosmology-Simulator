# Stage 8 restart/output closeout audit

## Scope

This closeout covers Stage 8 prompts 8.8 through 8.10: required restart equivalence coverage, restart/output diagnostics and provenance, and final integration audit. It assumes the previous Stage 8 work already established restart-safe boundary checks, schema-versioned restart validation, snapshot-vs-restart file-kind separation, output cadence persistence, stochastic module persistence contracts, and the reusable restart equivalence harness.

## Final Stage 8 audit table

| Requirement | Status | Evidence |
|---|---|---|
| Restart writes are boundary-safe | Pass | `writeRestartCheckpointHdf5()` validates `core::assertCanWriteCheckpointAtBoundary()` before serialization. |
| Unsafe restart attempts fail loudly | Pass | Boundary and payload validation throw exceptions with KDK/scheduler/PM context instead of warning-only behavior. |
| Scheduler truth is serialized | Pass | Restart schema persists `/scheduler/{bin_index,next_activation_tick,active_flag,pending_bin_index}` plus scheduler tick/max-bin attributes. |
| KDK phase state is serialized | Pass | `/integrator` attributes include current/last boundary kind, KDK half-step flag, restart-safe flag, and time-bin context. |
| PM sync state is serialized | Pass | `/integrator` and distributed TreePM state persist cadence, kick opportunity, field version, refresh opportunity, refresh step/scale factor, and force-validity state. |
| RNG/stochastic state is represented | Pass | `/stochastic_state` records deterministic module contracts, seeds, rank-local offsets, and committed step indices. Stateful RNG engines remain unsupported and are rejected rather than fabricated. |
| Output cadence survives restart | Pass | `/output_cadence` persists output enabled flag, restart policy, due flags, last completed step, next snapshot step, and stems. |
| Normal snapshots are semantically distinct from restart files | Pass | Root `cosmosim_file_kind` separates `science_snapshot` and `restart_checkpoint`; readers reject the wrong kind. |
| Direct-vs-restarted equivalence tests exist | Pass | Dedicated tests cover DM-only, TreePM cadence/PM metadata, hydro toy gas state, multirate bins, and output-enabled cadence. |
| Restart diagnostics/provenance are auditable | Pass | Schema v14 adds `/restart_diagnostics` with compact scheduler/PM/output/stochastic summaries and workflow restart events include schema/boundary/PM/output metadata. |
| No transient scratch serialized | Pass | Restart payload uses `RestartPersistentStateView` and authoritative runtime state; transient active views and scratch buffers remain outside the payload. |
| Unsupported half-step restart | Pass/explicitly unsupported | Intentionally represented half-step restart is not implemented; such states are rejected. |

## Tests added in this closeout

- `integration_restart_equivalence_dm_only`: direct 100 steps versus 40 + restart + 60, comparing particle and scheduler state.
- `integration_restart_equivalence_treepm`: same equivalence pattern with nontrivial PM cadence and PM sync metadata comparison.
- `integration_restart_equivalence_hydro_toy`: same equivalence pattern with gas cells and thermodynamic sidecar comparison.
- `integration_restart_equivalence_multirate_bins`: same equivalence pattern with multirate bin assignments, pending transitions, next activation ticks, and current tick comparison.
- `integration_restart_equivalence_output_enabled`: same equivalence pattern with output cadence/counter state comparison.

## Diagnostics metadata

Restart schema v14 writes `/restart_diagnostics` as an audit-only group. It includes schema identity, boundary labels, restart-safe decision, scheduler occupancy summary, PM cadence/field-version summary, output cadence summary, and stochastic module count. These fields are deliberately not used as continuation authority. The authoritative runtime truth remains the serialized state/integrator/scheduler/output/stochastic/distributed-gravity records plus the payload integrity hash.

## Remaining risks

- The new TreePM equivalence test exercises PM cadence and persistent PM sync metadata through the restart harness, not a full expensive force-solve trajectory. A future heavier validation should combine this harness with an actual small TreePM force step when CI runtime allows it.
- Stateful RNG engines remain unsupported by design. Current stochastic modules are treated as deterministic from serialized inputs; any future stateful RNG must add explicit engine serialization and validation.
- Half-step/local-substep restart remains unsupported. This is safer than pretending such restarts are valid, but a future schema would be required for intentionally represented half-step continuation.
