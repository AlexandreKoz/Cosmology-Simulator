# Stage 8 restart/output closeout audit

## Scope

This closeout covers the repaired Stage 8 restart/output maturity work after the post-audit upgrade pass. It includes restart-safe boundary checks, schema-versioned restart validation, snapshot-vs-restart file-kind separation, output cadence persistence, deterministic stochastic module contracts, a reusable direct-vs-restarted equivalence harness, stronger production-path equivalence tests, and restart read/write diagnostics.

## Final Stage 8 audit table

| Requirement | Status | Evidence |
|---|---|---|
| Restart writes are boundary-safe | Pass | `writeRestartCheckpointHdf5()` validates `core::assertCanWriteCheckpointAtBoundary()` before serialization. |
| Unsafe restart attempts fail loudly | Pass | Boundary and payload validation throw exceptions with KDK/scheduler/PM context instead of warning-only behavior. |
| Scheduler truth is serialized | Pass | Restart schema persists `/scheduler/{bin_index,next_activation_tick,active_flag,pending_bin_index}` plus scheduler tick/max-bin attributes. |
| KDK phase state is serialized | Pass | `/integrator` attributes include current/last boundary kind, KDK half-step flag, restart-safe flag, and time-bin context. |
| PM sync state is serialized | Pass | `/integrator` and distributed TreePM state persist cadence, kick opportunity, field version, refresh opportunity, refresh step/scale factor, and force-validity state. |
| RNG/stochastic state is represented | Pass for current modules | `/stochastic_state` records deterministic module contracts, seeds, rank-local offsets, and committed step indices. Stateful RNG engines remain unsupported and are rejected rather than fabricated. |
| Output cadence survives restart | Pass | `/output_cadence` persists output enabled flag, restart policy, due flags, last completed step, next snapshot step, and stems; equivalence tests compare due event sequences. |
| Normal snapshots are semantically distinct from restart files | Pass | Root `cosmosim_file_kind` separates `science_snapshot` and `restart_checkpoint`; readers reject the wrong kind. |
| Direct-vs-restarted equivalence tests exist | Pass | Dedicated tests cover DM-only, production TreePM kick continuation, production hydro solver continuation, multirate bins, output event sequence, and stochastic star-formation source terms. |
| Restart diagnostics/provenance are auditable | Pass | Schema v14 writes `/restart_diagnostics`; workflow emits both `restart.write.complete` and `restart.read.complete` operational events with schema/boundary/scheduler/PM/output/stochastic/hash metadata. |
| No transient scratch serialized | Pass | Restart payload uses `RestartPersistentStateView` and authoritative runtime state; transient active views and scratch buffers remain outside the payload. |
| Unsupported half-step restart | Pass/explicitly unsupported | Intentionally represented half-step restart is not implemented; such states are rejected. |

## Equivalence coverage

- `integration_restart_equivalence_dm_only`: direct 100 steps versus 40 + restart + 60, comparing particle, scheduler, integrator, PM, output, and stochastic contract state.
- `integration_restart_equivalence_treepm`: direct/restarted continuation where every harness step uses `gravity::TreePmCoordinator::solveActiveSet()` and compares the resulting trajectory plus PM cadence metadata.
- `integration_restart_equivalence_hydro_toy`: direct/restarted continuation where every harness step uses `hydro::HydroCoreSolver::advancePatch()` on a small periodic Sod-like problem and compares gas thermodynamics.
- `integration_restart_equivalence_multirate_bins`: direct/restarted continuation with multirate bin assignments, pending transitions, next activation ticks, and current tick comparison.
- `integration_restart_equivalence_output_enabled`: direct/restarted continuation with output cadence state and actual due snapshot/checkpoint event-step sequences.
- `integration_restart_equivalence_stochastic_sources`: direct/restarted continuation through production `physics::StarFormationModel` stochastic spawning with star sidecar and particle-growth comparison.

## Diagnostics metadata

Restart schema v14 writes `/restart_diagnostics` as an audit-only group. It includes schema identity, boundary labels, restart-safe decision, scheduler occupancy summary, PM cadence/field-version summary, output cadence summary, and stochastic module count. These fields are deliberately not used as continuation authority. The authoritative runtime truth remains the serialized state/integrator/scheduler/output/stochastic/distributed-gravity records plus the payload integrity hash.

The reference workflow now records both write-side and read-side restart operational events. `restart.write.complete` records the state being saved. `restart.read.complete` records the validated state after load, including the payload hash and compact scheduler/PM/output/stochastic summaries.

## Remaining risks

- FFTW-enabled PM coverage was not exercised in this environment because the local preset lacks FFTW development libraries. The TreePM equivalence test does exercise the production TreePM coordinator through the non-FFTW HDF5 preset.
- Stateful RNG engines remain unsupported by design. Current stochastic modules are deterministic from serialized inputs; any future stateful RNG must add explicit engine serialization and validation.
- Half-step/local-substep restart remains unsupported. This is safer than pretending such restarts are valid, but a future schema would be required for intentionally represented half-step continuation.
- The current output equivalence test compares due event sequences inside the restart harness. Full production resume-from-file workflow entry points remain a future extension.
