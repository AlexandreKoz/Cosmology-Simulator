# Stage 2 timestep hardening closeout

This note records the focused Stage 2 hardening pass that closes the remaining
scheduler-authority escape hatches found during adversarial review.

## Closed escape hatches

1. **Scheduler active-set provenance is mandatory for scheduler-derived work.**
   `StepOrchestrator::executeSingleStep` now rejects scheduler-derived active-set
descriptors unless the caller supplies the scheduler tick. Hierarchical workflows
should prefer `executeSchedulerSubstep`, which builds the descriptor from scheduler
authority and passes the tick internally.

2. **Time-bin particle reordering no longer consumes derived mirrors.**
   `ParticleReorderMode::kByTimeBin` is intentionally rejected by the legacy
`buildParticleReorderMap` entry point because `ParticleSoa::time_bin` is a derived
mirror. Production time-bin/rung reordering must call
`buildParticleReorderMapByScheduler(state, scheduler)`, which consumes scheduler
hot metadata directly.

3. **PM cadence now has a restartable authority DTO.**
   `PmSynchronizationPersistentState` captures cadence, kick opportunity,
field-version, last refresh, and pending refresh metadata. Import validates pending
refresh consistency so a skipped TreePM refresh cannot be hidden by restart.

4. **Scheduled element identity is explicit.**
   `ScheduledElementKey` documents the identity namespace. Current production
scheduling is particle-backed, with gas cells treated as particle-bound carriers.
Future AMR patches or independent mesh cells must enter through an explicit
`ScheduledElementKind` instead of overloading local row indices.

## Tests added or hardened

- scheduler-derived active-set without explicit tick is rejected;
- scheduler-backed time-bin reorder ignores deliberately stale mirrors;
- legacy mirror-based time-bin reorder throws;
- PM synchronization persistent state round-trips through a pending refresh;
- hot/cold sidecar layout test now uses scheduler authority for time-bin ordering.

## Remaining policy

`ParticleSoa::time_bin` and `CellSoa::time_bin` remain compatibility/diagnostic
mirrors for I/O, tests, and migration packets. They are not timestep authority.
New production code must not sort, schedule, or commit timestep decisions from
those mirrors. Scheduler-owned state remains the only authority for bin/rung,
next activation tick, active-set construction, and PM refresh legality.
