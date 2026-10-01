# M48-08 — streaming restart verification and borrowed force-cache persistence

_Date: 2026-09-30_

## Scope

M48-08 separates checkpoint write-time verification from actual restart restoration. The genuine `readRestartCheckpointHdf5()` path remains the complete owning restoration path. Ordinary checkpoint publication now verifies the committed current-schema HDF5 file directly against the live restart payload using bounded hyperslab reads.

The restart schema, schema version, dataset/group names, integrity algorithms, and restore semantics are unchanged.

## Previous checkpoint peak

Before this repair, `OutputRestartRuntime` wrote a checkpoint, exported an owning `GravityForceCachePersistentState`, then called `readRestartCheckpointHdf5()` and retained the resulting `RestartReadResult` while the live simulation remained resident. The verification event therefore admitted a second `SimulationState`, scheduler persistent arrays, a second force-cache population, distributed metadata, and other restore objects solely to prove a file that had just been written.

The force cache also paid a writer-side 32 B/DM-particle replica (ID plus acceleration XYZ) before serialization.

## New verification architecture

`verifyRestartCheckpointHdf5()` is a write-time verifier, not a restore reader. It:

- requires the current restart schema and runs the existing schema validator;
- reserves a fixed 16 MiB/rank application workspace through the existing `MemoryGovernor`;
- reads population datasets using bounded HDF5 hyperslabs;
- compares floating-point restart truth representation-exactly;
- obtains compact homogeneous-DMO metadata through logical `SimulationState` accessors rather than materializing sidecar lanes;
- obtains scheduler truth through scheduler logical accessors rather than `exportPersistentState()` for particle verification;
- streams `/distributed_gravity/state` against `DistributedRestartState::serializeTo()` through a comparing stream buffer;
- streams gravity force-cache identity and acceleration lanes directly against the borrowed live view;
- checks current v23 legacy/FNV and SHA-256 integrity metadata against the expected live payload;
- returns only small diagnostics/telemetry state.

The verifier result contains no population-sized state.

## Borrowed gravity force cache

`GravityRestartStateProvider` now exposes `restartForceCacheView()`. `GravityRuntime` validates cache generation, extent, and identity freshness, then returns spans borrowing:

- authoritative particle IDs from `SimulationState`;
- authoritative gas-cell IDs from `SimulationState`;
- live particle acceleration XYZ from `GravityRuntime`;
- live cell acceleration XYZ from `GravityRuntime`.

`RestartWritePayload` retains the historical owning force-cache pointer for direct library/test compatibility, but adds the non-owning `gravity_force_cache_view`. The production workflow uses the view. The writer and integrity traversal consume whichever representation is supplied without changing the on-disk schema.

`GravityForceCachePersistentState` remains the owning type used by genuine restart read/import.

## Compact DMO and uniform scheduler

For homogeneous DMO, explicit HDF5 `time_bin`, SFC/species/flags/owner/drift datasets are checked against logical accessors such as `particleTimeBin()`, `particleSfcKey()`, and `particleOwningRank()`. No compact metadata materialization is introduced.

Particle scheduler datasets are compared chunk-by-chunk through `binIndex()`, `nextActivationTick()`, `isElementActive()`, and `pendingBinIndex()`. The write-time verifier does not construct a second `TimeBinPersistentState` merely to compare the on-disk particle scheduler.

## Distributed state

The serialized `/distributed_gravity/state` byte dataset is not read into a complete string. The expected `DistributedRestartState` is serialized into a bounded comparing stream buffer which fetches corresponding HDF5 byte chunks, compares, and discards them.

The optional writer-side `owning_rank_by_item` vector remains intentionally retained in this campaign. Removing that approximately 4-B/particle metadata owner is deferred rather than broadening the mandatory repair.

## Integrity behavior

The writer still emits the same legacy FNV-1a identity and the same `sha256-canonical-le-v1` digest. The verifier checks those stored attributes against the digest of the expected live restart payload in addition to direct persisted-field comparisons. The file-format integrity contract is therefore unchanged.

## Memory model

Normal compact-DMO checkpoint verification no longer budgets or constructs:

- a second `SimulationState`;
- a readback particle scheduler population;
- a readback gas scheduler population solely for verification;
- a writer-side owning force-cache replica;
- a readback force-cache replica.

The production verification workspace is fixed at 16 MiB/rank of CHUÍ-owned diagnostic scratch. At eight ranks that is 128 MiB aggregate, independent of particle count. HDF5/libc internal allocations remain external-runtime memory.

The live gravity force cache remains owned once by `GravityRuntime`. At 512^3 DMO, eliminating the old exported force-cache replica removes the source-derived ~32 B/particle (~4 GiB aggregate) writer copy. Eliminating the complete write-time restore removes the historical dominant restart readback duplication event; this note does not claim measured RSS.

Actual resume still uses `restartReadCandidateStagingBytes()` and the full reader's MemoryGovernor admission because restoration legitimately constructs a candidate state.

## Observability

The runtime now reports `restart_verification_mode=streaming_field_exact` together with verification workspace high-water, bytes read, dataset count, and failed field. The old `restart.read.complete` write-time event has become `restart.verify.complete`; it no longer claims a full restore read occurred.

## Reproducibility

M48-08 changes staging/verification ownership only. The restart schema, logical values, scheduler truth, force-cache truth, distributed restart serialization, integrity algorithms, and actual restore behavior are intended to remain unchanged.

## Validation status

Tests: NOT RUN — intentionally deferred to the later M48 qualification campaign.

Build: NOT RUN — intentionally deferred by campaign instruction.

Static source review only is part of this campaign.

## Deferred work

- optional elimination of the writer-side distributed `owning_rank_by_item` vector;
- M48-09 final stretch-memory package;
- M48-Q whole-run 512^3 / 48-GiB qualification.
