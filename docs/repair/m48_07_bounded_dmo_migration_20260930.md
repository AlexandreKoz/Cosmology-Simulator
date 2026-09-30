# M48-07 — bounded homogeneous-DMO migration and exact-ownership transaction closure

_Date: 2026-09-30_

## Scope

M48-07 adds a compact production migration transaction for homogeneous dark-matter-only state while retaining the existing `ParticleMigrationRecord` / AMR migration path as the full-physics compatibility path. The repair does not redesign decomposition, restart verification, AMR migration, or the MPI fragmentation protocol.

## Transaction selection

`compact_homogeneous_dmo` is admitted only when every rank agrees that the live state is homogeneous DMO, the particle scheduler is `kUniformRungZero` with `max_bin == 0`, no gas/AMR/star/BH/tracer migration state is active, no particle-indexed module sidecar requires transport, and no heterogeneous per-particle softening state requires the generic representation. Rank-local compact flags, drift epoch, scheduler tick, and scheduler max-bin are checked for communicator consensus before compact admission. Any unsupported exception selects `generic_full_physics` before memory admission.

## Compact wire truth

The compact DMO wire record is explicitly encoded and versioned; it is not an MPI transfer of a padded native C++ struct. Per migrated particle it carries only:

- particle ID;
- exact logical SFC key;
- three comoving position doubles;
- three peculiar-velocity doubles;
- mass.

Species, destination owner, timestep bin, scheduler rung metadata, common drift epoch, uniform flags, and full-physics sidecars are implicit or transaction-level truth and are therefore not repeated per particle.

## Streaming exchange

The existing bounded fragment header and `exchangeBoundedAlltoallBytes` transport remain authoritative. Compact migration builds only per-peer vectors of authoritative local row indices, sorted deterministically by row. The producer encodes the next DMO record directly from canonical state into the bounded packet. The receiver assembles at most the bounded current record for each peer, validates it, decodes it directly into its deterministic candidate row, and retains only inbound particle IDs for O(M_in) duplicate checking.

The compact path does not construct outbound or inbound populations of `ParticleMigrationRecord`.

## Candidate transaction and commit

The old authoritative state remains untouched while the candidate is prepared and transferred. The compact candidate owns only:

- position x/y/z;
- velocity x/y/z;
- mass;
- particle ID;
- SFC lane only when exact post-migration SFC truth cannot remain implicit.

Kept rows are copied by one scan of canonical state against a sorted O(M_out) outbound-row selection. No population `preserved_indices` or removal mask is required. Inbound rows are placed after the kept rows in peer-rank/record-sequence order.

Before commit, inbound IDs are sorted to detect inbound duplicates, and the kept candidate prefix is scanned against that sorted O(M_in) set to reject an inbound ID that aliases a kept particle. When the configured exact ownership audit is enabled, the complete candidate logical ID range is also validated before publication. Candidate preparation/validation failures are coordinated collectively; the old state remains authoritative.

`SimulationState::commitCompactHomogeneousDmoCandidate()` swaps the validated candidate into authority without compact-to-generic metadata materialization. It restores compact scalar metadata and applies the existing particle-index/gravity-source generation invalidation semantics.

## Uniform scheduler reconstruction

The compact path does not export persistent scheduler identity records and does not build destination particle IDs. A replacement `kUniformRungZero` scheduler is prepared before commit, under the migration reservation, with the final element count and the existing `current_tick`. Its only population lane is the scheduler's required uint32 active-identity lane. Move publication after state commit is statically required to be non-throwing.

## Bounded exact ownership validation

`validateExactGlobalOwnershipPartition()` no longer materializes complete hash-owner partitions. It uses deterministic ID hashing plus a secondary high-bit radix prefix. For each prefix it first counts current/expected IDs by hash owner, globally reduces those counts, and deterministically refines any prefix whose largest owner bucket would violate the workspace contract. An accepted bucket is redistributed with the existing bounded byte exchange, sorted, compared exactly, released, and followed by the next bucket.

The hard auxiliary workspace contract is:

```text
k_exact_ownership_validation_workspace_limit_bytes = 64 MiB per rank
```

The bucket admission threshold reserves four ID-width shares inside that limit for current IDs, expected IDs, and bounded transport coexistence. World-size one follows the same radix-bucket architecture. Duplicate/missing/extra totals are exact; diagnostic ID vectors are capped at 16 samples per class. `ExactOwnershipPartitionReport::valid()` depends on exact counts, not sample-vector sizes.

Compact-DMO callers pass the canonical particle-ID lane directly. They no longer copy it into `local_owned_particle_ids[N]` solely to satisfy the validator API.

## Memory admission and telemetry

The existing `MemoryGovernor` remains the only memory-admission authority. The compact plan reserves before candidate allocation and charges physical owners that actually exist:

- compact final candidate;
- O(M_out) peer selections plus sorted outbound membership selection;
- bounded packet and record-assembly staging;
- O(M_in) inbound duplicate-ID staging;
- uniform scheduler active-identity candidate;
- 64 MiB/rank exact-ownership workspace when the audit is enabled.

It does not charge generic `ParticleMigrationRecord` capacities, dynamic module-record heaps, compact metadata materialization, population kept-index/removal arrays, scheduler identity records, destination-ID arrays, or population local-ID maps because those objects do not exist on this path.

Telemetry now identifies `migration_representation` and reports compact candidate, selection, packet/record-assembly, exact-workspace, reservation, wire-byte, and physical-round quantities separately. Generic telemetry retains its generic owners.

## Source-derived 512^3 memory model

For 512^3 = 134,217,728 particles globally, with balanced final population and no explicit SFC lane, the compact candidate is 64 B/particle = 8.00 GiB aggregate. The uniform scheduler active-identity lane is 4 B/particle = 0.50 GiB. The exact validator bound is 64 MiB/rank = 0.50 GiB aggregate at eight ranks. Packet staging is approximately 64 MiB/rank from four 16 MiB bounded transport buffers, plus a sub-KiB-per-rank fixed-record assembly term at eight ranks, for approximately 0.50 GiB aggregate. These fixed owners total approximately 9.50 GiB aggregate.

Selection costs 8 B per outbound particle across the two O(M_out) row-index views, and inbound duplicate staging costs 8 B per inbound particle. In the deliberately conservative all-particles-migrate case, those two terms contribute about 2.00 GiB aggregate, giving approximately 11.50 GiB plus negligible peer-vector headers. If exact SFC keys must be materialized, add 8 B/final particle = 1.00 GiB, giving approximately 12.50 GiB in that conservative case. Normal rebalances scale those O(M) terms with the actual migrated population rather than N.

For comparison, the retained generic source model allocates 720-byte native `ParticleMigrationRecord` objects for both outbound and inbound capacities. If every 512^3 particle migrated once, that record term alone is approximately 180 GiB aggregate before dynamic heaps. The generic DMO path additionally models compact-metadata materialization (~45 B/current particle when SFC is implicit), population scheduler/remap arrays, N-sized index maps, packet staging, and old/new state coexistence. M48-07 therefore removes the dominant generic convenience representations rather than subtracting a fixed number from the estimator.

These figures are source-derived ownership estimates, not measured RSS.

## Generic fallback and reproducibility

The generic full-physics transaction, `ParticleMigrationRecord`, generic wire codec, and `AmrPatchMigrationRecord` remain intact for gas/full physics, sidecars, heterogeneous scheduler/canonical metadata, AMR-coupled migration, and unsupported DMO exceptions. Decomposition decisions, migrated IDs, phase-space values, exact SFC values, scheduler tick, ownership partition, deterministic candidate ordering, and generation invalidation semantics are intended to remain unchanged; empirical equivalence testing is deferred.

## Build status

Successful build-only validation:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4
```

The MPI configure/build-only attempt was blocked by the available environment: CMake could not find MPI C++ support (`MPI_CXX_FOUND` / `mpi-cxx`). No new toolchain was installed.

Tests: NOT RUN — intentionally deferred to the later M48 qualification campaign.

## Deferred work

- M48-08 — streaming restart verification.
- M48-09 — final stretch-memory package.
- M48-Q — whole-run 512^3 / 48-GiB qualification.
- Generic/full-physics migration intentionally retains its broader native records and conservative memory model.
