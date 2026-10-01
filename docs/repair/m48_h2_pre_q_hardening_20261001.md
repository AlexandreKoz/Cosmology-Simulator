# M48-H2 pre-Q hardening — 2026-10-01

## Scope

M48-H2 closes the remaining source-level validation, compact-state, TreePM
query, and restart-writer memory seams before M48-Q. It does not change PM or
TreePM force mathematics, decomposition policy, migration wire formats,
snapshot/restart schemas, or scientific precision.

## Particle identity certification

The production `SimulationState` particle-ID validators no longer construct a
population-sized `std::unordered_set<uint64_t>` or a population-sized sorted
copy. `ParticleIdentityValidationCertificate` records the identity generation,
particle count, local uniqueness proof, and nonzero-ID proof. A valid
certificate makes repeated uniqueness/persistent-ID checks allocation-free.

When a certificate is stale, exact local validation uses deterministic raw-ID
radix-prefix refinement and one bounded sort bucket. The CHUI-owned bucket is
limited to 32 MiB/rank independently of particle count. The proof remains exact:
zero persistent IDs and duplicate IDs are still rejected. The serial IC reader
uses this same validator instead of materializing `vector<uint64_t>[N]`.
Startup establishes a persistent-ID certificate under the process
`MemoryGovernor` before normal identity reporting. Snapshot writing consumes an
existing certificate and reserves the same 32 MiB diagnostic bound only when
revalidation is required.

Identity generation is invalidated at population/identity replacement
boundaries (`resizeParticles`, generic migration commit, compact-DMO migration
commit). Generic migration re-stamps the certificate only after its exact local
transaction checks prove nonzero uniqueness. Compact-DMO migration transfers a
pre-existing persistent certificate only when every rank entered the migration
with one and the transaction's existing inbound/kept-ID checks preserve the
set contract. Pure particle reorder does not invalidate identity certification
because it is a permutation of the same ID set.

Arithmetic identity summaries (`count`, `sum`, square-sum, xor) no longer hide
a node-based hash. Workflow callers supply certified uniqueness. The legacy
standalone summary overload reuses the same exact bounded validator (32 MiB maximum bucket)
for callers without a `SimulationState` certificate.

## Ownership-invariant validation

`OwnershipValidationWorkspace` no longer owns `particle_ids[N]`. Star, black
hole, and tracer marker lanes are allocated only when those sidecars exist, and
cell-owner scratch is allocated only when patch/cell topology requires it.
Absent marker lanes mean logically false membership without storage.

For homogeneous DMO with no gas cells, stars, black holes, or tracers, a state
with a valid identity certificate has zero population-scale ownership-invariant
scratch. If identity certification is stale, the only particle-scale validation
workspace is the fixed <=32 MiB exact-ID bucket.

At the 512^3 / 8-rank reference population (16,777,216 particles/rank), the old
generic DMO ownership workspace modeled 11 B/particle (8-byte ID copy plus three
1-byte cold markers), or 176 MiB/rank. After H2 those lanes are absent for
certified homogeneous DMO.

## Exact global ownership validation

`validateExactGlobalOwnershipPartition()` no longer constructs complete logical
per-peer send/receive payload vectors. Each prefix is counted globally before
payload movement. Prefixes are refined until every destination's current plus
expected comparison set fits the planned bounded comparison budget.

MPI transport is streamed in deterministic bounded rounds: IDs are scanned into
one send stage, exchanged with `MPI_Alltoallv`, appended directly into the
bounded destination comparison partition, and the stages are reused. The
serial path uses the same prefix/count model and allocates only the counted
comparison partition rather than reserving the complete input population.

The declared 64 MiB/rank CHUI-owned workspace is physically bounded with
fixed storage rather than capacity-growth vectors:

- one fixed comparison store: 24 MiB;
- one fixed send stage: 15 MiB;
- one fixed receive stage: 15 MiB;
- exact rank metadata arrays plus prefix/diagnostic allowance: <=2 MiB;
- reserved headroom inside the contract: 8 MiB.

The MPI path also places collective preflight gates around rank-local round
planning, payload preparation, and payload consumption so a local preparation
failure is propagated before peers enter the next incompatible collective.

The validator reserves the complete 64 MiB diagnostic workspace through the
`MemoryGovernor` when one is supplied. Final-run validation supplies it
explicitly. Compact-DMO migration auditing remains covered by the enclosing
migration transaction reservation, whose admission plan already charges the
64 MiB exact-ownership bound. Duplicate, missing, and extra counts remain exact;
only diagnostic example-ID vectors are capped.

## Compact homogeneous-DMO core APIs

`buildParticleReorderMap()` now reads SFC/species keys through logical
`SimulationState` accessors. `reorderParticles()` reorders only physically
materialized metadata lanes; implicit time-bin/species/flags/owner/drift lanes
remain implicit. Homogeneous species/SFC ordering therefore degenerates to the
stable identity order when the logical key is uniform, without metadata
materialization. The homogeneous DMO species index remains its identity
representation after reorder.

Generic core transfer helpers (`packSpeciesTransferPacket()` and particle
migration record packing) now read time-bin, owner, drift, species, flags, and
SFC metadata through logical accessors, so a valid compact state cannot index
an intentionally absent dense lane. The optimized compact-DMO migration path
is unchanged and is not redirected through the generic full-physics record
path.

## TreePM top-domain query scratch

`TopLevelDomainHierarchy::ownersWithinCutoff()` no longer constructs an owner
vector and traversal stack per target. The caller prepares one function-local
query scratch object containing a node-count-bounded traversal stack, a
world-size-bounded unique-owner vector, and world-size owner markers. Multiple
top-domain leaves from the same owner are deduplicated before insertion, so the
owner vector cannot grow beyond its prepared world-size capacity. After
`prepare()`, one target query performs no dynamic heap allocation. Scratch is
local to the current TreePM transaction/batch path rather than global shared
state, preserving future thread-local/batch-local ownership and deterministic
sorted owner order.

## Restart distributed-owner write view

The checkpoint writer no longer materializes
`distributed_gravity_state.owning_rank_by_item[N]`. `DistributedRestartState`
can borrow the live `SimulationState` as its write-side ownership source and
exposes logical `owningRankItemCount()` / `owningRankAt()` accessors. The same
canonical serializer is used for borrowed writes and owning readback state, so
`item_count=N` and every `rank[i]=...` line retain the existing byte format.
Counting-stream sizing, streaming HDF5 writing, and integrity hashing all use
that same serializer. Deserialization/restoration continues to own
`owning_rank_by_item` exactly as before.

At 16,777,216 particles/rank, this removes the former 4 B/particle writer
staging vector (64 MiB/rank) from the checkpoint event. The existing bounded
streaming verification workspace remains the modeled restart-side temporary
owner; no second full-state ownership vector is charged.

## Memory-accounting changes

Directly affected estimators now charge the physical H2 architecture:

- local particle-ID validation: certificate hit = 0 allocation; stale proof <=32 MiB/rank;
- homogeneous-DMO ownership invariant validation: 0 population lanes with a valid certificate, otherwise the fixed ID proof bound;
- exact global ownership audit: <=64 MiB/rank for all CHUI-owned validator buffers together;
- distributed restart owner staging: 0 B/particle on the write path.

These are source-derived ownership limits, not measured RSS claims. M48-Q owns
whole-process RSS/high-water qualification.

## Source acceptance notes

The target production validators and identity-summary path no longer contain a
population-sized `unordered_set<uint64_t>`. Snapshot readiness still has lazy
`unordered_map` lookups for real gas parent/patch resolution; they are not
created for homogeneous DMO. Transaction-scoped `unordered_set`/`unordered_map`
uses remain in generic full-physics particle/gas/AMR migration validation. H2
does not redirect the accepted compact-DMO migration architecture through that
generic path.

The exact global validator no longer calls the complete-payload
`exchangeBoundedAlltoallBytes()` API and contains no complete per-peer
`send_payloads`/`recv_payloads`. The restart writer contains no population
`owning_rank_by_item.reserve/push_back` construction.

## Build-only validation

Executed successfully:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug --target cosmosim_harness -j4
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug --target cosmosim_harness -j4
```

The HDF5 build emitted the pre-existing `src/amr/amr_ghost_fill.cpp` warning for
an ignored `[[nodiscard]]` return value; it did not fail the build.

MPI configuration was attempted with:

```bash
cmake --preset mpi-release
```

and is blocked in this environment because CMake cannot find the MPI C++
development package (`Could NOT find MPI (missing: MPI_CXX_FOUND CXX)`; the
`mpi-cxx` pkg-config module is absent). No toolchain was installed for H2.

Tests: NOT RUN — intentionally deferred to M48-Q.

## Reproducibility and remaining qualification

H2 changes representation and validation allocation strategy only. Exact ID
semantics, deterministic owner ordering, canonical restart serialization, and
science/restart schemas are preserved. No numerical force, integration, or
precision policy is changed.

The intended next step is M48-Q: whole-run 512^3 / 48-GiB qualification. The
known compact decomposition event remains approximately 13.5 GiB aggregate in
the source model and is a qualification target, not an H2 defect.
