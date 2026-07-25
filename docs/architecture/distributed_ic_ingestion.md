# Distributed initial-condition ingestion

## Decision

External cosmological IC ingestion is a workflow-owned, explicitly configured,
multifile, bounded-staging pipeline. No MPI rank receives an implicit license to
read or materialize the complete universe. Distributed ingestion remains a
**provisional runtime capability** until the registered MPI acceptance matrix is
executed successfully on a machine with MPI, HDF5, and FFTW-MPI.

## Layering

- `core`: typed IC convention, complete direct-bridge values, Type2/Type3
  policies, validation, and normalized provenance. It has no dependency on
  `io`.
- `io`: strict audit manifest v3, file-set inspection, the authoritative
  dimensional conversion engine, sidecar construction, explicit wire format,
  serial/distributed readers, coverage proofs, and scientific validation.
- `parallel`: borrowed `MpiContext` plus lower-layer communication primitives.
  It does not initialize MPI or force initial-state replication.
- `workflows`: convention dispatch, config-relative path resolution, borrowed
  runtime services, manifest output, profiler publication, and the
  `already_partitioned` startup boundary.

## No-guess scientific boundary

`gadget_arepo_bridge_v1` is not a shorthand for kpc, solar masses, km/s, or a
particular GADGET velocity convention. A direct bridge requires explicit source
SI scales, coordinate frame, velocity convention, and all length/mass/velocity
h and scale-factor exponents. `manifest_v1` instead takes those values from a
strict authoritative manifest. Canonical CHUÍ input validates its canonical
header and does not use bridge defaults.

Runtime and `cosmosim_convert_ic` call the same C++ conversion implementation.
The Python scripts are launchers only and contain no independent scientific
conversion formulas.

## Distributed phases

1. Rank zero reads only enough of the requested first member to determine the
   deterministic ordered file list.
2. File indices are assigned by `file_index mod world_size`. Each member is
   opened, schema-inspected, and SHA-256 hashed exactly once by its assigned
   rank. More ranks than files are valid.
3. Bounded manifest fragments are gathered. Rank zero assembles and validates
   the canonical ordered manifest without rereading source files.
4. The canonical manifest is broadcast in bounded string chunks. Every rank
   deserializes it and agrees on its SHA-256 digest.
5. A globally deterministic chunk index assigns each `(file, PartType, begin,
   count)` token to exactly one reader rank. An exact collective coverage check
   requires one completion for each token.
6. The reader reads and converts one bounded chunk, validates all coordinates,
   velocities, masses, and required sidecars, then serializes explicit
   little-endian records into per-owner buckets.
7. Send layout, receive layout, allocation, exchange, wire validation, and
   deserialization are separate collective phases. Local exceptions are voted
   before the next MPI collective.
8. Source IDs and received IDs are hash-partitioned and balanced exactly for the
   active chunk. Each ID must contribute one source record and one final record.
9. Each rank appends only received owner-local records and rebuilds local
   species and sidecar indices.
10. Exact global duplicate-ID validation, count/species/mass reductions,
    ownership validation, and state invariants run before the workflow accepts
    the imported state.

## Collective failure protocol

Every potentially throwing local phase is wrapped by one phase guard:

1. catch the local exception without leaving the phase;
2. `MPI_Allreduce` a failure vote;
3. select the lowest failing rank;
4. broadcast a bounded diagnostic from that rank;
5. throw the same collective error on all ranks;
6. enter the next collective only when no rank failed.

| Phase | Local work that may fail | Failure vote before next collective | Next collective | Empty-rank behavior |
|---|---|---|---|---|
| File discovery | open first member, read validated file count, derive paths | yes | path-length/path broadcast | non-root contributes no local discovery work |
| Assigned inspection | open assigned files, validate HDF5, hash, serialize fragment | yes | metadata gather | ranks with no file emit an empty fragment |
| Manifest assembly | decode fragments, cross-file checks, supplied-manifest checks | yes | manifest broadcast | non-root returns an empty local assembly object |
| Manifest decode | allocation, JSON parse, digest calculation | yes | digest reductions | all ranks decode the same bounded payload |
| Chunk read/conversion | hyperslab reads, dimensional conversion, scientific checks | yes | coverage reductions | non-reader ranks hold an empty record vector |
| Owner/serialization | owner calculation and explicit wire encoding | yes | send-layout phase | zero-record ranks create empty per-peer buckets |
| Send layout | count/displacement overflow checks | yes | `MPI_Alltoall` counts | zero-count peers remain represented |
| Receive layout | received count/displacement checks | yes | allocation phase | zero receive counts are valid |
| Exchange allocation | flattened buffers and copies | yes | `MPI_Alltoallv` payload | zero-byte buffers use null/ignored pointers |
| Wire decode | length, species, field, owner, and sidecar validation | yes | source/final audit exchange | empty inbound payload decodes to no records |
| Reconciliation | exact source/final balance bucket preparation and decode | yes | audit exchanges/reductions | empty source or final sets are represented exactly |
| State append | resize and populate authoritative SoA/sidecars | yes | next chunk or finalization | empty append is valid |
| Finalization | species index, gas/tracer identity remap, ownership invariants | yes | global duplicate audit | empty owner ranks remain valid local states |
| Global duplicate audit | bounded sorted-run creation and external merge of validator partitions | yes | global count/species/mass reductions | empty validator partitions create no runs and still participate |
| Global totals | compensated local mass accumulation and global count/species/mass reductions | yes | workflow acceptance boundary | zero-particle ranks contribute exact neutral values |

MPI calls themselves are checked through their established runtime contract; the
phase protocol prevents C++ exceptions before a later collective from stranding
peer ranks.

## Exact coverage and identity

Coverage is not inferred from totals alone:

- every expected source file must contribute exactly one manifest fragment;
- every deterministic chunk token must be identical on all ranks and have
  exactly one reader;
- for every chunk, source and final IDs are hash-partitioned with signed source
  and final counts; each unique ID must balance as `(1,1)`;
- the final authoritative ID set is hash-partitioned in bounded rounds; each
  validator writes sorted temporary runs and performs a bounded-memory external
  merge, so duplicates across rounds are detected exactly without retaining the
  full validator partition in RAM.

This detects missing, duplicated, substituted, corrupted, unowned, and multiply
owned records without all-gathering all IDs onto every rank.

## Scientific and sidecar validation

Before routing and again after wire decode, each record must satisfy:

- finite x/y/z in the periodic box;
- finite velocity components;
- nonzero ID and finite positive mass;
- finite non-negative gas density and specific internal energy;
- valid stellar formation epoch, birth mass, and metallicity;
- finite positive black-hole subgrid mass and finite non-negative accretion rate;
- valid tracer parent, host, fraction, and mass-history fields.

The final rank must equal the deterministic x-slab owner. Required gas, star,
black-hole, and tracer sidecars are transmitted in the same wire record as the
particle and constructed exactly once.

## Bounded-memory and accounting contract

`mode.ic_chunk_particle_count` bounds HDF5 hyperslabs.
`mode.ic_staging_particle_count` bounds routing and exact validation rounds.
The global duplicate-ID audit spills sorted runs to rank-local temporary storage
and merges them with bounded in-memory state; temporary files are removed on
success or failure. The final owner-local authoritative state is not transient
staging.

The report separates:

- logical metadata bytes read;
- complete source bytes read for SHA-256;
- particle payload bytes actually read;
- converted payload bytes;
- serialized bytes;
- sent and received bytes;
- manifest metadata bytes communicated;
- peak rank-local transient staging bytes.

Peak staging uses vector capacities and includes coordinates, velocities,
masses, IDs, every species field, `ParticleRecord` storage, nested per-peer
buckets, flattened exchange buffers, count/displacement arrays, decoded records,
coverage data, and exact ID-reconciliation workspaces.

Prohibited designs include complete source hashing on rank zero, root/global
particle gathers, global-state broadcasts, all-rank global ID materialization,
unbounded outgoing-record construction, and legacy replicated-state compaction
of an already partitioned import.

## Wire safety

The distributed record is encoded field-by-field with fixed-width integers,
IEEE-754 values, and explicit little-endian order. Padding-sensitive raw C++
objects are never transmitted. Count and displacement arithmetic fails before
an `int` overflow can reach MPI.

## Acceptance evidence

Registered MPI tests cover one, two, and four ranks with direct-bridge,
canonical, and supplied-manifest ingestion. The matrix includes a two-file mixed
Gas/DM/Star/PartType5-BH fixture, alternate h/a/frame/velocity conversions,
PartType2-to-tracer and PartType3-to-star policy paths, exact duplicate
rejection, more ranks than files, growing-global-count staging checks,
local-authoritative-state capacity checks, deliberate lost and duplicated route
mutations, and one-rank-only fault injection at owner/serialization, send
layout, payload validation, deserialization, sidecar append, finalization, and
source/final reconciliation. Each fault test has a strict timeout.

Those tests are source-complete but were not executable in the repair
environment because `MPI_CXX` was unavailable. The capability therefore remains
provisional until the `mpi-hdf5-fftw-debug` preset passes on the user system.
