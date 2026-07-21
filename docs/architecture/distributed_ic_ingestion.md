# Distributed initial-condition ingestion

## Decision

External cosmological IC ingestion is a workflow-owned, explicitly configured,
multifile, bounded-staging pipeline. No MPI rank receives an implicit license to
read or materialize the complete universe.

## Layering

- `core`: typed `InitialConditionConvention` and Type2/Type3 species policies,
  validation, and normalized provenance. It has no dependency on `io`.
- `io`: strict audit manifest, file-set inspection, canonical conversions,
  sidecar construction, explicit wire format, serial/distributed readers, and
  scientific validation.
- `workflows`: convention dispatch, config-relative path resolution, borrowed
  MPI/profiler services, manifest output, and the `already_partitioned` startup
  boundary.
- `parallel`: process MPI context and later runtime migration/rebalancing. It
  does not force initial IC replication.

## Serial phases

1. Discover the deterministic source-member sequence.
2. Inspect every header and selected dataset.
3. Validate file-set cosmology, epoch, schema, counts, and 64-bit totals.
4. Build/validate the strict manifest.
5. Allocate only the final serial authoritative state.
6. Read and convert bounded chunks.
7. Construct all required species sidecars.
8. Sort IDs for exact duplicate detection and validate state invariants.

## Distributed phases

1. Rank zero inspects the file set and produces the canonical manifest.
2. The strict manifest is broadcast; every rank reconstructs and validates the
   same file schema and agrees on its SHA-256 digest.
3. A global ordered chunk index assigns each source chunk to exactly one reader
   rank (`chunk_index mod world_size`). More ranks than files are valid.
4. The reader converts a bounded chunk once and packs generic particle state and
   required sidecars into a fixed, explicit little-endian wire record.
5. The initial owner is evaluated directly from the canonical comoving x
   coordinate and runtime slab decomposition. Bounded `Alltoallv` rounds route
   records to final owners; no root/global state exists.
6. Each rank appends only received local records and rebuilds local species and
   sidecar indices.
7. IDs are hash-partitioned in bounded rounds for exact duplicate detection over
   the full unsigned 64-bit domain.
8. Collective reductions validate source/chunk coverage, global/species counts,
   global/species mass totals, ownership completeness/exclusivity, manifest and
   conversion agreement, and state invariants.

A rank-local failure is reduced before the next collective phase. This prevents
one rank throwing while peers wait in a later collective.

## Bounded-memory contract

`mode.ic_chunk_particle_count` bounds source hyperslabs.
`mode.ic_staging_particle_count` bounds routing and distributed validation
rounds. Import reports measure peak staging bytes. Final local state is not
counted as transient staging, but it is required to contain only records owned
by that rank.

Prohibited designs include:

- every rank reading all files;
- a root rank gathering all particles or sidecars;
- broadcasting a complete simulation state;
- building all outgoing records for an unbounded universe before communication;
- collecting all global IDs on every rank;
- calling the legacy replicated-state ownership compactor on already partitioned
  external state.

## Wire safety

The distributed record is encoded field-by-field with explicit widths and
little-endian byte order. It does not send padding-sensitive raw C++ objects.
The record bundles species tags and all required sidecar values so routing cannot
separate a particle from its authoritative auxiliary state.

## Validation and evidence

Required evidence comprises:

- serial single-file versus multifile equivalence;
- converter versus direct bridge equivalence;
- one-, two-, and four-rank reconstruction equivalence;
- h/scale-factor/frame/velocity conversion cases;
- gas, dark matter, star, black-hole, and explicit Type2/Type3 policy cases;
- malformed/inconsistent file sets;
- duplicate IDs within a file, across files, across reader ranks, and across
  final owner ranks;
- instrumentation showing exactly-once chunk reads, exactly-once ownership, and
  bounded peak staging.

MPI tests are compiled and run only when an MPI implementation and launcher are
available. A non-MPI build must not report distributed ingestion as supported.

## Canonical conversion boundary

`tools/convert_ic.py` uses the same manifest vocabulary and conversion equations
as runtime ingestion. It streams a validated external file set to one canonical
CHUÍ v1 HDF5 artifact plus a strict audit manifest, hashes source members,
detects duplicate unsigned 64-bit IDs exactly, and atomically finalizes output.
It can consume either an explicit bridge convention plus conversion parameters
or an existing version-2 audit manifest; manifest mode revalidates source order,
SHA-256 hashes, byte sizes, scientific parameters, species policies, and actual
HDF5 schemas before writing. The canonical output remains an interchange IC,
not a restart checkpoint.
