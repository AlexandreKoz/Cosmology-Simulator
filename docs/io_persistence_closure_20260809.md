# I/O Persistence Closure — 2026-08-09

This note records the persistent-format and runtime I/O changes introduced by the
research-grade I/O closure campaign. It is a migration note, not a replacement for
`docs/snapshot_hdf5_io.md` or `docs/restart_checkpointing.md`.

## Active persistent contracts

| Product | Active identity | Compatibility stance |
|---|---|---|
| Science snapshot | `chui_science_snapshot_v6` / version 6 | CHUÍ-authored identity with explicit AREPO Format-3 or GADGET-4 HDF5 export semantics; legacy CHUÍ snapshots remain readable where their semantics are known. |
| Restart checkpoint | `cosmosim_restart_v23` / version 23 | v23 adds canonical SHA-256 integrity and transactional publication; documented legacy restart schemas remain readable. |
| Distributed IC routing wire | wire v2 | Versioned common record plus species-specific payload; legacy v1 sizing is retained only for historical report interpretation. |

## Science snapshot semantic boundary

Science snapshots no longer use a single ambiguous "GADGET/AREPO" semantic label.
The writer selects an explicit dialect. Cosmological AREPO/GADGET-family output uses
one centralized conversion layer for scale-factor, Hubble-parameter, frame, and unit
semantics. In particular, the stored `Velocities` representation is converted from
CHUÍ's internal peculiar-velocity state and inverted on import rather than being copied
byte-for-byte.

The snapshot schema records the active dialect and unit/frame policy. Ambiguous
populated external particle types are rejected unless the caller supplies an explicit
species mapping. A partial external snapshot cannot silently become evolution-ready:
missing-field handling is a caller policy and field availability is reported.

## Logical multifile snapshots

MPI science output is one logical snapshot set rather than unrelated rank-local science
files. A typical set is:

```text
snapdir_042/
  snap_042.0.hdf5
  snap_042.1.hdf5
  ...
  snap_042.complete
```

Each member stores its local `NumPart_ThisFile`, common 64-bit logical totals using the
GADGET-family low/high-word convention where required, a common `NumFilesPerSnapshot`,
and common generation/schema/dialect/epoch/box/cosmology/unit/config/governance identity.
The completion product is now a versioned `chui_snapshot_set_v2` manifest binding exact member
filenames, contiguous indices, local/global counts, file sizes and per-member SHA-256 digests
plus a root manifest digest. A CHUÍ set is accepted only when all members and the manifest agree.

The current in-memory sidecar model still uses 32-bit dense local row indices. Therefore
**one member is intentionally limited to at most `UINT32_MAX` rows of a PartType**. This
is now an explicit implementation contract (`CHUILocalIndexWidthBits=32`) rather than
an implicit contradiction. Logical multifile totals remain 64-bit.

## Bounded-memory snapshot output

The science writer creates HDF5 datasets once and fills them through bounded hyperslab
chunks. Coordinates, velocities, IDs, masses, and sidecars are no longer materialized as
whole-species duplicate arrays. The write report records logical bytes, chunk writes, and
peak staging bytes. HDF5 chunk/compression reporting on read is obtained from dataset
creation properties rather than fabricated defaults.

## Full-physics science state

The v6 snapshot contract carries current gas, stellar, tracer, softening, and black-hole
science state. Model-specific BH quantities are deliberately named as CHUÍ extensions
instead of being mislabeled as canonical AREPO fields. Restart-only scheduler, RNG,
force-cache, decomposition, and AMR bookkeeping remain restart-only.

## Transactional publication

Snapshots, restarts, manifests, and snapshot completion markers now share one
transactional publication primitive:

1. create a unique same-directory temporary sibling;
2. write and finish library metadata;
3. close all HDF5/library handles;
4. optionally durable-sync the temporary file;
5. atomically publish/replace the destination;
6. durable-sync parent-directory metadata where the platform supports the contract.

There is no remove-old-before-rename fallback. Bare filenames are valid. Temporary
suffixes cannot resolve to the destination. POSIX durable publication uses file and
parent-directory synchronization. The guarded Windows implementation uses native
replace/flush APIs; it was not runtime-executed in the Linux campaign environment.

## Restart integrity migration

`cosmosim_restart_v23` adds canonical little-endian SHA-256 integrity identified as
`sha256-canonical-le-v1`. The legacy FNV-1a value is retained only for compatibility and
migration diagnostics. Both digests are produced during one logical payload traversal;
the checkpoint state is not walked a second time merely to obtain a second rendering of
the same digest.

## IC ingestion changes

- Default `kVerifiedIdentity` IC ingestion performs one authoritative inspection SHA-256 pass
  and stable source-file identity validation around ingestion. Explicit
  `kStrictFullRehash` adds a completion content rehash when the extra I/O is justified.
- The distributed wire is v2 and no longer transports all gas/star/BH/tracer sidecars in
  every dark-matter record.
- Successful routing bookkeeping was consolidated from the historical 30, then 23, to 20
  communicator-wide calls per routing batch while retaining exact reconciliation and
  duplicate-ID auditing.
- MPI wrappers request `MPI_ERRORS_RETURN`, check wrapped operation statuses, and treat a
  communicator-level MPI error as fatal instead of attempting another consensus on a
  possibly failed communicator.
- Manifest JSON parsing now has file/string/depth/container budgets, Unicode escape
  handling, control-character validation, and checked integer narrowing.
- Canonical-bundle metadata reads are bounded and path containment/special-file policy is
  explicit.
- Exact distributed duplicate-ID external-sort scratch can be rooted in a configured
  run-scoped scratch directory and is cleaned on normal/exceptional scope exit.

## Provenance/governance

Science output records schema/dialect identity, normalized configuration, configuration
hash, runtime provenance, naming-rule versions, member/set identity, and explicit local
index policy. Mandatory configuration truth is validated instead of silently recording an
empty value.

## Validation boundary

The campaign environment provided GCC and HDF5 but not an MPI runtime/development stack
or Windows. Serial/HDF5 behavior was built and exercised locally. MPI and Windows code
was implemented behind the existing feature/platform guards and is intended for capable
CI/hosts; no runtime claim is made for those unavailable environments.

## Post-campaign acceptance closure

The follow-up acceptance-gap campaign strengthens logical-set identity, cumulative read
budgets, full-species missing-field handling, runtime evolution-readiness classification,
Windows publication semantics, direct HDF5 validation, MPI snapshot-set CI coverage, and
I/O strict-warning CI. See `docs/io_post_campaign_acceptance_closure_20260809.md` and
`docs/repair/io_post_campaign_closure_matrix_20260809.md` for the current closure evidence.
