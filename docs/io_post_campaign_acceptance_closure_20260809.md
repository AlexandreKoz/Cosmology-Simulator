# I/O Post-Campaign Acceptance Closure — 2026-08-09

This note records the focused repair that followed the adversarial acceptance re-audit of the
research-grade I/O campaign. It is intentionally narrower than the preceding persistence
campaign: the v6 science-snapshot, v23 restart, transactional-file, conversion, streaming,
compact-wire, parser-hardening, and provenance architecture were preserved.

## Logical snapshot-set integrity

CHUÍ multifile snapshots are now treated as one scientific object rather than as a directory
of plausible HDF5 files. Set inspection compares common file kind, dialect, schema/version,
generation, epoch/redshift, box dimensions, cosmology, unit/frame/velocity conventions,
config hash, naming-governance versions, `NumFilesPerSnapshot`, and 64-bit logical type counts.
Member indices must form exactly `[0, NumFilesPerSnapshot)`. Discovery is generation/stem aware,
so unrelated HDF5 snapshots in one directory are not silently merged.

For CHUÍ-authored multifile output, `<generation>.complete` is a bounded versioned
`chui_snapshot_set_v2` manifest. It binds exact member filenames, indices, local counts, file
sizes and per-member SHA-256 digests to the common scientific identity. A SHA-256 over the
canonical manifest body binds the set metadata itself. The manifest is transactionally
published only after every member has been published and member-integrity evidence exists.
Normal CHUÍ reads verify member content hashes and fail closed on missing, stale, mixed,
reindexed or scientifically inconsistent members.

## Bounded multifile materialization

`SnapshotReadBudget` now applies to the logical snapshot rather than independently resetting
for each file. Preflight and merge accounting enforce maximum member count, particles, gas
cells, sidecar rows, datasets/attributes and cumulative materialized bytes with checked
arithmetic. The merged state is re-measured through the core memory-accounting API before it
is returned.

## Missing fields and readiness

The explicit missing-field policy is applied across gas, stars, tracers and black-hole CHUÍ
extensions. Values that were previously plausible silent inventions (for example star
metallicity, formation time, birth mass/identity and tracer ancestry/exchange state) are now
rejected, marked unavailable, or filled only under an explicit documented default policy. All
reconstructions/defaults remain visible in the report.

`analysis_ready` and `evolution_ready` are distinct. Evolution readiness now requires:

- no unavailable scientific fields;
- nonzero, unique persistent particle IDs;
- valid state ownership/sidecar invariants;
- every gas cell to have a resolvable authoritative patch owner or a resolvable local parent
  particle.

A scientifically readable snapshot can therefore be analysis-ready without being falsely
certified as safe to resume evolution. This does not turn science snapshots into restart
checkpoints.

## Persistent IDs

`SimulationState::validatePersistentParticleIds()` is the persistence-boundary predicate. It
requires every persistent particle ID to be nonzero and globally unique in the local state.
Algorithms that only need uniqueness may continue using the narrower uniqueness predicate.

## Independent HDF5 validation

`validateSnapshotSetHdf5()` independently validates HDF5 structure without reconstructing the
normal `SimulationState`. It checks logical-set completeness/identity, particle group presence,
required dataset rank/type/shape/counts, finite phase-space/science values, positive masses or
valid `MassTable` fallback, coordinate bounds for CHUÍ-authored files, required native
full-physics extension datasets, and global nonzero/unique particle IDs. Validation streams
bounded chunks. CHUÍ-native v6 storage widths remain exact; supported external
AREPO/GADGET-family floating/integer datasets may use 32- or 64-bit scalar storage and are
converted by HDF5 for validation.

Runtime science output preserves writer-to-reader round-trip verification and additionally
runs the independent validator over the completed logical set.

## Hydro IC materialization

`materializeRootHydroPatchIfMissing()` no longer exposes a partially initialized nonempty
`PatchSoa` to gas-cell identity synchronization. Cartesian layout and complete patch geometry,
dimensions and ownership are established before identity mirrors are committed and validated.
The `integration_initial_condition_manifest_runtime` regression therefore closes without
weakening `PatchSoa::isConsistent()` or suppressing ownership errors.

## Transactional publication on Windows

The Windows path no longer passes the unsupported `REPLACEFILE_WRITE_THROUGH` flag to
`ReplaceFileW`. Existing destinations use supported replacement semantics; new/racing
publication uses `MoveFileExW` with `MOVEFILE_WRITE_THROUGH` when durable publication is
requested. Durable mode explicitly flushes the temporary file before publication and the
published file afterward through `FlushFileBuffers`. The API/documentation does not claim a
POSIX-equivalent parent-directory `fsync` semantic on Windows. A dedicated Windows CI job
builds and runs the bare-path/existing-destination transactional test. Local Linux validation
cannot substitute for that Windows runtime result.

## IC integrity and MPI transport

`IcIntegrityMode` makes source-file integrity cost explicit:

- `kVerifiedIdentity` (default): one authoritative inspection SHA-256 pass plus stable file
  identity validation around ingestion;
- `kStrictFullRehash`: the above plus one completion content SHA-256 pass.

No metadata-only mode is mislabeled as cryptographic integrity.

The successful distributed routing contract is now version 3 at 20 communicator-wide calls per
routing batch, down from 23 in the preceding campaign and 30 historically. Post-exchange exact
reconciliation work is grouped into fewer consensus boundaries while exact count/ID/ownership
checks remain intact. Runtime MPI evidence is delegated to MPI-capable CI when the local
execution environment lacks MPI.

## MPI snapshot acceptance and warning governance

CMake now registers dedicated two- and four-rank logical snapshot-set integration tests. They
write one member per rank, publish the set manifest, validate global/local counts and IDs,
run the independent validator and logical reader, and prove a missing-member set fails closed.
The MPI HDF5 CI selection includes both tests.

`COSMOSIM_STRICT_IO_WARNINGS` applies a target-scoped first-party warning-as-error policy to
`cosmosim_io`; CI builds it with GCC and Clang. A separate Windows transactional CI job covers
the native publication branch. Third-party code is not pulled under the first-party warning
policy.

## Maintainability boundary

The closure campaign continues the responsibility split with dedicated snapshot-set,
independent-validator and readiness components. Required CHUÍ gas/star/tracer/BH native field
names, storage widths and finite/metallicity validation flags are centralized in
`src/io/internal/snapshot_field_contract.hpp`, shared by native dataset creation and the
independent validator. This reduces one of the highest-risk schema-drift paths without a broad
aesthetic rewrite of legacy restart/IC monoliths.

Large legacy implementation units still exist, but further splitting is ordinary maintenance
debt rather than a known persistence-contract acceptance defect. New snapshot validation and
set-integrity logic no longer accumulate inside `snapshot_hdf5.cpp`.

## Local validation boundary

The closure environment provided GCC 14.2 and HDF5 1.14.5 but no MPI development/runtime stack
and no Windows host. Locally validated paths include:

- HDF5-enabled `cosmosim_io` and `cosmosim_workflows` builds;
- HDF5-disabled `cosmosim_io` build;
- target-scoped strict GCC warning-as-error `cosmosim_io` build;
- focused snapshot/restart/IC/transaction/hydro-manifest tests.

MPI runtime and native Windows runtime behavior are not claimed as locally executed; dedicated
CI gates are part of the repository for those environments.
