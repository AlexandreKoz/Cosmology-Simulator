# CHUÍ science snapshots and AREPO/GADGET HDF5 interoperability

## Scope

Science snapshots are analysis/interchange products. They are deliberately separate from
restart checkpoints and from canonical IC bundles. The active native schema is
`chui_science_snapshot_v6` (version 6).

The public I/O contract separates four concepts:

- CHUÍ-authored science snapshots (`readCosmoSimScienceSnapshotHdf5`);
- AREPO Format-3 export semantics (`SnapshotDialect::kArepoFormat3`);
- GADGET-4 HDF5 export/import semantics where the same field contract is unambiguous
  (`SnapshotDialect::kGadget4Hdf5`);
- explicit external import (`importExternalSnapshotHdf5`), which requires a dialect and
  rejects ambiguous populated PartType2/PartType3 unless the caller supplies a species map.

The legacy `writeGadgetArepoSnapshotHdf5`/`readGadgetArepoSnapshotHdf5` names remain as
compatibility entry points, but they delegate to the versioned CHUÍ science-snapshot
implementation. They no longer define a fictional merged GADGET/AREPO dialect.

## HDF5 shape

Canonical groups are:

- `/Header`
- `/Config`
- `/Parameters`
- `/Provenance`
- `/Units`
- `/PartType0..5` when populated

Canonical phase-space datasets retain external spellings: `Coordinates`, `Velocities`,
`Masses`, and `ParticleIDs`. Optional `/ParticleTypeN` hard-link aliases may be emitted by
policy, but `/PartTypeN` is the native science layout.

The header records `NumPart_ThisFile`, 64-bit logical totals through
`NumPart_Total`/`NumPart_Total_HighWord`, `NumFilesPerSnapshot`, epoch/cosmology, schema,
dialect, generation/member identity, naming-rule versions, and the explicit local-index
policy. CHUÍ currently uses 32-bit dense local particle/cell indices, so one member is capped
at `UINT32_MAX` rows per PartType; logical multifile totals remain 64-bit.

## Cosmological storage semantics

For comoving production runs the default export dialect is AREPO Format 3. Conversion is
centralized in `src/io/snapshot_conversion.cpp`; solver code does not scatter `a`/`h`
conversions.

For the current code-unit convention the stored external values are:

| Quantity | CHUÍ internal semantic | AREPO/GADGET-family stored semantic |
|---|---|---|
| Coordinates | comoving code coordinate | internal coordinate multiplied by `h` |
| Velocities | physical peculiar code velocity | peculiar velocity divided by `sqrt(a)` |
| Masses | code mass | internal mass multiplied by `h` |
| Density | CHUÍ comoving code density | internal density divided by `h^2` |
| Pressure | CHUÍ comoving code pressure | internal pressure divided by `h^2` |
| InternalEnergy | specific internal energy in configured velocity-squared code units | same numerical code-unit value |
| GravitySofteningComoving | comoving code length | same conversion as Coordinates |

The inverse conversion is applied by the reader. The `/Units` and `/Header` metadata record
the selected dialect and configured unit names so a CHUÍ-authored file is self-describing.
Isolated physical-coordinate runs default to `chui_native` rather than pretending to be an
AREPO cosmological snapshot.

## Distributed logical snapshots

MPI science output is one logical snapshot set, not independent rank-local science files.
The current scalable topology writes one member per MPI rank without gathering global state
to rank zero:

```text
<shared_run>/snapdir_042/
  <stem>_042.0.hdf5
  <stem>_042.1.hdf5
  ...
  <generation_id>.complete
```

Each member has its own `NumPart_ThisFile`, the same global totals and
`NumFilesPerSnapshot`, and common schema/epoch/generation metadata. After all members have
been transactionally published and collectively accepted, root writes a bounded completion
marker binding generation ID, member count, and global counts. CHUÍ multifile reads reject
sets whose marker is absent or disagrees with the HDF5 members.

`inspectSnapshotSet` performs shallow discovery/consistency checks without materializing the
full payload. Standard external multifile sets may be imported when the explicit dialect and
species mapping are sufficient.

## Scientific fields

PartType0 is emitted from authoritative dense gas-cell state. Stable gas identity is preserved
through `GasCellIDs`, optional parent lineage, validity masks, and owning-patch identity.
Current science fields include InternalEnergy, Density, Pressure, Metallicity,
StarFormationRate/effective-ISM diagnostics where applicable, plus CHUÍ temperature and
sound-speed diagnostics.

PartType4 includes formation time, metallicity, birth mass, immutable birth identity, and
current CHUÍ stellar-evolution cumulative diagnostics.

PartType5 now preserves the current black-hole science sidecar. Established phase-space names
remain canonical; model-specific quantities use explicit `CHUI_` extension names such as
`CHUI_BHSubgridMass`, accretion rate, Eddington ratio, feedback/cumulative energy, duty-cycle
accounting, and host-cell relation. They are not mislabeled as canonical AREPO fields when
CHUÍ semantics are project-specific.

PartType3 is CHUÍ tracer state only for CHUÍ-authored files. Arbitrary external PartType2/3
semantics are not guessed.

## Missing-field and hostile-input policy

External imports use `SnapshotMissingFieldPolicy`. The default is fail-closed `Reject`.
`MarkUnavailable` returns analysis-only state with availability reporting; generated/default
values require explicit policy and cannot silently masquerade as evolution-ready state.
`SnapshotReadResult::requireEvolutionReady()` enforces this boundary.

Read budgets bound particle count, materialized bytes, individual dataset bytes, and attribute
bytes. Dataset rank/dimensions/types are checked before allocation. Imported state is checked
for finite phase space, nonzero IDs, mass/species invariants, and CHUÍ gas constraints.

## Transaction and storage policy

Snapshot members use the shared transactional-file layer: unique same-directory temporary,
HDF5 completion/close, optional durable file sync, atomic publish/replace, and directory sync
where the platform provides it. The writer never deletes the previous valid final file before
publishing the replacement.

Producer-side staging is bounded by `SnapshotIoPolicy::chunk_particle_count` and HDF5
hyperslab streaming; HDF5 chunking is not confused with full-species in-memory staging.
Reports record peak staging bytes, logical bytes written, chunk writes, and storage properties
actually inspected from the HDF5 datasets. Compression/chunk values are never fabricated as
zero merely because the reader did not look.

## Validation boundary

The production workflow preserves its CHUÍ writer-to-reader scientific roundtrip, but this is
not treated as proof of external semantics. Snapshot-set inspection independently checks HDF5
structure/counts/epoch/set completeness without reconstructing the full simulation state.
MPI and Windows runtime acceptance remain dependency/platform-specific and must not be
claimed when those environments were not executed.
