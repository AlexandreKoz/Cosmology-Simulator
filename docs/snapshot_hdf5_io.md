# GADGET/AREPO-Compatible HDF5 Snapshot I/O

## Scope

CosmoSim now provides a dedicated snapshot HDF5 reader/writer for a GADGET/AREPO-compatible
schema while keeping internal naming and ownership rules authoritative. External names are
centralized in `GadgetArepoSchemaMap` and are not scattered through solver modules.

## Schema and provenance

- Canonical groups are centralized in `GadgetArepoSchemaMap`:
  - `/Header`
  - `/Config`
  - `/Provenance`
  - `/PartType0..5` with optional `/ParticleType0..5` aliases
- Header group `/Header` stores canonical GADGET/AREPO attributes, including:
  - `NumPart_ThisFile`, `NumPart_Total`, `NumPart_Total_HighWord`, `MassTable`
  - `Time`, `Redshift`, `BoxSize`, `Omega0`, `OmegaLambda`, `OmegaBaryon`, `HubbleParam`
- Snapshot schema metadata:
  - `CosmoSimSchemaName="gadget_arepo_v1"`
  - `CosmoSimSchemaVersion=1`
- Config metadata group:
  - `/Config` attribute `normalized` containing normalized text config dump
- Provenance metadata group:
  - `/Provenance` attributes for compiler/build/git/config hash/timestamp/hardware

## Dataset aliasing and compatibility

- Canonical export writes `/PartTypeX` groups.
- Optional alias hard links `/ParticleTypeX` are written for tolerant downstream readers.
- Reader accepts aliases for common fields:
  - `Coordinates|Position|POS`
  - `Velocities|Velocity|VEL`
  - `Masses|Mass`
  - `ParticleIDs|ParticleID|ID`

## Safety and determinism

- Writer uses a `*.tmp` path and atomic rename finalize in the same directory.
- Dataset chunking is explicit through `SnapshotIoPolicy::chunk_particle_count`.
- Compression is optional and surfaced through policy.

## Conservative assumptions

- Current writer serializes particles by species mapping:
  - gas -> PartType0
  - dark matter -> PartType1
  - tracers -> PartType3
  - stars -> PartType4
  - black holes -> PartType5
- Tracer-sidecar provenance fields are stored in `/PartType3` when present:
  - `TracerParentParticleID`
  - `TracerInjectionStep`
  - `TracerHostCellIndex`
  - `TracerMassFractionOfHost`
  - `TracerCumulativeExchangedMassCode`
- `MassTable` is currently written as zeros and per-particle masses are emitted.
- Parallel HDF5 is not yet implemented, but APIs are structured to enable backend extension.

## Feature-flag behavior and test coverage

- `COSMOSIM_ENABLE_HDF5=ON`:
  - `test_unit_snapshot_hdf5_schema` validates centralized canonical names.
  - `test_integration_snapshot_hdf5_roundtrip` performs a real write+read round trip and asserts:
    - state payload integrity,
    - schema name/version,
    - normalized config text preservation,
    - provenance field preservation (`schema_version`, `git_sha`, `enabled_features`,
      `config_hash_hex`),
    - dataset alias resolution reporting.
- `COSMOSIM_ENABLE_HDF5=OFF`:
  - snapshot writer and reader are intentionally unavailable and throw
    `COSMOSIM_ENABLE_HDF5=OFF` runtime errors.
  - OFF-mode tests validate graceful unavailability only; they are not treated as I/O
    correctness validation for the HDF5 path.
