# Output schema guide

This guide documents snapshot, restart, and provenance payload contracts currently implemented in CosmoSim.

Authoritative interfaces:

- `include/cosmosim/io/snapshot_hdf5.hpp`
- `include/cosmosim/io/restart_checkpoint.hpp`
- `include/cosmosim/core/provenance.hpp`

## 1) Snapshot schema (GADGET/AREPO interoperability)

Current schema identity:

- `schema_name = gadget_arepo_v1`
- `schema_version = 1`

Logical groups:

- `/Header`
- `/Config`
- `/Provenance`
- `/PartType0` ... `/PartType5`

Canonical fields and accepted read aliases:

- `Coordinates` (`Coordinates`, `Position`, `POS`)
- `Velocities` (`Velocities`, `Velocity`, `VEL`)
- `Masses` (`Masses`, `Mass`)
- `ParticleIDs` (`ParticleIDs`, `ParticleID`, `ID`)

`SnapshotIoPolicy` controls compression, chunk size, and whether alias groups are emitted.

## 2) Restart schema

Current restart identity:

- `name = cosmosim_restart_v2`
- `version = 2`

Restart payload includes:

- `SimulationState`
- `IntegratorState`
- `HierarchicalTimeBinScheduler` persistent state (`TimeBinPersistentState`)
- normalized config text and normalized config hash
- provenance record
- payload integrity hash

Compatibility rule is explicit through `isRestartSchemaCompatible(version)`.

## 3) Provenance payload

`ProvenanceRecord` persists:

- schema tag (`provenance_v1`)
- build identity (`git_sha`, compiler id/version, build preset, feature flags)
- deterministic config hash
- UTC timestamp and hardware summary
- rank attribution (`author_rank`)

## 4) Naming and stability conventions

- Output stems are restricted to `[A-Za-z0-9_-]` for deterministic naming.
- Normalized config snapshots are written with run outputs.
- Schema names/versions are part of compatibility behavior and must be updated intentionally.

## 5) Change procedure for schema-affecting work

When changing snapshot/restart/provenance fields:

1. Update the corresponding interface headers and implementation.
2. Update `docs/configuration.md` if config keys or normalized text semantics changed.
3. Update validation expectations in `docs/validation_plan.md`.
4. Add/update tests in `tests/unit` + `tests/integration` + `tests/validation` as applicable.
5. Record rationale in `docs/architecture/decision_log.md`.
