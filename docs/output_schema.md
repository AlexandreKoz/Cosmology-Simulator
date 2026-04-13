# Output schema guide

This guide documents snapshot, restart, and provenance payload contracts currently implemented in CosmoSim.

Authoritative interfaces:

- `include/cosmosim/io/snapshot_hdf5.hpp`
- `include/cosmosim/io/restart_checkpoint.hpp`
- `include/cosmosim/core/provenance.hpp`
- `include/cosmosim/core/profiling.hpp` (operational run-event report)

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

## 2.1) Field ownership table (snapshot vs restart)

| Ownership | Fields |
|---|---|
| Shared metadata contract | normalized config text/hash, provenance payload, schema identity |
| Snapshot-only (interoperable science output) | `/Header` cosmology attrs, `/PartTypeN` particle datasets, read aliases (`Position`, `VEL`, `ID`, etc.) |
| Restart-only (exact continuation state) | full `SimulationState` hot/cold lanes, `StateMetadata`, module sidecars + schema versions, `IntegratorState`, scheduler persistent state (`bin_index`, `next_activation_tick`, `active_flag`, `pending_bin_index`), payload integrity hashes |

Snapshot and restart intentionally remain separate contracts: snapshot is analysis/interchange oriented;
restart is execution-resume oriented.

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
- Reference workflow operational report (`reference_operational_events.json`) is schema-tagged (`schema_version = 1`) and includes:
  - run label,
  - provenance config hash linkage (`provenance_config_hash_hex`),
  - severity-count summary and status,
  - structured events (`event_kind`, `severity`, `subsystem`, optional step/time/scale context, message, payload map).

## 5) Change procedure for schema-affecting work

When changing snapshot/restart/provenance fields:

1. Update the corresponding interface headers and implementation.
2. Update `docs/configuration.md` if config keys or normalized text semantics changed.
3. Update validation expectations in `docs/validation_plan.md`.
4. Add/update tests in `tests/unit` + `tests/integration` + `tests/validation` as applicable.
5. Record rationale in `docs/architecture/decision_log.md`.

## Compatibility notes (2026-04-13)

- No external snapshot dataset names were changed.
- Restart schema version/name were not changed (`cosmosim_restart_v2`, version `2`).
- Restart contract enforcement was tightened: missing continuation-critical metadata
  now fails fast with explicit errors instead of producing weak checkpoints.
- New operational diagnostics output is additive only; snapshot/restart/provenance schemas and compatibility contracts are unchanged.
