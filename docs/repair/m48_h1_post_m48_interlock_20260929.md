# M48-H1 — Post-M48 interlock, compact-state compatibility, and ghost fallback closure

Date: 2026-09-29

## Scope

M48-H1 is a repair/hardening pass over the accepted M48-01 through M48-06 memory architecture. It does not change TreePM force equations, KDK ordering, timestep criteria, snapshot/restart schemas, or the compact-state authority model.

## Compile-regression repairs

- `attachSchedulerFieldsToParticleMigrationRecords()` now reads scheduler truth through `binIndex()`, `nextActivationTick()`, and `pendingBinIndex()` and no longer depends on a removed persistent-state materialization.
- the gas-cell time-bin mirror assertion reads scheduler truth through `binIndex()` rather than removed raw hot metadata.
- `PmSlabHaloExchangeResult` now carries `exchange_sequence`, preserving protocol-generation diagnostics after direct-to-cache halo staging.

## PM halo result and telemetry semantics

TreePM production halo receives remain direct-to-cache. `PmSlabHaloExchangeResult::left_halo/right_halo` remain legacy materialized lanes and may be empty in production. `exchange_sequence`, peers, depth, and byte counts describe the committed logical exchange. `TreePmDiagnostics::pm_halo_value_count` is defined as the total logical left-plus-right halo values for one scalar force component, derived from committed halo geometry rather than legacy temporary vectors.

## Compact DMO restart compatibility

Distributed restart ownership publication now expands logical ownership with `particleOwningRank(row)` for exactly `particles.size()` rows. Restart topology validation compares against logical particle count and validates each logical owner without requiring dense owner metadata.

Restart exact-equivalence compares logical time-bin, SFC key, species, flags, owner, drift time, and drift scale factor through representation-aware accessors. Compact and materialized states therefore compare by logical value rather than physical vector shape.

## Compact DMO analysis compatibility

`ParticleDiagnosticsView` and `HaloParticleView` now expose a homogeneous-DMO species policy. Diagnostics and FOF treat an empty dense species span as valid only when that explicit policy is set. Distributed FOF emits `DarkMatter` tags directly into its existing wire format and does not construct a full local species vector.

## Compact-state eligibility hardening

`compactHomogeneousDmoMetadata()` now fails closed unless canonical particle storage and sidecar extents are consistent, particle IDs cover the full population, and SFC storage is either absent or exactly population-sized. Nonuniform imported SFC values remain materialized and exact.

## Non-empty particle-ghost closure

The generic particle-ghost fallback now stores `LocalGhostDescriptor` records only for actual remote-owned ghost rows. Each demand record carries its canonical `local_index`, owner, stable particle ID, and epoch. The authoritative particle payload remains a borrowed read-only view.

Incoming request resolution no longer builds a hash for every locally owned particle. Only after actual incoming IDs arrive does the receiver allocate request-sized lookup state, then scan authoritative local IDs once while excluding the demand-scaled ghost-row set. Duplicate requests, multiply-resolved requested IDs, missing IDs, stale epochs, and ownership mismatches fail closed. Legacy dense descriptor callers remain supported through the sentinel `local_index` convention.

Thus metadata retained for non-empty ghost demand scales with actual ghost/request traffic rather than the complete local population; the state scan itself remains O(N) with O(demand) auxiliary storage.

## Borrowed homogeneous-DMO eligibility

Compact metadata is used as authority-level proof for homogeneous dark-matter species and local ownership instead of re-proving those facts row-by-row. Finite/non-negative mass and periodic-coordinate validation remains explicit. A generation-stamped eligibility certificate is reused while both `particleIndexGeneration()` and `gravitySourceGeneration()` remain unchanged; compact states with a materialized softening-override mask are not cached by this certificate.

## Deferred item

`TopLevelDomainHierarchy::ownersWithinCutoff()` still has the pre-existing per-target vector/stack allocation behavior. It is the intentionally deferred P2 performance item from M48-H1; changing its OpenMP scratch ownership safely is not required for correctness of this interlock repair.

M48-07 migration-packet redesign, M48-08 streaming restart verification, and M48-09 tree/PM stretch compaction remain out of scope.

## Validation

Commands executed:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4
```

The configure completed successfully. The build was resumed after command-window interruptions and ultimately completed all configured CPU-debug targets successfully. No CTest, MPI simulation, HDF5 roundtrip campaign, numerical comparison, or runtime RSS qualification was run.

Runtime/numerical/MPI qualification remains deferred to the dedicated qualification campaign.
