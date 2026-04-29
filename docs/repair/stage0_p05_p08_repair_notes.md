# Stage 0 P0-05..P0-08 repair notes

Status: implemented as a targeted runtime-truth hardening pass.

Scope: P0-05 gas cell identity, P0-06 softening ownership, P0-07 config-derived runtime values, and P0-08 active-set construction ownership. No new physics was added.

## P0-05 gas cell identity

Gas cells now have explicit identity lanes in `GasCellSidecar`:

- `gas_cell_id`
- `parent_particle_id`

For the current particle-bound gas-cell lane, these IDs are derived from canonical gas-particle ordering by `refreshGasCellIdentityFromParticleOrder()`. `debugAssertGasCellIdentityContract()` now checks that the gas-cell identity lanes match gas particle IDs and detects duplicate gas-cell IDs. `SimulationState::rebuildSpeciesIndex()` refreshes identity automatically when the current state has a one-to-one gas-particle/gas-cell layout. States that still use standalone hydro-cell tests without particle-bound gas are not silently promoted to particle-bound identity.

Restart/checkpoint state now serializes and hashes the gas-cell identity lanes. A restart missing these datasets is rejected rather than defaulting identity silently.

## P0-06 softening ownership

Softening now distinguishes between a numeric softening value and true per-particle override authority.

- `ParticleSidecar::gravity_softening_comoving` stores materialized values.
- `ParticleSidecar::has_gravity_softening_override` records whether a row owns a real override.
- `TreeSofteningView` now accepts source/target override masks.
- Softening resolution uses per-particle values only when no mask is supplied for legacy callers or when the mask marks the row as an explicit override.
- Species defaults and global fallback remain authoritative for unmasked rows.

Migration records preserve both numeric softening values and true override state, so species/default materialization is no longer confused with override authority.

## P0-07 config-derived runtime values

`DerivedRuntimeConfig` and `deriveRuntimeConfig()` define the typed owner for values derived from normalized config:

- time begin/end in code units
- axis-aware comoving box size
- TreePM grid shape
- per-species resolved gravity softening defaults
- normalized config hash and hash hex

This keeps raw config, normalized config, and derived runtime constants separate. Runtime modules should consume this derived object rather than independently recomputing equivalent constants.

## P0-08 active-set construction ownership

`ActiveSetDescriptor` now carries source authority/generation metadata:

- scheduler provenance flags
- state particle/cell generation at construction
- scheduler tick at construction

`makeSchedulerActiveSetDescriptor()` builds descriptors derived from scheduler output, and `debugAssertActiveSetDescriptorFresh()` catches stale descriptors after state generation changes or index invalidation. The orchestrator validates descriptors before stage callbacks consume them.

## Test coverage touched

Targeted tests were updated for:

- gas-cell identity refresh and stale hydro view detection
- softening override mask semantics and restart preservation
- typed derived runtime config ownership
- active-set descriptor provenance and generation invalidation

## Build note

The repair was configured with the CPU debug preset. In this sandbox, full C++ builds timed out before completion, so local verification should run the normal presets:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4
ctest --preset test-cpu-debug --output-on-failure -R "gas|softening|config|active|scheduler|time_integration"
ctest --preset test-cpu-debug --output-on-failure
```

If HDF5 is available, also run:

```bash
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug -j4
ctest --preset test-hdf5-debug --output-on-failure -R "restart|checkpoint|softening|gas|provenance"
```
