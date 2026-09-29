# M48-06 — Uniform rung-zero scheduler and homogeneous canonical DMO state

Date: 2026-09-29

## Scope

M48-06 specializes two physical representations without creating new authorities or changing the numerical model. The generic hierarchical scheduler and generic/full-physics particle metadata remain available as exact fallbacks. The production `max_bin == 0` path uses a uniform rung-zero scheduler representation, and validated all-owned homogeneous DMO state may use compact canonical metadata.

No snapshot schema or restart logical schema is changed by this campaign.

## M48-06A — uniform scheduler representation

`HierarchicalTimeBinScheduler` now has an explicit representation policy:

- `kGenericHierarchical`: the historical per-element bin, activation, pending-transition, membership, active-sort and candidate lanes remain authoritative.
- `kUniformRungZero`: valid only for `max_bin == 0`; scheduler authority is the current tick, element count, substep state and low-cardinality diagnostics plus one explicit `uint32` identity active-index lane.

For the uniform representation the following population lanes are physically absent and their retained capacities are released:

- `bin_index`
- `next_activation_tick`
- `active_flag`
- `pending_bin_index`
- bin-membership storage
- position-in-bin storage
- active-sort scratch
- candidate-bin storage
- candidate-source storage

The identity active-index lane is rebuilt only when the population changes. `beginSubstep()` exposes that existing range and opens the logical substep; `endSubstep()` performs the rung-zero no-transition reconciliation and advances `current_tick` exactly once. Candidate criteria still execute, but their scheduler diagnostics are aggregate-only because no nonzero bin is representable.

Representation-aware accessors expose logical bin, next-activation, active and pending state without materializing generic arrays. Explicit `exportPersistentState()` remains a compatibility API and may materialize complete logical vectors when a caller deliberately asks for them; production restart serialization does not require that materialization for the uniform particle scheduler.

### Scheduler persistence

The restart schema remains explicit. In uniform mode `/scheduler/bin_index`, `/scheduler/next_activation_tick`, `/scheduler/active_flag`, and `/scheduler/pending_bin_index` are written using bounded constant-fill HDF5 staging rather than four N-sized vectors. Restart integrity is computed from the same logical field sequence and values as the materialized representation.

Historical restart arrays are still accepted. Import compresses a scheduler only when `max_bin == 0`, committed bins are zero, activation ticks equal the saved current tick, active flags describe a closed checkpoint boundary, and pending bins are unset/rung-zero compatible. Otherwise the generic representation is retained.

## M48-06B — homogeneous canonical DMO metadata

`SimulationState` now has an explicit particle metadata representation policy:

- `kMaterializedGeneric`
- `kHomogeneousDmo`

The compact policy is selected only after validation proves a gas/star/BH/tracer-free local DMO population with uniform flags, local ownership, rung-zero time-bin mirrors, and a common drift epoch.

### Dense canonical DMO truth

The stable compact DMO population keeps these mandatory dense lanes:

- `position_x_comoving`, `position_y_comoving`, `position_z_comoving`
- `velocity_x_peculiar`, `velocity_y_peculiar`, `velocity_z_peculiar`
- `mass_code`
- `particle_id`

That is the 64-byte-per-particle physical/identity floor before optional exception sidecars.

### Implicit/uniform logical metadata

For eligible compact DMO these values are represented as validated scalar/policy state rather than N independent values:

- particle scheduler/time-bin mirror = 0
- species = dark matter
- particle flags = one uniform scalar
- owning rank = local rank
- last drift time = one common scalar
- last drift scale factor = one common scalar

The KDK all-active drift updates the common drift epoch directly. A subset drift cannot silently make compact metadata heterogeneous; such a path must materialize explicitly before mutation.

### SFC policy

If imported/current SFC keys are uniform, their value is stored once and the dense SFC capacity is released. If SFC keys are nonuniform, the exact imported lane remains materialized while the rest of the homogeneous DMO metadata may still be compact. No Morton/SFC regeneration replaces imported truth.

### Particle species index

`ParticleSpeciesIndex` has an identity representation for homogeneous DMO. In that mode:

- dark-matter count = N
- all other species counts = 0
- `localIndex(global) = global`
- `globalIndex(DarkMatter, local) = local`

No `uint32[N]` global/local identity maps are retained. APIs that would imply hidden materialization of a dense identity vector fail explicitly; production consumers iterate the mathematical identity range or use logical accessors.

## Materialization and mutation boundaries

Population-scale metadata expansion is explicit. `SimulationState::materializeParticleMetadata()` is not triggered by a normal read accessor, and changing a compact population through generic resize fails closed until the caller establishes a governed materialization boundary.

Current generic migration remains a compatibility boundary. If a rebalance is actually required, the migration admission plan includes the worst-case compact-DMO metadata/species-index expansion before the central `MemoryGovernor` reservation is committed. Only then is generic metadata materialized. After a successful migration commit, compact DMO is restored only if the homogeneous invariants still hold.

Species/full-physics mutations likewise require explicit materialization before a heterogeneous row can become authoritative. Stable DMO operation does not oscillate between representations every step.

## Gravity, decomposition and hot-path consumers

The M48-05 borrowed homogeneous-DMO gravity path now reads logical metadata and does not require dense species, owner or drift lanes. Generic gravity fallback and prediction use representation-aware accessors where compact metadata may be present.

Runtime decomposition source views carry an explicit uniform-particle-metadata policy. DMO decomposition therefore does not create replacement species/owner identity vectors merely to satisfy span-shaped legacy interfaces.

## Snapshot and restart compatibility

Science snapshots retain the existing GADGET-style logical particle fields and PartType routing. Homogeneous DMO is routed directly to PartType1 without materializing `species_tag[N]`.

Restart files retain all historical logical datasets. Uniform canonical metadata and uniform scheduler lanes are expanded during serialization using bounded constant-value chunks. Logical integrity hashes are based on the values and historical field order, not on physical in-memory representation bytes.

The current full restart reader still reads explicit datasets. Complete integrity verification occurs before eligible metadata is compacted and redundant capacity released. This intentionally leaves the broader M48-08 readback/verification redesign deferred.

## Memory accounting result

Current ownership reporting remains capacity-based, so released vectors disappear naturally from `collectSimulationMemoryReport()` and scheduler `ownedCapacityBytes()`; nonuniform/optional exception lanes remain visible.

For the stable homogeneous rung-zero DMO model represented by the current source:

| Owner | Generic model | M48-06 compact model | Source-derived saving |
|---|---:|---:|---:|
| Scheduler population state | ~33 B/particle | 4 B/particle | ~29 B/particle |
| Canonical particle + metadata + identity species maps | ~109 B/particle | 64 B/particle | ~45 B/particle |
| Combined | ~142 B/particle | ~68 B/particle | ~74 B/particle |

At 512^3 particles, 74 B/particle is approximately 9.25 GiB aggregate source/model-derived ownership reduction: about 3.625 GiB from scheduler storage and 5.625 GiB from canonical metadata/species-index storage. A nonuniform materialized SFC exception adds 8 B/particle to the compact canonical state and correspondingly reduces that saving.

These numbers describe source/model ownership, not measured process RSS. Allocator behavior, MPI/FFTW/HDF5 internals, page residency and runtime overlap remain to be qualified separately.

## Reproducibility contract

The representation change preserves the intended logical contract:

- global timestep criteria and accepted physical timestep are unchanged
- KDK stage ordering is unchanged
- rung-zero `current_tick` advances once per completed substep
- all production particles remain active each global step
- physical phase-space/mass and stable particle IDs are unchanged
- species, owner, flags, time-bin and drift-epoch logical values are unchanged
- TreePM borrowed-source inputs remain canonical
- science snapshot logical fields are unchanged
- restart logical datasets and integrity ordering are unchanged

No empirical bitwise/runtime equivalence is claimed in this implementation pass.

## Deferred work

- Generic migration packet architecture remains; M48-07 is still responsible for bounded DMO-native migration redesign.
- Full restart restoration still materializes the explicit on-disk logical metadata before post-verification compaction.
- Complete restart verification/readback duplication remains M48-08.
- Tree resolved epsilon and later PM/FFT storage campaigns remain M48-09 scope.
- Runtime RSS/scaling and numerical qualification are deferred to the dedicated validation campaign.

## Validation status

Build/tests: NOT RUN — intentionally deferred by campaign instruction.

Only static patch review/hygiene is permitted for this campaign.
