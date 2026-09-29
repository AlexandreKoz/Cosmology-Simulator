# M48-05 — Homogeneous DMO Borrowed Gravity Sources and Target/Map Collapse

Date: 2026-09-29. Mode: focused production memory-architecture repair.

## Scope and representation

M48-05 adds one explicit runtime representation policy shared by the gravity
source snapshot and memory estimator:

```text
materialized_generic
borrowed_homogeneous_dmo
```

The policy is selected from authoritative runtime state. It is not a user-facing
configuration mode and it does not create a second gravity-source authority.
`SimulationState` remains the canonical owner; `GravitySourceSnapshot` remains a
non-owning immutable view.

## Borrowed homogeneous-DMO eligibility

The borrowed path fails closed. The state-only predicate requires all of the
following before the runtime can credit the memory reduction:

- zoom long-range correction is disabled;
- there are no authoritative gas cells;
- the local particle count fits the 32-bit TreePM local-index contract;
- particle and particle-sidecar extents are consistent;
- every local particle is owned by the current MPI rank;
- every local particle has the canonical dark-matter species tag;
- no particle carries an explicit gravity-softening override;
- canonical coordinates and masses satisfy the existing finite/non-negative
  source invariant;
- periodic coordinates are already canonical inside the configured box;
- the exact effective dark-matter softening is finite, non-negative, and inside
  the currently certified TreePM epsilon/r_s envelope.

The force-evaluation predicate additionally requires the force target span to be
the dense identity population and rejects borrowing whenever the current
scheduler set proves that any source needs independent inactive-coordinate
prediction. Eligibility is an O(N) read-only scan with O(1) owned scratch; it
does not build a source-row vector merely to discover the identity case.

## Borrowed source ownership

For `borrowed_homogeneous_dmo`, `GravitySourceSnapshot` aliases the canonical
particle SoA lanes directly:

```text
SimulationState::particles.position_x_comoving
SimulationState::particles.position_y_comoving
SimulationState::particles.position_z_comoving
SimulationState::particles.mass_code
```

The runtime therefore releases retained capacity for the redundant materialized
source owners:

```text
gravity_runtime.source_x
gravity_runtime.source_y
gravity_runtime.source_z
gravity_runtime.source_mass
gravity_runtime.source_species
gravity_runtime.source_softening
gravity_runtime.source_softening_mask
gravity_runtime.source_particle_row
gravity_runtime.source_cell_row
```

Homogeneous species is represented implicitly as dark matter. The canonical
`ParticleSidecar::species_tag` and optional canonical softening sidecars remain
untouched; canonical-state compaction belongs to a later campaign.

The snapshot is stamped with gravity-source generation, force-evaluation epoch,
and particle-index generation. `GravityRuntime` also keeps the particle-index
and source-generation stamps while the borrowed representation is active.
Migration/decomposition commit invalidates those stamps and drops the borrowed
active span. The step-workspace span is also cleared from `GravityRuntime` after
the solve so a raw borrowed pointer cannot survive a later scratch reset.

## Prediction-aware fallback

The generic prediction path remains authoritative when an inactive source is not
already represented at the requested force epoch. It keeps the existing
cosmological predictor:

```text
x_pred = x + v_x * drift_factor
y_pred = y + v_y * drift_factor
z_pred = z + v_z * drift_factor
```

including drift-time/scale validation, future-epoch rejection, cosmological
drift-factor evaluation, finite-factor checks, periodic wrapping, and source
finite/negative-mass validation.

The prediction membership mask now follows the scheduler-current particle set,
not the possibly widened force-cache rebuild target set. Consequently, a dense
force-cache rebuild may reuse the identity target-index scratch while still
falling back to materialized predicted source coordinates when the scheduler
has inactive rows. Borrowed canonical coordinates are never mixed with stale
inactive sources.

Top-domain geometry refit consumes the selected source view. If predicted
inactive coordinates are present, it refits even when the persistent gravity
source generation has not changed, so routing geometry follows the coordinates
actually passed to TreePM.

## Target and mapping collapse

For the borrowed path:

```text
source index == canonical particle row == target slot == active force slot
```

The runtime releases retained capacity for the redundant identity owners:

```text
gravity_runtime.owned_local_index_by_particle
gravity_runtime.owned_local_index_by_cell
gravity_runtime.owned_leaf_cell_mask
gravity_runtime.target_particle_row
gravity_runtime.target_cell_row
gravity_runtime.local_active_source_index
gravity_runtime.local_active_particle_row
gravity_runtime.active_slot_by_particle
gravity_runtime.active_slot_by_cell
gravity_runtime.force_refresh_particle_rows
```

The single explicit target/source index lane is
`TransientStepWorkspace::gravity_particle_index_scratch`. It is populated with
`0..N-1`, used as the TreePM active/source-index span, and can then be reused by
the existing all-particle direct drift view. `GravityRuntime` does not allocate a
second population-sized `uint32` vector beside it.

Particle-to-active-slot lookup is centralized in the representation policy:
borrowed DMO maps row directly to slot; generic operation uses the retained
`m_active_slot_by_particle` map. Zero-cell DMO owns no cell mapping capacity.

## Species and softening semantics

The borrowed snapshot exposes an empty species span to express a proven
homogeneous population. The exact dark-matter scalar softening is selected from
the configured species policy when it is enabled, otherwise from the configured
fixed gravity softening, and is passed as the existing Tree softening fallback.
The tree's resolved `m_source_softening_epsilon_comoving` remains unchanged.

Any explicit per-particle softening override rejects the borrowed representation
and uses the generic source staging. No softening approximation is introduced.

## TreePM and PM data flow

Top-domain refit, health checks, previous-acceleration lookup, force-cache
scatter, and the TreePM call consume representation-aware source/target views.
TreePM receives the selected source XYZ/mass spans directly. The M48-02 direct
indexed PM interpolation path remains unchanged:

```text
governed uint32 target/source index
  -> direct indexed source-coordinate access
  -> PM interpolation
```

No compact target XYZ array is introduced.

## Fallback admission and memory governor

The state-only preflight may credit the borrowed representation only when its
base invariants are provable. If the force-evaluation shape later requires the
generic representation (for example, inactive prediction or a subset target),
`GravityRuntime` computes the generic-versus-borrowed incremental difference,
collectively admits that additional `PhaseResident` amount through the existing
`MemoryGovernor`, and commits it before materialized O(N) source/map allocation.
A low-memory rank therefore cannot silently allocate the generic live set after
a borrowed-sized reservation.

The one all-particle `uint32` lane remains a first-class governed scratch-arena
owner. The gravity estimator charges that lane once in its borrowed model; the
gravity phase lease excludes those same bytes because the
`TransientStepWorkspace` arena admits and commits the physical allocation
separately.

## Memory model

The exact physical reduction is representation- and feature-dependent, so this
campaign does not convert the old audit coefficient into an RSS claim. The
source-derived ownership changes for ordinary zero-cell homogeneous DMO are:

```text
removed from GravityRuntime fast path
  source XYZ + mass                    32 B / particle
  copied species                        4 B / particle
  source particle/cell row lanes        8 B / particle
  target source/global/particle/cell
    identity lanes                     16 B / particle
  particle -> source map                4 B / particle
  particle -> active-slot map           4 B / particle
  duplicate force-refresh identity      4 B / particle when live
  copied source softening               8 B / particle when materialized
  copied softening override mask        1 B / particle when materialized

retained explicit identity owner
  governed uint32 scratch               4 B / particle
```

Thus the historical approximately 72 B/particle (~9 GiB aggregate at 512^3,
~1.125 GiB/rank at eight balanced ranks) opportunity remains the correct order
of magnitude for the duplicated identity architecture. The current source also
has an optional materialized per-source softening lane, so the direct owner
reduction can be larger when that staging was populated. Conversely, the exact
instantaneous delta depends on whether the force-refresh vector is live and on
retained allocator capacity.

The memory estimator now models `materialized_generic` and
`borrowed_homogeneous_dmo` explicitly. Borrowed mode charges no copied source
staging or runtime identity maps and charges exactly one `uint32` target/source
index lane. When the relative-force MAC is configured, the estimator also
retains a conservative charge for the lazy previous-acceleration magnitude
lane; that scientific criterion is not folded into the identity-map saving.
Persistent force cache, active acceleration output, tree construction,
periodic tree staging, PM resources, communication arena, halo cache, and other
real owners remain charged.

These are **source/model-derived ownership reductions**, not measured RSS.
Allocator RSS return and whole-process runtime qualification are deferred.

## Diagnostics and accounting

Runtime telemetry now reports:

```text
gravity_source_representation
borrowed_source_count
materialized_source_count
predicted_inactive_source_count
```

The existing predicted-inactive diagnostic is preserved for compatibility.
`GravityRuntime::memoryReport()` continues to report every physical vector
capacity; in borrowed mode the obsolete generic owners have zero retained
capacity rather than being hidden from accounting.

## Reproducibility and numerical scope

M48-05 changes representation, ownership, and mapping only. It is intended to
preserve TreePM split mathematics, PM assignment/interpolation, source order,
force accumulation order, KDK kick factors, Hubble drag, force-cache validity,
source-generation semantics, prediction arithmetic, softening semantics,
zoom/full-physics fallback behavior, restart schema, and migration wire format.
No empirical numerical-equivalence claim is made in this implementation pass.

## Deferred work

This campaign intentionally leaves the following for later work:

- dense homogeneous canonical metadata in `SimulationState`;
- generic scheduler/time-bin storage and uniform-rung compaction;
- the tree's population-sized resolved epsilon lane;
- materialized source prediction for subset/hierarchical force epochs;
- generic migration records and bounded migration redesign;
- implicit all-active TreePM index ranges;
- runtime RSS, rank-count equivalence, force-error, restart, and long-run
  qualification.

Build/tests: NOT RUN — intentionally deferred by campaign instruction.
