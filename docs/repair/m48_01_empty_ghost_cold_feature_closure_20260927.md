# M48-01 empty-ghost and cold-feature memory closure

Date: 2026-09-27. Mode: focused production repair.

## Architectural closure

The generic particle-ghost refresh now determines local demand by scanning the
authoritative `particle_sidecar.owning_rank` lane without allocating a
population-sized helper. Preparation failures are coordinated first, then all
ranks use one collective maximum reduction to decide whether any rank has ghost
demand. When demand is globally empty, the lifecycle is still invalidated,
committed at the current `GhostLayerEpoch`, and required fresh; profiling records
a successful zero-byte refresh. No local descriptor table, complete gravity
payload replica, owned-ID hash, per-peer request vectors, or wire payload is
materialized on that path.

When demand exists, the current descriptor/ID/epoch validation remains active,
but payload packing now borrows a `ReadOnlyGhostExchangeView` over canonical
particle lanes instead of cloning IDs, positions, masses, and velocities into a
second full SoA. Only requested wire rows and received ghost rows are owned by
the protocol. Received generic particle ghosts are validated against the
exchange plan and immutable particle IDs before their gravity/kinematic lanes
are committed; hydro lanes remain forbidden on this generic path.

Black-hole timestep lookup storage is now cold. The population-sized
`bh_row_by_particle` vector is materialized only when AGN physics is enabled and
black-hole sidecar rows actually exist. Existing enabled-path extent, bounds,
and duplicate-row checks are unchanged.

TreePM high-resolution source/target masks and spatial/file membership
classification are now materialized only when zoom long-range correction is
enabled. Non-zoom solves pass empty high-resolution spans and retain zero mask
capacity in a fresh runtime. The gravity preflight model likewise charges the
high-resolution classification bytes only for zoom-enabled estimates.

The production all-active KDK path no longer constructs an all-particle direct
view before the pre-kick gravity runtime can admit governed scratch. The
pre-kick admits `TransientStepWorkspace::gravity_particle_index_scratch` before
the subsequent drift view, including the cached-force path that intentionally
skips the full force-refresh process preflight. Fresh-force paths preserve the
existing post-preflight admission point. Subset and standalone compatibility
paths remain available.

## Ownership, accounting, and scientific scope

This patch changes storage ownership and lifetime only. It does not change
TreePM force equations, TSC/deconvolution policy, softening, integration order,
timestep criteria, zoom classification semantics when zoom is enabled,
snapshot/restart schemas, or configuration keys. Existing actual-capacity
reporting continues to expose the workspace scratch arena and gravity-runtime
vectors; after this closure the compatibility all-particle vector and non-zoom
mask capacities no longer appear as mandatory production owners. No RSS saving
is claimed here: runtime qualification remains pending.

The later M48 work on PM indexed interpolation/spectral compaction, tree
workspace/topology changes, communication-arena unification, homogeneous-state
policy, scheduler redesign, migration redesign, restart streaming, in-place
FFT, mixed precision, and topology compression is intentionally deferred.

Build/tests: NOT RUN — intentionally deferred by campaign instruction.
