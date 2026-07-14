# Tree gravity solver (Barnes-Hut)

## Scope and ownership

`TreeGravitySolver` is the standalone, nonperiodic Barnes-Hut implementation in
`include/cosmosim/gravity/tree_gravity.hpp` and `src/gravity/tree_gravity.cpp`.
It owns transient tree topology, Morton ordering, node multipoles, softening
envelopes, and traversal scratch. Caller-owned particle arrays remain the
authoritative source state. `TreePmCoordinator` reuses the tree implementation
for the short-range branch but supplies its own periodic geometry and screened
kernel; periodic behavior is therefore documented separately in
`docs/tree_pm_coupling.md`.

The tree uses:

- SoA node storage (`TreeNodeSoa`);
- stable Morton ordering;
- explicit build, multipole, and active-target traversal stages;
- monopole or second-order multipole force evaluation;
- per-source and per-target softening sidecars;
- geometric, COM-distance, or relative-force opening criteria.

Standalone tree construction remains Euclidean and nonperiodic. Do not call it
a periodic solver merely because TreePM applies minimum-image distances around
the same traversal machinery.

## Multipole definition and force convention

Let

\[
\mathbf d = \mathbf x_{\rm com}-\mathbf x_{\rm target},\qquad
I_{ij}=\sum_a m_a s_{a,i}s_{a,j},\qquad
\mathbf s_a=\mathbf x_a-\mathbf x_{\rm com}.
\]

The stored traceless quadrupole and the separately stored trace of the raw
second central moment are

\[
Q_{ij}=3I_{ij}-{\rm Tr}(I)\delta_{ij},\qquad
T={\rm Tr}(I).
\]

For the unsoftened Newtonian kernel and this COM-to-target displacement
convention, the accepted-node acceleration through quadrupole order is

\[
\mathbf a = GM\frac{\mathbf d}{r^3}
-G\frac{Q\mathbf d}{r^5}
+\frac{5G}{2}\frac{(\mathbf d^TQ\mathbf d)\mathbf d}{r^7}.
\]

The signs above are intentional. The former opposite-sign expression was not
consistent with `d = x_com - x_target` and is not the current implementation.

The production tree also supports Plummer-equivalent pair softening. A softened
radial kernel is not harmonic, so `Q` alone is insufficient. With

\[
f(r)=(r^2+\epsilon_{\rm pair}^2)^{-3/2},
\]

the code evaluates the true second-order central-moment expansion

\[
a_i = G\left[M d_i f
+I_{ij}d_j\frac{f'}{r}
+\frac{T}{2}d_i\frac{f'}{r}
+\frac{d_i}{2}(\mathbf d^TI\mathbf d)
\left(\frac{f''}{r^2}-\frac{f'}{r^3}\right)\right].
\]

`TreeNodeSoa::second_moment_trace` preserves `T`, while `I` is reconstructed
from `Q` and `T`. Child-to-parent accumulation applies the parallel-axis terms
to both quantities. This is also the moment basis used by screened TreePM.

## Softening contract

Source epsilon resolution order is:

1. per-particle source sidecar, when present and selected by its override mask;
2. source species table;
3. `TreeSofteningPolicy::epsilon_comoving`.

Target epsilon resolution order is:

1. compact active-target sidecar and override mask;
2. compact active-target species table;
3. source-indexed sidecar/species fallback when the active target is a local
   source row;
4. the scalar policy fallback.

Every leaf pair uses

\[
\epsilon_{ij}=\max(\epsilon_i,\epsilon_j).
\]

An accepted node uses the conservative aggregate
`max(target_epsilon, node_softening_max)`. A heterogeneous-softening node is
forced open while the target is inside its size-plus-softening envelope. This
is deliberately conservative; a single node maximum can still bias an
accepted mixed-softening node and is not a substitute for a calibrated
adaptive-softening model.

## Multipole-acceptance criteria

For node width `l = 2h`, COM distance `r`, and COM displacement from the
geometric node center `delta_com`, the supported criteria are:

1. Geometric Barnes-Hut:

   \[
   l/r < \theta.
   \]

2. COM-distance Barnes-Hut (default):

   \[
   (l+\delta_{\rm com})/r < \theta.
   \]

3. Relative force-error control:

   \[
   G M l^2 \le
   \alpha\max(|\mathbf a_{\rm previous}|,a_{\rm floor})r^4.
   \]

The relative criterion uses `relative_force_tolerance` for `alpha` and
`relative_force_acceleration_floor_code` for `a_floor`. Its safety behavior is:

- an internal node containing the target is always opened, preventing an
  accepted multipole from including the target's own mass;
- a missing or non-finite previous acceleration selects the deterministic
  COM-distance criterion, using `opening_theta`;
- a finite zero history uses `a_floor` and therefore remains well-defined;
- leaves remain exact P2P evaluations;
- the same softening-envelope guard applies after any MAC accepts a node.

`TreeGravityTargetView::previous_acceleration_magnitude_code` and the equivalent
legacy overload argument are optional compact active-target lanes. They do not
own force history. In the production workflow, a valid committed force cache
with matching dense-row generation supplies the previous magnitude; otherwise
the deterministic fallback is used.

## Configuration, provenance, and migration notes

The workflow maps the normalized config keys directly to `TreeGravityOptions`:

- `numerics.treepm_tree_opening_criterion` (`geometric`, `com_distance`, or
  `relative_force_error`);
- `numerics.treepm_tree_opening_theta`;
- `numerics.treepm_tree_relative_force_tolerance`;
- `numerics.treepm_tree_relative_force_acceleration_floor`.

The same values are recorded in `provenance_v6`. The default remains
COM-distance opening. The public default
`TreeGravityOptions::opening_theta` in
`include/cosmosim/gravity/tree_gravity.hpp` is now `0.7`, aligned with the typed
production config and the Ewald-certified runtime-default profile. Callers that
relied on aggregate/default construction at the former public-header value
`0.6` must set `opening_theta=0.6` explicitly to preserve that behavior. Adding
the relative criterion itself did not silently select it for existing runs.

`TreeGravityOptions::gravitational_constant_code` remains an explicit caller
input. Production `ReferenceWorkflow` no longer assumes dimensionless `G=1`:
it constructs `core::UnitSystem` from the frozen `units.*` config and calls the
additive public helper
`core::newtonGravitationalConstantCode(const UnitSystem&)`, declared in
`include/cosmosim/core/units.hpp`. This converts the
physical SI Newton constant into `L_code^3/(M_code T_code^2)`, with
`T_code=L_code/V_code`. Standalone dimensionless tests may continue to pass an
explicit `G`; external workflow integrations should migrate away from hard-coded
unity when their configured units are physical.

Public-interface changes that callers must account for are:

- `TreeOpeningCriterion::kRelativeForceError`;
- `TreeGravityOptions::{relative_force_tolerance,
  relative_force_acceleration_floor_code}`;
- the optional previous-acceleration span on `TreeGravityTargetView` and the
  active-set overload;
- `TreeNodeSoa::second_moment_trace`, which is required whenever node storage is
  inspected or serialized by non-production tooling;
- `TreeGravityProfile::{build_count, multipole_refresh_count, opened_nodes}`;
  callers mirroring or aggregate-initializing the public profile must account
  for the appended observer counters.

The tree itself is transient and is rebuilt after drift, migration, or
decomposition changes. No tree topology or multipole cache is restart truth.

## Validation and profiling

- `tests/unit/test_tree_gravity.cpp` covers all three MACs, fallback/floor
  behavior, target-containing guards, multipole signs, raw second moments, and
  force/work diagnostics.
- `tests/integration/test_tree_gravity_vs_direct.cpp` compares quasi-uniform and
  clustered fixtures against direct summation and checks the expected
  error/work trend across opening thresholds.
- `tests/integration/test_tree_pm_coupling_periodic.cpp` covers the periodic
  TreePM use of monopoles, quadrupoles, all three MAC paths, seam geometry, and
  translation invariance.
- `bench/bench_tree_gravity.cpp` exposes build/traversal timing and interaction
  counters.

`TreeGravityProfile` distinguishes tree builds and multipole refreshes from
traversal work through `build_count`, `multipole_refresh_count`,
`visited_nodes`, `accepted_nodes`, `opened_nodes`, and
`particle_particle_interactions`. `opened_nodes` counts internal-node descent,
so MAC changes can be audited without inferring work from timing alone. These
are observer counters; they do not own or invalidate tree state.

The isolated integration thresholds remain `max_relative_error < 0.02` and
`mean_relative_error < 0.01` for the tight opening in the covered fixtures.
These are regression limits for those tests, not a universal astrophysical
force-error guarantee.
