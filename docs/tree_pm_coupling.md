# TreePM coupling

## Residual traversal locality — October 2026

Leaves address contiguous ranges of `TreeMortonOrdering::sorted_particle_index`.
Each residual evaluation now gathers the caller-owned wrapped source positions
and masses into tree order (`treepm.tree_ordered_source_*` lanes) and visits
local active targets in tree order (`treepm.target_visit_order`). Leaf pairs
therefore read sequential memory instead of gathering through the sorted index,
and consecutive targets in a block share node cache lines.

- The lanes are derived scratch regathered on every evaluation (O(N), a few
  milliseconds at 262k sources); they are never authoritative, persisted or
  reused across source generations. SimulationState rows are not permuted, so
  ownership, restart, migration and force-cache contracts are unaffected.
- Per-target walk order and per-target summation order are unchanged; each
  active slot keeps exactly one writer. Forces are bitwise identical to the
  previous traversal (verified on 4096 targets of a clustered 64^3 box).
  Only block-ordered floating diagnostics (sum of squares) change order.
- Memory: 4 doubles + 2 uint32 per local source (the second is the optional
  group-walk leaf map) and 1 uint32 per target,
  reported in the TreePM memory report and in
  `gravity.estimate.treepm_tree_ordered_source_lanes`.
- Measured (clustered 64^3 box, 64^3 PM, adaptive alpha=0.005, 4 threads):
  residual traversal 11.9 s -> 5.9 s with identical pair counts.

## Adaptive estimator and cutoff-node policy — October 2026

Two adaptive-policy options (typed config keys in `docs/configuration.md`):

- `adaptive_error_estimator = relative_monopole` accepts a node when
  `G M l^2 <= alpha |a_previous| r^4`, the leading monopole truncation scale,
  without the derivative proxy. The derivative proxy charges second-order terms
  that accepted quadrupoles have already removed (about 9 times the monopole
  term for `f ~ r^-3`) and costs one `erf`/`exp` per tested node. With valid
  history the COM-center offset is not computed. Neither estimator is a bound;
  alpha is calibrated by the measured errors below.
- `cutoff_node_policy = com_distance` lets a node that straddles `r_cut` reach
  the MAC; an accepted straddling node is evaluated at its COM distance and is
  skipped when that distance exceeds `r_cut`. Mass inside an accepted node
  beyond `r_cut` then feeds the residual kernel untruncated, where the residual
  is the Gaussian tail (`F_SR(r_cut=5 r_s) ~ 6e-3` of Newtonian). `contain`
  keeps the previous opening-to-leaves behavior.
- `adaptive_maximum_opening_angle` may now be up to 1.0.

The strict reference policy ignores both options. Defaults reproduce the
previous behavior exactly.

Measured on the synthetic clustered 64^3 box
(`tools/generate_clustered_dmo_box.py`, 10 length units, epsilon 0.005,
random row order) with `bench_tree_pm_sweep`, Release `-O3 -march=native`,
4 OpenMP threads, same-snapshot bootstrap history, Gaussian pair lookup,
leaf 16, block 64. Errors are |a - a_ref|/|a_ref| of the TOTAL Tree+PM force
over 4096 targets, against the strict reference policy on the same mesh
(so split/mesh error is excluded). Times are one force evaluation.

| PM | Estimator | Cutoff nodes | alpha | max angle | Tree s | PM s | pairs/target | rms err | p99 err | max err |
|---|---|---|---|---|---|---|---|---|---|---|
| 64 | residual_derivative | contain | 0.005 | 0.5 | 6.05 | 0.10 | 2740 | 3.1e-5 | 1.0e-4 | 1.3e-4 |
| 64 | relative_monopole | contain | 0.005 | 0.5 | 4.72 | 0.10 | 2026 | 1.6e-4 | 4.8e-4 | 7.0e-4 |
| 64 | relative_monopole | com_distance | 0.005 | 1.0 | 1.89 | 0.10 | 769 | 5.3e-4 | 1.3e-3 | 4.5e-3 |
| 64 | relative_monopole | com_distance | 0.0025 | 1.0 | 2.56 | 0.14 | 1039 | 2.7e-4 | 7.2e-4 | 2.5e-3 |
| 64 | relative_monopole | com_distance | 0.01 | 1.0 | 2.07* | 0.29* | 566 | 1.0e-3 | 2.7e-3 | 7.1e-3 |
| 128 | residual_derivative | contain | 0.005 | 0.5 | 3.92 | 0.18 | 1897 | 2.8e-5 | 9.8e-5 | 1.5e-4 |
| 128 | relative_monopole | contain | 0.005 | 0.5 | 3.12 | 0.19 | 1481 | 1.4e-4 | 4.7e-4 | 7.1e-4 |
| 128 | relative_monopole | com_distance | 0.005 | 1.0 | 1.32 | 0.20 | 535 | 4.0e-4 | 1.2e-3 | 4.6e-3 |
| 128 | relative_monopole | com_distance | 0.0025 | 1.0 | 1.70 | 0.19 | 737 | 2.1e-4 | 6.0e-4 | 2.5e-3 |
| 128 | relative_monopole | com_distance | 0.01 | 1.0 | 1.11 | 0.22 | 389 | 7.3e-4 | 2.3e-3 | 6.3e-3 |

`*` measured while a compiler job shared the cores; pair counts are exact.
Before the locality change the first row took 11.9 s. The production deck
`configs/treepm_evolved_snapshot_fast.param.txt` uses the 128^3 /
relative_monopole / com_distance / alpha 0.0025 / angle 1.0 row. These are
force-accuracy results on one synthetic workload; trajectory, growth and
power-spectrum qualification of the selected settings remain separate gates.

## Group walk — October 2026 (opt-in)

`TreePmOptions::group_walk_enabled` / `numerics.treepm_group_walk_enabled`
(default false) walks the tree once per group instead of once per target. A
group is the run of active targets, in tree-order visit sequence and within one
residual block, that lie in the same tree leaf (at most `max_leaf_size`).

Group-conservative acceptance (`src/gravity/internal/tree_pm_group_walk.hpp`):

- prune a node when the minimum distance between its AABB and the group AABB
  exceeds `r_cut`;
- open a node whose AABB touches the group AABB (no target may lie inside an
  accepted node), and under `contain` a node not entirely within `r_cut` of
  every group point;
- accept a node only if `G M l^2 <= alpha a_min r_min^4` and
  `l < theta_max r_min`, with `r_min` the distance from the node COM to the
  group AABB (a lower bound on every target's COM distance) and `a_min` the
  smallest previous-acceleration magnitude in the group.

Every target therefore meets the per-target `relative_monopole` criterion with
margin, and the measured error decreases. Accepted multipoles are evaluated
per target exactly as in the per-target walk (including the `com_distance`
zero beyond `r_cut`). Accepted leaves become tree-order source ranges, merged
when contiguous, and are evaluated by a branch-free pair loop over the
tree-ordered lanes with the Gaussian Hermite table (AVX2 gathers); the
compiler vectorizes it (`#pragma omp simd`). A Chebyshev/Estrin polynomial
alternative was measured slower on the reference CPU and was not kept.

Eligibility: adaptive policy, `relative_monopole`, Gaussian pair lookup,
uniform softening, local (non-distributed) residual path, leaf SIMD and
full-residual lookup off, no spatial work history, and
`(r_cut / 2 r_s)^2 < 63` (table range). Typed config rejects enabling it
without those settings. At run time a group falls back to the per-target walk
(before any accumulation) when a target has non-uniform softening, missing or
sub-floor force history (first evaluation), or periodic coordinates outside
the certified bounds; isolated/open boundaries therefore always fall back.

Differences from the per-target walk: summation order (not bitwise
reproducible against it), coincident distinct pairs are not counted in
`direct_pair_evaluations`, and pairs/target rise while node visits fall.

Measured (same clustered 64^3 box and settings as the table above, 128^3 PM,
alpha 0.0025, angle 1.0, `com_distance`, 4 threads):

| Walk | Tree s | pairs/target | nodes/target | rms err | p99 err | max err |
|---|---|---|---|---|---|---|
| per target | 1.44 | 737 | 351 | 2.1e-4 | 6.0e-4 | 2.5e-3 |
| group (leaf 16) | 1.06–1.11 | 947 | 102 | 1.5e-4 | 4.4e-4 | 2.5e-3 |

Leaf 32 groups cost more (1.5 s, 1605 pairs/target); leaf 8 matches leaf 16.
The recovery Ewald qualification gates the group-walk fast profile against the
unchanged Ewald targets on all five fixtures.

## Attempt 1 recovery interfaces — 20261012

The [implementation handoff](repair/treepm_deep_hotpath_attempt1_closure_20261012.md)
documents unvalidated source changes and exact manual gates. All new arithmetic
switches default false. The original analytic multipole implementation remains
available through `radial_coefficients_enabled=false`; the new analytic shared
derivatives retain the raw second-moment trace and do not use a pair table.
Nonfinite shared coefficients fall back to that original implementation.
Adaptive MAC derivatives remain on the original analytic path in both modes.

Periodic metric selection uses actual original-source and unwrapped-node
center/COM bounds per axis and target. Certified deltas in [-L,L] use one
correction; arbitrary multi-box inputs use nearbyint. Half-box rounding bands
and non-nearest rounding retain the general helper. Canonical source bounds
also bypass redundant floor normalization during unwrap/refit preparation.
Imported targets are independently certified against the receiving local tree;
remote top-domain/packet geometry keeps its existing general arithmetic.

COM-center offsets are lazy and measured in the common unwrapped tree frame.
They are omitted for leaves, independently mandatory fast descent and a selected
MAC that does not use them. No per-node cache is retained. Full-forensic mode
keeps overlapping rejection accounting. In fast mode geometry/MAC skip counters
refer to unevaluated tests, never implicit passing guards.

`TreePmResidualCoefficientTable` is a new immutable, homogeneous-pair-only
optional API. Its domain and error envelope are described in the handoff;
fallback remains the analytic residual. `TreePmLeafScratch` and the internal
leaf kernel stage bounded blocks in existing tree order, dispatch AVX2 at
runtime where supported and sum in the original source order. These additive
public options/records require recompilation; existing callers retain defaults.
`TreeGravityProfile` adds explicit refit attempted/accepted/rejected and geometry
preparation time. Historical `reuse_accepted` includes exact reuse or refit;
`rebuild_reason` and the new refit fields distinguish them. No wire, restart or
particle storage format changes accompany these APIs.

## Scope

`TreePmCoordinator` combines the periodic or isolated PM long-range field with
a Barnes-Hut short-range residual. The production periodic split is

\[
r_s={\tt asmth\_cells}\,\Delta_{\rm mesh},\qquad
r_{\rm cut}={\tt rcut\_cells}\,\Delta_{\rm mesh},
\]

where the current runtime derives

\[
\Delta_{\rm mesh}=\sqrt[3]{(L_x/N_x)(L_y/N_y)(L_z/N_z)}.
\]

The Fourier-space long-range filter is

\[
F_{\rm LR}(k)=\exp(-k^2r_s^2),
\]

and the complementary real-space force factor is

\[
F_{\rm SR}(r)=\operatorname{erfc}(q)
+\frac{2q}{\sqrt{\pi}}e^{-q^2},\qquad q=\frac{r}{2r_s}.
\]

`r_cut` is an explicit truncation of the residual, not part of the analytic
Gaussian identity. Before that truncation, the PM and short-range kernels are
the complementary TreePM split. Both branches return the same scale-free
comoving Newtonian kernel and use the same physical Newton constant converted
to configured code units, `G_code`. Neither branch multiplies by a scale-factor
power.

Periodic operation requires
`r_cut < 0.5 min(Lx,Ly,Lz)`. The residual tree evaluates one minimum image per
source; a cutoff at or beyond half the shortest axis would make that
single-image contract ambiguous. Typed config rejects invalid derived geometry,
and `TreePmCoordinator` repeats the check for direct API callers. The remedy is
to increase PM resolution or reduce `rcut_cells`, followed by validation of the
changed force profile.

`pm_options.gravitational_constant_code` and
`tree_options.gravitational_constant_code` must therefore be identical at the
TreePM API boundary. A mismatch is rejected; the caller-owned standalone-tree
option is not silently repurposed. Production `ReferenceWorkflow` derives this
value once with `core::newtonGravitationalConstantCode(UnitSystem)` from the
frozen length, mass, and velocity units. For cosmological physical peculiar
velocity `u=a dx/dt`, TreePM does not apply Hubble drag or scale factors. The
collisionless KDK and gas `ComovingGravityExpansionSource` consume the common
kernel and apply `du/dt + H u = A/a^2` to their respective state.
`scale_factor` remains force-build validity/source-time metadata and does not
rescale `A`.

Window deconvolution is a periodic-PM transfer correction. Once mode policy
resolves to isolated/open gravity, typed config validation requires
`numerics.treepm_enable_window_deconvolution=false`; the coordinator does not
silently reinterpret the periodic default. The internally constructed focused
zoom correction likewise uses the isolated/open operator with deconvolution
disabled.

## Periodic tree geometry

Periodic traversal cannot be made correct by applying a minimum-image delta
only at the final force evaluation. Tree topology, COMs, moments, node bounds,
the MAC, cutoff pruning, and remote summaries must use compatible geometry.

For periodic TreePM, the coordinator therefore constructs a transient derived
source frame independently on each axis:

1. wrap every finite source coordinate into `[0,L_axis)`;
2. stable-radix order the wrapped lane by its IEEE-754 key;
3. locate the largest circular gap, with a deterministic anchor tie-break;
4. unwrap values below the post-gap anchor by adding `L_axis`.

The resulting interval is the shortest contiguous source interval for that
axis. The tree is built in this per-axis unwrapped frame, so topology, root
extent, COMs, quadrupoles, and raw second moments all remain compact when a
physical cluster crosses `0/L`. The caller-owned wrapped particle lanes remain
authoritative and are still used for PM deposition.

Tree construction owns one reusable population-scale workspace: two
`uint64_t` key lanes and one 32-bit local-index lane. Periodic preprocessing
stores each axis's original wrapped values directly in its final unwrapped
coordinate lane, encodes the non-negative IEEE-754 representation into the two
key lanes for the stable byte-radix order, selects the same largest-gap anchor,
and then transforms that final coordinate lane in place. X, Y, and Z reuse the
same key storage sequentially; TreePM does not retain separate wrapped/ordered
axis vectors.

During short-range traversal, node-center, COM, and particle deltas are reduced
with the axis-specific minimum image. Periodic AABB distance functions also
accept unwrapped intervals, so cutoff pruning cannot discard a nearby image of
a seam-crossing node. Rectangular `Lx`, `Ly`, and `Lz` are supported; no cubic
alias is used when axis lengths are supplied.

Hierarchy packets mark this representation with `geometry_frame=1` and carry
the unwrapped node bounds. `geometry_frame=0` means ordinary Euclidean bounds.
Mixed frames in one exchange are rejected. Deterministic x/y/z seam, edge,
corner, rectangular-box, leaf/internal-node, monopole/quadrupole, geometric/COM
MAC, and independent integer-image translation tests are in
`tests/integration/test_tree_pm_coupling_periodic.cpp`.

## Short-range force and cutoff contract

For leaf pairs, the tree first applies the same pair-softening law as standalone
tree gravity,

\[
\epsilon_{ij}=\max(\epsilon_i,\epsilon_j),
\]

then multiplies the softened force by `F_SR(r)`. The periodic source-target
distance is minimum-image on each axis.

Cutoff behavior is enforced at three levels:

1. a node is pruned when its minimum periodic AABB distance exceeds `r_cut`;
2. an accepted internal node must lie completely within `r_cut`, otherwise it
   is opened;
3. a leaf pair with `r > r_cut` is skipped.

Accepted quadrupoles use a second-order expansion of the complete radial scalar

\[
f(r)=F_{\rm SR}(r)(r^2+\epsilon_{\rm pair}^2)^{-3/2}.
\]

The implementation contracts `f'` and `f''` with the raw central second moment
`I`, reconstructed from the traceless quadrupole and
`TreeNodeSoa::second_moment_trace`. Merely multiplying the Newtonian quadrupole
by `F_SR` omits derivatives of the screen and is not the implemented model.
The standalone Newtonian quadrupole sign correction described in
`docs/tree_gravity_solver.md` applies to the unsoftened limit.

To keep rank-local forest topology from degrading strict rank-count
equivalence, TreePM monopole mode opens internal nodes to exact leaves. A
screened quadrupole may be accepted only after the configured MAC, cutoff, and
softening guards pass and `l/r < 0.08`. This extra decomposition-stability
envelope is separate from the user-selected MAC.

## Relative force-error MAC in TreePM

TreePM supports geometric, COM-distance, and relative-force opening. The
relative criterion is

\[
G M l^2 \le \alpha\max(|\mathbf a_{\rm previous}|,a_{\rm floor})r^4.
\]

`TreePmForceAccumulatorView::previous_acceleration_magnitude_code` is an
optional compact target lane. The production workflow derives it only from a
valid committed gravity-force cache whose particle-row generation matches the
current state. Missing or non-finite history selects the deterministic
COM-distance fallback; a finite zero uses the configured floor. The history
value and its presence flag are transported with remote target requests, while
the request's force epoch prevents mixing it with another distributed solve.

## Target and source ownership

Source arrays are rank-owned source truth. An active target need not be a local
source:

- when `target_pos_x/y/z_comoving` are absent, every active index addresses a
  local source row and supplies its self-interaction identity;
- when all three target-position lanes are present, they are authoritative;
  `UINT32_MAX` denotes a target with no local source/self identity.

This permits a zero-source rank to own targets without fabricating a mass or a
dummy source. A source-only rank may have an empty active set. An all-empty
collective round is also legal. Empty ranks contribute zero density and force.
Partial target-position triplets, mismatched compact lanes, and out-of-range
source-indexed targets are rejected.

## Distributed local-essential-tree routing and sparse request/response protocol

The distributed short-range path is locality-driven. It keeps a compact global
top-level domain map while retaining detailed tree topology only on the source
rank:

1. each rank builds a tree from rank-owned authoritative sources;
2. each rank publishes exactly one compact top-level source-domain leaf with
   owner rank, bounds, source count, source generation, tree-build generation,
   decomposition epoch, force epoch, and exchange sequence;
3. target owners traverse those top-level bounds against `r_cut` and construct
   the actual peer set required by local targets;
4. an MPI distributed-graph communicator is built from that sparse peer set,
   including reverse edges so source-only and target-only ranks remain legal;
5. selected target work is exported only to those peers;
6. local tree traversal proceeds while the nonblocking neighborhood request
   payload is in flight;
7. each destination evaluates received targets against its detailed local tree
   and returns validated partial accelerations through the same sparse graph.

This is a target-export LET model: remote detailed tree nodes are not globally
replicated. The globally known structure is intentionally limited to one compact
domain ownership leaf per rank. Rank count therefore controls only the compact
top-level routing map; short-range payload communication is determined by
actual domain/target overlap.

An empty local tree emits one explicit zero-source top-level sentinel. It carries
the current geometry frame and semantic identities but never intersects a cutoff
query. This preserves collective participation for zero-source ranks without
creating fake gravity work.

Top-level domain records use the existing fixed, versioned little-endian tree
wire schema. Short-range target requests and acceleration responses remain
versioned little-endian records rather than raw C++ object representations.
Protocol identity covers source/destination rank, exchange sequence, strong
decomposition/source/tree/force identities at the serialization boundary, batch
token, request/target identity, target geometry/softening, previous-force-scale
presence, returned acceleration components and (wire v2) the nearest
direct-pair separation squared.

The compact top-level domain cache is reusable only when every rank agrees that
its decomposition epoch, gravity source generation, and local tree-build
generation still match. A collective cache-validity vote is performed before
any rank skips the top-level exchange, so cache reuse cannot make ranks diverge
into different collective paths. Under the current rung-zero implementation the
tree build generation normally advances every solve, so reuse is conservative
rather than optimistic.

Sparse request counts use `MPI_Neighbor_alltoall`; request payloads use
`MPI_Ineighbor_alltoallv`, allowing useful local traversal to overlap transport;
responses use `MPI_Neighbor_alltoallv`. The short-range path no longer uses
communicator-wide `MPI_Alltoall/Alltoallv` payload phases. World collectives are
still used where they represent true global consensus/failure choreography, and
the compact one-leaf-per-rank domain map is still globally exchanged.

Before exchange, ranks agree on protocol identities and batch policy. Decoders
reject wrong versions, stale/mixed epochs, wrong peers, malformed payloads,
non-finite fields, duplicate/unexpected responses, and incomplete response
coverage. Count multiplication, cumulative byte displacements, and `MPI int`
limits remain checked. Rank-local failures are coordinated before peers enter
the next matching communication phase.

The reusable exchange workspace owns world-sized compact count/displacement
metadata plus payload buffers that grow only with actual exchanged work. Memory
reporting exposes current/capacity/high-water bytes. LET diagnostics expose:

- candidate and communicating peers;
- exported and imported targets;
- sent/received wire bytes;
- LET workspace high-water;
- discovery, communication, and remote traversal time;
- local work overlapped with request transport;
- communication wait time and overlap efficiency;
- local/remote interaction work and imbalance metrics.

The current top-level routing scan is intentionally simple and compact. A future
hierarchical rank-domain tree can reduce peer-discovery CPU work at very large
rank counts without changing the target-export protocol or solver ownership.

## PM cadence and cache validity

Production configuration retains `numerics.treepm_update_cadence_steps=1`.
Default `numerics.hierarchical_max_rung=0` refreshes PM on every authorized
force-refresh surface. Optional P6 DMO block KDK instead refreshes at explicit
coarse endpoints and evaluates Tree-only forces between them, with no PM
interpolation there. It is source-implemented and unqualified. Cadence greater
than one still lacks a specified predictor and fails config validation. The integrator owns the PM synchronization event, kick opportunity,
field version, last refresh opportunity, build step, and build scale factor,
and the workflow requires rank consensus before collective PM work.

The coordinator also validates its transient long-range field before honoring
an explicit lower-level reuse request. Its compatibility signature contains the distinct gravity source generation,
PM-field version, force-evaluation epoch/scale factor, code gravitational
constant, split scale, x/y/z box lengths, assignment scheme, boundary condition,
PM decomposition mode, and window-deconvolution flag. A missing or incompatible
signature makes reuse fail coherently; it is not silently changed into a solve.
Explicit refresh/reuse votes are reduced first, and a mixed vote throws before
any rank enters PM density or FFT collectives. This signature is an invalidation
guard, not a predictor for cadence greater than one.

Tree topology is rebuilt for each default force call; P6 reuse/refit is opt-in. Its explicit
`TreeBuildGeneration`, together with `GravitySourceGeneration` and
`DecompositionEpoch`, prevents stale top-level-domain/LET reuse after source
mutation, migration, or rebuild. The workflow decomposition epoch advances only after
an actual globally committed particle-ownership transition and is restored from
restart; rank-local dense-row generations are not substituted for it. Because
PM field ownership remains the fixed FFT slab map, particle decomposition epoch
alone is intentionally not a PM-field invalidator. Dense-row acceleration
history is invalidated immediately on ownership change. The lower-level
`solveActiveSetWithPmCadence` refresh flag remains test/future-integration
surface only; production rung zero does not exercise local-bin PM reuse.

The tree workspace borrow is also the explicit build-validity boundary. A warm
rebuild invalidates the prior nodes/ordering before overwriting retained Morton
storage, builds the new key/permutation in place, uses the workspace index lane
as radix scratch and then recursive partition scratch, and publishes the next
`TreeBuildGeneration` only after topology and multipoles complete. Failed
construction therefore leaves the prior generation number unpublished and the
partially rebuilt tree inaccessible through the guarded node/ordering views.

TreePM acceleration lanes store scale-free `A`. The post-step adaptive
timestep criterion combines those lanes with comoving softening through
`computeComovingGravityTimeStep`. Its public input carries `eps_com`, unscaled
`|A|`, and the committed `IntegratorState.current_scale_factor`; the helper
validates finite positive `a` and converts to comoving-coordinate acceleration
`|A|/a^3`. The resulting production criterion is
`dt_grav=eta sqrt(a^3 epsilon_com/|A|)`. The `A/a^2` peculiar-velocity response
used by kicks/hydro is not the coordinate acceleration for this length-based
criterion, while the generic `computeGravityTimeStep` remains
coordinate-neutral.

## Accuracy posture

The independent reference in `tests/support/periodic_ewald_reference.*` is a
double-precision, unsoftened, small-N rectangular Ewald sum. It includes real
images, reciprocal modes, central-self removal, periodic self images, the
omitted `k=0` homogeneous-background convention, independent `G`
normalization, compensated summation, configurable alpha, and independent real
and reciprocal truncation limits. It does not call production PM or TreePM
kernels.

`tests/validation/test_periodic_ewald_reference.cpp` proves convergence under
stricter image/mode limits, alpha independence, integer-box translation
invariance, zero self force, symmetric cancellation, sign, and linear scaling
with `G`.

`tests/validation/test_tree_pm_ewald_accuracy.cpp` reports absolute RMS,
relative L2, median/p90/p95/p99/maximum normalized errors, a documented
denominator floor, mass-weighted net force, seam translation drift, and
separation classes around the Gaussian transition and cutoff. Its diagnostic
matrix also crosses CIC/TSC with deconvolution off/on, monopole/quadrupole with
all three MACs, mesh-axis/diagonal directions, and separations on both sides of
the split transition and cutoff. The certified runtime-default profile obtains
`opening_theta=0.7` from the typed config default rather than duplicating a test
literal. The relative-MAC diagnostic receives the converged Ewald force
magnitude as a deterministic stand-in for compatible previous-epoch history.
Only the exact default profile below is a release accuracy gate; diagnostic
variants print the same metrics without weakening that gate. For the current
FFTW-backed, unsoftened, `G=1`, `a=1` validity metadata, `asmth=1.25`,
`rcut=6.25`, quadrupole, COM-distance profile:

- TSC with matched window deconvolution is the certified configuration and must
  satisfy `relative_L2 <= 1e-2` and `p99 <= 5e-2` on the covered fixtures;
- observed worst certified metrics in the 2026-07-13 run were
  `relative_L2 = 8.059667017e-3` and `p99 = 8.157551429e-3`, both from
  the cancellation-dominated `4^3` DMO initial lattice. That fixture's
  mass-weighted net-force fraction was `5.985844233e-17`;
- on that DMO lattice, Ewald/PM/tree/total force L2 norms were
  `10.16780756` / `8.296819318` / `1.789203432` / `10.08586460`, and its
  integer-image translation drift was `1.523007584e-13`;
- default-setting cutoff-classification diagnostics included
  `relative_L2=6.209e-3` on the
  cutoff-aligned axis fixture, `8.072e-3` on the cutoff diagonal,
  `6.152e-3` immediately beyond cutoff on-axis, and `2.397e-3` at 1.1 cutoff
  on the diagonal;
- the separate rectangular seam fixture had `relative_L2=2.003945987e-4`,
  `p99=2.383921190e-3`, net-force fraction `4.861775315e-17`, and translation
  drift `3.031142721e-14`;
- CIC remains a compatibility/diagnostic profile for this mesh and split: its
  observed uniform result (`relative_L2 = 1.0368e-2`, `p99 = 7.5286e-2`) and
  the earlier split-transition diagnostic (`relative_L2/p99 = 1.2138e-2`)
  exceed the certified targets.

An additional asserted MAC-comparison fixture uses independent targets so no
self/P2P leaf work can hide the opening-policy behavior. It holds quadrupoles,
assignment, split, and force reference fixed and calibrates the geometric and
COM thresholds to `theta=sqrt(alpha)=0.07071067811865475` for
`alpha=0.005`, below the TreePM `l/r<0.08` decomposition-stability envelope.
Its FFTW results are:

| MAC | relative L2 | p99 | visited | accepted | opened | particle pairs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| geometric | `4.541113485e-4` | `6.721724831e-4` | 48 | 44 | 4 | 0 |
| COM-distance | `4.541439201e-4` | `6.721724831e-4` | 56 | 51 | 5 | 0 |
| relative force-error | `4.541347197e-4` | `6.721724831e-4` | 64 | 58 | 6 | 0 |

The test requires every MAC to accept and open internal nodes, beat all-pairs
work, satisfy fixed Ewald error gates, and keep the three relative-L2 values
within 25 percent. This is a matched-policy diagnostic, not a replacement for
the runtime-default `theta=0.7` certification.

The former `rcut=4.5` profile remains a legal explicit compatibility setting
but is diagnostic, not certified: separation fixtures around its hard cutoff
showed roughly nine-percent error. The default `rcut=6.25` gives
`r_cut/r_s=5`; it improves cutoff-transition accuracy at the cost of a larger
short-range search volume, more tree work, and potentially more distributed
target traffic.

These are deterministic small-N validation envelopes, not a universal
high-dynamic-range cosmological accuracy certification. The older
minimum-image direct test remains useful for short-range and split regression
but is not described as an Ewald reference.

## Top-level domain geometry freshness and fallback

The compact one-leaf-per-rank top-domain map is **derived routing geometry**,
not restart truth. It has two independent identities:

- `DecompositionEpoch` advances only on a committed ownership change (migration
  / rebalance / restart restore). Drift does **not** advance it.
- `GravitySourceGeneration` (owned by authoritative `SimulationState`) advances
  whenever the physical source set mutates, including drift that moves sources.

Workflow lifecycle:

1. Seed leaves are installed once at segment/restart boundaries and again only
   after an ownership decomposition change, stamped with the current
   `GravitySourceGeneration`. An install replaces two distinct sets inside
   `GravityRuntime`: the stable **seed leaves**
   (`m_authoritative_top_domain_seed_leaves` — decomposition-local leaf
   identities and SFC intervals, allowed to include currently-empty groups)
   and the **current published leaves** (`m_authoritative_top_domain_leaves`).
2. Before each force solve (after the compact source view is rebuilt, and only
   on paths that will actually solve), `GravityRuntime` performs an O(N_local)
   bound refit that consumes the **seed set**, never the previously published
   result: seed ownership/epoch/SFC intervals are retained; only per-leaf AABB
   bounds, entity counts, and the published generation stamp advance. The
   published set may omit currently-empty seed groups (finite bounds, no
   zero-entity routing packets), while the seed set keeps their SFC intervals
   so repeated empty-group disappearance cannot progressively coarsen the
   routing partition. A current source whose SFC key falls outside every seed
   interval is still assigned to the nearest seed leaf (ownership authority
   exceeds old SFC range membership; `out_of_seed_range_source_count` records
   it) and never dropped. Drift changes the source generation but not the
   decomposition epoch and never rewrites the seed set.
3. `commitParticleDecompositionChange()` explicitly invalidates geometry
   freshness when an ownership epoch advances: between commit and the
   subsequent install, `authoritativeDomainGeometryMatches(...)` is false and
   the published routing view is empty, so stale geometry is unusable. This
   does not clear ownership data, migration state, or the epoch.
4. `TreePmCoordinator` selects authoritative leaves only when the generation
   stamp equals `options.source_generation` **and** the existing conservative
   coverage check still covers every local source. Any failure falls back to
   the local gravity-tree root packet and records a specific
   `TreePmDomainGeometryFallbackReason`:
   `kNoGeometryInstalled`, `kStaleSourceGeneration`,
   `kDecompositionEpochMismatch`, `kSourceCoverageFailure`, or
   `kGeometryPreparationFailure`.

A graph-cache hit is never geometry validation. Geometry is not serialized in
restart/snapshot schema; after restart the seed path reinstalls it. Single-rank
and MPI-disabled builds remain valid; there is no hard-coded world size in the
refit.

P6 hierarchical freshness rule: all canonical DMO source positions are
advanced to every fine evaluation epoch. `GravitySourceGeneration` advances
once after each all-source drift, so P2 geometry refresh and collective fallback
continue to consume canonical sources, with no predicted-source token alias.
Lazy inactive prediction would require a separate future contract; it is not
part of this implementation.

## Residual traversal counters and timing truth

Residual traversal work is recorded in two non-overlapping counter bundles:

- `local_owned_targets` — serial local solve and local work overlapped with
  request transport;
- `incoming_remote_targets` — peer-evaluated incoming remote targets.

Identity (exact, same solve):

```text
residual_pair_evaluations
  == local_pair_evaluations + incoming_remote_pair_evaluations
```

`tree_profile.particle_particle_interactions` and the combined
`tree_profile` visited/accepted/opened/cutoff counters remain the **sum** of
both bundles (compatibility with prior aggregate semantics).
`remote_pairs_pruned_by_bounds` is unchanged and still counts only the
`!skip_self` path. Standalone `TreeGravitySolver` PPI contributions are not on
the TreePM residual path.

Timing split (all additive into existing profile totals):

- `PmProfileEvent.total_ms` / `profile.pm_profile.total_ms` spans the PM phase
  entry (including long-range refresh when taken) through the tree short-range
  start — PM total alone, not including short-range traversal.
- `tree_wall_ms_recent` in the workflow event is
  `tree_short_range_ms` alone; `pm_wall_ms_recent` is `pm_profile.total_ms`
  alone.
- Remote-phase split:
  `incoming_request_decode_validation_ms` (wire decode, record-count,
  peer/epoch/identity validation, duplicate-identity hashing, finite-value
  checks, and response structure preparation), then
  `incoming_remote_target_compute_ms` (validated incoming target force
  evaluation against this rank's local tree only, plus the inseparable direct
  acceleration write into the prepared response slot — no decode, hash,
  validation, serialization, or communication), then
  `incoming_response_encode_pack_ms` (response byte encoding, size check, and
  direct pack into the response send payload);
  `protocol_validation_ms` (response count/displacement layout and response
  payload buffer sizing), `protocol_consensus_ms`
  (`coordinate_protocol_failure` duration), `response_exchange_ms` (response
  `MPI_Neighbor_alltoallv` call only).
  `let_remote_traversal_ms` remains a compatibility alias equal to
  `incoming_remote_target_compute_ms`. `let_communication_ms` is the transport
  subset (count exchange, request wait, and response exchange); protocol
  consensus is reported separately in `protocol_consensus_ms`.

## OpenMP residual execution (P3)

Short-range residual evaluation parallelizes **between targets only**, in
validated blocks selected by `numerics.treepm_residual_block_size` (default 64). Three
regions use the same worker-safe evaluator:

1. non-distributed local residual;
2. local work overlapped with sparse request transport;
3. pure incoming remote-target compute (after main-thread decode/validation).

Contract:

- one immutable local tree shared by all workers; no MPI calls from workers
  (`MPI_THREAD_FUNNELED` only);
- shared mutation is limited to unique active-slot writes and one integer
  counter bundle per planned worker, plus separate optional spatial scratch;
  counters merge exactly after join. The deterministic floating diagnostic
  remains one `sum_sq` value per logical
  configured target block and is reduced in fixed block order; no atomics appear in
  hot node/pair loops;
- each worker owns one bounded DFS stack slot of `S = 1 + 7 * D` entries
  (`D = TreeGravitySolver::maxDepth()`), backed by one contiguous
  `m_worker_stack_storage` of `T * S` `TreeLocalIndex` slots. The tree builder
  enforces `kMaximumTreeDepth`; preflight uses that same bound and runtime
  rejects a recorded depth above it;
- source softening is reused from the immutable resolved build lane. Target
  softening spans, species tags, finite/non-negative values, and independent-
  target sidecar requirements are validated on the main thread before OpenMP;
  workers use the allocation-free, non-throwing unchecked resolver against that
  validated view. There is no full-active `double` softening lane;
- `local_short_range_sum_sq` is reduced deterministically (block partials in
  block order); integer counters are `O(worker_count)` and are prepared before
  the distributed request is posted;
- `openmp_observed_workers` is the maximum actual team size observed in the
  local, distributed local-overlap, and incoming-target OpenMP regions;
  serial execution reports one;
- worker exceptions are captured under named OpenMP critical sections and
  rethrown on the main thread after join (they never escape the region);
- builds without OpenMP (`COSMOSIM_HAVE_OPENMP=0`) keep a serial path with the
  same kernel and the same numerical order.

P6.1 times each existing logical block with two `steady_clock::now()` reads,
outside its target loop, in all three regions (including serial execution).
The executing worker adds block duration to `elapsed_work_ms` and keeps
`block_work_ms_max`; join-time merges sum the former and take the maximum of
the latter separately for local-owned and incoming targets. Timing includes
target preparation, traversal, force-slot writes and local block diagnostics;
it excludes scheduling waits and MPI. Worker-region time summaries retain
their existing min/max/sum semantics over these block-duration sums. Clock
reads are O(blocks), and the target evaluator performs no wall-clock reads.

Hot `TreePmTraversalCounters` remain `alignas(64)` and contain only scalar
work counters and timing summaries. Separate `alignas(64)` worker scratch owns
the fixed 64-bin work/target arrays. It is allocated/grown and zeroed once at
the outer workspace-preparation scope only when spatial feedback is enabled,
before distributed request posting. Owned-target integer work accumulates
across local batches; incoming evaluation receives no histogram pointer.
After all joins, fixed worker/bin-order reduction precedes the existing MPI
integer sums and unchanged 0.5-decayed work/target history update. With P5-W
disabled, execution performs no histogram allocation, clearing or updates.
Capacity retained from an earlier enabled solve remains separately reported.

Diagnostics provenance: `openmp_compiled`, `openmp_configured_workers`,
`openmp_observed_workers`, `residual_local_target_count`,
`residual_incoming_target_count`, and
`residual_worker_scratch_high_water_bytes`. The current preflight terms are:

```text
M_worker_stack = T * (1 + 7 * kMaximumTreeDepth) * sizeof(TreeLocalIndex)
M_worker_counters = T * sizeof(TreePmTraversalCounters)
M_optional_spatial = T * kTreePmResidualSpatialCounterBytes
M_block_diagnostics = ceil(A / 64) * sizeof(double)
```

Preflight conservatively includes `M_optional_spatial` even when disabled,
preserving coverage of the prior worker-scratch envelope without a new input
policy. The retained runtime report exposes actual worker-stack, hot-counter,
separate optional spatial-counter and block-diagnostic capacities. Distributed
target metadata is bounded by the current communication batch; no `double[A]`
target-softening allocation is retained or admitted.

## LET exchange memory estimate

`gravity::estimateTreePmExchangeMemory` (single auditable arithmetic source in
`gravity_memory.cpp`) derives the short-range exchange peak from
`tree_exchange_batch_bytes` with checked arithmetic and the existing
`planSparseTreePmRound` clamp:

```text
b = floor(B / 96)                    # targets per peer per batch
d = max(R - 1, 0)                    # peer rounds
M_wire = 2 * d * b * (96 + 88)       # request + response payloads
known  = wire
       + structured request storage
       + count/mask/accumulator/metadata terms
       + transient codec workspace b * (2*(96+88)+24)
```

Wire record widths are the shared public constants
`kTreePmShortRangeRequestWireBytes = 96` and
`kTreePmShortRangeResponseWireBytes = 88` (wire v2: the response appends the
serving rank's nearest direct-pair separation squared, `+inf` when none, for the
hierarchical encounter timestep criterion; it was 80 in wire v1). Host `sizeof` of the request/response
packet structs equals these widths, and that identity is compile-time enforced
by `static_assert` next to the packet definitions in `tree_pm_coupling.cpp`;
if a supported ABI ever pads the structs, the build fails rather than letting
wire record bytes silently stand in for host object bytes. Unordered-set
bucket/allocator overhead in the codec is intentionally excluded and documented
as an uncertainty rather than double-counted. Runtime events split high-water
into `let_wire_buffer_high_water_bytes` (four reusable payload buffers) and
`let_known_workspace_high_water_bytes` (wire + structured + metadata +
transient codec, read from retained `vector.capacity()` for CHUÍ-owned
storage where accessible; the transient codec term remains a modeled
conservative upper envelope). `let_high_water_bytes` remains a wire-only
compatibility alias. Preflight (a conservative "might require" estimate) and
runtime retained-capacity high-water are conceptually distinct and need not be
numerically identical. `MemoryGovernor` remains the sole admission authority;
the estimate feeds the existing preflight, never a second governor.

## Diagnostics and validation entry points

`TreePmDiagnostics` reports local source/active-target/tree-node counts, global
empty-source and empty-target rank counts, remote hierarchy packets, unique
communicating peers, PM solve/reuse counts, cached halo-value count, and local
FFT slab dimensions. It also reports split/cutoff scales, split composition
error, cutoff pruning, residual pair work with the local/incoming-remote split
above, request/response packets and bytes, batch/peer participation,
zero-request targets, peer pressure imbalance, zoom gather bytes, local/remote
residual norms, remote-phase timer splits, LET wire/workspace high-water
split, and domain-geometry freshness/fallback fields
(`domain_geometry_source_generation`, `current_gravity_source_generation`,
`domain_geometry_fresh`, `domain_geometry_fallback_used/reason`,
`domain_geometry_uncovered_source_count`). Tree profiling separately counts
builds, multipole refreshes, visited/accepted/opened nodes, and
particle-particle interactions. The coordinator memory report includes
reusable PM, tree, active-set, unwrapped-coordinate, and exchange workspaces.

Relevant gates are:

- `integration_tree_pm_coupling_periodic` and its MPI np2/np3/np4 entries;
- `integration_pm_slab_halo_exchange_mpi_two_rank`, `_three_rank`, and
  `_four_rank`;
- `validation_periodic_ewald_reference`;
- `validation_tree_pm_ewald_accuracy`;
- `validation_phase2_mpi_gravity_{single,two,three,four}_rank`;
- `validation_dmo_zeldovich_workflow_single_rank` and MPI np2/np3/np4.

See `docs/gravity_production_readiness.md` for current pass/limitation status.

## Public-interface migration notes

- P6.1 removes the raw `spatial_work` / `spatial_targets` arrays from
  `TreePmTraversalCounters`; planner consumers continue to use
  `TreePmDiagnostics::spatial_work_per_target` and
  `spatial_work_history_solves`. Native consumers must rebuild for the smaller
  counter layout. `elapsed_work_ms` now sums logical-block durations and
  `block_work_ms_max` gives the largest block duration per target family.
  `kTreePmResidualCounterBytes` remains the actual hot `sizeof`;
  `kTreePmResidualSpatialCounterBytes` is the separate optional scratch width.
  No config, force, snapshot/restart or serialized event-field schema changes
  accompany this cleanup. Existing worker timing values cover block work
  rather than individually timed traversals; qualification remains pending.
- `TreePmForceAccumulatorView` adds optional previous-acceleration and explicit
  target-position spans. Existing source-indexed callers retain prior behavior.
- `TreePmOptions` adds `decomposition_epoch` and `force_epoch`; distributed
  workflow callers must supply coherent runtime values.
- `TreePmDiagnostics` adds explicit local/global occupancy, hierarchy/peer,
  PM solve/reuse/halo, local slab-dimension, split residual-pair
  (local/incoming-remote), remote-phase timer, LET wire/workspace high-water,
  and domain-geometry freshness/fallback counters. Callers using aggregate
  initialization or mirroring this public type must account for the appended
  fields.
- `TreePmOptions` adds `authoritative_geometry_source_generation` (geometry
  freshness stamp) alongside the existing `source_generation`; both must be
  supplied coherently or authoritative routing falls back with
  `kStaleSourceGeneration`.
- `TreePmDomainGeometryFallbackReason` and its name helper are additive public
  enum/string APIs for the fallback contract above.
- `parallel::refitAuthoritativeTopDomainLeaves` and
  `TopDomainGeometryRefitDiagnostics` are additive public refit APIs; seed
  leaves retain owner/epoch; empty/non-finite inputs fail closed or omit empty
  leaves rather than publishing non-finite bounds. The refit's seed input is
  the stable decomposition-local seed set, not the previously published
  result, so empty-group omission cannot erode SFC partition identity across
  refreshes.
- `workflows::GravityRuntime::installAuthoritativeTopDomainLeaves` takes the
  source generation the leaves cover and replaces both the stable seed set and
  the current published set;
  `commitParticleDecompositionChange()` invalidates geometry freshness so
  `authoritativeDomainGeometryMatches` is false until reinstall;
  `authoritativeDomainGeometryMatches` queries freshness. No restart/snapshot
  schema change accompanies this interface; geometry remains derived state.
- Short-range wire width constants `kTreePmShortRangeRequestWireBytes` and
  `kTreePmShortRangeResponseWireBytes` are public additive constants.
- Hierarchy packet wire version, exchange sequence, force epoch, and geometry
  frame are now part of the distributed contract. External test tools that
  constructed packet structs without them must use the version-1 defaults and
  current epochs.
- `core::newtonGravitationalConstantCode(const UnitSystem&)` is an additive
  public units helper declared in `include/cosmosim/core/units.hpp`. Workflow
  integrations that formerly hard-coded
  `G_code=1` must construct the frozen-config `UnitSystem` and call this helper;
  standalone gravity tests may still choose an explicit dimensionless `G`.
  This changes no snapshot/restart dataset and adds no second config lane.

## P6 acceptance, source identity and tree retention

**SOURCE-IMPLEMENTED / VALIDATION PENDING.** `kStrictReference` remains default:
the selected existing MAC, all safety guards, quadrupole-only internal
acceptance and `width/r < 0.08` remain. Internal monopoles descend in this policy.
Adaptive mode may accept internal monopoles or quadrupoles, retaining the
existing complete second moment/trace force expansion, source softening,
self/inside geometry and complete-node cutoff containment.

For CHUI's screened Plummer residual,

\[
 F_i=GM d_i f(r),\quad
 f=(r^2+\epsilon_{pair}^2)^{-3/2}+(S-1)/r^3,
\quad S=\operatorname{erfc}(q)+2q e^{-q^2}/\sqrt\pi,
\quad q=r/(2r_s).
\]

The adaptive accuracy inequality is

\[
 GM\max\left(l^2/r^4,\frac{\rho^2}{2}
 [3|f'|+|r f''-f'|]\right)\le\alpha |A_{previous}|,
\quad l=2h,\quad\rho=\sqrt3h+|COM-center|.
\]

Derivatives are evaluated for this exact residual kernel. The second-moment
scale is charged even for quadrupoles. This is an accuracy **proxy**, not a
proven bound on the omitted quadrupole remainder or a qualified alpha. `A` is
the prior/reference unscaled total acceleration. It is never replaced by a
fake floor acceleration. Missing/incompatible, nonfinite or <= floor history
uses the existing COM-distance geometric criterion. The independent adaptive
maximum angle defaults to 0.25 and is validated <=0.5; it is an unqualified
safeguard, not a replacement of strict 0.08 or a claimed scientific constant.
Rejection reasons are separately observable and may overlap.

Workflow source generation now reaches tree build directly. Trustworthy
nonzero generation skips the O(N) content hash; zero retains the old hashing
fallback. Stale-source checks are retained. Exact-generation reuse also
requires count, softening/build options, boundary/box, ownership epoch and
nonzero stable source-row layout generation. Saturated generation is not reuse
authority. Generic views with unknown row layout rebuild.

Motion refit keeps the original periodic unwrap anchor and checks every source
strictly inside its original octant leaf cell. Original cells are reconstructed
from four root scalars and immutable octant child slots, without new per-node
or per-source metadata. Failure rebuilds using the existing largest-gap unwrap
and Morton/radix path. Success preserves child/source membership, refits
bounding cubes bottom-up with outward rounding, resets COM/moments and
recomputes all multipoles, then publishes a new tree-build generation. Logical
Morton keys are cleared after refit; the leaf permutation remains legal but is
not advertised as a fresh coordinate Morton order. Refitted bounds may differ
from a fresh build, so force-error/equivalence and restart roundoff need explicit
qualification. No bitwise claim is made for optional refit across restart.

`TreeGravityProfile::rebuild_reason` codes: 0 disabled, 1 invalid tree, 2 source
identity changed/unavailable, 3 frame/ownership change, 4 options/softening,
5 identical reuse, 6 certified refit, 7 original leaf escape, 8 unknown/changed
row layout. Rejected membership time survives a subsequent build. Reference
full rebuild is retained with both switches false.

Accumulator API migration: optional PM-only outputs are active-slot indexed;
optional Tree-only outputs are source-row indexed and require local indexed
targets. Ordinary callers leave these spans empty and retain total-force
behavior. `short_range_only` is an integrator split surface, forbids PM
refresh/capture and zoom correction, and returns short-range-only diagnostics.
Existing explicit PM cache-reuse compatibility rules are unchanged.

## Unified recovery implementation (2026-10-10)

Production DFS pushes children in reverse octant order and visits increasing
octants. There is no child-distance sorting or child-order geometry. This changes
summation order, preserving identities, pair softening, MACs and inclusive cutoff.
Worker stacks keep the existing `1+7*max_depth` bound and overflow exception.

Node-center nearest-image deltas are evaluated once per visited node with all
three axis lengths. AABB cutoff/containment uses squared distances, with the
former square-root comparison retained within eight machine epsilons of the
cutoff and for extreme/nonfinite squared scales. Leaves omit COM/MAC geometry;
internal nodes reuse the same center deltas for containment. General `nearbyint`
minimum-image behavior, including half-box ties and multiple wraps, is retained.

The analytic pair law is still
`(r*r+epsilon*epsilon)^(-3/2) - L(r,a)/r^3`, with `epsilon=max(epsilon_i,epsilon_j)`.
`S(0)=1`, `L(0)=0`. For `q=r/(2a)<1/8`, evaluate
`L/q^3 = (4/sqrt(pi))*sum(n=0..7, (-q*q)^n/(n!*(2n+3)))`.
The long-range coefficient limit is `1/(6*sqrt(pi)*a^3)`. At exactly zero
radius with positive finite softening the residual coefficient is
`epsilon^(-3) - 1/(6*sqrt(pi)*a^3)`, evaluated by a scaled subtraction to avoid
`inf-inf`. At zero softening the established coincident coefficient remains
zero. Production counts coincident direct interactions and returns their exact
zero vector before multiplying an overflowing coefficient by zero. Analytic and
lookup selections share this treatment; nonzero residuals and cutoff are unchanged.
Checked APIs reject nonfinite/negative inputs; unchecked functions require
prevalidated inputs. Representational overflow at extreme finite physical
scales is not an additional supported arbitrary-precision contract.

`TreePmGaussianCoefficientTable` is an additive public diagnostic/kernel API.
It tabulates `H(t)=L(sqrt(t))/t^(3/2)`, `t=(r/(2a))^2`, using 4096 Hermite cells
on `[0,64]`. Near zero the same analytic series is used. Endpoints/outside values
use analytic fallback. The standalone target is **absolute H error <=5e-11**;
this is an unexecuted regression budget, not a relative residual-force bound.
The enabled coordinator lazily reserves and owns one immutable 65,552-byte
array object (`sizeof` is authoritative), retained until coordinator destruction.
Initialization precedes OpenMP traversal and participates in collective failure
handling. Split changes use new dimensionless coordinates/normalization without
rebuilding the table. Uniform and variable pair softenings share this table.

`TreePmOptions::gaussian_pair_lookup_enabled=false` and
`full_mac_diagnostics=true` preserve analytic/full accounting defaults. The
additive `residual_block_size` direct-API option is bounded to `[1,4096]`, defaults
to 64. Production exposes the narrower validated set 16/32/64/128/256 through
the typed configuration path; memory preflight uses that actual block size.
Nondefault embedded block sizes require the caller to pass the matching
`GravityMemoryEstimateInput::residual_block_size`. Forces retain per-target
source order; changing block size can change floating diagnostic reduction
grouping and scheduling, so it belongs in normalized experiment provenance.

In adaptive acceptance, the Newtonian operand of the residual-error `max` is
checked before expensive derivatives. A proven rejection preserves its existing
relative-MAC counter. Fast diagnostics additionally skip the entire selected MAC
after independent rejection and explicitly count skipped MACs. Before COM
geometry, target-containing or cutoff-intersecting nodes can descend directly;
`skipped_geometry_evaluations` distinguishes unevaluated softening/envelope/MAC
tests. Full diagnostics retain overlapping failed-guard counts.
No force-history floor is fabricated. Existing `alpha=0.005`, maximum angle
`0.25`, strict envelope, softening, and cutoff containment remain enforced.
Quadrupoles retain analytic derivative-consistent radial terms and the full
second-moment trace; the pair table never supplies multipole derivatives.

## Group-walk leaf culling, target sizing and telemetry (2026-10-25, untested)

Group interaction lists now keep one entry per accepted source leaf with its
cell center (minimum-image offset from the group center) and half width. Each
target rejects a leaf when `sum_a max(0, |minimg(c_a - p_a)| - h)^2 >
r_cut^2 (1 + 1e-9)`, a conservative lower bound on every pair separation in
that leaf; surviving tree-contiguous leaves are re-coalesced for the
vectorized pair kernel. `numerics.treepm_group_walk_target_size` (0 = legacy
one-leaf groups) groups up to N tree-order targets within one residual block
and an AABB extent of `0.5 r_cut`. Owned targets use the group walk in both
the local-only and the distributed (LET-overlapped) residual paths. Incoming
remote LET targets are grouped per peer after a Morton sort (21 bits per axis
on a lattice of the group extent) and use the same group evaluator through a
sample/emit adapter; ineligible groups fall back to the per-target walk
(`incoming_group_walk_*` counters).

Round 2 (2026-10-25, unqualified): the tree exposes `topologyGeneration()`
(changes only on full builds, unchanged by certified refits), so the
row → tree-position map and the group-walk leaf map are rebuilt only on
topology change (`topology_cache_rebuilds/reuses`). Large target visit orders
use a stable two-pass 16-bit radix sort. Uniform-mass source populations
release the tree-order mass lane and pass the species mass as a scalar
(same per-pair expression, bitwise identical). PM particle routing reserves
are sized from the measured remote x-plane workload (per-peer per-round
record maxima, rank-independent round model) via `planDensityRouting` /
`planForceInterpolationRouting`, planned before the communication arena is
configured; the PM-only path measures and reduces locally. The PM backend
works in place (`real` aliases the padded `fourier` array, rows of
`2(nz/2+1)` doubles). `gravity.treepm_let` events report
`group_walk_status`, groups, grouped/fallback targets by reason, candidate and
culled leaves and candidate particles. The PM long-range force may use the
opt-in `numerics.treepm_pm_force_gradient = finite_difference_4` backend
(see `docs/configuration.md`). Record:
[campaign 2026-10-25](repair/treepm_optimization_campaign_20261025.md).
