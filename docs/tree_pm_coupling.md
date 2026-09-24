# TreePM coupling

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
presence, and returned acceleration components.

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

Production configuration currently requires both
`numerics.treepm_update_cadence_steps = 1` and
`numerics.hierarchical_max_rung = 0`. Every integrator-issued,
rank-coordinated production force-refresh surface rebuilds the long-range PM
field. Cadence greater than one lacks a validated predictor/interpolator;
mixed-rung KDK lacks per-element kick/drift epochs. Both unsupported semantics
therefore fail at config validation rather than being presented as production
maturity. The integrator owns the PM synchronization event, kick opportunity,
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

Tree topology is rebuilt for each current production force call. Its explicit
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

Freshness-token design note (hierarchical KDK, documentation only): current
P2 uses `GravitySourceGeneration` as the routing-geometry physical-state
freshness token. That is valid for the current production all-active/rung-zero
workflow, where every source moves with the single global timestep. A future
mixed-rung KDK implementation may predict inactive sources to a
force-evaluation epoch without advancing canonical source generation the same
way; hierarchical KDK must then revisit routing-geometry freshness to include
the prediction/evaluation epoch. No scheduler state for that exists today.

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
compile-time blocks of 64 (`k_residual_block_size`; not a config key). Three
regions use the same worker-safe evaluator:

1. non-distributed local residual;
2. local work overlapped with sparse request transport;
3. pure incoming remote-target compute (after main-thread decode/validation).

Contract:

- one immutable local tree shared by all workers; no MPI calls from workers
  (`MPI_THREAD_FUNNELED` only);
- shared mutation is limited to unique active-slot writes and one integer
  counter bundle per planned worker; counters merge exactly after join. The
  deterministic floating diagnostic remains one `sum_sq` value per logical
  64-target block and is reduced in fixed block order; no atomics appear in
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

Diagnostics provenance: `openmp_compiled`, `openmp_configured_workers`,
`openmp_observed_workers`, `residual_local_target_count`,
`residual_incoming_target_count`, and
`residual_worker_scratch_high_water_bytes`. The current preflight terms are:

```text
M_worker_stack = T * (1 + 7 * kMaximumTreeDepth) * sizeof(TreeLocalIndex)
M_worker_counters = T * 7 * sizeof(uint64_t)
M_block_diagnostics = ceil(A / 64) * sizeof(double)
```

The retained runtime report exposes the actual worker-stack, worker-counter,
and block-diagnostic capacities. Distributed target metadata is bounded by the
current communication batch; no `double[A]` target-softening allocation is
retained or admitted.

## LET exchange memory estimate

`gravity::estimateTreePmExchangeMemory` (single auditable arithmetic source in
`gravity_memory.cpp`) derives the short-range exchange peak from
`tree_exchange_batch_bytes` with checked arithmetic and the existing
`planSparseTreePmRound` clamp:

```text
b = floor(B / 96)                    # targets per peer per batch
d = max(R - 1, 0)                    # peer rounds
M_wire = 2 * d * b * (96 + 80)       # request + response payloads
known  = wire
       + structured request storage
       + count/mask/accumulator/metadata terms
       + transient codec workspace b * (2*(96+80)+24)
```

Wire record widths are the shared public constants
`kTreePmShortRangeRequestWireBytes = 96` and
`kTreePmShortRangeResponseWireBytes = 80`. Host `sizeof` of the request/response
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
