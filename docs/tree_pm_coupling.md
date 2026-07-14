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
2. sort the wrapped lane;
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

## Distributed hierarchy and request/response protocol

The current distributed short-range algorithm is correctness-first and aimed
at workstations and small clusters:

1. each rank builds a tree from rank-owned sources;
2. every rank contributes a bounded breadth-first hierarchy summary;
3. target owners use the periodic bounds and `r_cut` to select peer ranks;
4. selected targets are exported in bounded batches;
5. each destination evaluates the target against its local source tree;
6. validated partial accelerations return to the target owner.

An empty local tree emits one explicit zero-source root sentinel. It carries
`source_count=0`, zero mass, the current geometry frame, decomposition epoch,
force epoch, and exchange sequence. This distinguishes "participated with no
sources" from missing rank coverage. The sentinel never intersects a cutoff
query and therefore creates no request or force.

Hierarchy records use a version-1, fixed 152-byte, explicit little-endian wire
encoding. Short-range requests and responses use version-1 explicit
little-endian records of 96 and 80 bytes respectively. They do not transmit raw
C++ object representations or uninitialized padding. Protocol identity covers:

- origin, destination, and source ranks as applicable;
- exchange sequence;
- decomposition and force epochs;
- batch token and request ID;
- target identity;
- previous-force-scale presence;
- target coordinates and softening;
- response acceleration components.

Before an exchange, ranks reach consensus on exchange sequence, epochs, and
batch-byte policy. Decoders reject wrong versions, misaligned/truncated payloads,
wrong peers, stale/mixed epochs, non-finite numerical fields, invalid flags,
duplicate hierarchy IDs, missing/multiple rank roots, out-of-range batch slots,
unexpected responses, duplicate responses, and incomplete expected-response
coverage. Count multiplication, cumulative displacements, byte alignment, and
`MPI int` limits are checked. Every global batch enters request and response
collectives even when a rank has zero records, using a safe non-null empty
buffer.

At the public solve boundary the actual active MPI world also fingerprints the
PM mesh `nx/ny/nz`, physical/operator options, protocol epochs, layout mode, and
softening metadata before PM collectives. This prevents different per-rank mesh
shapes or controls from selecting incompatible FFT and exchange paths. An
MPI-enabled process without an active MPI session remains a serial library
caller. A collective full-domain serial reference call inside a multi-rank MPI
job evaluates its short-range residual locally on every rank and does not enter
the distributed hierarchy exchange.

Failure coordination starts before that wire protocol. A rank-local exception
while constructing the periodic unwrapped source lanes or building the local
tree is caught and reduced across the actual MPI world before any rank enters
the residual hierarchy/request exchange. Constructor-local PM-layout metadata
is likewise not rejected on one rank in isolation; layout/context coherence is
checked at solve entry and voted before collective work. Thus a bad rank makes
every peer leave the same preflight phase instead of stranding peers in the
next collective.

The same solve-level vote covers post-exchange PM halo-cache commit, allocation
and filling of compact active-position/PM-force/zoom-correction target buffers,
and worst-case residual traversal-stack reservation. The traversal stack
reserves the complete local-node bound before visiting a target, so traversal
cannot introduce an unvoted growth allocation between later collectives.

Per-peer count/displacement arrays, request/response payloads, remote partial
accelerations, and expected/received response-count lanes are owned by the
reusable `TreeExchangeWorkspace`. World-size-dependent allocations occur
inside the coordinated exchange preflight, and batch-dependent growth occurs
inside the already coordinated request-preparation phase. These capacities are
reported by `TreePmCoordinator::memoryReport()` as transient MPI buffers. This
removes unvoted per-batch response bookkeeping allocations without changing
force values, wire bytes, or deterministic ordering.

The hierarchy control plane is still an `MPI_Allgatherv`, and request/response
data still use communicator-wide `MPI_Alltoallv` count/data phases. Although
cutoff pruning avoids many payload records, this is not a locally essential
tree (LET), sparse neighborhood transport, or large-cluster scaling claim.

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
an explicit lower-level reuse request. Its compatibility signature contains
the force epoch, force-evaluation scale factor, code gravitational constant,
split scale, x/y/z box lengths, assignment scheme, boundary condition, PM
decomposition mode, and window-deconvolution flag. A missing or incompatible
signature makes reuse fail coherently; it is not silently changed into a solve.
Explicit refresh/reuse votes are reduced first, and a mixed vote throws before
any rank enters PM density or FFT collectives. This signature is an invalidation
guard, not a predictor for cadence greater than one.

Tree topology and hierarchy packets are rebuilt for each force call. Their
decomposition and force epochs prevent stale packet reuse after migration or a
field-generation change. The workflow decomposition epoch advances only after
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

## Diagnostics and validation entry points

`TreePmDiagnostics` reports local source/active-target/tree-node counts, global
empty-source and empty-target rank counts, remote hierarchy packets, unique
communicating peers, PM solve/reuse counts, cached halo-value count, and local
FFT slab dimensions. It also reports split/cutoff scales, split composition
error, cutoff pruning, local/remote pair work, request/response packets and
bytes, batch/peer participation, zero-request targets, peer pressure
imbalance, zoom gather bytes, and local/remote residual norms. Tree profiling
separately counts builds, multipole refreshes, visited/accepted/opened nodes,
and particle-particle interactions. The coordinator memory report includes
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
  PM solve/reuse/halo, and local slab-dimension counters. Callers using
  aggregate initialization or mirroring this public type must account for the
  appended fields.
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
