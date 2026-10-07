# P6 pre-production performance architecture implementation

**SOURCE IMPLEMENTATION ONLY**

**NO BUILD EXECUTED**

**NO TESTS EXECUTED**

**NO SCIENTIFIC VALIDATION EXECUTED**

**VALIDATION REQUIRED IN FOLLOW-UP CAMPAIGN**

Every new implementation here is **SOURCE-IMPLEMENTED / VALIDATION PENDING**.
The requested filename uses campaign label 20261006; the workspace session date
is 2026-10-05. This report contains no new runtime/performance/science evidence.

P6.1 follow-up (2026-10-07): per-target clock reads were removed; the existing
64-target blocks now own work timing. Optional spatial-worker histograms are
separate from aligned hot counters. No force mathematics, P5-W weighting,
configuration or default changed. **SOURCE IMPLEMENTED; BUILD NOT RUN;
TESTS NOT RUN; BENCHMARKS NOT RUN; VALIDATION PENDING.** Details:
[P6.1 cleanup note](p6_1_treepm_hotpath_instrumentation_cleanup_20261007.md).

## Input and execution scope

Input revision: `a88de7ec1098b89c8d8e6679431a8466067ebc62`.
Initial tracked diff/status was clean. Unrelated untracked `opencode.jsonc` and
`session-ses_efe3.md` were preserved and excluded from delivery. Existing source
edits survived the accidental interruption; implementation resumed in place.
Mode: feature implementation with directly related repairs. Repository/module
ownership, AGENTS, current status, contribution/build workflow, runtime truth,
TreePM, scheduler, restart, memory, P1–P4 and M48 authority were inspected before
source changes. No Git reset or namespace migration occurred.

No CMake configuration, compiler/linker, CTest, unit/integration tests, MPI
launch, CHUI run, force comparison, benchmark or cosmological validation ran.
No test suites were added. Existing validation floors are deferred by the
explicit campaign instruction, not waived for release.

## Retained reference and preservation contracts

Default selection remains strict TreePM (including quadrupole-only internal
acceptance under historical width/distance <0.08), global/all-active KDK,
full tree rebuild, and the existing compact decomposition policy. Original
relative/geometric MAC selections remain available within strict mode.

Local/incoming pair accounting remains disjoint with exact combined identity.
P2 source-fresh geometry remains independent of ownership epoch, with
collective conservative fallback and reusable QueryScratch. Both existing
OpenMP target traversals and worker-private bounded DFS scratch remain.
MPI calls stay on the funneled owner. Compact decomposition/startup records,
streamed work recomputation, dense local indices, bounded communication arena,
borrowed DMO sources, compact rung-zero scheduler, bounded migration, streaming
restart verification, compact PM/spectral storage and homogeneous cold-state
elimination are retained. No per-worker tree/mesh, N-by-rank history, full-peer
normal routing, fast-math, new dependency, parser, or memory governor was added.

## A — observability

Before: TreePM had established phase/LET/P1 accounting but ambiguous accepted
nodes, sparse preparation timing and insufficient acceptance/worker causality.
Now: trusted identity/hash, source/periodic preprocess, softening, Morton,
topology, multipoles, node/depth, rebuild/refit reason and time are exposed.
Local/incoming targets, node visits/opens, internal multipoles, leaves, actual
direct pairs, cutoff and individual rejection families are separate. Existing
aggregate counters remain. Rejection counters overlap; opens count unique
descents. Timer parents/children are explicitly documented.

Aligned worker bundles contain current-force targets, visits, pairs, multipoles,
summed block work time and maximum block duration. P6.1 stores the two optional
64-bin count arrays separately, prepared only when spatial feedback is enabled.
Summaries expose bounded nonempty worker regions, extrema and sums/derivable
means. There are no node/pair atomics or retained run-wide block tables.
PM axis factor/inverse/normalization,
density-routing wait, halo and total inclusive timing augment existing fields.
LET adds arena high-water and export multiplicity alongside existing geometry,
candidate-peer, graph/cache and transport metrics. Fixed-name lifetime counters
preserve report totals when event detail is evicted.

Files: profiling headers/source, TreePM/Tree/PM headers/source, gravity-memory
model, gravity runtime, analysis runtime, reference workflow/time coordinator.
Memory: bounded current-force/worker owners; real counter sizeof enters the
existing worker budget. Risks: schema consumers, overlapping timers and MPI
worker-region interpretation. Follow-up: matrix C/D/H/P and event semantics.

## B — cheap defects and bounded history

Before: workflow generation did not reach tree build; a zero-generation content
hash could run unnecessarily. Adaptive-global metadata incorrectly said
unsupported. Events/cadence/stage traces retained run-length detail.
Now: authoritative generation reaches build, zero retains hashing/stale-source
fallback, and identity path is reported. Capability metadata separates supported
adaptive global KDK, default rung-zero, and provisional opt-in hierarchy.

Profiler events use a global recent-256 ring with lifetime severity counters and
bounded ordered report view; event recording has a control-path mutex, while
scope/counter shards remain. Cadence detail retains 256 decisions, stage audit
128 names (16 sequences), and every stage is checked using a scalar lifetime
position. Operational JSON is version 2; ordinary profiler JSON stays version 1.
Payload/string capacities are reported; cadence warm capacity is included in
the gravity phase model. Prospective per-event admission remains incomplete.
The single-drift adaptive DMO history path keeps valid old row flags as an
accuracy scale while invalidating the old kick cache.

Files: core profiling/config, gravity runtime, runtime capabilities, analysis,
reference/time workflow and public report types. Reference numerical behavior
is preserved. Risks: existing report/config-hash fixtures require migration and
qualification; count retention differs intentionally from full history.
Follow-up: long-run bounded retention, exact totals, invalid-history cases and
finite-budget accounting; no such run occurred here.

## C — adaptive TreePM acceptance

Before: selected MAC plus all guards and hard strict 0.08 quadrupole envelope;
internal monopoles descended. Now: `kStrictReference` preserves that path.
Opt-in adaptive mode separates an accuracy proxy from an independent maximum
geometric angle and safety guards.

CHUI residual is `F_i=GM d_i f(r)`,
`f=(r^2+eps_pair^2)^(-3/2)+(S-1)/r^3`,
`S=erfc(q)+2q exp(-q^2)/sqrt(pi)`, `q=r/(2 r_s)`.
The implemented inequality is

`GM max(l^2/r^4, rho^2 [3|f'|+|r f''-f'|]/2) <= alpha |A_previous|`,

with `l=2h`, `rho=sqrt(3)h+|COM-center|`. Derivatives use this screened,
Plummer-softened residual. Complete second moment and trace are retained.
The proxy charges second-order scale even with quadrupoles; it is **not a
certified bound on the quadrupole remainder**. No foreign-code alpha is assumed
correct. Finite compatible previous/reference scale-free total A > configured
floor is used. Missing, incompatible, zero/tiny/nonfinite history selects
controlled COM-distance geometry without a fake acceleration. Maximum angle
(default 0.25, disabled mode, validated <=0.5) is a safeguard, not calibrated
accuracy truth. Self/inside, source softening compatibility/heterogeneity,
complete-node cutoff, finite and distributed checks remain.

Files: new internal acceptance header, TreePM/public policy diagnostics, config,
gravity runtime. No new N-sized permanent history is introduced in reference
mode; existing total-force cache supplies history. Risks: proxy calibration,
monopole acceptance, tiny forces and screened cutoff. Follow-up: E–J, rank/thread
comparison and rejection-reason accounting before considering defaults.

## D — P5-W spatial decomposition

Before: measured rank totals could not identify expensive spatial regions.
Now: fixed 64 regions from existing 30-bit decomposition SFC keys, separate from
tree Morton keys; actual owned-target visits + direct pairs and activations
are globally reduced and decayed with factor 0.5. Rate is decayed work /
decayed targets, zero before history. Existing measured-tree-pair coefficient
consumes regional rate. Local pair rank feedback is replaced; incoming service,
communication and rank timing remain, protecting P1 double-count correction.
Hierarchy aggregates completed-block rank work before legal rebalance.

Files: distributed-memory header/source, TreePM coordinator, gravity runtime,
time coordinator, config. Fixed arrays, compact records unchanged; components
remain streamed/recomputed. Reference switch false retains prior weighting.
No permanent particle map/history or O(NR) allocation. History is derived and
cold on restart; subsequent cuts can differ. Work units and coefficients need
calibration, not blind signal addition. Follow-up: K/N/P, active-frequency
balance, tied keys, seam motion, incoming service and memory legal phases.

## E — shared-memory completion

Before: target traversal already threaded; independent source/key, multipole,
PM local sweeps remained serial. Now: static OpenMP source/key work preserves
radix order; root multipole subtrees run independently with child completion
and canonical parent combination. PM spectral operator/gradient sweeps, mean
subtraction, normalization and independent local interpolation are threaded.
No MPI worker calls, deposition atomics, private meshes or duplicate trees.

Files: tree_ordering.cpp, tree_gravity.cpp, tree_pm_coupling.cpp, pm_solver.cpp,
profiling/memory headers. Existing serial paths remain when OpenMP is off or
small-loop conditions apply. Deposition and FFT execution are unchanged:
many-to-one deposition needs a bounded deterministic tile strategy and the
build has no clean FFTW-thread linkage to plumb without dependency changes.
Risks: OpenMP-off guards, nested teams, empty slabs and exact arithmetic promises.
Follow-up: A–D/H/P with source/binary provenance.

## F — conservative topology reuse/refit

Before: every solve rebuilt. Now: two independent false-default policies retain
full rebuild. Identical reuse requires valid nonsaturated source generation,
count, stable nonzero dense-row token, ownership epoch, frame/box, build options
and uniform softening. Motion reuse additionally certifies every assigned row
strictly inside its original octant leaf cell, recursively reconstructed from
original root four doubles and immutable child slots. This derives a concrete
membership invariant from current topology; it is not a displacement heuristic.

Periodic refit uses original unwrap anchors in existing coordinate staging.
Unknown/reordered rows, leaf escape/boundary ambiguity, regressed generation,
heterogeneous softening or option/frame/ownership changes rebuild. A rejected
membership check leaves the tree untouched. Successful refit invalidates the
tree while reverse preorder refits enclosing cubes with outward rounding,
resets COM/moments and recomputes multipoles. Only completion publishes a new
TreeBuildGeneration. Logical Morton keys are cleared; retained permutation is
legal membership, not a fresh Morton ordering. Retained capacities do not gain
an artificial ratchet; no second tree or new population-scale scratch.

Files: tree_gravity header/source, TreePM header/source, config, gravity runtime.
Memory: O(1) certificate/frame identities; existing 3N coordinate and node lanes
reused. Reasons 0–8 and rejected-validity/refit/rebuild timings are observable.
Risks: bounding-cube conservatism, finite rounding, periodic seams, restart
roundoff and source-generation caller authority. Full-build topology can differ,
so refit equivalence needs numerical qualification; no bitwise promise across
restart. Follow-up: L/E/N/P, explicit leaf escapes, mass changes, zero mass,
row reorder, migration, frame/softening changes and empty ranks.

## G — hierarchical power-of-two KDK

Before: scheduler machinery existed but production required global rung zero;
all target forces were refreshed each step. Now: optional `1..12` collisionless
DMO blocks extend the same scheduler/orchestrator and typed stage path.

From existing `dx/dt=u/a`, `du/dt+Hu=A/a^2`, use `p=au`:
`dx/dt=p/a^2`, `dp/dt=A/a`. Drift holds p fixed,
`x1=x0+a0*u0*integral(dt/a^2)`, `u1=u0*a0/a1`.
Endpoint half-kicks add `A_endpoint*integral(dt/a)/a_endpoint`.
Existing background/time conventions and 64-sample quadrature are retained;
the optional canonical operator differs from reference factors and must be
convergence-qualified.

At full sync, mandatory fresh total/split force establishes history and criteria.
Coarse dt <= `2^M * min(gravity_dt)`, further limited by cosmology, explicit dt,
ordered output and endpoint. Quantum = coarse /2^M. Each row receives the largest
power-of-two interval <= its criterion. Bins/quantum freeze within block. Only
scheduler subsets receive Tree force/kicks; all sources drift every fine tick,
updating common scalar drift epoch/source generation once. Tree row generations
must exactly match current source state for kicks. MAC history is separate;
within block it may be fine Tree + coarse PM and cannot authorize a stale kick.

PM kicks occur only at two coarse endpoints. Fine solves are explicitly
short-range-only: zero PM solve/reuse counts, no stale interpolation. Tree
components are accumulated directly, avoiding total-PM cancellation. PM source
and committed version authorize nonzero PM kicks. Every block bootstraps in
both uninterrupted and restarted execution; it costs two PM solves/block but
makes split caches reconstructible. All eight stages are audited; zero-gas
hydro/source work and fine analysis/output are skipped. Local preparation is
collectively guarded through existing FailureCoordinator.

Final all-active closing kick completes before scheduler closure, integrator
coarse commit, cosmological checks, bounded migration, analysis and output.
`step_index` counts coarse blocks; ticks count fine drifts. Mirror bin state
may be dense, but homogeneous species/owner/drift metadata remain compact.
Rung-zero default and its optimized representation remain untouched.

Files: core scheduler/integration/state, time/gravity/runtime composition,
capabilities/config, workflow initialization, restart writer/resume validation.
New memory: 56 bytes/local row split/gen cache + optional byte mirror, existing
~33 bytes/row scheduler on 64-bit ABI, bounded 65544-byte table. Initial scheduler
and later capacities, gravity phase and table use existing governor; actual
capacities reconcile. Generic bounded migration carries current scheduler state.
Restart v23 unchanged: full-sync bins/tick/next activation, common epoch, PM
cadence/source identity and total force history suffice. Derived split caches
rebuild. Fine snapshots/checkpoints are forbidden.

Risks/follow-up: M/N/O/P plus empty/uneven ranks, bin changes, expansion/near-zero
forces, PM split convergence, failure paths and output events. Exact continuation
and scientific equivalence are **not validated**. Lazy drift is deferred (O(N)
source drift remains). Gas/source multirate evolution needs a coupled operator,
not guessed clocks; it is explicitly rejected. The smallest prerequisite is a
specified hydro/source synchronization and force-coupling law plus payload
ownership for those intervals, outside this collisionless TreePM implementation.

## H/I — distributed control and memory completion

The graph must discover whole-active-set peer adjacency before bounded request
layout. Removing the second owner query with cached variable peer lists would
restore population-scale peer storage; batching graph mutation needs protocol
qualification. Existing double query, reusable scratch, bounded wire arena,
compile-time layouts and collective validation remain; no LET rewrite.

All optional large owners above enter existing process/phase admission and
memory reports. Bounded cadence capacity has a gravity phase model and actual
string accounting. Event payloads have bounded count and capacity reporting,
but prospective per-event admission is still missing. Compact planner transient
fit is not a resulting-rank process-fit proof: full future-rank live-set
prediction in cuts remains incomplete; runtime admission can reject such ranks.
No new max-memory option. Follow-up P/Q/R must exercise coexistence with planner,
migration, PM, tree/refit, workers, arena and output/restart.

Optional output cleanup: existing step/code-time schedules preserved. In
hierarchy they clip/coincide with coarse sync. No explicit cosmological target
list, HDF5 redesign or publication-topology change was added (lower priority).

## Configuration, persistence and interface migration

New false-default keys: `numerics.treepm_adaptive_acceptance_enabled`,
`numerics.treepm_identical_source_tree_reuse_enabled`,
`numerics.treepm_topology_refit_enabled`,
`parallel.decomposition_spatial_work_enabled`.
New angle `numerics.treepm_adaptive_maximum_opening_angle=0.25`, finite `(0,0.5]`.
Existing `numerics.hierarchical_max_rung=0` now accepts DMO-only opt-in 1..12.
All are typed/validated/normalized; new keys change config hashes. Strict hash
checks are retained, so old workflow checkpoints need matching provenance or an
explicit future migration, not silent acceptance.

No science snapshot or restart schema change (restart v23). Operational JSON
version 2 has lifetime totals and bounded recent detail. Public tree refit,
TreePM policy/counters/component outputs, scheduler block closure/capacity and
orchestrator failure-seam APIs have same-patch module migration notes. Ordinary
callers leave optional components/directives empty and retain reference behavior.

## Static evidence and handoff

Executed inspection includes `git rev-parse HEAD`, `git status --short`,
`git diff --stat`, targeted `git diff -- <paths>`, `git diff --check`, `rg`,
`rg --files`, `sed`/`cat` documentation/source reads, and Python standard-library
file editing/archive inventory/hash operations. `git diff --check` passed during
source review; final hygiene/ZIP inventory is recorded with delivery.
Final tracked `git diff --check` returned 0 with no diagnostics. New-file
`git diff --no-index --check /dev/null <path>` inspection produced no whitespace
diagnostics (exit 1 denotes the expected new-file difference in no-index mode).
No repository test/hygiene script, build/configuration or scientific command ran.
Static inspection is not C++ compilation or validation.

Suggested branch: `perf/p6-preproduction-mega-closure`.
Suggested commit: `perf: implement pre-production gravity and scaling architecture`.
Suggested PR: `P6: implement pre-production gravity, scaling, and hierarchical runtime architecture`.
No branch/commit/PR was created; Git metadata is read-only in this workspace.

Delivery: changed/new files only, original repository-relative paths, including
this report and `p6_validation_campaign_plan_20261006.md`. Build trees, binaries,
runtime outputs/caches and unrelated local config/history are excluded.

## Exact changed/new file inventory

```text
CURRENT_STATUS.md
README.md
docs/architecture/decision_log.md
docs/architecture/runtime_truth_map.md
docs/build_instructions.md
docs/configuration.md
docs/memory_governance.md
docs/output_schema.md
docs/parallel_distributed_memory_contracts.md
docs/profiling.md
docs/repair/p6_preproduction_performance_architecture_implementation_20261006.md
docs/repair/p6_validation_campaign_plan_20261006.md
docs/repair_open_issues.md
docs/repair_state_recap.md
docs/restart_checkpointing.md
docs/state_model_memory_layout.md
docs/time_integration.md
docs/tree_gravity_solver.md
docs/tree_pm_coupling.md
docs/validation_plan.md
include/cosmosim/core/config.hpp
include/cosmosim/core/profiling.hpp
include/cosmosim/core/time_integration.hpp
include/cosmosim/core/time_scheduler.hpp
include/cosmosim/gravity/gravity_memory.hpp
include/cosmosim/gravity/pm_solver.hpp
include/cosmosim/gravity/tree_gravity.hpp
include/cosmosim/gravity/tree_pm_coupling.hpp
include/cosmosim/parallel/distributed_memory.hpp
include/cosmosim/workflows/reference_workflow.hpp
include/cosmosim/workflows/time_coordinator.hpp
src/core/config.cpp
src/core/profiling.cpp
src/core/simulation_state_ownership.cpp
src/core/simulation_state_species.cpp
src/core/time_integration.cpp
src/gravity/gravity_memory.cpp
src/gravity/internal/tree_pm_acceptance.hpp
src/gravity/pm_solver.cpp
src/gravity/tree_gravity.cpp
src/gravity/tree_ordering.cpp
src/gravity/tree_pm_coupling.cpp
src/io/restart_checkpoint.cpp
src/parallel/distributed_memory.cpp
src/workflows/analysis_runtime.cpp
src/workflows/gravity_runtime.cpp
src/workflows/output_restart_runtime.cpp
src/workflows/reference_runtime_composition.cpp
src/workflows/reference_workflow.cpp
src/workflows/runtime_capabilities.cpp
src/workflows/time_coordinator.cpp
```

Inventory: 51 modified/new repository files; no build/test artifacts.
