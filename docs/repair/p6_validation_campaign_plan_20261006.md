# P6 second campaign: implementation qualification plan

**NO PART OF THIS MATRIX WAS EXECUTED IN THE P6 SOURCE CAMPAIGN.**

**NO TEST IMPLEMENTATION WAS ADDED BY THAT CAMPAIGN.**

**SOURCE-IMPLEMENTED / VALIDATION PENDING.**

This is the required handoff, not a qualification report. Input source revision
was `a88de7ec1098b89c8d8e6679431a8466067ebc62`; qualify the exact subsequent patch/
commit identified by the changed-files artifact and manifest. Do not replace
existing validated P1–P4/M48 reference paths or make optional modes defaults
before numerical and memory gates are accepted.

## Preparation and acceptance authority

Apply the changed-files ZIP onto its recorded revision without unrelated local
files. Freeze the resulting source commit/diff hash, compiler, dependency/provider
families, CPU/rank/thread affinity, normalized config and input-state hashes.
First repair real compile/fixture failures with auditable source changes; do not
hide them by weakening assertions or errors. The new normalized keys change
config hashes and operational JSON advances to v2; update fixtures/consumers
explicitly. Old restart config hashes must not be bypassed.

Implement the smallest missing focused qualification cases in this second
campaign using existing test/validation architecture. Existing starting points
in CMakeLists include `unit_tree_gravity`, `unit_time_integration`,
`unit_tree_pm_split_kernel`, `integration_tree_gravity_vs_direct`,
`integration_tree_pm_coupling_periodic` and MPI variants,
`integration_reference_workflow_distributed_treepm_mpi_two_rank`,
`integration_hierarchical_time_bins`, `integration_hierarchical_timestep_regression`,
restart/snapshot suites, `validation_convergence`,
`validation_tree_pm_ewald_accuracy`, `validation_dmo_zeldovich_workflow_*`,
and `validation_power_spectrum_mesh_consistency`. Confirm registrations on the
qualification checkout; existing scheduler tests alone do not certify the new
workflow block operator.

Freeze numerical acceptance thresholds before performance tuning. Use current
reference tolerances (`validation/reference/validation_tolerances_v1.txt`) where
applicable, and document any new force-percentile, PM cadence, evolution/growth
or restart tolerance with scientific justification. No alpha/angle is promoted
because it runs faster. Near-zero absolute errors need separate thresholds;
relative-error division by a fabricated acceleration floor is not evidence.

## A — compilation/build matrix

Run the existing configure/build/test triples on a dependency-capable host:

| Configure preset | Build preset | Test preset |
| --- | --- | --- |
| `cpu-only-debug` | `build-cpu-debug` | `test-cpu-debug` |
| `hdf5-debug` | `build-hdf5-debug` | `test-hdf5-debug` |
| `pm-hdf5-fftw-debug` | `build-pm-hdf5-fftw-debug` | `test-pm-hdf5-fftw-debug` |
| `mpi-serial-hdf5-fftw-debug` | `build-mpi-serial-hdf5-fftw-debug` | `test-mpi-serial-hdf5-fftw-debug` |
| `mpi-hdf5-fftw-debug` | `build-mpi-hdf5-fftw-debug` | `test-mpi-hdf5-fftw-debug` |
| `mpi-hdf5-fftw-release` | `build-mpi-hdf5-fftw-release` | `test-mpi-hdf5-fftw-release` |

For each row run `cmake --preset <configure>`, `cmake --build --preset <build>`,
`ctest --preset <test> --output-on-failure`. Include `cpu-only-release` with
`build-cpu-release` for performance comparisons and `asan-debug` with
`cmake --build build/asan-debug` / `ctest --test-dir build/asan-debug` for the
small focused source/memory cases. These commands are **future work**.
Confirm C++20, MPI-only/no-FFTW diagnostic policy, OpenMP-off compilation, FFTW
slab/transposed layouts and provider-family compatibility. MPI science-output
acceptance requires real Parallel HDF5; serial HDF5 is a failure-policy gate,
not an alternative production output topology. No new FFTW-thread dependency.

## B — CPU Debug floor

Run `./scripts/ci/check_repo_hygiene.sh`, CPU Debug configure/build/full CTest
and applicable focused floors. Require preserved default config behavior,
normalized dump reparsing, invalid angle/rung/source-physics rejection,
capability eligibility/maturity, canonical audit, ring ordering/eviction and
lifetime counters. Check empty populations, saturated identities, incompatible
force histories and source fingerprints. Attribute failures to baseline versus
patch with exact source/binary evidence; no blanket pre-existing classification.

## C — OpenMP OFF / 1 / 2 / 4

Configure an isolated CPU/PM checkout with `-DCOSMOSIM_ENABLE_OPENMP=OFF`; do
not overwrite the main preset build. Run default and opt-in force cases with
OpenMP OFF, then enabled `OMP_NUM_THREADS=1,2,4`, fixed binding/places. Compare
source generation, topology/order, moments, vectors, total/worker counts and
memory. Require exact integer counters and pair identity; freeze appropriate
floating equality/tolerance. Exercise >=1024 key/interpolation loops, >=256-node
multipole root branches, single-leaf/deep/skewed trees, no work workers and
empty slabs. Verify no MPI calls from workers and no worker-hot false sharing.

## D — MPI 1 / 2 / 4 fixed-core cases

Use four fixed CPU cores: rank/thread pairs (1,4), (2,2), (4,1), plus (1,1),
(2,1), (4,1) and thread variations when needed to separate rank/thread effects.
Record launcher binding and `MPI_THREAD_FUNNELED` level. Same stable global IDs,
same physical saved source/target state, deterministic decomposition order.
Include empty and strongly uneven ranks, tied SFC keys, periodic seam domains,
collective stale-geometry fallback/refit recovery and failed admission on one
rank. Check sparse peers, graph cache hits, packets/identity/version validation
and no deadlock. Per-rank combined pairs must equal local + incoming.

## E — strict versus adaptive force vectors

Use frozen positions/masses/softening and identical PM mesh/operator. Compare
strict reference to adaptive angle/tolerance settings, monopole/quadrupole
separately, with direct/Ewald oracle on feasible subsets. Cover valid prior A,
missing/zero/tiny/nonfinite A, changed row/ownership/source generation, initial
bootstrap, compatible single drift, hierarchy fine histories, heterogeneous
softening (guarded descent), self/near geometry and cutoff-overlapping nodes.
Verify the exact proxy/derivatives against the implemented residual, including
trace terms. Inspect rejected guard families and accepted internal nodes;
leaves must not be counted as multipoles. A proxy is not an a priori error bound.

## F — representative saved cosmological states

Freeze at least linear/high-redshift, moderately clustered and low-redshift
nonlinear states, with both compact halos and diffuse/void regions, periodic
seam clusters, cancellation/near-zero acceleration and unequal softening where
supported. Record snapshot/input hashes, cosmology, force normalization,
softening/split/cutoff, meshes, source generation and target sampling policy.
Include current supplied-run states if available; absence must be reported,
not replaced by a claim that a uniform lattice represents clustered evolution.

## G — error distributions

For each E/F setting record per-component and vector errors, RMS, maximum,
p50/p90/p95/p99/p99.9, and absolute errors for near-zero reference force.
Distinguish PM, Tree residual and total vectors; also report accepted-node
subsets and difficult geometry/softening/cutoff categories. Publish full sample
counts and bounded reproducible samples/histograms. Rank/thread comparison
must use the same physical target IDs. Reject outlier-tail regressions even if
median improves. Freeze tolerances before selecting recommended alpha/angle.

## H — pair/node/multipole work and causal timing

For each identical force case record preparation/hash/order/topology/moments,
node count/depth, local/incoming targets, visits/opens, internal multipoles,
leaves, direct pairs, cutoff and all rejection counters. Verify aggregate
identity and immutable per-target accumulation order. Record worker-region
min/max/mean and summed work versus wall time, PM inclusive/nested timing,
routing/FFT/axis/halo/interpolation, peers/export multiplicity, arena high-water
and graph/geometry reuse. Run enough matched warm repetitions for variability;
do not add nested timers. Verify ring eviction retains lifetime totals.

## I — strict/adaptive short 64^3 evolution

Same IC phases, cosmology, softening, PM mesh, global KDK/output schedule and
fixed resources. Strict default versus individually selected adaptive settings,
short matched redshift/time interval, several timestep resolutions. Track
trajectory/velocity differences, forces, work, runtime, conservation/health,
step/source generation and publication epochs. Add refit/work/hierarchy only in
separate cases after their gates pass; isolate causes rather than toggling all.

## J — power spectrum / growth

For I and later hierarchy cases compare matched growth histories, Zel'dovich
linear behavior, power spectrum ratios/residuals over documented trustworthy k,
and low/high-resolution consistency. Same analysis windows/deconvolution and
particle selection. Separate sampling/mesh noise, timestep, MAC and PM-cadence
errors. Use current scalable FFT diagnostics; small direct DFT is oracle only.
Qualification thresholds and uncertainty must be explicitly accepted.

## K — work-aware decomposition balance

Reference rank totals versus 64-bin spatial mode on F/I at 1/2/4 ranks and
hierarchical active fractions. Inspect regional integer work/EMA/targets/rates,
work coefficient units, incoming service and communication contribution,
per-rank residual wall/pairs/visits, cut stability, migrated volume and memory.
Require compact records/dense local indices and no duplicate local pair feedback.
Test zero history, tied keys, shifted hot regions, cold restart, disabled mode
and legal coarse-only rebalance. Quantify balance versus migration/PM cost.

## L — topology reuse/refit equivalence

Identical source reuse against full rebuild; then controlled motions strictly
inside original cells, boundary ties, deliberate escape, periodic anchor seam
crossing, source count/row permutation, ownership, box/frame, options/softening,
mass/zero-mass changes and saturated/regressed/unknown generations. Verify
rejection reasons, original membership proof, enclosing refitted bounds,
child completion, multipoles, generation publication, no partial valid tree,
cleared logical Morton keys and no capacity ratchet. Compare forces/work to
fresh build with documented tolerance; bitwise equality across different legal
topologies is not assumed. Include adaptive policy and hierarchy only after
isolated refit cases. Confirm zero extra N-node/source certificate arrays.

## M — hierarchical/global KDK convergence

Test collisionless free expansion (p constant), two-body/simple gravity in the
core operator, and periodic TreePM DMO cosmology in workflow. Use M=1/2/4 and
maximum admitted rung 12 on a tiny case; equal bins and mixed occupancies,
empty bins/ranks and rung reassignment at successive full sync. Compare matched
physical endpoints to refined global reference, halving quantum/coarse PM
interval separately. Check 64-sample cosmological integral error, opening and
closing interval factors, common source epochs, exact Tree cache generation,
PM source/version, force-only bootstrap, fine zero PM solve/reuse counts,
two PM solves/coarse block and source/tick overflow rejection. Verify stale
kick and unsynchronized migration fail closed. All-source fine drift overhead
is part of measured cost; do not claim lazy prediction. Gas/source modes reject.

## N — restart equivalence

Default rung-zero unchanged, then hierarchical strict/full-build split at
multiple closed coarse boundaries before/after rung changes and migration.
Compare uninterrupted versus restart at identical physical endpoints: state
vectors keyed by ID, integer ticks/bins/next activation, common drift time/a,
PM opportunity/version, source identity, future output names/events and cache
invalidations. Derived split forces bootstrap in both paths. Require exact
continuation where the accepted contract promises it; diagnose any difference
before relaxing it. Separately qualify refit and cold spatial-history restart
with explicit tolerance/ownership semantics. Reuse existing v23 checksum,
streaming verification, roundtrip, legacy reader and corruption tests. No fine
checkpoint. Reject malformed scheduler extents and incompatible config hashes.

## O — coherent snapshots

At 1/2/4 ranks using Parallel HDF5, snapshot/restart only after closing all coarse
kicks and legal migration. All particles share header time/a and common drift
epoch, unique global IDs/count/mass and periodic coordinates. Exercise time-event
clipping, step cadence, simultaneous checkpoint/snapshot, final endpoint
publication/de-duplication and overwrite refusal. Confirm schema unchanged and
one shared transactional science file; rejected fine output must not create a
partial committed artifact.

## P — finite process-memory budgets

Use existing process/gravity budgets and reserves, zero-gas reference plus
hierarchy. Cover success just above measured peaks and collective rejection
below required admission, on uneven ranks. Inventory capacities and simultaneous
live sets: canonical state, generic scheduler (~33N ABI dependent) and bin
mirror, new 56N split/gen caches, total history, compact planner/migration,
periodic/tree/refit lanes, workers, PM spectra/halo, arena, 65544-byte timeline,
event/cadence/stage retention, output/restart staging and opaque reserves.
Verify scratch commitments are not double-counted, capacity remains governed,
no geometric ratchet, no full per-worker mesh/tree and zero-gas cold owners
stay absent. Long event runs must stay bounded and lifetime totals exact.
Prospective event-payload admission and predictive resulting-rank fit are
known incomplete areas; distinguish runtime rejection from startup fit proof.

## Q — 128^3 scaling

After A–P acceptable: real matched 128^3 strict reference, then individual
adaptive/spatial/refit/hierarchy modes, 1/2/4 ranks and fixed-core configurations.
Report force/evolution error alongside causal phase work, migration and peak
process memory. Select production candidate settings only from qualified error
and measured work; weak scaling comparisons must keep physics/mesh/cadence fixed.

## R — real 256^3 smoke

Actual 16,777,216 particles, not extrapolation or a smaller renamed case. Use a
host/rank layout with finite admitted process budgets and record complete owner
model, provider family, source/binary identity, IC hash and resource affinity.
Run a small real synchronized interval, save/read coherent output/restart,
inspect peak RSS/arena/scheduler/Tree/PM and force health. Test each accepted
optional mode in a defensible order. A startup-only dry plan is not a smoke run.
If resource admission fails, report exact phase/rank/budget; do not weaken it.

## S — production output cadence

Use a deliberate sparse step/code-time cadence appropriate to the scientific
interval; avoid hundreds of snapshots from a short step modulus. Compare I/Q/R
with and without output to separate HDF5 cost. Hierarchical cadence counts
coarse blocks, and time events align via clipping. Explicit redshift/scale-factor
lists were not implemented; do not imply they exist. Verify final endpoint,
future restart cadence, shared paths and unchanged HDF5 schema/publication.

## T — exact provenance and final decision

Every table/plot/run records input revision, resulting qualified source commit
and patch hash, binary checksum/build metadata, normalized config/hash,
compiler/options, dependency versions/provider family, MPI threading/binding,
OpenMP threads, IC/state hashes, mesh/softening/split/cutoff, history policy,
refit/decomposition/hierarchy switches, cadence and process budgets. Record
actual commands and return codes plus environment limitations. Keep default
reference evidence separate from opt-in evidence and historical campaigns.

Deliver a qualification report mapping A–T to artifacts/pass/fail/blocked and
named residual defects, with force/error tails, convergence, restart/output and
finite-memory evidence. Defaults require an explicit accepted qualification
change; implementation availability and faster work counts do not certify
science or production readiness.
