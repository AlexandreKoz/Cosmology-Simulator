# Unified TreePM performance recovery and hierarchical qualification — 20261010

Mode: feature implementation with targeted correctness repairs, using the current
checkout. The starting working tree was clean. The older supplied archive was
not used. The external P0 profile is accepted as supplied evidence; none of its
measurements was repeated here.

Compilation: **NOT RUN — deferred by user**  
Tests: **NOT RUN — deferred by user**  
Performance qualification: **NOT RUN — deferred by user**

Static checks performed: `git diff --check` returned exit code 0. `ast.parse`
accepted both new Python tool files without importing or executing them. The
38-entry changed-file inventory was checked against the Git working tree.

No configure, compiler, numerical test, simulation, benchmark, sweep, profiler,
dependency installation, or CI hygiene script was executed. Python AST parsing
and `git diff --check` are static checks only. Written assertions and budgets
below are acceptance targets, not passing evidence. No 2x speedup, memory fit,
validated nonlinear scale range, or capability promotion is claimed.

## Implementation and evidence

### P1 — bounded deterministic child traversal

`src/gravity/tree_pm_coupling.cpp` previously computed periodic child-center
distances and sorted up to eight children at each descent. P0 attributed 6.9%
of sampled self cycles to that helper. Production now inserts existing children
in reverse octant order into the same bounded DFS stack, visiting octants 0–7.
There are no sorting-only distances or heap allocations. Stack capacity checks,
child identity, MAC and cutoff contracts remain. Incoming MPI target evaluation
uses the same evaluator.

This changes force accumulation order. Target traversal order is independent of
OpenMP scheduling; legacy force checksums can change through roundoff. The new
direct-reference fixture checks actual pair counts and force agreement, then
compares 1/4-worker results. Locality might regress on some trees; measure the
complete traversal before accepting a performance conclusion.

### P2 — stable analytic residual and gated coefficient lookup

The pair law remains, for positive separation,

```text
K_SR = (r^2 + epsilon_pair^2)^(-3/2) - L(r,a_split)/r^3
L = erf(q) - 2 q exp(-q^2)/sqrt(pi), q = r/(2 a_split)
epsilon_pair = max(epsilon_source, epsilon_target)
```

The Gaussian factor limits are now S(0)=1 and L(0)=0. For q<1/8, the analytic
path evaluates the n=0..7 series

```text
H(t) = L(q)/q^3 = 4/sqrt(pi) sum_n (-t)^n / [n! (2n+3)], t=q^2
L/r^3 = H(t)/(8 a_split^3)
lim_(r->0) L/r^3 = 1/(6 sqrt(pi) a_split^3).
```

Sequential scale divisions avoid premature a_split^3 overflow/underflow. This
avoids subtracting S from one and avoids dividing an underflowed L by r^3.
The omitted dimensionless series term is below 2e-20 on that interval. At
exact r^2=0 the existing finite-softening coefficient convention is retained;
the displacement vector is zero. Checked public pair/factor APIs still reject
nonfinite/negative separations and invalid softening/split values. Extreme
inputs whose mathematical coefficient exceeds floating-point range are not a
promise of finite representability.

`TreePmGaussianCoefficientTable` stores H and its analytic derivative at 4097
knots, with cubic Hermite interpolation on t in [0,64]. Near zero it uses the
same stable series; outside the interval or unrepresentable fast normalization
it uses the analytic coefficient. Table creation is outside every worker/pair
loop. Ordinary enabled leaf pairs evaluate the interpolation directly, without
erf/erfc/exp/sqrt. The independent Plummer term remains per pair, so variable
softening needs no additional table. Absolute H error <=5e-11 is a **written
qualification target**, not an established bound on relative residual force:
softened and long-range terms can nearly cancel.

Accepted-node monopoles/quadrupoles retain the established analytic radial
derivatives and full second moment including its trace. The lookup is never
used for multipoles. The MAC proxy uses stable -L in place of S-1; its derivative
formula and second-order interpretation are retained. No spectral PM gradient,
assignment/deconvolution, cosmological normalization or split policy is changed.

Two typed `.param.txt` booleans, with canonical validation, defaults, normalized
dump and tests, are added:

| Key under numerics | Default | Meaning |
|---|---|---|
| treepm_gaussian_pair_lookup_enabled | false | Experimental direct leaf-pair interpolation |
| treepm_full_mac_diagnostics | true | Preserve overlapping forensic MAC counters |

The normalized configuration/provenance hash includes both keys. Old input
files receive these defaults; normalization changes their hash. Existing strict
restart compatibility checks are retained. There is no restart/schema migration
or bypass in this patch. Existing adaptive alpha=0.005 and maximum angle=0.25
semantics remain; no defaults for mesh, leaf size, rungs or OpenMP are changed.

### P3 — reused periodic geometry and staged MAC

The evaluator computes the general per-axis nearest-image node-center deltas
once. It reuses them for AABB pruning and target containment. Squared AABB
distances avoid ordinary cutoff square roots; an 8-epsilon neighborhood of the
boundary and extreme/nonfinite squared values use the old sqrt comparison to
retain rounding-sensitive decisions. Leaves skip COM, offset and MAC geometry.
Internal nodes compute the maximum envelope only when needed. The established
general nearbyint convention remains for rectangular boxes, multiple wraps,
noncanonical coordinates and half-box ties. COM-center offset remains the
established unwrapped offset. No approximate periodic pruning is introduced.

Cheap independent guards are evaluated before expensive MAC derivatives.
Default full accounting retains overlapping rejection counters. Opt-in fast
accounting counts all cheap guards, then records `skipped_mac_evaluations` when
those guards already force descent; skipped selected/relative MAC checks are
not counted as passes or failures. Profile/event metadata explicitly identifies
the accounting mode. Unique opened-node and actual pair counters remain actual
work. In either mode, the first operand of the adaptive proxy's max can reject
before derivative evaluation, using the same arithmetic and rejection counter.
Invalid/unavailable acceleration history still uses legitimate geometric
fallback; no artificial history floor is introduced.

### P4 — frozen-source split-aware experiments

The historical `bench/bench_tree_pm_snapshot.cpp`,
`tools/chui_treepm_oracle.py` and `tools/chui_treepm_profile.py` are absent from
this checkout. Their historical implementation or CLI was not invented.
`bench_tree_pm_sweep` instead uses existing production gravity interfaces and
the existing independent periodic Ewald support. It accepts immutable sorted
particle CSV and fixed target IDs, mesh, leaf size, bounded residual block size,
policy, kernel, counter mode, repeated solves and optional explicit same-snapshot
strict TOTAL bootstrap history. That history is not original simulation KDK
history. Every timed solve evaluates all sources as targets and refreshes PM.
Accuracy/reference generation is outside measured solve time.

JSONL records include PM phase totals, tree build/multipoles/traversal, pair/node/
opening/multipole counts, observed workers, solver physical retained capacity,
process RSS and peak RSS where available. Process peak includes benchmark-owned
sources, outputs and optional accuracy/history arrays; it is not a whole-runtime
production memory certification. The optional small-N Ewald comparison requires
image convergence and uniform softening. Unsupported naive-DFT large meshes are
rejected. This driver deliberately uses one MPI owner even in MPI builds.

`tools/treepm_sweep.py` exports one HDF5 DMO group in explicit solver units and
runs mesh 64/96/128/192, leaf 4/8/16/32 and optional block/thread matrices. It saves
exact commands, input/ID/binary/reference SHA-256 values, per-candidate JSONL and
stderr, and timing/error CSV/Markdown reports. Production still has one
authoritative `.param.txt` configuration system. `TreePmOptions::residual_block_size`
is an API/benchmark option validated in [1,4096]; default 64 and dynamic,1 OpenMP
scheduling remain. Nondefault callers must budget retained per-block scratch.

Each candidate derives physical split/cutoff from its own mesh spacing. Direct
long-double short references are regenerated at that split with frozen masses,
positions, softening and target IDs. Total-force comparison accepts an externally
qualified converged CSV with exactly matching IDs. Wall time alone cannot
qualify a split. The older force-error-map benchmark reference was also repaired:
it multiplied softened inverse-r^3 by Gaussian screening, which did not match
the production softened-minus-unsoftened-long-range law.

### P5 — existing hierarchical KDK safeguards and qualification

The existing dispatcher and P6 ownership implementation are extended, not
replaced. Preflight now verifies the exact all-source identity row set, positive
quantum, synchronized common source drift epoch, and supplies scheduler pointers
at the core dispatch seam. The inside-KDK marker opens only after all-active
preparation succeeds. Existing all-source drift, generation increments, PM
cache/version checks, fixed-within-block rung assignment and synchronized
endpoint PM half-kicks remain.

Source inspection found that the previous coarse interval scaled the minimum
softening/gravity criterion by the largest rung period without an independent
velocity/PM-resolution crossing restriction. The new additional restriction is
derived at the synchronized epoch:

```text
ell = min(dx,dy,dz,a_split)
v_com = |u_peculiar|/a_cosmology
g_com = |A_scale_free|/a_cosmology^3
v_com dt + g_com dt^2/2 <= ell
dt = 2 ell / [v_com + hypot(v_com,sqrt(2 g_com ell))]
```

Global maxima of speed and acceleration are obtained through two unconditional
collectives, including empty ranks. The root is evaluated in long double. This
is a local frozen-force resolution estimate; later nonlinear force growth is
not bounded by it. Existing softening/gravity, expansion, requested timestep,
scheduled output and endpoint restrictions remain and can only tighten the
accepted interval. No empirical fraction is added. The cosmological law remains
p=a_cosmology u, dp/dt=A/a_cosmology, dx/dt=p/a_cosmology^2. No opening/closing
kick law was speculatively rewritten. Mandatory synchronized bootstrap remains:
equivalent restart/total-history initialization without it has not been proved.

Actual workflow qualification now covers a small cosmological DMO mode plus a
close pair, global fine KDK versus hierarchical coarse/refined runs, deterministic
replay, synchronized restart continuation, source/PM epochs and MPI 1/2/4 ranks.
The existing Zel'dovich validation gains hierarchical fundamental-mode growth,
production power-spectrum refinement and saved 1/2/4-rank state equivalence.
Its scale target is **only k=2 pi/L, the six fundamental Cartesian modes of that
fixture**. There is no already validated hierarchical scale range. The offline
`tools/hierarchical_dmo_qualification.py` adds assertions for matched larger-run
global/coarse/refined snapshots over an explicitly requested k interval, with
identical CIC estimation, a quarter-Nyquist ceiling and at least eight modes per
shell. Those are estimator safeguards, not scientific certification.

## Ownership, memory, MPI and reproducibility

- The sole new production retained allocation is the coordinator-owned immutable
  coefficient table: `sizeof(TreePmGaussianCoefficientTable)` (65552 bytes of
  arrays on ordinary double platforms), independent of N, scale and softening.
  Governed coordinators reserve persistent-cache capacity before allocation and
  commit only after construction. Preparation failure uses existing collective
  consensus. Memory reports expose physical size and its governed commitment.
  Whole-runtime estimates include it; overlapping gravity phase leases exclude
  this separately committed allocation. No population-scale cache is added.
- Dimensionless coefficients need no rebuild when split changes. Disabling the
  option retains the bounded allocation until coordinator destruction, which
  frees storage before releasing its commitment. The governor must outlive the
  coordinator. Readers are const after initialization; no process-global mutable
  table/scratch or cross-instance race is introduced. No table is persisted.
- Existing owner-local and incoming-target evaluation share the same changes.
  No LET wire fields, source ownership, target identity, rank topology or
  collective ordering is changed by traversal/lookup. Coarse restriction
  reductions are unconditional; existing collective failure propagation remains.
- Default analytic force arithmetic changes at small separation and force sums
  change with traversal order. Deterministic replay must use equivalent versions,
  normalized configurations and rank/thread contracts. Old numerical checksums
  should be compared by scientific errors, not relabeled as new bitwise checks.
- No 256^3/512^3 fit claim is made. Allocator overhead, opaque libraries and
  workflow overlap still require existing whole-runtime measurements.

## Tests written, not executed

| Source | Assertions/target |
|---|---|
| tests/unit/test_tree_pm_split_kernel.cpp | Limits, tiny r, softening ratios, extreme split scales, series transition, 131073-point interpolation grids, edges/fallback, invalid checked inputs, squared cutoff and full/fast MAC equivalence |
| tests/unit/test_config_parser.cpp | Typed defaults, normalized roundtrip and invalid booleans |
| tests/integration/test_tree_pm_coupling_periodic.cpp | Rectangular seams, wraps, half-box/cutoff neighbors, exact direct-pair count, force tolerance, lookup decision equivalence, bounded governed table lifetime, 1/4-thread replay |
| tests/integration/test_hierarchical_timestep_regression.cpp | Derived displacement root, malformed criteria, actual core hierarchical dispatcher free cosmological drift/momentum, sparse/all-active transitions, generation and endpoint checks, invalid row identity |
| tests/integration/test_hierarchical_treepm_workflow.cpp | Actual global/hierarchical workflow refinement, replay, endpoint safety, restart continuation, PM cadence, mass, optional 2/4 ranks |
| tests/validation/test_tree_pm_ewald_accuracy.cpp | Existing Ewald budgets applied to strict/adaptive analytic/lookup with declared same-fixture reference history |
| tests/validation/test_dmo_zeldovich_workflow.cpp | Hierarchical fundamental growth/power refinement and 1/2/4-rank physical-state comparison |

Existing hierarchical-bin, restart-equivalence, memory-governor, distributed
TreePM and global Zel'dovich tests remain required regression gates. Unit tests
using `assert` must be run in Debug as well as any Release scientific gates.

## Manual commands — user execution only

Run from repository root. These commands have **not** been executed here.

### Build and baseline

```bash
./scripts/ci/check_repo_hygiene.sh
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure

cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R 'unit_tree_pm_split_kernel|unit_config_parser|unit_memory_governor|integration_tree_pm_coupling_periodic|integration_hierarchical|integration_restart_equivalence_treepm|validation_tree_pm|validation_hierarchical_dmo|validation_dmo_zeldovich'

cmake --preset mpi-hdf5-fftw-release
cmake --build --preset build-mpi-hdf5-fftw-release
ctest --preset test-mpi-hdf5-fftw-release --output-on-failure
```

MPI Release requires the existing Parallel-HDF5/FFTW-MPI provider contract; no
dependency availability was probed in this campaign. Focused named gates:

```bash
ctest --test-dir build/pm-hdf5-fftw-debug --output-on-failure -R '^validation_tree_pm_recovery_ewald_accuracy$'
ctest --test-dir build/pm-hdf5-fftw-debug --output-on-failure -R '^integration_hierarchical_treepm_workflow$|^validation_hierarchical_dmo_growth_single_rank$'
ctest --test-dir build/mpi-hdf5-fftw-release --output-on-failure -R 'integration_tree_pm_coupling_periodic|integration_reference_workflow_distributed_treepm_mpi_two_rank|integration_restart_equivalence_treepm'
ctest --test-dir build/mpi-hdf5-fftw-release --output-on-failure -R 'integration_hierarchical_treepm_workflow|validation_hierarchical_dmo'
```

CTest hierarchical rank-equivalence fixtures automatically arrange the required
np1/np2/np4 artifact producers. The Ewald recovery test uses the existing
first-light softening ratio and inherited scientific limits; tolerances were
not relaxed. Compare Ewald image levels separately for new fixtures.

### Frozen snap_1698 force checks and full-target measurements

Supply actual P0 normalized physical metadata in shell variables first:
`SNAPSHOT`, `POSITION_SCALE`, `MASS_SCALE`, `EPS_COMOVING`, `BOX_X`, `BOX_Y`,
`BOX_Z`, `G_CODE`, `SCALE_FACTOR`, `ASMTH_CELLS`, `RCUT_CELLS`. Positions, boxes
and epsilon must use the same comoving solver length units; masses/G must share
the same code units. The CLI does not infer unit/h factors. The driver reports
scale-free acceleration A, before cosmological kick factors. Choose a fresh
output directory for each experiment.

```bash
python3 tools/treepm_sweep.py export "$SNAPSHOT" \
  --output /tmp/snap_1698_particles.csv --target-ids /tmp/snap_1698_targets_4096.txt \
  --targets 4096 --epsilon "$EPS_COMOVING" \
  --position-scale "$POSITION_SCALE" --mass-scale "$MASS_SCALE"

TREEPM_BENCH="$PWD/build/mpi-hdf5-fftw-release/bench_tree_pm_sweep"
TREEPM_IDS=/tmp/snap_1698_targets_4096.txt
# To reproduce the historical target population, set TREEPM_IDS to its saved
# one-integer-ID-per-line file instead. The generated even-ID set is different.
TREEPM_PHYSICAL=(--box-x "$BOX_X" --box-y "$BOX_Y" --box-z "$BOX_Z" \
  --asmth "$ASMTH_CELLS" --rcut "$RCUT_CELLS" --epsilon "$EPS_COMOVING" \
  --g "$G_CODE" --scale-factor "$SCALE_FACTOR")

python3 tools/treepm_sweep.py run --exe "$TREEPM_BENCH" \
  --particles /tmp/snap_1698_particles.csv --target-ids "$TREEPM_IDS" \
  --output /tmp/treepm_1698_baseline_compare --meshes 64 --leaves 16 --blocks 64 \
  --threads 1,2,4 --policies strict,adaptive --kernels analytic,lookup \
  --history fallback --accounting full --repeats 5 --warmups 1 --oracle \
  "${TREEPM_PHYSICAL[@]}"

python3 tools/treepm_sweep.py run --exe "$TREEPM_BENCH" \
  --particles /tmp/snap_1698_particles.csv --target-ids "$TREEPM_IDS" \
  --output /tmp/treepm_1698_bootstrap_compare --meshes 64 --leaves 16 --blocks 64 \
  --threads 4 --policies strict,adaptive --kernels analytic,lookup \
  --history bootstrap --accounting fast --repeats 5 --warmups 1 --oracle \
  "${TREEPM_PHYSICAL[@]}"
```

With 262144 CSV sources every timed solve is a complete 262144-target Tree+PM
refresh; `--targets` affects accuracy selection only. Verify count, source/ID
hashes, physical split, observed workers and phase/counter records. The second
command labels reconstructed same-snapshot history explicitly. Saved historical
oracle output cannot be reused after a changed split. The absent P0 scripts'
native oracle-file format/commands cannot be reproduced from this checkout;
use its exact saved IDs and this driver's regenerated matching direct reference,
or adapt the known external format outside this patch.

### Mesh/leaf/block tuning with total-force validation

```bash
python3 tools/treepm_sweep.py run --exe "$TREEPM_BENCH" \
  --particles /tmp/snap_1698_particles.csv --target-ids "$TREEPM_IDS" \
  --output /tmp/treepm_1698_split_sweep --meshes 64,96,128,192 --leaves 4,8,16,32 \
  --blocks 32,64,128 --threads 4 --policies strict,adaptive --kernels analytic,lookup \
  --history fallback --accounting fast --repeats 5 --warmups 1 --oracle \
  --total-reference "$CONVERGED_TOTAL_FORCE_CSV" "${TREEPM_PHYSICAL[@]}"
```

`CONVERGED_TOTAL_FORCE_CSV` requires columns id,total_x,total_y,total_z in the
same A units, identical source population and target IDs, masses/softening,
cosmology and box, with independently established periodic convergence. A CSV
hash alone does not prove that identity. Review wall/PM/build/traversal tradeoffs,
short and total L2, P95/P99/max, absolute and small-force tails. No winner is
certified automatically. The matrix is intentionally explicit and may take
substantial user-side time; start with fewer axes if appropriate.

For a small independent periodic check (diagnostic synthetic source population):

```bash
for level in 2 4; do
  OMP_NUM_THREADS=4 "$TREEPM_BENCH" --count 64 --mesh 64 --targets 64 \
    --policy adaptive --history bootstrap --kernel lookup --accounting fast \
    --repeats 5 --warmups 1 --ewald-level "$level" \
    > "/tmp/treepm_small_ewald_level${level}.jsonl"
done
```

Check image convergence and unchanged existing Ewald budgets before interpreting
total accuracy. For a timing-only perf attribution, omit force/oracle output:

```bash
OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE perf record -g -o /tmp/treepm_1698_perf.data -- \
  "$TREEPM_BENCH" --particles /tmp/snap_1698_particles.csv --mesh 64 --leaf 16 \
  --policy strict --kernel analytic --accounting full --repeats 5 --warmups 1 \
  "${TREEPM_PHYSICAL[@]}"
perf report -i /tmp/treepm_1698_perf.data
```

P0's supplied 52.958s strict and 40.616s fallback medians are comparable only
when exact frozen state, split, softening, acceptance, targets, build and hardware
match. A new traversal changes sums; do not demand old bitwise force checksums.

### Hierarchical scientific convergence

The named integration/validation commands above exercise actual workflow
dispatch, restart and MPI equivalence. Keep `hierarchical_max_rung=0` for the
reference; compare nonzero rung input decks with the same physics/initial state
and successively halved coarse limits. For larger user-run matched snapshots,
with an explicitly chosen scale target in inverse snapshot-coordinate units:

```bash
python3 tools/hierarchical_dmo_qualification.py \
  --initial "$DMO_INITIAL_SNAPSHOT" --global-reference "$DMO_GLOBAL_FINAL" \
  --hierarchical "$DMO_HIERARCHICAL_FINAL" --refined "$DMO_HIERARCHICAL_REFINED_FINAL" \
  --mesh 64 --bins 6 --k-min "$QUALIFICATION_K_MIN" --k-max "$QUALIFICATION_K_MAX" \
  --max-power-relative-error 0.02 --max-growth-relative-error 0.01 \
  --output /tmp/hierarchical_dmo_convergence.json
```

This offline tool uses numpy/h5py supplied by the user, not installed here. It
does not run simulations, choose a scientifically validated k interval, infer
unit conversions or replace production spectrum diagnostics. The current
written in-repo linear test only targets the fundamental six modes.

## Limitations and remaining acceptance work

- All C++ and new scientific assertions await compilation and execution; source
  inspection is not compile evidence. Actual speedups and lookup error targets
  are unmeasured. Analytic remains default and lookup is experimental.
- Dedicated long-time isolated two-body orbit, nonlinear cosmological scale-range
  convergence, large-rung efficiency, broad restart topologies, and hierarchical
  3/8-rank qualification are not added here. The short close-pair workflow is
  not a substitute for an isolated-orbit error study. Existing all-rank gravity
  fixtures remain required. Gas/source multirate KDK is outside this DMO path.
- New lookup/governor tests cover bounded allocation retention/release, while
  distributed allocation-failure injection and full-runtime high-water acceptance
  still need the existing failure/memory qualification matrix.
- The new benchmark is single-owner. Its optional history and full accuracy
  arrays are diagnostic-only O(N) storage, explicitly excluded from production
  hot paths. It exports one single-file HDF5 group; multi-file snapshot sets and
  historical external oracle adapters remain user-side integration extensions.
- No unnecessary synchronized bootstrap was established; none was removed.
  Coarse displacement restricts a local frozen-force estimate, not future-force
  growth. Current hierarchical runtime capability remains provisional.

## Suggested pull request and artifact

Branch: `perf/treepm-unified-hotpath-and-hierarchical-qualification`  
Title: **TreePM Performance Recovery: Fast Residual Kernel, Traversal Optimization, Split Tuning, and Hierarchical KDK Qualification**

The source-only patch ZIP contains only the modified/new files listed below,
with repository-relative paths. It excludes all build/runtime products, snapshots,
oracle outputs and unrelated files. The modified working tree is authoritative.

## Exact changed-file manifest

<!-- CHUI_CHANGED_MANIFEST -->

38 files: 32 modified, 6 new.

| Repository-relative path | Kind | Change |
|---|---|---|
| [CMakeLists.txt](../../CMakeLists.txt) | Build registration | Register the sweep benchmark and hierarchical/Ewald qualification tests. |
| [README.md](../../README.md) | Documentation | Describe opt-in provisional hierarchy and link the unqualified source handoff. |
| [bench/bench_tree_pm_force_error_map.cpp](../../bench/bench_tree_pm_force_error_map.cpp) | Benchmark | Repair the direct reference to use the implemented softened residual law. |
| [bench/bench_tree_pm_sweep.cpp](../../bench/bench_tree_pm_sweep.cpp) | New benchmark | Add full-target frozen-source split-aware timing and force diagnostics. |
| [docs/build_instructions.md](../../docs/build_instructions.md) | Documentation | Provide actual deferred build and test presets. |
| [docs/configuration.md](../../docs/configuration.md) | Documentation | Document typed flags, defaults and normalized hash implications. |
| [docs/memory_governance.md](../../docs/memory_governance.md) | Documentation | Document table allocation, commitment and lifetime. |
| [docs/profiling.md](../../docs/profiling.md) | Documentation | Document sweep inputs, counters and measurement conventions. |
| [docs/repair/treepm_unified_performance_recovery_20261010.md](../../docs/repair/treepm_unified_performance_recovery_20261010.md) | New documentation | Record implementation, invariants, commands, limitations and this manifest. |
| [docs/repair_open_issues.md](../../docs/repair_open_issues.md) | Documentation | Record pending scientific, distributed and performance acceptance gates. |
| [docs/repair_state_recap.md](../../docs/repair_state_recap.md) | Documentation | Record source progress without inventing validation evidence. |
| [docs/restart_checkpointing.md](../../docs/restart_checkpointing.md) | Documentation | Document unchanged schema and config/continuation boundaries. |
| [docs/time_integration.md](../../docs/time_integration.md) | Documentation | Document coarse crossing restriction and dispatcher safeguards. |
| [docs/tree_pm_coupling.md](../../docs/tree_pm_coupling.md) | Documentation | Document new API, interpolation, geometry, ordering and counter contracts. |
| [docs/validation_plan.md](../../docs/validation_plan.md) | Documentation | Document written assertions and pending qualification scope. |
| [include/cosmosim/core/config.hpp](../../include/cosmosim/core/config.hpp) | Production configuration | Add the two conservative typed TreePM control defaults. |
| [include/cosmosim/core/time_scheduler.hpp](../../include/cosmosim/core/time_scheduler.hpp) | Public interface | Expose the unit-aware displacement timestep input/helper. |
| [include/cosmosim/gravity/gravity_memory.hpp](../../include/cosmosim/gravity/gravity_memory.hpp) | Public interface | Expose optional coefficient capacity in gravity estimates. |
| [include/cosmosim/gravity/tree_gravity.hpp](../../include/cosmosim/gravity/tree_gravity.hpp) | Public interface | Expose skipped MAC work in aggregate profiles. |
| [include/cosmosim/gravity/tree_pm_coupling.hpp](../../include/cosmosim/gravity/tree_pm_coupling.hpp) | Public interface | Define flags, block control, counters and governed table ownership. |
| [include/cosmosim/gravity/tree_pm_split_kernel.hpp](../../include/cosmosim/gravity/tree_pm_split_kernel.hpp) | Production/public interface | Stabilize analytic limits and add the bounded direct-pair coefficient table. |
| [src/core/config.cpp](../../src/core/config.cpp) | Production configuration | Parse/validate/default/normalize both flags through the existing pipeline. |
| [src/core/time_integration.cpp](../../src/core/time_integration.cpp) | Production | Enforce block preparation invariants and derive the crossing timestep root. |
| [src/gravity/gravity_memory.cpp](../../src/gravity/gravity_memory.cpp) | Production | Include physical coefficient capacity in preflight arithmetic. |
| [src/gravity/internal/tree_pm_acceptance.hpp](../../src/gravity/internal/tree_pm_acceptance.hpp) | Production internal | Stage MAC work while distinguishing forensic and skipped accounting. |
| [src/gravity/internal/tree_pm_geometry.hpp](../../src/gravity/internal/tree_pm_geometry.hpp) | New production internal | Reuse axis distances and preserve cutoff rounding boundaries. |
| [src/gravity/tree_pm_coupling.cpp](../../src/gravity/tree_pm_coupling.cpp) | Production | Use octant traversal, reused geometry, optional interpolation and governed ownership. |
| [src/workflows/gravity_runtime.cpp](../../src/workflows/gravity_runtime.cpp) | Production | Wire typed flags, distinct commitments and accounting metadata. |
| [src/workflows/time_coordinator.cpp](../../src/workflows/time_coordinator.cpp) | Production | Apply collective independent coarse PM displacement restriction. |
| [tests/integration/test_hierarchical_timestep_regression.cpp](../../tests/integration/test_hierarchical_timestep_regression.cpp) | Tests | Assert displacement roots and actual dispatcher drift/activation/epoch invariants. |
| [tests/integration/test_hierarchical_treepm_workflow.cpp](../../tests/integration/test_hierarchical_treepm_workflow.cpp) | New tests | Assert actual cosmological workflow refinement, replay, restart and PM boundaries. |
| [tests/integration/test_tree_pm_coupling_periodic.cpp](../../tests/integration/test_tree_pm_coupling_periodic.cpp) | Tests | Assert pair identity, boundaries, table governance and thread reproducibility. |
| [tests/unit/test_config_parser.cpp](../../tests/unit/test_config_parser.cpp) | Tests | Assert new typed defaults, normalization and invalid boolean rejection. |
| [tests/unit/test_tree_pm_split_kernel.cpp](../../tests/unit/test_tree_pm_split_kernel.cpp) | Tests | Assert analytic limits, interpolation budgets, input rejection and geometry/MAC equivalence. |
| [tests/validation/test_dmo_zeldovich_workflow.cpp](../../tests/validation/test_dmo_zeldovich_workflow.cpp) | Scientific tests | Add hierarchical fundamental growth/power refinement and MPI artifact equivalence. |
| [tests/validation/test_tree_pm_ewald_accuracy.cpp](../../tests/validation/test_tree_pm_ewald_accuracy.cpp) | Scientific tests | Apply existing Ewald budgets to strict/adaptive analytic/lookup candidates. |
| [tools/hierarchical_dmo_qualification.py](../../tools/hierarchical_dmo_qualification.py) | New diagnostic tool | Assert matched offline global/hierarchical spectrum and position convergence. |
| [tools/treepm_sweep.py](../../tools/treepm_sweep.py) | New diagnostic tool | Export frozen sources and save reproducible split/leaf/block/thread sweep reports. |
