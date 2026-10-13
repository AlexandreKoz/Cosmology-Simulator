# TreePM deep-hot-path test regression closure — 20261012

Repair mode; authoritative checkout: `perf/treepm-deep-hotpath-recovery`,
starting commit `a2a70ea`, initially clean working tree. The reference ZIP does
not override the checkout. All changes remain in the working tree; no commit,
reset, merge or push was performed.

**Every new or changed test: NOT RUN — deferred to maintainer.**
Compilation, CPU tests, MPI tests, Ewald accuracy and performance benchmarks:
**NOT RUN — deferred to maintainer.** Static inspection is not execution evidence.

## Defect dispositions

| Defect | Original behavior | Root cause | Repair | Validation |
|---|---|---|---|---|
| A: split-only force fixtures | Threw `short-range split surface cannot request a PM refresh/capture or zoom correction` | `short_range_only=true` passed to ordinary refresh-requesting `solveActiveSet`, including helper-mediated calls | Explicit full/split test dispatch; split uses `solveActiveSetWithPmCadence(..., false)` for coincidence, mixed-softening/cutoff/thread tests, SIMD/table/fallback/reuse/refit fixtures | FIXED IN SOURCE — NOT YET TESTED; periodic integration |
| Split safety controls | No dedicated positive/contradictory-request pair | Fixture failure prevented intended numerical checks | Accepted finite, nonzero pair with zero PM solves/reuses; reject ordinary and explicit refresh, valid zoom request and PM capture using `invalid_argument` plus split-guard context | FIXED IN SOURCE — NOT YET TESTED; periodic integration |
| B: experiment base | Loader failed at missing `mode.mode` | Incomplete authoritative input | Small explicit cosmo-cube, generated-IC, TreePM, source-free base; all ten variants checked, reloaded through generator/typed loader, compared after only named changes | FIXED IN SOURCE — NOT YET TESTED; experiment integration |
| Experiment invalid/overwrite controls | Existing-directory failure only checked nonzero exit | Failure could have occurred before output-directory protection | Check exit 1 and output-directory diagnostic, preserve sentinel/manifest/all deck hashes; missing mode and invalid tuning must fail in their own contexts without publication | FIXED IN SOURCE — NOT YET TESTED; experiment integration |
| Experiment workspace collision | Fixed directory removed at test start | Concurrent invocations shared destructive scratch ownership | Random per-invocation workspace; cleanup only its own successful run | FIXED IN SOURCE — NOT YET TESTED; experiment integration |
| C: tuning/flag round trips | Failed before positive assertions; empty defaults fixture also invalid | Missing required `[mode]` / `mode.mode` | Complete typed baselines; block/leaf boundaries, all performance booleans, normalized reload, FNV/SHA provenance equality | FIXED IN SOURCE — NOT YET TESTED; config unit |
| False-positive TreePM negatives | Missing-mode failure could satisfy generic rejection | Invalid baseline and catches lacking target context | Admit canonical baseline first, replace one value, catch `ConfigError`, require intended key or validation context; legal controls retained | FIXED IN SOURCE — NOT YET TESTED; config unit |
| D: PM expected storage | Old three-complex/four-full-scalar accounting and obsolete 512³ lower bound | Tests predated retained `PlanResources` layout | Independently derived real/two-complex/global-axis expectations; rank-local slab/transposed geometry, odd sizes, zero-width slabs, serial/distributed 512³ planning and overflow checks | FIXED IN SOURCE — NOT YET TESTED; gravity unit |
| PM actual retained storage | Formula alone could miss capacity/report drift | Needed comparison with live allocation owner | Tiny odd-grid serial report check; periodic integration checks slab and transposed plans on active FFTW-MPI ranks, including backend query flag and physical capacity | FIXED IN SOURCE — NOT YET TESTED; gravity unit and periodic MPI matrix |
| PM estimator/owner | Current source retains real, Fourier, potential and axis vectors | Old tests describe obsolete allocations | Estimator and production MemoryGovernor path preserved after tracing allocation/report sites | ALREADY CORRECT — SOURCE VERIFIED; runtime/backend qualification pending |
| Hierarchical rung diagnostic | Expected `does not yet carry per-element kick/drift epochs` | Test predates opt-in collisionless hierarchy | Check current cosmo-cube/source-free admission category; valid rung 1 and 12 controls, invalid range/type cases; default remains zero | FIXED IN SOURCE — NOT YET TESTED; config unit |
| Hierarchical integration/restart guards | Audit suspected a stale integration rejection | Actual stale string was in config unit; inspected integration strings match current guard sites | Preserve identity-row, missing-event, stale-field, open-KDK and active-interval negative tests | ALREADY CORRECT — SOURCE VERIFIED; all runtime cases pending |
| Release assertions | Gravity numerical/memory checks used standard `assert` | Manual NDEBUG builds could erase checks | Convert gravity checks to established `CHUI_REQUIRE`; existing `unit_assertions_enabled` now deliberately defines NDEBUG after its target-policy check and probes evaluation/rejection with an independent exit condition | FIXED IN SOURCE — NOT YET TESTED; Release assertion gate |
| Additional numerical/MPI failures | No additional execution evidence in this session | Execution explicitly prohibited | Preserve references/tolerances and request exact failing command/log if later qualification exposes another issue | UNCONFIRMED — NEEDS RUNTIME EVIDENCE |

## Owners and unchanged production contracts

- `gravity::TreePmCoordinator`, public `tree_pm_coupling.hpp`: force split and
  PM cache; ordinary `solveActiveSet` requests refresh. The explicit cadence
  surface controls refresh. `short_range_only` forbids refresh, capture and
  zoom; the production guard is retained. A split-only request bypasses PM
  rather than asking to reuse an uninitialized PM cache.
- `core::FrozenConfig` / `loadFrozenConfigFromString/File`, `src/core/config.cpp`:
  typed admission, normalization and provenance authority. Required mode,
  numeric ranges, booleans and enums remain production rules.
- `PmSolver::Impl::PlanResources`, `src/gravity/pm_solver.cpp`: retained FFT
  allocations and actual vector-capacity reports. `src/gravity/gravity_memory.cpp`
  composes these with separate mesh/tree/cache/communication/unknown reserves.
- `parallel::PmSlabLayout`: rank-local real-space ownership; the current
  `pencil` option selects transposed x/y slabs, not a true 2-D FFT backend.
- `core::HierarchicalTimeBinScheduler` and `StepOrchestrator`: scheduler truth,
  active intervals and synchronized PM event commit. Workflow restart admission
  is in `src/workflows/output_restart_runtime.cpp`; IO owns persistent schemas.
- `tools/treepm_experiment_configs.cpp`: typed, normalized experiment generation;
  tests exercise this executable without changing its physics configuration.

No production implementation, public interface, configuration default/key,
restart/snapshot schema, MPI ownership, provenance algorithm, force law or
scientific tolerance changes. Reproducibility of production simulations is
unchanged. Random scratch identifiers affect only test artifact isolation.
No new benchmark hooks are needed for this test-only repair.

## Typed fixture and negative-test details

The tuning/performance tests now include `mode=cosmo_cube`, `ic_file=generated`,
`gravity_solver=treepm` and a legal 16³ mesh. The former empty defaults input is
replaced with the same valid baseline. All emitted canonical text is reloaded;
round trips retain normalized text, stable config hash and normalized SHA-256.

`requireNumericsRejected` first reloads the complete normalized baseline, then
replaces exactly one existing value and checks a `ConfigError` with the intended
key/context. Coverage includes zero, negative, unsupported block choices, leaf
1025, overflow text, fractional integral text, NaN/infinities, invalid and legacy
boolean spellings, invalid enums, unsafe PM cadence, cutoff ambiguity,
coarsest-cell split compatibility and excessive softening. Legal block sizes are
16/32/64/128/256; leaves include 1, 8 and 1024; adaptive angle controls include
0.25 and its supported 0.5 boundary. No artificial opening-theta upper bound or
forbidden optimization-flag combination is introduced.

The experiment test verifies `reference`, `optimized_matched`,
`matched_gaussian`, `matched_full_residual`, `matched_simd`, `matched_radial`,
`matched_reuse`, `matched_refit`, `adaptive_025`, and
`adaptive_050_experimental`. It checks TreePM, global rung-zero KDK, block/leaf
64/16, shared mesh/softening/physical inputs, distinct run IDs and exact intended
flag differences. Only the last explicitly experimental variant uses angle 0.5.
After removing only each variant's advertised changes, the entire normalized
deck must equal reference. Each deck is also supplied back to the actual typed
loader through the generator and its same named variant must be byte-identical.
The supplied base explicitly uses the generator's matched tolerance 0.005;
reference and optimized-matched retain identical acceptance policy.

## PM retained allocation derivation

Source evidence: `PlanResources`, `planForLayout`,
`determinePmPlanLocalStorageLayout`, `ensureSpectralAxisMetadata`,
`Impl::appendMemoryReport`, `estimatePmPlanResourcesMemory` in
[`pm_solver.cpp`](../../src/gravity/pm_solver.cpp), and the plan-resource entry
in [`gravity_memory.cpp`](../../src/gravity/gravity_memory.cpp).

Let `Kz = floor(nz/2)+1`, `A` be allocated local complex coefficients and
`M = 2*(nx+ny+Kz)` global-axis doubles (three wave-number and three window lanes).
One plan retains:

| Storage | Serial elements | Distributed elements |
|---|---|---|
| `real` doubles | `nx*ny*nz` (out-of-place, no serial padding) | `2*A`, with real z stride `2*Kz` |
| `fourier` complex doubles | `nx*ny*Kz` | `A` |
| `potential_k` complex doubles | `nx*ny*Kz` | `A` |
| Spectral-axis metadata doubles | `M` | `M`, replicated per rank |
| Full scalar spectral arrays | zero | zero |

Thus bytes are `R*sizeof(double) + 2*A*sizeof(complex<double>) + M*sizeof(double)`.
Serial `A=nx*ny*Kz`. Distributed preflight without active MPI uses the fixed
ceil-block slab partition, at least one allocation coefficient, and for
transposed slabs the maximum of input x-slab and output y-slab requirements.
Logical output and allocated capacity can differ. Active FFTW-MPI uses the
backend allocation query and checks ownership against the rank's layout.

For a 512³ mesh on eight equal slabs, each rank has 64 real x planes and
`A=64*512*257`; the formula above includes padded real storage and global-axis
metadata. This is a planning test and allocates no 512³ arrays.

Live tests check logical bytes against independently counted buffers and require
reported retained capacity >= logical size and high-water >= capacity. They do
not pretend allocator excess capacity is zero. Fresh-plan preflight predicts
one plan's array requirements; live reports sum all retained cache entries.
FFTW plan internals/library workspace remain unknown external reserves, separate
from known CHUÍ arrays. Grid force/potential storage, cached plans, output overlap,
canonical state and communication are accounted in their existing owner reports;
these tests do not establish a complete process fit or reduce a 48 GiB estimate.

Overflow controls reach valid layout geometry first, then exercise real byte,
complex byte, metadata byte, padded-stride and total-sum overflow rejection with
`overflow_error` and matching context. Positive odd/non-cubic, truncated,
transposed and zero-slab geometries precede these controls.

## Hierarchical rejection evidence

The obsolete text occurs in
`testProductionHierarchicalRungsFailClosed` in `test_config_parser.cpp`.
`validateFrozenConfig` admits rungs 1–12 only for cosmo-cube DMO with the listed
source physics disabled. Rung zero remains default. The test still rejects the
original unsupported zoom/source-enabled state and adds supported controls.

Inspected integration contracts remain distinct:

- `test_hierarchical_time_bins.cpp`: malformed persistent vectors are rejected
  by `importPersistentState` with `invalid_argument`.
- `test_hierarchical_timestep_regression.cpp`: exact all-source identity rows
  are required before block mutation; synchronized endpoints need a refresh
  event and current step/scale field epoch before cadence commit.
- `test_hierarchical_treepm_workflow.cpp`: real dispatcher/restart path preserves
  opt-in force/refinement and provenance checks. Restart requires matching rung
  schemas/ticks, no active interval or pending transition, closed synchronized
  collisionless KDK and consistent drift epochs. Existing negative checks use
  these contexts rather than accepting arbitrary exception text.

No scheduler acceptance repair was justified by source evidence. Additional
hierarchical failures need their exact runtime log before changing expectations.

## Focused harness review

Reviewed the named integration/config/memory files, split-kernel unit, Ewald
qualification, PM plan implementation/reports, generator and CMake registrations.
Full-force helper calls remain ordinary refresh calls; retained PM-cadence tests
still use the explicit API and reject incompatible cache epochs. Independent
pair references, identities, mixed softening, cutoff/seam cases, scalar/SIMD,
lookup, reuse/refit and governed fallback assertions retain their thresholds.
The MPI routing fixtures retain rank-specific qualification gates. The new live
PM checks run on all ranks when FFTW-MPI is available; the serial memory unit
intentionally does not initialize MPI and cannot qualify active MPI allocations.

CTest already registers the repaired executables and generator script, so no
CMake registration change was necessary. Hierarchical workflow tests require
HDF5+FFTW; np2/np4 additionally require MPI and Parallel HDF5. Ewald acceptance is
FFTW-gated. `validation_tree_pm_recovery_ewald_accuracy` gates strict/adaptive
analytic/Gaussian lookup profiles; it does not independently certify every newer
full-residual/SIMD/radial mode. Their integration equivalence checks and broader
scientific qualification remain pending. No test is disabled or marked expected
failure by this patch.

## Manual validation order and exact registrations

**All commands below: NOT RUN — deferred to maintainer.** Run from repository
root. Stop at the first meaningful failure and retain the command and full log.
`--no-tests=error` prevents a missing registration from looking like success.

1. Static hygiene, then CPU-only Debug configuration/build:

   ```bash
   bash scripts/ci/check_repo_hygiene.sh
   cmake --preset cpu-only-debug
   cmake --build --preset build-cpu-debug
   ```

2. Fixture admission and cheap local invariants first:

   ```bash
   ctest --preset test-cpu-debug -R '^unit_config_parser$' --no-tests=error --output-on-failure
   ctest --preset test-cpu-debug -R '^integration_treepm_experiment_configs$' --no-tests=error --output-on-failure
   ctest --preset test-cpu-debug -R '^unit_tree_pm_split_kernel$' --no-tests=error --output-on-failure
   ctest --preset test-cpu-debug -R '^unit_tree_gravity$' --no-tests=error --output-on-failure
   ```

3. Periodic force, scheduler and composed reference workflow:

   ```bash
   ctest --preset test-cpu-debug -R '^integration_tree_pm_coupling_periodic$' --no-tests=error --output-on-failure
   ctest --preset test-cpu-debug -R '^integration_hierarchical_(time_bins|timestep_regression)$' --no-tests=error --output-on-failure
   ctest --preset test-cpu-debug -R '^integration_reference_workflow$' --no-tests=error --output-on-failure
   ```

4. Complete registered CPU suite after focused checks:

   ```bash
   ctest --preset test-cpu-debug --no-tests=error --output-on-failure
   ```

5. Release assertion effectiveness, including the registered NDEBUG probe:

   ```bash
   cmake --preset cpu-only-release
   cmake --build --preset build-cpu-release
   ctest --test-dir build/cpu-only-release -R '^(unit_assertions_enabled|unit_config_parser|unit_tree_gravity|unit_tree_pm_split_kernel|integration_tree_pm_coupling_periodic)$' --no-tests=error --output-on-failure
   ```

   There is no plain `test-cpu-release` preset. Existing test-target `-UNDEBUG`
   policy remains; the assertion probe deliberately defines NDEBUG inside its
   translation unit after verifying that policy. Its independent return code
   requires exactly-once evaluation and rejection by `CHUI_REQUIRE`.

6. Dependency-complete MPI/HDF5/FFTW Release build (C++20 toolchain, Ninja,
   MPI C/C++, Parallel HDF5 C, FFTW and FFTW-MPI development files):

   ```bash
   cmake --preset mpi-hdf5-fftw-release
   cmake --build --preset build-mpi-hdf5-fftw-release
   ctest --preset test-mpi-hdf5-fftw-release -R '^(unit_assertions_enabled|unit_config_parser|unit_tree_gravity|integration_treepm_experiment_configs|integration_tree_pm_coupling_periodic)$' --no-tests=error --output-on-failure
   ```

7. Independent Ewald reference and accelerated Gaussian recovery accuracy:

   ```bash
   ctest --preset test-mpi-hdf5-fftw-release -R '^validation_(periodic_ewald_reference|tree_pm_ewald_accuracy|tree_pm_ewald_accuracy_first_light_softening|tree_pm_recovery_ewald_accuracy)$' --no-tests=error --output-on-failure
   ```

   CTest supplies the registered finite-softening environment and recovery
   executable argument. These are single-process Ewald checks on the FFTW
   backend; they do not replace the distributed rank matrix.

8. Relevant registered MPI force, active PM geometry, scheduler/restart and
   workflow qualification, serialized to avoid rank oversubscription:

   ```bash
   ctest --preset test-mpi-hdf5-fftw-release -R '^integration_tree_pm_coupling_periodic_mpi_(two|three|four|eight)_rank$' --parallel 1 --no-tests=error --output-on-failure
   ctest --preset test-mpi-hdf5-fftw-release -R '^integration_(pm_periodic_mode_mpi_two_rank|pm_zero_width_slab_mpi_four_rank)$' --parallel 1 --no-tests=error --output-on-failure
   ctest --preset test-mpi-hdf5-fftw-release -R '^validation_phase2_mpi_gravity_(single|two|three|four|eight)_rank$' --parallel 1 --no-tests=error --output-on-failure
   ctest --preset test-mpi-hdf5-fftw-release -R '^integration_(hierarchical_treepm_workflow(_np[24])?|reference_workflow_distributed_treepm_mpi_two_rank|runtime_app_mpi_treepm_smoke_two_rank|restart_equivalence_treepm)$' --parallel 1 --no-tests=error --output-on-failure
   ```

Classify failures before changing anything: configure/dependency/toolchain
failure; compile/link failure; fixture admission/construction failure (before
force/scheduler assertions); real numerical failure after target path; MPI or
backend capability/launcher failure; unsupported environment/missing conditional
registration; scientific qualification still pending. An MPI-disabled run or a
CPU diagnostic DFT result is not production FFTW/MPI qualification. No benchmark
or cosmological simulation command was executed or is needed to claim this
source-only repair; performance and long-interval science remain pending.

## Changed-file manifest

- `tests/integration/test_tree_pm_coupling_periodic.cpp`
- `tests/integration/test_treepm_experiment_configs.cmake.in`
- `tests/unit/test_config_parser.cpp`
- `tests/unit/test_tree_gravity.cpp`
- `tests/unit/test_assertions_enabled.cpp`
- `docs/repair/treepm_deep_hotpath_test_regression_closure_20261012.md`
- `docs/repair_state_recap.md`
- `docs/repair_open_issues.md`
- `docs/validation_plan.md`

## Static checks and handoff

Source inspection used `rg`, `sed`, `cat`, Git status/branch/diff and allocation/
registration tracing. `git diff --check`: passed (no whitespace errors).
`./scripts/ci/check_repo_hygiene.sh`: exit 126, permission denied; the script is
tracked and present as mode 100644/0644. `bash scripts/ci/check_repo_hygiene.sh`:
passed, including artifact/naming, config aliases, preset inventory, FFTW-MPI
wiring and shell syntax checks. This invokes only static CMake script/preset
inspection, not CMake configuration. No unrelated script mode was changed.
No configure, build, CTest, Python test suite, MPI executable, scientific
benchmark or simulation was run. Python was used only to edit/package text.

Suggested branch: `fix/treepm-deep-hotpath-test-regressions` (the existing feature
branch may retain the repair if still unmerged). Suggested commit:
`Fix TreePM regression fixtures and scientific test reliability`.
Suggested PR title: **Repair TreePM Hot-Path Qualification Tests, Configuration
Fixtures, and PM Memory Assertions**. Leave review/merge to the maintainer.
