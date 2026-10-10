# P1–P5 qualification, snapshot integrity and test reliability repair — 20261010

Mode: repair. The authoritative checkout started clean on
`perf/treepm-unified-optimization`; neither archive was used. The branch and
accepted P1–P5 implementation were preserved. No commits, resets or branch
switches were performed.

Compilation: **NOT RUN — deferred by user**  
Tests: **NOT RUN — deferred by user**  
MPI qualification: **NOT RUN — deferred by user**  
Scientific convergence: **NOT RUN — deferred by user**  
Performance benchmarks: **NOT RUN — deferred by user**

Only source inspection, searches, editing, static Python grammar parsing,
`git diff --check`, and source packaging were performed. No CMake configuration,
build, test, benchmark, simulation, dependency installation or profiling was run.
Written checks below are acceptance criteria, not passing evidence.

## Defect disposition and source evidence

| Issue | Disposition | Repair or evidence |
|---|---|---|
| Python 3.9 `removeprefix` in MassTable path | FIXED IN SOURCE — NOT YET TESTED | `tools/treepm_sweep.py` delegates to `chui_snapshot.species_index`: exact `PartType[0-5]` validation works with Python 3.8, including explicit Masses. |
| Single HDF5 member mistaken for complete population | FIXED IN SOURCE — NOT YET TESTED | Both offline consumers use `tools/chui_snapshot.py`, following production `src/io/snapshot_set.cpp` v2/v3 closure and `snapshot_hdf5.cpp` v6 metadata. |
| Incompatible sources or opaque external total reference | FIXED IN SOURCE — NOT YET TESTED | Export/checksum sidecars, scientific normalized provenance comparison, explicit unit conversions and bound external-reference identity/evidence. |
| Missing Gaussian coefficient at r=0, epsilon>0 | FIXED IN SOURCE — NOT YET TESTED | `tree_pm_split_kernel.hpp` subtracts the finite Gaussian limit using scaled arithmetic; nonzero-radius formula unchanged. |
| Coincident vector multiplication can form 0*infinity | FIXED IN SOURCE — NOT YET TESTED | Production direct-pair traversal counts exact coincident interactions and contributes zero before coefficient multiplication. |
| Two new default keys invalidate legacy normalized hash | FIXED IN SOURCE — NOT YET TESTED | Explicit `treepm_normalized_defaults_v1` migration, verified v23/v7 only, literal projection and authoritative roundtrip. |
| Opaque-hash/older/partial checkpoint migration | BLOCKED — equivalence cannot be proven by the implemented contract | Refused; no warning-based admission. Exact matching older checkpoints retain their existing reader/admission rules. |
| Hierarchical endpoint accepts missing/stale PM directive | FIXED IN SOURCE — NOT YET TESTED | Core dispatcher requires a solver-executed coarse PM refresh with exact step/scale epoch before PM cadence commit. |
| Recovery CTest can depend on an environment switch | FIXED IN SOURCE — NOT YET TESTED | Mandatory `--recovery-qualification` argument, four strict/adaptive x analytic/lookup gates, direct-pair and accepted-multipole assertions. |
| NaN passes `abs(error)>tolerance` comparators | FIXED IN SOURCE — NOT YET TESTED | Restart harness and periodic vector comparator require finite values and positive-form bounds; negative controls added. |
| Official Release CMake assertions disappear | ALREADY CORRECT — NO CHANGE REQUIRED to build flags | Existing test-only `-UNDEBUG`/`/UNDEBUG`, applied to all `test_` executables, plus compile-time `unit_assertions_enabled` sentinel. Priority sources additionally use always-active `CHUI_REQUIRE`. |
| CTest defaults overwrite explicit timeout/directory | FIXED IN SOURCE — NOT YET TESTED | Default properties now apply only if unset. Numeric and `_npN` rank names receive processor counts/labels. |
| CI misses recovery gate / lacks required Parallel HDF5 | FIXED IN SOURCE — NOT YET TESTED | Existing PM job selects recovery, kernel, reader, hierarchy and restart cases; MPI jobs install the parallel provider already required by presets. |
| Growth/power tests are merely smoke tests | ALREADY CORRECT — NO CHANGE REQUIRED to scientific bounds | `test_dmo_zeldovich_workflow.cpp` asserts independent particle-mode growth/phase and production spectrum growth; its hierarchical fixture checks refinement and six fundamental modes. This is narrow linear coverage, not nonlinear certification. |
| Historical unregistered analysis/AMR sources | ALREADY CORRECT — NO CHANGE REQUIRED to current qualification registration | See inventory below: superseded sources are not counted as executed tests. |
| Other numerical tolerance/order failures without logs | UNCONFIRMED — NEEDS RUNTIME EVIDENCE | No golden values or tolerances changed to match output. |
| Native historical P0 oracle recomparison | BLOCKED — external tools/dataset format absent | See manual force handoff; no historical implementation or command invented. |
| Long-time orbit/nonlinear/many-rung scientific qualification | NOT COMPLETED — requires numerical campaigns | Existing provisional boundary remains; no capability promotion. |

## Numerical and runtime contracts

The direct coefficient remains

```text
K_SR(r,epsilon,a_split) = (r^2+epsilon^2)^(-3/2) - L(r,a_split)/r^3
L = erf(r/(2a_split)) - r/(sqrt(pi)*a_split)*exp(-r^2/(4a_split^2))
```

For positive finite epsilon the exact coincident coefficient is
`epsilon^-3 - 1/(6*sqrt(pi)*a_split^3)`. The previous branch returned only
`epsilon^-3`. Subtraction is scaled by the smaller length before sequential
cubed division, avoiding `inf-inf` when the true result exceeds double range.
Representational overflow keeps the correct sign; arbitrary precision is not
promised. For epsilon=0 at exact coincidence the established zero coefficient
convention remains. Both kernel choices share the coincident path. Production
explicitly handles an exact zero displacement vector; it does not drop the
interaction from pair counters. The positive-radius series/transition, table,
Gaussian normalization, Plummer pair epsilon, inclusive cutoff and analytic
multipole derivatives are unchanged. No Gaussian term is subtracted twice.

P1 child ordering, P3 geometric/MAC decisions, full second-moment quadrupoles,
TSC/window deconvolution, cosmological normalization and periodic axis conventions
are retained. This repair adds no general force summation reorder. Historical
pre-P1 numerical checksums can already differ by legitimate accumulation roundoff.
Lookup stays opt-in, full MAC diagnostics stay default, and
`hierarchical_max_rung=0` stays the safe default.

The hierarchical guard runs in the existing preparation/failure-coordination
surface. It does not change kick or drift equations or introduce rank-local
collective branching. Missing PM events cannot be confused with a passing epoch
check. The workflow's existing ownership, active-set, source-generation and PM
cache checks remain authoritative.

## Complete logical science snapshot contract

The shared reader is offline tooling, not a production snapshot/restart reader
or an alternative `.param.txt` parser. It consumes:

- Modern CHUI science v6 with a valid production completion marker, whether
  selected by marker, one member, or a directory containing exactly one marker.
- Manifest v2 with validated body SHA-256, exact ordered members, member SHA-256,
  sizes, local/global counts and scientific identity.
- Manifest v3 for one distributed science file with
  `distributed_science_readback_v1`; current member SHA-256 is recorded as
  consumed provenance, without claiming the marker stores a byte hash.
- Explicit external legacy or CHUI `gadget_arepo_v1` through v5 input using a
  declared dialect and complete known unit metadata/override. Header member
  counts/indices still must establish the complete population. Unsupported
  revisions and malformed modern inputs cannot bypass closure with a legacy flag.

It rejects missing/duplicate/inconsistent members, bad manifests, inconsistent
epoch/geometry/generation/config identifiers, bad counts, duplicate persistent
IDs across members and species, invalid shapes/types, nonfinite positions or
velocities, missing/invalid MassTable fallback, and nonpositive masses. IDs are
sorted independently of row order. Checkpoint file kinds/state groups are refused.
All populated types are checked, even when only PartType1 is retained. Empty
selected populations fail explicitly. Optional source identity/generation fields
and normalized-config representations/fingerprints are checked where present.

Rectangular box axes are retained. Stored units, optional CGS/a/h metadata and
velocity conventions must agree. CHUI-native values retain declared code units;
interoperable stored coordinates/masses divide by h to recover the declared
internal units, and stored velocity multiplies by sqrt(a). The reader returns
stored arrays and records conversion factors; export applies its explicit chosen
length/mass factors in chunks. Epsilon is supplied in solver comoving length
units and stays fixed across mesh sweeps. Hierarchical k is explicitly in inverse
stored-coordinate units, with identical units/dialect required for all states.

Provenance contains exact marker/member paths and SHA-256 checksums, closure
schema, global/selected counts, sorted selected/all-species uint64 ID fingerprints,
epoch, geometry, cosmology, dialect, units/conversions and available config/source
identifiers. Files changing size/mtime during consumption are refused. Export
sidecars additionally bind frozen CSV and target checksums, epsilon and solver box.
`--fixed-target-ids` preserves original oracle ID order and rejects duplicates,
missing IDs and malformed rows. Sweep admission checks the sidecar rather than
silently accepting altered targets or epoch/units/softening.

External total-force comparison requires a JSON manifest with these exact fields:

```text
reference_total_sha256  particles_sha256  target_ids_sha256
scale_factor  box_solver_units:[Lx,Ly,Lz]  epsilon_comoving  g_code
reference_kind: independent_periodic_ewald | converged_total
qualification_evidence: nonempty description/path of independent convergence evidence
```

This binds scientific identity and explicit evidence; a declaration is not itself
a proof of Ewald convergence. Do not substitute the TreePM solver as its own total
reference. The benchmark's separate long-double direct residual is regenerated
for each changed split. Wall time alone cannot certify a configuration.

Hierarchical snapshot comparison requires matching normalized scientific
provenance. Its literal comparison ignores output/parallel sections and only
named rung/time-resolution/end-condition keys, allowing matched convergence
runs while retaining mesh, split, softening, cosmology and initial-condition
settings. It is deliberately a comparison of emitted provenance, not another
configuration loader. Epochs, IDs, masses, geometry and unit conventions must
also agree as appropriate to initial versus final states.

Python staging is 65,536 rows by default. Persistent arrays contain only selected
XYZ/mass/optional velocity and IDs; one all-species ID lane enables global
uniqueness and is released before selected sorting. No retained production cache,
particle layout or MemoryGovernor accounting changes are introduced by these
external allocations. Full selected state and FFT estimator still consume O(N)
plus mesh memory; large-problem memory fit is not claimed.

CMake previously declared no numeric global Python minimum. These two workflows,
the shared reader and their new regression explicitly support Python >=3.8.
NumPy/h5py must support the chosen interpreter. No unrelated repository Python
support policy was redefined; no P0 tools are installed in this checkout to scan.

## Restart migration and reproducibility

The reader sets a transient verification receipt only after complete logical
payload validation and integrity checks. Mismatched hashes can migrate only
verified `cosmosim_restart_v23` / `provenance_v7` checkpoints with full original
normalized text, matching stable/strong identities, a closed restart-safe global
KDK boundary and rung-zero schedulers. The current typed config must retain lookup
false/full diagnostics true. Deleting exactly those two canonical lines must
reproduce the original normalized text byte for byte; authoritative config loading
must regenerate the exact current normalization. Other missing keys, changed
physics/output/topology or unsupported historical defaults fail.

Raw read config/provenance are preserved. New outputs retain a distinct versioned
compatibility block in existing `derived_runtime_state`: original stable/strong
config hashes, original payload SHA, hex-encoded original text and serialized
provenance, plus distinct derived fingerprints. Later exact resumes preserve it.
No persistent schema is added. Normalized migration text is capped at 1 MiB per
representation and original serialized provenance at 4 MiB. This is bounded
metadata owned by the workflow/report and existing provenance outputs; it adds
no per-particle storage. Same-topology, source epochs, cache authority and all
existing restart checks still run under collective failure propagation.

Exact admission compares literal normalized text in addition to the stable hash.
Old opaque hashes, incomplete provenance, unverified synthetic payloads, nondefault
new settings, partial default-key transitions and hierarchical migration fail
closed. The implemented compatibility path intentionally does not promise all
historical checkpoints can resume.

## Test reliability inspection and regression matrix

Static inventory found 164 `test_*.cpp` sources. 132 still contain ordinary
`assert()`, protected for official targets by the existing test-only CMake flags.
There was no confirmed disabled-assert defect in those configured targets.
The new helper retains the original conditions/tolerances and evaluates each
expression once; it does not change production NDEBUG or add a framework.

Three source names lack current literal CMake registrations:
`test_validation_amr_cross_code_contract.cpp` and
`test_validation_amr_hydro_convergence.cpp` contain the same test logic as their
registered unit successors (only helper include paths differ). The old
`test_validation_phase3_power_spectrum.cpp` has an invalid indexed comparison of
8/16 meshes with equal bin counts and shared artifact cleanup. It is superseded
by registered `validation_power_spectrum_mesh_consistency`, which explicitly
matches physical bin edges/centers/mode counts below 0.75 coarse axis Nyquist,
uses the frozen tolerance table and tests incompatible-k refusal. The historical
source is preserved, not resurrected as a qualification gate. Eight restart
sources are registered via CMake foreach, not missing tests. This inventory is
not a claim every test in the repository has been scientifically rederived.

| Test | Previous problem | Repair | Validation status |
|---|---|---|---|
| `unit_tree_pm_split_kernel` | Coincident residual expectation omitted Gaussian coefficient; new edge/invalid/history coverage needed | Independent long-double zero/approach/lookup assertions, extrema, immutable lookup, invalid parameters and absent/invalid adaptive-history fallback; always-active requirements | NOT RUN — deferred by user |
| `integration_tree_pm_coupling_periodic` | NaN could pass vector comparator | Finite positive-form comparison and NaN/Inf negative controls; exact coincident vector/counter regression; MPI failure diagnostics; existing physical bounds retained | NOT RUN — deferred by user |
| `validation_tree_pm_recovery_ewald_accuracy` | Recovery selection depended on environment; profiles were not all mandatory accuracy gates | Mandatory CLI selection, four unchanged Ewald accuracy gates, selected direct-pair path and separate accepted-multipole requirements; explicit experimental output | NOT RUN — deferred by user |
| Baseline/finite-softening Ewald cases | Already independent and gated by FFTW | Certified baseline force/conservation bounds retained; no fallback certification | NOT RUN — deferred by user |
| `integration_hierarchical_timestep_regression` | Missing/stale PM endpoint not rejected in core commit | Actual dispatcher fault injection; requires precise refusal, no cadence/step commit; always-active existing drift/momentum tests | NOT RUN — deferred by user |
| `integration_hierarchical_treepm_workflow`, `_np2`, `_np4` | Legacy default-key migration untested; invalid resumed synchronization needed negative checks | Real workflow checkpoint write/read/migration/resume and phase equivalence; incompatible physics/flags/provenance/integrity/open-state refusals | NOT RUN — deferred by user |
| `integration_hierarchical_time_bins` | Relied on build-level assertion protection | Always-active original scheduler assertions | NOT RUN — deferred by user |
| `unit_config_parser` | New boolean invalid-value cases incomplete | Invalid strings rejected by typed parser with named key; defaults and normalized roundtrip preserved; always-active assertions | NOT RUN — deferred by user |
| `integration_restart_checkpoint_roundtrip`, `integration_restart_equivalence_treepm` | Relied on build-level assertion protection | Always-active original integrity/continuation requirements | NOT RUN — deferred by user |
| `integration_restart_equivalence_harness` | NaN/Inf could pass scalar comparator | Finite operands/tolerance and `abs<=tol`, explicit nonfinite refusal cases; unchanged tolerances | NOT RUN — deferred by user |
| `integration_snapshot_hdf5_roundtrip`, `integration_snapshot_set_mpi_*_rank` | Relied on build-level assertions; MPI errors could lose diagnostics | Always-active checks, named missing-member refusal, root diagnostic/abort handling; empty rank and collective closure preserved | NOT RUN — deferred by user |
| `python_treepm_snapshot_tools` | No standalone complete-snapshot regression independent of bindings | New production-shaped v2/v3 positive/negative fixtures, both consumer ingestion, units/IDs/provenance/reference identity/fixed-target checks | NOT RUN — deferred by user |
| `validation_dmo_zeldovich_workflow*`, `validation_hierarchical_dmo*` | Suspected smoke-only coverage | Existing independent growth/phase/spectrum/refinement and rank artifact fixtures inspected and retained; no new physical certification claimed | NOT RUN — deferred by user |
| `unit_assertions_enabled` | Suspected Release false-green behavior | Existing compile-time sentinel and test-only UNDEBUG verified in source; no global flag changes | NOT RUN — deferred by user |

The recovery profiles use the existing relative-L2/P99 accuracy targets and
periodic/reference invariants. The established certified baseline still enforces
its net-force threshold. Experimental adaptive tests are not relabeled certified
momentum-conservation tests: the baseline symmetric 1e-12 conservation bound is
not silently asserted to certify adaptive approximations. Diagnostic matrix
alternatives remain clearly diagnostic. Multipole acceptance is mandatory on a
separate independent-target fixture for each of the four profiles.

CTest preserves explicit directory/timeout properties, sets numeric/_npN MPI
processor accounting, and registers the offline Python gate independently of
bindings. Python <3.8/no interpreter produces an explicit configure status;
missing NumPy/h5py gives exit 77 and CTest SKIP. Malformed data do not skip. A
CPU-only run does not register FFTW scientific certification; direct executable
recovery selection without FFTW fails. Existing runtime smoke dependency guards
and real MPI launchers are retained. CI now selects the new gate and installs the
Parallel-HDF5 provider required by its existing MPI presets. The PM job also
runs the reader test with `/usr/bin/python3`, whose NumPy/h5py are installed by
its apt step, so an unrelated CMake-discovered Python cannot silently turn this
CI dependency gate into a skipped qualification.

No force/reference golden values, scientific budgets, particle counts or difficult
geometries were weakened. Remaining unconfirmed staleness needs actual failure
logs, source configurations and independent reference evidence before edits.

## Manual validation commands — DO NOT infer prior execution

Run from the repository root. The MPI preset requires MPI, Parallel HDF5, FFTW
and FFTW-MPI as documented by existing provider checks. Python 3.8 requires
compatible NumPy/h5py. Dependencies were not probed or installed in this campaign.

### Build, registration and complete CPU suite

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug -N
ctest --preset test-cpu-debug --output-on-failure

cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug -N

cmake --preset mpi-hdf5-fftw-release
cmake --build --preset build-mpi-hdf5-fftw-release
ctest --preset test-mpi-hdf5-fftw-release -N
```

### Focused numerical, restart, I/O and Python gates

```bash
ctest --preset test-cpu-debug --output-on-failure -R '^(unit_assertions_enabled|unit_tree_pm_split_kernel|unit_config_parser|unit_memory_governor|integration_hierarchical_time_bins|integration_hierarchical_timestep_regression)$'
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R '^integration_tree_pm_coupling_periodic$'
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R '^validation_tree_pm_(ewald_accuracy|recovery_ewald_accuracy|ewald_accuracy_first_light_softening)$'
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R '^(integration_hierarchical_treepm_workflow|validation_hierarchical_dmo_growth_single_rank|validation_dmo_zeldovich_workflow_single_rank)$'
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R '^(integration_restart_checkpoint_roundtrip|integration_restart_equivalence_harness|integration_restart_equivalence_treepm|integration_provenance_roundtrip|integration_hierarchical_treepm_workflow)$'
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R '^(unit_snapshot_hdf5_schema|integration_snapshot_hdf5_roundtrip|python_treepm_snapshot_tools)$'
python3.8 tests/python/test_treepm_snapshot_tools.py
ctest --preset test-mpi-hdf5-fftw-release --output-on-failure -R '^unit_assertions_enabled$'
```

The always-active Gaussian source can also be built as the actual target
`test_unit_tree_pm_split_kernel`; the commands above build it as part of the
preset. Restart compatibility is exercised inside the real
`integration_hierarchical_treepm_workflow`, not a fictional separate target.
HDF5 fixtures use private temporary workspaces and are generated by their tests.

### Real distributed qualifications and convergence fixtures

```bash
ctest --preset test-mpi-hdf5-fftw-release --output-on-failure -R '^integration_tree_pm_coupling_periodic_mpi_(two|three|four|eight)_rank$'
ctest --preset test-mpi-hdf5-fftw-release --output-on-failure -R '^integration_snapshot_set_mpi_[2348]_rank$'
ctest --preset test-mpi-hdf5-fftw-release --output-on-failure -R '^integration_hierarchical_treepm_workflow(_np[24])?$'
ctest --preset test-mpi-hdf5-fftw-release --output-on-failure -R '^validation_hierarchical_dmo_'
ctest --preset test-mpi-hdf5-fftw-release --output-on-failure -R '^validation_dmo_zeldovich_workflow_(single_rank|mpi_np[2348])$'
ctest --preset test-mpi-hdf5-fftw-release --output-on-failure -R '^integration_reference_workflow_distributed_treepm_mpi_two_rank$'
ctest --preset test-mpi-hdf5-fftw-release --output-on-failure -R '^integration_runtime_app_mpi_treepm_smoke_two_rank$'
./scripts/ci/check_repo_hygiene.sh
```

CTest hierarchical artifact fixtures arrange np1/np2/np4 producers for rank
comparison. Avoid running multiple CTest invocations against the same build's
shared qualification artifact outputs concurrently. A serial CTest pass is not
MPI evidence. The small linear fixture's six fundamental modes do not establish
a nonlinear validated k range.

For larger matched user-run snapshots, supply complete logical sets at matching
final epochs and explicit k bounds in inverse stored coordinate units:

```bash
python3.8 tools/hierarchical_dmo_qualification.py \
  --initial "$DMO_INITIAL_COMPLETE" --global-reference "$DMO_GLOBAL_COMPLETE" \
  --hierarchical "$DMO_HIERARCHICAL_COMPLETE" --refined "$DMO_REFINED_COMPLETE" \
  --mesh 64 --bins 6 --k-min "$QUALIFICATION_K_MIN" --k-max "$QUALIFICATION_K_MAX" \
  --max-power-relative-error 0.02 --max-growth-relative-error 0.01 \
  --output /tmp/chui_hierarchical_qualification.json
```

This requires matched scientific normalized provenance and compatible masses,
IDs and geometry. Keep the global rung-zero reference and halve coarse limits
for refinement. Requested shells must have >=8 modes and lie below quarter
Nyquist. Choosing and validating a physical scale range remains a user campaign.

### Frozen snap_1698 force and full-population P0 comparison

Set these from the saved P0 configuration/provenance, not this report's defaults:
`SNAPSHOT_COMPLETE`, `P0_TARGET_IDS` (the original 4,096 persistent IDs),
`EPS_COMOVING`, `BOX_X/Y/Z`, `ASMTH_CELLS`, `RCUT_CELLS`, `G_CODE`, `SCALE_FACTOR`.
Use the exact epoch, not the rounded redshift. Default export conversions must
produce the saved solver units; if another explicit conversion is needed, supply
`--position-scale` and `--mass-scale` and adjust box/epsilon/G consistently.

```bash
python3.8 tools/treepm_sweep.py export "$SNAPSHOT_COMPLETE" \
  --output /tmp/chui_snap_1698_particles.csv \
  --target-ids /tmp/chui_snap_1698_ids.txt --fixed-target-ids "$P0_TARGET_IDS" \
  --targets 4096 --epsilon "$EPS_COMOVING"

TREEPM_BENCH="$PWD/build/mpi-hdf5-fftw-release/bench_tree_pm_sweep"
TREEPM_PHYSICAL=(--box-x "$BOX_X" --box-y "$BOX_Y" --box-z "$BOX_Z" \
  --asmth "$ASMTH_CELLS" --rcut "$RCUT_CELLS" --epsilon "$EPS_COMOVING" \
  --g "$G_CODE" --scale-factor "$SCALE_FACTOR")

python3.8 tools/treepm_sweep.py run --exe "$TREEPM_BENCH" \
  --particles /tmp/chui_snap_1698_particles.csv --target-ids /tmp/chui_snap_1698_ids.txt \
  --output /tmp/chui_treepm_1698_repair_reference --meshes 64 --leaves 16 --blocks 64 \
  --threads 1,2,4 --policies strict,adaptive --kernels analytic,lookup \
  --history fallback --accounting full --warmups 1 --repeats 5 --oracle \
  "${TREEPM_PHYSICAL[@]}"
```

The leaf/block values above are explicit experiment choices: replace the leaf
with the saved P0 value if different. Verify output configuration `particles`
is 262,144, all-source hashes, exact target IDs, mesh/split/cutoff/softening,
FFTW flag, workers and history policy. Every timed call covers all targets and
includes PM; the fixed IDs select only force comparisons. Reports contain
wall/PM/tree-build/traversal medians, actual pairs/nodes/opened/multipoles,
retained capacities/RSS, matching-split L2/P95/P99/max, absolute tails and
small-force absolute maxima. Bootstrap history, if requested with a fresh output
directory and `--history bootstrap`, is labeled same-snapshot reconstruction,
not original KDK history.

**Native saved P0 oracle invocation is blocked:** this checkout lacks
`bench/bench_tree_pm_snapshot.cpp`, `tools/chui_treepm_oracle.py` and
`tools/chui_treepm_profile.py`, plus their historical numeric output format.
The command above regenerates an independent long-double matching-split direct
reference for the exact saved IDs using the actual available driver. Recompare
the saved oracle with its original external tool in the environment that contains
it; do not call the regenerated reference a verification of the saved file.
Saved baseline checksums/config hashes must be supplied externally.

### Split-aware mesh/leaf experiments and optional perf

With an independently qualified total-force CSV and identity/evidence JSON as
specified above, run:

```bash
python3.8 tools/treepm_sweep.py run --exe "$TREEPM_BENCH" \
  --particles /tmp/chui_snap_1698_particles.csv --target-ids /tmp/chui_snap_1698_ids.txt \
  --output /tmp/chui_treepm_1698_repair_sweep --meshes 64,96,128,192 \
  --leaves 4,8,16,32 --blocks 32,64,128 --threads 4 \
  --policies strict,adaptive --kernels analytic,lookup --history fallback \
  --accounting fast --warmups 1 --repeats 5 --oracle \
  --total-reference "$CONVERGED_TOTAL_FORCE_CSV" \
  --total-reference-manifest "$CONVERGED_TOTAL_REFERENCE_MANIFEST" \
  "${TREEPM_PHYSICAL[@]}"

OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE perf record -g -o /tmp/chui_treepm_repair_perf.data -- \
  "$TREEPM_BENCH" --particles /tmp/chui_snap_1698_particles.csv --mesh 64 --leaf 16 \
  --policy strict --kernel analytic --history fallback --accounting full \
  --warmups 1 --repeats 5 "${TREEPM_PHYSICAL[@]}"
perf report -i /tmp/chui_treepm_repair_perf.data
```

Use fresh output directories. Review each candidate's own physical split and
matching short reference, total Ewald/converged error and tails alongside timing.
P0's supplied strict/fallback timing medians are historical reference evidence,
not a speedup measured here. Do not extrapolate memory fit or scaling to 512^3.

## Static checks and packaging

`git diff --check` returned exit 0. Python `ast.parse(...,
feature_version=(3,8))` accepted the shared reader, both consumers and the new
Python regression source without importing them or executing tests. Static test
inventory and diff inspection are not build, CTest registration or physics proof.
Repository hygiene is supplied as a manual command and was not run.

The source-only ZIP contains exactly the changed/new files listed below, preserving
repository-relative paths. No build trees, outputs, snapshots, profiler files,
executables or caches are included. The modified working tree is authoritative.

Suggested branch: `fix/treepm-p1p5-qualification-and-test-reliability` (leave the
current P1–P5 draft branch/history intact if repairs stay on that PR).

Suggested PR title: **Repair TreePM P1–P5 Qualification, Snapshot Integrity,
Restart Compatibility, and Test Reliability**.

## Exact changed-file manifest

Generated from the working tree after the source/documentation repairs. Files
marked new did not exist at the start. All other entries are modified.

Total: **38 files**.

### Production C++

- `include/cosmosim/gravity/tree_pm_split_kernel.hpp`
- `src/core/time_integration.cpp`
- `src/gravity/tree_pm_coupling.cpp`

### Restart/compatibility C++

- `include/cosmosim/io/restart_checkpoint.hpp`
- `include/cosmosim/workflows/reference_workflow.hpp`
- `src/io/restart_checkpoint.cpp`
- `src/workflows/output_restart_runtime.cpp`
- `src/workflows/reference_workflow.cpp`

### Python tooling

- `tools/chui_snapshot.py` — new
- `tools/hierarchical_dmo_qualification.py`
- `tools/treepm_sweep.py`

### Scientific tests and support

- `tests/integration/restart_equivalence_harness.hpp`
- `tests/integration/test_hierarchical_time_bins.cpp`
- `tests/integration/test_hierarchical_timestep_regression.cpp`
- `tests/integration/test_hierarchical_treepm_workflow.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/integration/test_restart_equivalence_harness.cpp`
- `tests/integration/test_restart_equivalence_treepm.cpp`
- `tests/integration/test_snapshot_hdf5_roundtrip.cpp`
- `tests/integration/test_snapshot_set_mpi.cpp`
- `tests/integration/test_tree_pm_coupling_periodic.cpp`
- `tests/python/test_treepm_snapshot_tools.py` — new
- `tests/support/test_require.hpp` — new
- `tests/unit/test_config_parser.cpp`
- `tests/unit/test_tree_pm_split_kernel.cpp`
- `tests/validation/test_tree_pm_ewald_accuracy.cpp`

### CMake/CTest and CI

- `.github/workflows/ci.yml`
- `CMakeLists.txt`

### Documentation

- `docs/build_instructions.md`
- `docs/profiling.md`
- `docs/repair/treepm_p1p5_qualification_and_test_reliability_repair_20261010.md` — new
- `docs/repair/treepm_unified_performance_recovery_20261010.md`
- `docs/repair_open_issues.md`
- `docs/repair_state_recap.md`
- `docs/restart_checkpointing.md`
- `docs/tree_pm_coupling.md`
- `docs/validation_plan.md`
- `tools/README.md`
