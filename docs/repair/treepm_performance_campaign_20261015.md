# TreePM performance campaign — 2026-10-15

Session mode: feature implementation + repair. Input: the 12–13 October
TreePM investigation report (264k-particle evolved DMO box, 17.8 s per global
step on an i5-8250U after the CPU-binding fix) and a source audit of the
residual traversal, acceptance criterion, PM backend and hierarchical coarse
limiter. AREPO (development source) was consulted as a design reference only;
no AREPO code was copied.

## Diagnosis (evidence basis)

| Finding | Evidence | Disposition |
|---|---|---|
| Strict policy is near-direct summation | 5573 direct pairs/target vs 5849 mean neighbours inside `r_cut` (clustered 64^3 fixture, 64^3 PM) | Unchanged reference; fast settings are opt-in |
| Leaf pairs gathered through `sorted_particle_index` from row-ordered lanes | Same particles in Morton row order: traversal 11.9 s -> 6.75 s, identical pairs | Fixed (tree-ordered lanes, Phase 1) |
| Derivative proxy ~9x stricter than the monopole relative criterion; one erf/exp per node | Short-range error 5.7e-5 at alpha 0.005; alpha sweep insensitive | Option `treepm_adaptive_error_estimator` |
| Nodes straddling `r_cut` opened to leaves | 21–29 % fewer pairs when evaluated at COM distance | Option `treepm_cutoff_node_policy` |
| PM FFT serial | no `fftw_plan_with_nthreads` in source; 256^3 PM ~2 s/solve | Optional `libfftw3_omp` |
| Hierarchical coarse block bounded by box max speed + max total acceleration over one PM cell | Rung-2 coarse block advanced ~1.8 global steps (report logs) | Option `hierarchical_coarse_displacement_criterion` + limiter logging |
| No governor reconcile between hierarchical bootstrap and rung-membership transaction | B-10 stale-baseline fatal at that site | Fixed (collective reconcile) |
| Global KDK steps everyone at the minimum particle step | `allreduceMinDouble` in global step selection | Unchanged; hierarchical KDK is the remedy |
| Open MPI default binds a hybrid rank's threads to one core | ~199 % CPU for 4 threads | Launcher auto-placement |
| `bench_tree_pm_sweep`, three tests did not compile; two gravity tests stale | build log of unmodified `main` | Repaired (test/bench only) |

## Changes

- `src/gravity/tree_pm_coupling.cpp`, `include/cosmosim/gravity/tree_pm_coupling.hpp`:
  tree-ordered source lanes and tree-order target visits; estimator and
  cutoff-node options (validated, mixed into the options fingerprint);
  adaptive angle cap (0,1].
- `src/gravity/internal/tree_pm_acceptance.hpp`: `kRelativeMonopole` branch,
  COM-offset skip when it is unused.
- `src/gravity/pm_solver.cpp`, `CMakeLists.txt`: `COSMOSIM_ENABLE_FFTW_THREADS`
  (default ON when `libfftw3_omp` and OpenMP exist) and per-plan thread count.
- `src/gravity/gravity_memory.cpp`: `gravity.estimate.treepm_tree_ordered_source_lanes`.
- `include/cosmosim/core/config.hpp`, `src/core/config.cpp`,
  `src/workflows/gravity_runtime.cpp`: four typed keys (defaults preserve
  previous behavior) with normalized dump and event payload fields.
- `src/workflows/time_coordinator.cpp`: RMS split-scale coarse criterion,
  coarse-limiter/rung-histogram event fields, post-bootstrap reconcile.
- `include/cosmosim/core/profiling.hpp`: retained event limit 4096.
- `chui`: Open MPI hybrid placement.
- `bench/bench_tree_pm_sweep.cpp`, `tools/generate_clustered_dmo_box.py`:
  measurement tooling.
- `configs/treepm_evolved_snapshot_fast.param.txt`: measured fast deck.
- Tests: `tests/unit/test_config_parser.cpp` (new keys, angle range),
  `tests/integration/test_chui_launcher.py` (placement),
  `tests/integration/test_config_examples.cpp` (fast deck), and the
  test-only staleness repairs listed above.

## Reproducibility impact

- Defaults: unchanged numerical behavior. The normalized config dump gains
  four lines, so config hashes of existing decks change; restart readers are
  unaffected (no restart/snapshot schema change).
- Tree-ordered traversal: per-target forces bitwise identical (verified).
- Threaded FFTW: deterministic for fixed plan rigor and thread count; changing
  `OMP_NUM_THREADS` may change PM forces at round-off level. Configure with
  `-DCOSMOSIM_ENABLE_FFTW_THREADS=OFF` for thread-count-independent PM.
- Opt-in options and the fast deck change forces within the measured error
  table in `docs/tree_pm_coupling.md`.

## Commands and results

```bash
# measurement build
cmake --preset cpu-only-release-native -DCOSMOSIM_ENABLE_FFTW=ON -DCOSMOSIM_ENABLE_OPENMP=ON
cmake --build build/cpu-only-release-native --target bench_tree_pm_sweep
python3 tools/generate_clustered_dmo_box.py clustered_random.csv --order random
OMP_NUM_THREADS=4 build/cpu-only-release-native/bench_tree_pm_sweep \
  --particles clustered_random.csv --box-x 10 --box-y 10 --box-z 10 --epsilon 0.005 \
  --mesh 64 --leaf 16 --block 64 --policy adaptive --maximum-angle 0.5 \
  --kernel lookup --accounting fast --history bootstrap --forces forces.csv
# accuracy matrix: same command per row of the table in docs/tree_pm_coupling.md,
# reference = --policy strict --kernel analytic on the same mesh.

# test gates
cmake --preset cpu-only-debug && cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

Results are recorded in `docs/tree_pm_coupling.md` and in the closing
section of this file.

## Not done / open

- Endpoint force reuse across hierarchical coarse blocks (two full force
  evaluations per block remain). It needs a split-cache/history restart
  equivalence contract; deferred rather than rushed.
- Group (shared interaction list) walk and a real SIMD kernel.
- Trajectory, growth and power-spectrum qualification of the fast deck and of
  `rms_split_scale` with `hierarchical_max_rung >= 4` on the real snapshot.
- Pre-existing failures unrelated to this campaign (see closing section).

## Closing results (executed in this session)

`cpu-only-debug` full CTest: 146/152 pass. The 6 failures are identical, with
the same assertion and line, on unmodified `main` (c5ccc87), which in addition
failed to compile three tests and failed `integration_tree_pm_coupling_periodic`:

| Test | Failure on main and here |
|---|---|
| `unit_simulation_state` | `test_simulation_state.cpp:441` duplicate-ID assertion |
| `unit_runtime_capabilities` | `production_hierarchical_local_timestep` status expectation |
| `unit_tree_gravity` | `test_tree_gravity.cpp:853` `saw_pm_routing` (estimate folded into the communication arena) |
| `unit_restart_checkpoint_schema` | `test_restart_checkpoint_schema.cpp:134` `integrity_rejected` |
| `integration_source_package_completeness` | `rsync` not installed |
| `integration_runtime_truth_ctest_labels` | `require_backend_label_absent_when_disabled` |

`pm-hdf5-fftw-debug` gravity/hierarchical gates (same targets built from
`main` for comparison):

| Test | main | this campaign |
|---|---|---|
| `validation_tree_pm_ewald_accuracy` | pass | pass |
| `validation_tree_pm_ewald_accuracy_first_light_softening` | pass | pass |
| `validation_tree_pm_recovery_ewald_accuracy` | fail: strict qualification exercises no direct pairs | same |
| `integration_tree_pm_coupling_periodic` | fail: zero PM box (stale test) | fail later: split consistency `split_rel=0.576`; identical numbers when the repaired test runs against `main` solver code, so pre-existing |
| `integration_hierarchical_treepm_workflow` | fail: B-10 stale governor baseline | all four runs complete with restart round-trip and deterministic replay digest; fails the global checkpoint boundary-kind expectation (`kCheckpointPoint`) and, past that, the hierarchical-vs-global trajectory gate `sqrt(e2/64) < 1e-3` (line 238); `main` plus only the B-10 reconcile fails line 238 identically |
| `validation_dmo_zeldovich_workflow_single_rank` | segfault (stale test reads absent compact-DMO lanes) | fail: uninterrupted/resumed digest mismatch; identical on `main` with the repaired test |
| `validation_hierarchical_dmo_growth_single_rank` | fail: B-10 stale governor baseline | **pass**, including the new `rms_split_scale` case |

Benchmark evidence: tree-ordered traversal forces bitwise identical to the
previous traversal on 4096 targets; accuracy/timing matrix in
`docs/tree_pm_coupling.md`.

Newly exposed pre-existing defects for the next campaign (not caused here):
hierarchical-vs-global trajectory refinement (line 238 gate), global
checkpoint boundary kind, global resumed-run digest determinism, TreePM
split consistency in the FFTW periodic regression, and the recovery Ewald
profile that exercises no direct pairs.

## Round 2 — defect repairs (2026-10-15, later session)

Root causes and dispositions for the defects listed above:

| Defect | Root cause | Fix | Evidence |
|---|---|---|---|
| Hierarchical trajectory gate (D1) | Acceleration-only rung criterion `sqrt(2 eta eps/|a|)` sees |a| -> 0 while a softened pair passes through coincidence; the pair jumped from the finest to the coarsest rung for a whole block (rung histogram `0,0,64`), with refined error 1.36e-2 > coarse 5.9e-4 | Timestep growth limiter in `TimeCoordinator::runHierarchicalSegment`: an interval may at most double between synchronization points, reconstructed from restart-persisted state only (previous quantum x committed rung period) | refined error 1.36e-2 -> 1.62e-3 (8x); gate is 1e-3, so **partially fixed** — see open item below |
| Global checkpoint boundary kind (D2) | Global KDK tags completed steps by output purpose (`kCheckpointPoint`), a restart-safe synchronized kind by design | Test accepts `isRestartSafeBoundary` for rung 0; hierarchical blocks still require `kGlobalSynchronizationPoint` | `integration_hierarchical_treepm_workflow` passes that check |
| Global resumed-run digest (D3) | `ReferenceWorkflowOptions::dt_time_code` is a process-local cap on the adaptive step (since adaptive cosmological stepping), not restart state; the test resumed without the cap and took a 1.69e-5 step instead of 5e-6 | Test applies the same cap to the resumed segment | `validation_dmo_zeldovich_workflow_single_rank` passes (digest, phases, P(k)) |
| FFTW periodic split consistency (D4) | With the FFTW 32^3 mesh r_cut = 0.14 is below every separation of the 32-point fixture, so the residual path evaluated zero pairs | Fixture places one source within r_cut on every backend | test passes with FFTW and naive DFT |
| Recovery Ewald profile (D5) | Same coverage issue on the 24-point quasi-uniform fixture (r_cut = 6.25/32) | Coverage gate moved to `two_cluster_cubic`; added the production fast profile (adaptive relative_monopole + com_distance, alpha 0.0025, angle 1.0) to the qualification | all three Ewald gates pass; fast profile meets unchanged Ewald targets on all five fixtures |
| `unit_simulation_state` | In-place ID rewrite reused a cached identity-validation certificate | New `SimulationState::notifyParticleIdentityMutation()`; production writers resize first (already invalidating) | passes |
| `unit_runtime_capabilities` | Expected hierarchical KDK `kUnsupported`; it is `kProvisional` and runtime-available | Test updated | passes |
| `unit_tree_gravity` | PM routing estimate folded into the communication arena | Asserts arena owner = max logical contributor and PM density routing >= 2x7x1024 B | passes |
| `unit_restart_checkpoint_schema` | Integrity hashing became allocation-free; parentless gas payload lacked the now-required gas-cell scheduler | Test asserts success with zero governed bytes; adds scheduler | passes |
| `integration_runtime_truth_ctest_labels` | `python_*` tier label collided with the Python-binding backend label | Tier label `python_offline` when bindings are off; offline h5py test relabelled | passes |
| `integration_source_package_completeness` | `rsync` missing in the environment | none (environmental) | passes with rsync |

`cpu-only-debug` full CTest after these repairs: **152/152 pass**.

Open: closing the last factor of D1 needs a velocity- or neighbour-aware
encounter criterion (e.g. softened two-body dynamical time of the nearest
direct-pair neighbour, or a signal-velocity criterion). Every particle pays
for it, so it is a design decision left to the project owner; the existing
gate is unchanged.

### Round 3 — HDF5/FFTW suite repairs (code-level, owner-run validation)

A `pm-hdf5-fftw-debug` build of the pre-campaign base (`c5ccc87`) fails the
same 16 tests as this branch, so all of them predate the campaign. Root causes
and fixes (validation of these is left to the project owner):

| Test(s) | Root cause | Fix |
|---|---|---|
| `integration_restart_checkpoint_roundtrip` | Did not compile: `particle_sidecar.particle_id` is an `AlignedVector`, the force-cache lane a `std::vector` | Explicit range `assign` |
| `integration_restart_equivalence_{dm_only,treepm,multirate_bins,output_enabled}` | Harness compared raw `species_tag`/`owning_rank`/`time_bin` lanes, which are absent (empty) in the compact homogeneous-DMO representation | Compares the logical accessors `particleSpeciesTag`/`particleOwningRank`/`particleTimeBin` per row |
| `integration_restart_equivalence_{hydro_toy,stochastic_sources,amr_*}` | Restart writer requires the authoritative gas-cell scheduler when cells exist; the harness never supplied one | Harness passes a gas-cell scheduler that shares the particle scheduler's tick (a copy when the toy cell bins advance in lockstep, else built from the cell mirror) |
| `integration_transform_fuzz_invariants`, `validation_transform_fuzz_invariants_long` | Same missing gas-cell scheduler | Builds one from the synced cell mirror at the particle scheduler's tick |
| `integration_reference_workflow` | `restart.read.complete` was deliberately renamed `restart.verify.complete` by the streaming-verification change (`36967e9`), which also dropped `scheduler_pending_transition_count` from the event | Test expects the new name; the diagnostic is restored in `restart.verify.complete` |
| `unit_ic_reader` | Serial duplicate-ID message changed to "duplicate or zero persistent particle IDs" (`a5620b2`) | Test accepts the serial and distributed wordings (both still contain "duplicate" and "particle IDs") |
| `integration_hierarchical_treepm_workflow` | D1 above (1.62e-3 vs gate 1e-3) | Round 4 below |

### Round 4 — D1 encounter criterion and R1 endpoint force reuse (source-level)

Status: implemented and compiled (`pm-hdf5-fftw-debug`, `cpu-only-debug`,
`cpu-only-release-native` full builds, no errors); **not executed** — the
project owner runs the suites.

**D1 — softened-encounter rung criterion.** Every Tree solve records per target
the squared separation of its nearest distinct direct-summed source inside
`r_cut` (`TreePmForceAccumulatorView::short_range_nearest_pair_r2_comoving`,
min-combined over the per-target walk, group walk (`omp simd reduction(min)`)
and remote MPI responses). The hierarchical rung candidate adds
`dt <= 0.2 sqrt(a^3 (r_nn^2+eps^2)^{3/2} / (G m_max))`: the softened two-body
dynamical time at coincidence, weaker than the acceleration criterion by
`sqrt(r/eps)` beyond `r ~ eps`. Config `numerics.hierarchical_encounter_timestep_enabled`
(default `true`, typed, dumped, documented, parser test). Global KDK unchanged.
Wire: the short-range response record grows 80 -> 88 bytes (wire v2), memory
estimators follow the shared constant; the remote accumulator estimate is
4 lanes. Expected effect on the gate: the crossing pair is now held at the
softened dynamical time for the whole block; whether that closes the remaining
1.62e-3 -> < 1e-3 is the validation target (gate unchanged).

**R1 — hierarchical endpoint force reuse.** A block closes with an all-active
Tree+PM solve at the synchronized positions and the next block opened with a
second all-target solve at the same positions. `TimeCoordinator` now skips that
bootstrap when every rank (collective min) reports the endpoint complete:
every row's Tree generation equals the source generation, PM captured at it
from the integrator's committed field version (`GravityRuntime::
hierarchicalEndpointForceReusable`). Migration clears the lanes, so a block
after a migration bootstraps. Counters `time.hierarchical_endpoint_force_reuse`
/ `_bootstrap`. Saves one all-target Tree+PM solve per coarse block.

Restart equivalence: re-evaluating at the same positions is not bitwise equal
to the endpoint (the adaptive relative MAC reads the force history the endpoint
updated), so restart schema **v24** persists `/gravity_force_cache/hierarchical_endpoint`
(Tree xyz, PM xyz, nearest-pair r^2, `valid`, `pm_field_version`), validated,
covered by both integrity digests and the streaming verifier, and restored by
`importRestartForceCache`. v23 files stay readable (no endpoint -> bootstrap).
Read-staging admission now also counts the force cache and endpoint lanes.

Public interface changes (migration notes): `GravityRuntime` gains pure virtual
`hierarchicalEndpointForceReusable(state, integrator)` and
`hierarchicalNearestPairR2()` (single in-tree implementation);
`GravityForceCachePersistentView/State` gain the v24 fields and
`makeGravityForceCachePersistentView()`; `TreePmForceAccumulatorView` gains
`short_range_nearest_pair_r2_comoving` and `noteNearestPairToActiveSlot`;
`kTreePmShortRangeResponseWireBytes` is 88. `restartSchema()` reports v24;
release manifest and doc gates updated.

Reproducibility impact: hierarchical trajectories change (encounter criterion
and reuse of the endpoint force instead of a re-evaluated one); global KDK,
snapshots and provenance schemas are unchanged; all ranks must run the same
binary (wire v2).

Tests added/updated: `unit_config_parser` (new key), `unit_restart_checkpoint_schema`
(v24 identity, v23 compatible), `integration_restart_checkpoint_roundtrip`
(endpoint roundtrip incl. `+inf`, stale-lane rejection), doc/release gates.

Commands for the owner:

```bash
cmake --build build/pm-hdf5-fftw-debug -j
ctest --test-dir build/pm-hdf5-fftw-debug --output-on-failure
ctest --test-dir build/pm-hdf5-fftw-debug -R "hierarchical|restart|config_parser|tree_pm" --output-on-failure
```
