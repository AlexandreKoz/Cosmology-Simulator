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
