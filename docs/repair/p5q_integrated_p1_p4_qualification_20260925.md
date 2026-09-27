# P5-Q integrated P1–P4 qualification — 2026-09-25

**Session mode:** Repair mode (qualification-driven defect repair only).
**Scope:** integrated qualification of the TreePM scalability repairs merged in
PR #106 (`1dbe004`, baseline `108abfe` = `5e6ca4c^1`), repair of the concrete
defects this qualification exposed, and an acceptance matrix with a readiness
verdict. P5-W, force-model changes, and new science campaigns are out of scope.

## 1. Environment

| Item | Value |
| --- | --- |
| Host | Linux, g++ 13.3.0, CMake 3.28.3, Ninja |
| MPI | OpenMPI 4.1.6 (`/usr/bin/mpiexec`) |
| HDF5 | 1.10.10 parallel (`h5pcc.openmpi`, `H5Pset_fapl_mpio` probe green) |
| FFTW | 3.3.10 + FFTW-MPI |
| OpenMP | 4.5 (GCC) |
| Presets | `cpu-only-debug`, `mpi-hdf5-fftw-debug`, ad-hoc `build/cpu-noomp-debug` |
| Not available | `gh`, network CI inspection, `h5dump`, `h5py` |

## 2. Repairs applied this qualification

| # | Defect | Class | Fix | Evidence |
| --- | --- | --- | --- | --- |
| R1 | Budget-reject assertion in end-to-end fixture only recognized `gravity.` ownership; startup planner rejects surface as `parallel.decomposition.startup_planner` (`src/workflows/migration_balance_runtime.cpp:2326`). | P1-era fixture gap | Assertion accepts both owners. | cpu floor green; noomp floor green. |
| R2 | Red-pressure fixture budget/overlap too small for the compressed streaming admission path. | P1-era fixture gap | Budget `1400000000`, overlap `1331000000` in fixture. | cpu floor green. |
| R3 | Catchup fixture used an unpinned `dt_time_code`, letting the adaptive CFL dominate the catchup step count. | P1-era fixture gap | Pin `catchup_options.dt_time_code = 1.0e-4`. | cpu floor green. |
| R4 | MPI-only compile break: `persistent_limit` declared inside `try` but captured by the post-validation lambda (`src/parallel/distributed_memory.cpp`, region `#if COSMOSIM_ENABLE_MPI` from :4104). Introduced by P1–P4 (`1dbe004`/`db33a75`); baseline kept it outside the `try`. | **P1–P4 regression (compile)** | Hoist `std::uint64_t persistent_limit = 0U;` before the `try`, assign after validation. | `cmake --build --preset build-mpi-hdf5-fftw-debug` → RC 0 (385 targets); without the fix the MPI preset cannot build. |
| R5 | `integration_reference_workflow_distributed_treepm_mpi_two_rank` deadlocked: per-rank output roots diverged under collective single-file publication (`H5Pset_fapl_mpio(..., MPI_COMM_WORLD)` at `src/io/snapshot_hdf5.cpp:1087` with `collective_single_file = single_file_snapshot && world_size > 1`, `src/workflows/output_restart_runtime.cpp:494`), while only rank 0 creates the shared snapshot directory → OpenMPI ompio barrier mismatch inside `PMPI_File_open`. | Pre-existing (baseline hangs identically, RC 124) | Shared-across-ranks roots (`direct`, `first`, `resumed`, `bad_rank_count`) in `tests/integration/test_reference_workflow_distributed_treepm_mpi.cpp`; restarts remain per-rank via `rankNNN`-suffixed run directories. | `mpiexec -n 2` test → RC 0; final full MPI suite: **#9 Passed (0.86 s)**. All P1/P2/P3 asserts inside the test pass at 2 ranks. |
| R6 | `integration_reference_workflow` time-cadence scenario failed at line 734 (`final_time_code - 0.0104`): the fixture's adaptive gravity CFL dt (≈1.49e-4, hardcoded factor 0.2 at `src/workflows/time_coordinator.cpp` particle path; option `dt_time_code` is only a min-cap at `:1061`) sits below the 2e-4 option cap, shifting the event phase so step 4 stopped at 0.0103; the `dt == 0.0002` restart assert was unreachable (checkpoint dt is the latched unclipped resume dt from the last clipping output boundary, `:1085-1086` + `:1203-1205`). Pre-existing (baseline fails the same assert line). | Pre-existing fixture drift | `max_global_steps` 4→5 (endpoint clip reaches 0.0104 on step 5); restart-dt assert now anchors on the last `time.output_event_clip` event's `unclipped_dt_time_code` (real resume-dt contract). | #5 **Passed** in cpu floor, noomp floor, and final MPI suite. |
| R7 | OpenMP-off qualification: end-to-end P3 thread matrix requested `omp_threads ∈ {1,2,4}`; non-OpenMP builds fail closed on `>1` (`src/core/config.cpp:1542`). | New (P3 fixture vs OFF build) | Thread sizes `#ifdef _OPENMP {1,2,4} #else {1}`. | `build/cpu-noomp-debug` floor → 150/150. |
| R8 | OpenMP-off qualification: `integration_config_examples` loads `configs/chui_firstlight_32_smoke.param.txt` (`omp_threads = 4`) as-is → `ConfigError` on non-OpenMP builds. | Pre-existing test/config tension | Under `#else`: assert the fail-closed rejection message, then reload with a serial team to keep every other firstlight invariant under test. | noomp floor → 150/150. |
| R9 | `chui` launcher documented as `./chui run ...` (docs/build_instructions.md) but tracked mode `100644`; `integration_chui_launcher_mpi_runtime_two_rank` failed with `Permission denied`. Pre-existing (baseline script fails identically). | Pre-existing repo defect | Mode `100644 → 100755` recorded via `git update-index --chmod=+x chui` (working tree has `core.fileMode=false`, so a plain `chmod` alone would not persist in git). | #281 rerun → **RC 0**. |

No physics, numerics, force-model, or schema behavior was changed; all repairs
are fixture/test/build-interface scoped except R4 (compile-correctness fix in
`src/parallel/distributed_memory.cpp`).

## 3. Acceptance matrix

### 3.1 P1–P4 deferred validation gates (from `docs/repair_open_issues.md` rows)

| Gate | Command | Outcome |
| --- | --- | --- |
| G1 CPU floor | `cmake --preset cpu-only-debug` → `cmake --build --preset build-cpu-debug` → `ctest --preset test-cpu-debug -j 4` | **150/150 passed** (re-run after all edits, RC 0). |
| G2 OpenMP OFF floor | `cmake -S . -B build/cpu-noomp-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCOSMOSIM_ENABLE_TESTS=ON -DCOSMOSIM_ENABLE_BENCHMARKS=ON -DCOSMOSIM_ENABLE_MPI=OFF -DCOSMOSIM_ENABLE_HDF5=OFF -DCOSMOSIM_ENABLE_FFTW=OFF -DCOSMOSIM_ENABLE_CUDA=OFF -DCOSMOSIM_ENABLE_PYTHON=OFF -DCOSMOSIM_ENABLE_LTO=OFF -DCOSMOSIM_ENABLE_OPENMP=OFF` → build → `ctest -j4` | Build RC 0 (453 targets); **150/150 passed** after R7/R8. `feature_openmp=false` confirmed in configure log and `CMakeCache.txt`. |
| G3 OpenMP worker telemetry (P3) | G1/G2 runs of `test_integration_reference_workflow` P3 block: `openmp_configured_workers == requested`, `openmp_observed_workers == configured`, `total == local + remote`, identical pair totals across thread sets | Passed (thread matrix `{1,2,4}` on OpenMP builds, `{1}` on non-OpenMP builds). |
| G4 P1 observability identity | Same block: 5 `gravity.treepm_let` events per 2-step run, event-by-event pair identity | Passed (cpu + noomp floors). |
| G5 P1/P2/P3 in distributed context | `ctest -R integration_reference_workflow_distributed_treepm_mpi_two_rank` (2 ranks, Parallel HDF5) | **Passed (RC 0)**: geometry freshness/fallback fields, pair identity, and startup-planner budget-reject (R1) all green. |
| G6 compact dense-index contract | `test_parallel_distributed_memory` — duplicate/missing/out-of-range `local_index` rejection with byte-exact messages across 256 vectors × 24 permutations | Passed (registered in G1 floor). |
| G7 startup transient admission | `test_parallel_distributed_memory` — transient byte model, overflow rejection, governor reject/accept | Passed (registered in G1 floor). |
| G8 SFC rebalance 2-rank | `ctest -R integration_distributed_sfc_rebalance_mpi_two_rank` | **Passed** in final MPI suite (3-/4-rank variants remain pre-existing failures, §3.3). |
| G9 MPI-only compile | `cmake --preset mpi-hdf5-fftw-debug` + build | **RC 0** after R4; configure reports `feature_hdf5_parallel=true`, MPIO probe, `fftw_mpi`, OpenMP 4.5. |
| G10 time-cadence + catchup integration scenario | `ctest -R integration_reference_workflow` in MPI/HDF5 build | **Passed** after R6 (and R3). |

### 3.2 MPI/HDF5/FFTW suite (integration reality)

`ctest --preset test-mpi-hdf5-fftw-debug -j4` → **231/308 passed (77 failed)**
before the R9 chmod; #281 passes on rerun after R9, so **76 residual failures**,
all classified below. Full logs: `build/mpi_ctest_final.log` (and
`build/mpi_ctest_full.log` for the pre-repair baseline of this session).

### 3.3 Residual failure inventory (all baseline-reproducible)

Every cluster below was re-executed against baseline worktree `108abfe`
(`build/baseline-wt/bld`, same preset/flags) and failed with the same message,
or (where noted) fails for reasons outside the P1–P4 diff (the P1–P4 diff is
25 files: gravity/parallel/workflows + docs + unit test only).

| Tests | HEAD symptom | Baseline | Class |
| --- | --- | --- | --- |
| #11 isolated_pm 2-rank | `report.completed_steps == 2` | fails, `tree pseudo hierarchy exchange returned mixed epochs or geometry frames` | Pre-existing failure; **HEAD symptom change is a risk hypothesis, not a confirmed regression** — needs a dedicated follow-up. |
| #12 hydro 2-rank | `duplicate node identity` | identical message | Pre-existing |
| #67–70 gravity_gas np2/np3 | `mpi authoritative gas gravity contract: GasCellIdentityMap ... diverged` | identical | Pre-existing (check in unchanged `src/core/simulation_state.cpp`) |
| #71–76 star_formation mpi + ism_eos | `ConfigError: periodic TreePM rectangular-grid ... split_scale=0.000125, coarsest_spacing=0.000125` | identical (2 args sampled) | Pre-existing (check in unchanged `src/core/config.cpp`) |
| #133 pm_periodic 2-rank | `Distributed PM solve did not report routed density records` | identical | Pre-existing |
| #142/143 sfc_rebalance 3/4-rank | `checkHardMemoryCollectivePreparation` assert (test line 341) | identical | Pre-existing (2-rank green = G8) |
| #151–190 ic_reader 1/2/3/4/8-rank matrix (40) | asserts at test lines 1082/1093 | identical (1-rank + 4-rank sampled) | Pre-existing |
| #230 gas_cell_migration 2-rank | `migration transport round limit is too small for one fragment header per rank` | identical | Pre-existing (check in unchanged `src/workflows/internal/migration_wire.cpp`) |
| #231 amr_patch_migration 2-rank | `AMR patch migration wire scheduler record count does not match gas-cell count` | identical | Pre-existing |
| #240–244 tree_pm_coupling_periodic serial + 2/3/4/8-rank | `TreePM consistency failure: split_rel=0.576497, pm_only_rel=0.566082, split_pruned_nodes=153, split_pair_skips=381` | **bit-identical numbers** | Pre-existing (passes only in FFTW-off cpu preset → FFTW-era drift) |
| #265 installed_package_consumer | `Installed CosmoSim package requires Parallel HDF5, but the consumer-resolved HDF5/MPI dependency family cannot compile/link H5Pset_fapl_mpio` | identical | Pre-existing (environment HDF5/MPI family resolution in probe) |
| #266 configure_presets | `serial-HDF5 rejection failed for an unexpected reason: ... H5Pset_fapl_mpio ... cannot compile/link` | identical | Pre-existing (same environment class) |
| #276/277 snapshot smoke 2/3-rank | `rerunning a committed single-file snapshot namespace unexpectedly succeeded` | identical (2-rank script run) | Pre-existing |
| #278 restart smoke | missing `restart_001_rank000.hdf5` | identical | Pre-existing |
| #281 launcher | `Permission denied` | identical | Pre-existing → **repaired by R9 (now green)** |
| #294/295 phase2 gravity 4/8-rank | `tree pseudo hierarchy exchange returned duplicate node identity` | identical (4-rank sampled) | Pre-existing |
| #301 dmo single-rank | `uninterrupted/resumed deterministic state digest mismatch; passing_ranks=0/1` | identical | Pre-existing (scientific validation failure, pre-dates P1–P4) |
| #304–307 dmo np2/3/4/8-rank | **timeout 300 s with zero output** (hang); reproduced manually: RC 124, empty log | np2 identical at baseline (RC 124, empty log) | Pre-existing hang |
| #308 rank_equivalence | Not Run (depends on #304) | cascade | Pre-existing cascade |
| #284 validation_convergence | Timeout at `-j4` | — | **Environment contention artifact**: passes standalone (`./build/mpi-hdf5-fftw-debug/test_validation_convergence` → RC 0), repeatedly. |

**P1–P4 regressions found in the MPI suite: none.** The only confirmed
P1–P4-induced defect was the R4 compile break (fixed). All other failures
reproduce on the untouched baseline in this environment.

## 4. Readiness verdict

- **P1–P4 deferred validation gates: GREEN** (G1–G10), with the R4 compile
  fix as the only source change required.
- **Integrated MPI/HDF5/FFTW suite: NOT release-ready on this host** — 76
  residual failures are pre-existing (baseline-reproducible), dominated by
  distributed IC reading (40), TreePM coupling consistency (5), star-formation
  MPI configs (6), validation hangs (7), and environment-probe tests (3).
  They are **not attributable to P1–P4** but block any "full suite green"
  claim.
- **CI correspondence is unverifiable here**: `.github/workflows/ci.yml`
  runs a restricted MPI regex (line 116) that includes several tests that
  fail on this host at baseline (`#11 #12 #133 #142 #230 #241 #294 …`);
  whether upstream CI is green cannot be checked without `gh`/network
  (named blocker).
- **Risk hypotheses (not confirmed defects):** #11 symptom change
  (completed_steps vs mixed-epochs); #276 "unexpectedly succeeded" rerun
  semantics may indicate a namespace-collision guard that no longer fires —
  both failed at baseline too, so they need dedicated follow-up rather than
  P5 repair.
- **No schema, config-key, restart, or force-model changes** were made;
  reproducibility impact is limited to fixture/failure-surface behavior.

## 5. Follow-ups (recommended, out of P5-Q scope)

1. Dedicated campaign for the distributed-IC reader assert family
   (`tests/integration/test_distributed_ic_reader_mpi.cpp:1082/1093`).
2. TreePM split-vs-pm-only consistency drift under FFTW
   (`split_rel 0.576497` vs `0.566082`) — numerics investigation.
3. DMO validation hang (zero output at np≥2, digest mismatch at np1) —
   scientific validation blocker, predates P1–P4.
4. Environment probes #265/#266 — container/host HDF5-MPI family coherence.
5. Baseline-vs-HEAD symptom deltas (#11) — reproduce on a dependency-complete
   CI host.
