# Round 4 — final pre-validation engineering completion (2026-10-26)

**Milestone claimed: ENGINEERING CODE-COMPLETE (source level) for the TreePM / MPI /
hierarchical / multirate / rollback / memory campaign.**
**Status of every item: IMPLEMENTED — NOT COMPILED — NOT TESTED — OWNER VALIDATION REQUIRED.**
**Not claimed: successful compilation, passing tests, or scientific validation.**

Owner instruction for this round: no build, configure, test, MPI, Python,
benchmark, lint or CI. The source and the diffs were inspected statically.
Static inspection cannot guarantee that the code compiles.

Earlier records:
* [round 3](multirate_stabilization_round3_20261026.md)
* [AREPO-grade campaign](arepo_grade_campaign_20261026.md)
* [TreePM optimization campaign](treepm_optimization_campaign_20261025.md)

## A. Implementation completion matrix

| Item | Files | Change | Production algorithm modified | Test source | Closed at source level | Numerical uncertainty |
| --- | --- | --- | --- | --- | --- | --- |
| P0 rollback restoration/invalidation | `src/workflows/time_coordinator.cpp`, `src/workflows/gravity_runtime.cpp`, `include/cosmosim/workflows/gravity_runtime.hpp`, `src/gravity/tree_pm_coupling.cpp`, `include/cosmosim/gravity/tree_pm_coupling.hpp` | Explicit transaction API (`begin`/`commit`/`rollbackCycleTransaction`). The rollback drops every source-history-dependent derived state: full tree rebuild via `invalidateSourceHistoryForRollback`, ghost layer, borrowed-source certificates. It re-imports the endpoint at the advanced generation. Fingerprint guard; collective failure on any rank-local error. | Yes (multirate only) | Workflow (3), (6), (8), (9), (10); rank-local failure injection | Yes | Bitwise replay with reference tree settings; roundoff with reuse/refit (documented policy) |
| P0 transaction boundaries / MPI agreement | `time_coordinator.cpp` | Admission by `allreduce(min)`. Begin, rollback and reconcile go through `FailureCoordinator`. Accept and reject are decided from globally reduced criteria. Segment start/end check that no transaction is open. | Yes | Workflow (7), (8) failure injection, np2/np4 | Yes | — |
| P1 memory-pressure fallback | `time_coordinator.cpp`, `include/cosmosim/core/config.hpp`, `src/core/config.cpp`, `docs/configuration.md` | Segment-lifetime governed store: capacity reused and the reservation grown, never duplicated. `bad_alloc` counts as a refusal. A refused rank releases its storage. New key `numerics.treepm_multirate_rollback_memory_policy` (`single_block_fallback` \| `fail`). Fallbacks are recorded as cycle limiter `rollback_memory`. | Yes (multirate only) | Workflow (7): all/one/several ranks, allocation failure, first-only, `fail` policy; config parser | Yes | Under memory pressure the schedule depends on the budget (explicit, recorded; `fail` forbids it) |
| P2 multirate numerical contract | `include/cosmosim/core/time_integration.hpp`, `include/cosmosim/core/config.hpp`, `docs/time_integration.md` | Stale "equal-duration M blocks" comments corrected to the variable-tiling contract. Fallback block legality proven (`coarse_dt < T ≤ Δt_disp`). | Comments only | Driver: cosmological kick coefficients (exact momentum), cosmological oscillator | Yes | Order unmeasured until the owner runs it |
| P2 soft-stop accounting | `time_coordinator.cpp` | `multirate_soft_stop_extension_blocks` is now derived from accepted steps. Previously rejected attempts inflated it. | Diagnostic only | Workflow `check_blocks` | Yes | — |
| P3 transactional diagnostics | `gravity_runtime.cpp`, `src/workflows/reference_workflow.cpp`, `include/cosmosim/workflows/reference_workflow.hpp` | Bounded journal (256) per cycle transaction: committed into the accepted window, or moved to execution cost on rollback. Accepted history can no longer be evicted by a rejected attempt. `commit` is `noexcept` (both windows reserved at `begin`). | Diagnostic only | Workflow (8): accepted cadence records equal the uninterrupted run elementwise | Yes | — |
| P4 tests | `tests/integration/test_hierarchical_multirate_convergence.cpp`, `tests/integration/test_hierarchical_treepm_workflow.cpp`, `tests/unit/test_config_parser.cpp`, `CMakeLists.txt` | Cosmological driver suites. Workflow sections (7)–(10), including the failure injection. Skewed/empty-rank fixture. Config key test. `COSMOSIM_CTEST_HIERARCHICAL_TIMEOUT_SECONDS` (1800 s) deadlock bound for the workflow tests. | — | Yes | Yes | Order gates `[1.8, 2.2]` are targets, not results |
| P5 earlier "Not implemented" items | docs only | Classified in §E. The multirate analytic two-body gate and the MPI decomposition fixture listed as missing in the AREPO-grade record now exist. | — | — | Yes (classified) | — |
| P6 static audit | all of the above | See §G. | — | — | Yes | — |

## B. State-ownership and rollback contract (inventory)

Classes:
* **A**: authoritative, restored.
* **I**: derived, invalidated.
* **R**: derived, safely reusable (validity condition given).
* **C**: execution-cost diagnostic (counts rejected work).
* **T**: accepted-trajectory diagnostic, rolled back.

| State element (mutated during a PM cycle by …) | Class | Rule |
| --- | --- | --- |
| Positions, velocities `u = p/a` (drift, tree kicks, opening PM kick) | A | Lane copy at opening, restored in place |
| SFC key lane (tree builds may write it) | A | Copied/restored |
| Common drift epoch (time, a) (every fine tick) | A | Restored via `updateAllParticleDriftEpoch` |
| `StateMetadata` (`step_index`, `scale_factor`) | A | Copied/restored |
| Gravity source generation (every fine tick) | A* | Advanced strictly beyond every value the attempt used (monotone; never restored backwards) |
| Identity, mass, species/flags/rank, softening and mask, particle index generation, metadata representation, cell/gas containers | frozen | Fingerprint at opening; a mismatch at rollback fails closed |
| `IntegratorState` (time, a, step, dt, PM sync state, previous quantum, boundary flags, `pm_source_generation`) | A | Copied/restored; `pm_source_generation` re-stamped |
| Particle and cell scheduler (bins, ticks, activation, candidates) | A | Persistent state exported/imported, then the time-bin mirrors are resynchronized |
| Active particle sets | I | Rebuilt by `beginSubstep` from the restored scheduler |
| Pending checkpoint/snapshot/output boundary | A | Copied/restored (interior requests stay latched) |
| Force cache, per-row validity, Tree/PM endpoint split, nearest-pair lane, PM field version | A (endpoint) | Re-imported from the opening endpoint at the new generation |
| Relative-MAC previous-acceleration history | I | Recomputed every solve from the re-imported Tree+PM lanes (generation-checked) |
| Built tree, identical-source reuse candidate, refit candidate, tree moments | I | `invalidateSourceHistoryForRollback` forces a full build |
| Tree-order / leaf-of-position maps | I | Topology-generation keys reset (rebuilt after the next full build) |
| Tree-ordered source lanes, target visit order, group-walk worker scratch | R | Regathered every residual evaluation |
| Source prediction memo | R | Recomputed per call; not used by hierarchical KDK (predictor fails closed) |
| Borrowed-source eligibility certificates | I | Cleared explicitly (also generation-keyed) |
| Particle ghost layer | I | Invalidated explicitly; it is refreshed every solve anyway |
| LET domain leaves and hierarchy | R | Keyed by a content fingerprint of the geometry |
| Sparse peer graph (communicator) | R | Communication topology only; payloads exchanged every solve; ownership frozen |
| LET payloads, remote request/response buffers, communication arena | R | Phase-local, rebuilt every solve |
| PM routing plans | R | Validated per use (plan validity keys from the AREPO-grade round) |
| PM mesh / long-range field validity | R | No PM refresh inside a cycle, so the mesh describes the restored opening sources; the closing surface always refreshes |
| Uniform-mass certificate | R | Mass frozen and fingerprinted; keyed by source-mass generation and pointer |
| Domain geometry (top-domain leaves) | I | Reinstalled (`install_geometry`) at the new generation |
| Decomposition epoch | frozen | No rebalance or migration inside a cycle |
| MemoryGovernor snapshot reservation | resource | Segment-lifetime, grown not duplicated; released at segment end, on failure (RAII) and on local refusal |
| Runtime memory baseline | I | `reconcile_runtime()` after rollback |
| PM refresh/reuse counts, TreePM cadence records | T | Journal per cycle; moved to `multirate_rejected_*` on rollback |
| `hierarchical_block_records`, `completed_steps`, soft-stop extension | T | Staged per cycle / restored step index / derived from accepted steps |
| Profiler timings, stage counters, soft-stop counter, `multirate_rejected_*`, `multirate_pm_cycle_rollbacks` | C | Never rolled back |

**Retry semantics.**
* A genuine rejection retries with `T ← min(T/2, Δt_disp(now))`.
* The test seam retries with the same `T`, up to the per-cycle count.
* A single-block cycle cannot be rejected. The retry bound is 64; above it the
  run fails closed.
* Each retry re-captures into the same store, reusing its capacity.

## C. Numerical integration contract

* **Reference policy (`exact`).** PM at both endpoints of every coarse block.
  It is unchanged by rounds 3 and 4, and the multirate code paths are
  inactive.
* **Experimental multirate policy.** The composition is
  `K_PM(T/2) ∘ B_1 ∘ … ∘ B_m ∘ K_PM(T/2)`, with `m` variable and the blocks
  tiling `T` exactly.
  * `T` is fixed at opening as the minimum of: the target blocks × the opening
    block interval, the PM displacement bound, the next output/end boundary
    (exact clip), and the retry limit. The memory fallback sets
    `T = one block`.
  * Each inner block is re-planned by the reference criteria (quantum ≤
    particle criterion), with the end-rounding nudge.
  * The displacement bound is re-checked at every interior boundary.
  * Long-range kick ownership:
    * the opening factor is `∫_{t0}^{t0+T/2} dt/a / a0`;
    * the closing factor is `∫_{a_mid}^{a_end} … / a_end`;
    * the opening uses the cycle-opening PM components;
    * the closing PM is refreshed at the closing positions;
    * the interior uses the Tree only.

  The driver's epoch checks reject a missing or double PM surface. The
  cosmological driver test makes this an exact-momentum gate.
* **Restart boundaries.** Cycle ends only: soft stop at the step budget, the
  `t_end`/output clip, interior checkpoint/output/analysis deferral, and
  refusal of non-restart-safe resumes. There is no interior-cycle
  serialization, and none is claimed; no stop source can end inside a cycle.
* **Resource fallback.** See B and `docs/time_integration.md` §5.

## D. Validation manifest (owner; in order; NOT run)

```bash
# 1. CPU-only compilation + tests (includes the CPU-only driver suite)
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
ctest --test-dir build/cpu-only-debug -R 'integration_hierarchical_multirate_convergence|integration_hierarchical_timestep_regression|unit_config_parser' --output-on-failure

# 2. FFTW/HDF5 compilation (serial)
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure

# 3. MPI compilation
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug

# 4. Gravity unit tests
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'unit_tree_gravity|unit_tree_pm_split_kernel|unit_gravity_communication_arena|integration_tree_gravity_vs_direct|integration_tree_pm_coupling_periodic|integration_tree_pm_fine_tick_caches' --output-on-failure

# 5. FFT layout regression
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'integration_pm_in_place_fft_layout' --output-on-failure

# 6. MPI routing / deadlock tests (ctest timeouts are the deadlock bound)
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'integration_pm_routing_admission_mpi|integration_tree_pm_remote_group_walk_mpi|integration_mpi_collective_safety|integration_dmo_empty_rank_smoke_mpi|integration_pm_zero_width_slab_mpi' --output-on-failure

# 7. TreePM force accuracy
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'validation_gravity_invariant_suite|validation_tree_pm_ewald_accuracy|validation_tree_pm_recovery_ewald_accuracy|validation_phase2_mpi_gravity_single_rank' --output-on-failure

# 8. Hierarchical integration tests (1/2/4 ranks)
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'integration_hierarchical_' --output-on-failure

# 9. Rollback and resource-fallback tests: sections (3), (6)-(10) of
#    integration_hierarchical_treepm_workflow{,_np2,_np4} (step 8). Raise the
#    deadlock bound if needed:
#    cmake -DCOSMOSIM_CTEST_HIERARCHICAL_TIMEOUT_SECONDS=3600 build/mpi-hdf5-fftw-debug

# 10. Restart equivalence
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'integration_restart_equivalence_harness|integration_restart_checkpoint_roundtrip|integration_parallel_two_rank_restart|integration_runtime_app_mpi_restart_smoke' --output-on-failure

# 11. Cosmological DMO validation
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'validation_dmo_zeldovich_workflow_single_rank|validation_hierarchical_dmo_growth|validation_hierarchical_dmo_rank_equivalence' --output-on-failure

# 12. Matched AREPO comparison (needs numpy/h5py and matched AREPO output)
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'python_compare_dmo_arepo_identity|python_treepm_snapshot_tools' --output-on-failure
python3 tools/compare_dmo_arepo.py --help   # then run on the matched pair with --arepo-id-offset as appropriate

# 13. Memory and performance benchmarks (release build)
cmake --preset mpi-hdf5-fftw-release && cmake --build --preset build-mpi-hdf5-fftw-release
ctest --test-dir build/mpi-hdf5-fftw-release -R 'unit_memory_governor|unit_memory_accounting|unit_process_memory|unit_benchmark_harness_smoke' --output-on-failure
# bench_tree_pm_coupling / bench_pm_solver / bench_tree_pm_sweep per docs/profiling.md
```

Build-tree paths assume the presets' default binary directories
(`build/<preset>`). Check `CMakePresets.json` if yours differ.

## E. Known remaining limitations (classified)

| Limitation | Classification |
| --- | --- |
| Multirate accuracy/order, force error, long-run P(k) and halo statistics, AREPO agreement | Scientific validation pending |
| Order gates `[1.8, 2.2]` in the driver suites may need asymptotic-regime step sizes | Scientific validation pending |
| Interior-cycle restart | Deliberately unsupported configuration (no stop source can end inside a cycle; enforced) |
| Memory-budget-dependent schedule under `single_block_fallback` | Deliberately unsupported for reproducibility (use `fail`) |
| Tree reuse/refit runs agree to roundoff after a rollback or restart | Deliberately unsupported (bitwise) configuration; documented policy |
| Gas/source terms, AMR in hierarchical/multirate KDK | Future physics feature outside current scope (DMO-only enforced at segment start) |
| Lossy precision compaction, ID narrowing | Optional optimization not justified by evidence |
| LET payload delta exchange, asynchronous neighbour collectives, per-particle lazy drift, node-level analytic tree prediction, zero-copy periodic unwrap | Optional optimization not justified by evidence |
| Unresolved source-level defects | **None known.** Compilation errors are possible and are the owner's first validation step. |

## F. Commit history

See the final handoff message for the commit identifiers.
Branch: `claude/vigilant-fermat-i8jvqw`.

## G. Static integration audit (source inspection only)

* **New pure virtuals.** All five `GravityRuntime` additions are implemented in
  `GravityRuntimeImpl` (the only implementation in `src`, `tests` and
  `include`). The removed `diagnosticsMark`/`discardDiagnosticsSince` have no
  remaining references.
* **Config.** The new key has an enum, parse/name functions, a default-table
  entry, normalized dump, consume, docs row and a parser test. Normalized-text
  round trip is tested.
* **Includes.** `time_coordinator.cpp` adds `<cstring>` and `<new>`; the
  driver test adds `cosmology.hpp`, `<algorithm>`, `<exception>` and
  `<utility>`; the workflow test uses `<functional>`, `<iomanip>` and
  `<sstream>` (all present).
* **Designated initializers.** New `ReferenceWorkflowOptions` fields are
  appended after the existing ones, so existing in-order initializer lists stay
  valid. The `HierarchicalBlockRecord` initializer follows declaration order.
* **Persisted state.** No snapshot, restart or provenance schema change. The
  normalized config text gains one key, so config hashes change.
* **Exceptions.** The `noexcept` methods (`commitCycleTransaction`,
  `invalidateSourceHistoryForRollback`, `releaseStorage`) do not allocate.
* **Exact reference.** No change: all new code is guarded by `multirate_pm`
  or plan enablement.
* **Corrected contradictions.**
  * Stale equal-duration comments in `config.hpp` and `time_integration.hpp`.
  * The old `discardDiagnosticsSince` reference in `docs/time_integration.md`.
  * Inflated soft-stop accounting.
  * An unaccounted retained capacity after a governor refusal.
  * Byte-wise fingerprint cost (now word-wise).
