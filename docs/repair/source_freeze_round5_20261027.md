# Round 5 — final TreePM optimization completion and source freeze (2026-10-27)

**Milestone: ENGINEERING CODE-COMPLETE — VALIDATION PENDING.**

* **Not done this round:** no compilation, configuration, test, MPI, Python,
  lint or benchmark was run (owner instruction).
* **Not claimed:** compilation, passing tests, or scientific validation.
* **Static inspection only:** it cannot guarantee that the code compiles.

**Scope.** Four priorities, with no other change:
* P0: a row-aligned rollback force restore;
* P0: rollback memory governance;
* P1: hierarchical preparation caches;
* P1: regression coverage and a static audit.

Earlier records:
* [round 4](multirate_round4_completion_20261026.md)
* [round 3](multirate_stabilization_round3_20261026.md)

## A. Implementation summary

| Priority | Files | Change | Algorithm changed | Tests |
| --- | --- | --- | --- | --- |
| P0 row-aligned force restore | `include/cosmosim/workflows/gravity_runtime.hpp`, `src/workflows/internal/gravity_force_endpoint.{hpp,cpp}` (new), `src/workflows/gravity_runtime.cpp`, `src/workflows/time_coordinator.cpp`, `CMakeLists.txt` | <ul><li>New `TransactionalForceEndpoint`: row-aligned, ID-free, capacity-reused.</li><li>New `GravityRuntime::{capture,validate}TransactionalForceEndpoint`.</li><li>`rollbackCycleTransaction` restores row-aligned instead of calling the restart importer: no `unordered_map`, no allocation, every precondition checked before the first write.</li><li>The ordinary importer keeps its ID remap, extracted to `internal::remapForceEndpointById` with its messages and semantics unchanged.</li><li>The coordinator rollback is two collective phases: read-only validation on every rank, then the restore.</li></ul> | No (multirate rollback plumbing only) | `unit_gravity_force_endpoint`; workflow sections (3), (8), and the certificate-refusal failure injection |
| P0 memory governance | `time_coordinator.cpp`, `src/core/time_integration.cpp`, `include/cosmosim/core/time_scheduler.hpp`, `reference_workflow.hpp` | <ul><li>The reservation is derived from buffer types and retained capacities; the fixed 160 B estimate is removed.</li><li>Release-before-grow, so old and new allocations never overlap.</li><li>The reservation is reconciled with the actual capacity (grow before releasing).</li><li>`HierarchicalTimeBinScheduler::exportPersistentStateInto` removes the per-cycle temporary scheduler copies.</li><li>New report fields: `…_peak_capacity_bytes`, `…_allocations`.</li></ul> | No | Workflow (3), (7), (8): first allocation, reuse across cycles and rejections, refusal on all/one/several ranks, allocation failure after a granted reservation, capacity ≤ reservation, 0 retained at end |
| P1 identity row cache | `time_coordinator.cpp`, `reference_workflow.hpp` | <ul><li>The identity row list is reused while (row count, `particleIndexGeneration`, scratch span) is unchanged.</li><li>The arena is reset on a miss. **Defect fixed:** every block used to allocate a new 4N-byte list in a never-reset arena, unbounded growth.</li></ul> | No (bitwise) | Workflow (3), (8), (10), (11): cached vs uncached bitwise, one rebuild per segment, survives rollback, empty/skewed ranks |
| P1 max-mass cache | `include/cosmosim/core/simulation_state.hpp`, `src/core/simulation_state_ownership.cpp`, mass writers (`black_hole_agn`, `tracer_support`, `star_formation`, `gas_cell_ownership`, `particle_ghost_runtime`, `simulation_state_active_views`, `ic_stream_ingestion`, `python_bindings`), `time_coordinator.cpp` | <ul><li>New `SimulationState::particleMassGeneration()`: advanced by every dense-row change and by `notifyParticleMassMutation()`, which every in-place mass writer now calls.</li><li>The local maximum mass is cached on (row count, mass generation, lane address).</li><li>The global reduction is kept every block, so the collective count is unchanged.</li></ul> | No (bitwise) | `unit_particle_mass_generation`; workflow (11) homogeneous and mixed-mass, cached vs uncached |

**Commits.** See the handoff message.

## B. Rollback correctness contract

* **Row identity.** It is proven, not assumed from the particle count. The
  `TransactionalForceEndpoint` carries:
  * the particle count;
  * `particleIndexGeneration()`, which every reorder, resize, compaction,
    migration and population replacement advances;
  * the cycle-opening identity/layout certificate: a word-wise FNV-1a over the
    ID lane, mass, species, flags, rank, softening and its mask, the index
    generation, the metadata representation and the cell/gas extents.

  All three must match at restore.
* **Direct force restoration.** `restoreForceEndpointRowAligned` copies ten
  lanes by row:
  * total x/y/z (the force cache);
  * Tree x/y/z;
  * PM x/y/z;
  * the nearest-pair r².

  Then `restoreTransactionalForceEndpoint` re-stamps the force-cache
  generations (particle index, cell index, gas identity, source), the per-row
  validity, the per-row hierarchical force generation and the PM source
  generation at the new source generation, and the PM field version from the
  endpoint.
* **Allocation and exception safety.** The restore performs no allocation (the
  unit test counts global `operator new`), and no write happens before every
  check has passed.
* **Source generation.** It is advanced strictly beyond every value the
  rejected attempt used, before the restore stamps it.
* **Cache invalidation.**
  * Forced full tree build: built tree, reuse and refit candidates, topology
    maps.
  * Ghost layer and borrowed-source certificates.
  * Content-keyed caches are retained, with the validity conditions given in
    the round-4 inventory.
* **MPI coordination.**
  * Admission: `allreduce(min)` of the local verdict. A governor refusal and
    `std::bad_alloc` both count as refusal.
  * Rollback phase 1: a read-only validation on every rank, then
    `FailureCoordinator::rethrowCollectiveFailure`, so nothing is restored
    anywhere if any rank refuses.
  * Rollback phase 2: the restore, followed by a collective failure check.
  * A failure on any rank fails every rank, and no rank stays inside a
    collective.
* **Resource ownership.** The coordinator's segment-lifetime
  `MultirateRollbackStore` owns all snapshot buffers and the governor
  reservation. It is released at segment end, on failure (RAII) and on a local
  refusal.

## C. Memory accounting (derived model)

| Buffer | Element | Logical size | Retained capacity | Lifetime / owner | Reallocates |
| --- | --- | --- | --- | --- | --- |
| position x/y/z, velocity x/y/z | `double` | N | ≥ N (exact after `reserve` on growth) | segment / store | only on growth (release first) |
| SFC key | `uint64` | N_sfc (N or 0) | ≥ N_sfc | segment / store | only on growth |
| endpoint total/tree/PM xyz + nearest r² | `double` | N | ≥ N | segment / store | only on growth |
| particle scheduler bin/next/active/pending | `u8/u64/u8/u8` | N | ≥ N | segment / store | only on growth (`exportPersistentStateInto`) |
| cell scheduler (same lanes) | — | N_cells (0 DMO) | ≥ N_cells | segment / store | only on growth |
| `StateMetadata` strings | `char` | short | `capacity()+1` each | segment / store | small, reconciled after capture |
| `IntegratorState`, `PendingOutputBoundary`, epoch, fingerprint | POD | fixed | in the store object | segment | never |

* **Per-particle rollback storage:** 48 + 8 + 80 + 11 = **147 B** (139 B
  without a materialized SFC key), plus 11 B per gas cell.
* **Fixed overhead:** the metadata string capacities plus 25 × 64 B of
  allocator allowance.
* **Reusable storage:** in steady state (no growth, reservation held) a capture
  allocates nothing and the reservation is unchanged.
* **Temporary peaks:**
  * Growth releases the old storage before allocating, so old and new never
    overlap.
  * The reservation can grow during reconciliation; that is a transient
    over-count of one grant, never an under-count.
  * The restore and the rollback allocate nothing. The previous generic
    importer built an O(N) hash map, which is removed from this path.
* **Governor reservation policy:**
  * The exact fresh-storage bytes are reserved before allocating.
  * After the capture, `retainedCapacityBytes()` is compared with the
    reservation, which is grown if needed. A refusal at that point refuses
    admission.
  * Reported: `peak_capacity_bytes ≤ peak_bytes`, and `retained_bytes = 0` at
    the end.
* **Known limits.**
  * The model covers buffer capacity, not process RSS.
  * Allocator headers beyond the 64 B allowance, alignment padding inside the
    allocator, and fragmentation are not bounded exactly.
  * The scheduler's own import (`importPersistentState`) uses the scheduler's
    existing capacity, which the scheduler accounts for.

## D. Scheduler caching

| Cache | Contents | Validity key | Invalidated by | Avoided work |
| --- | --- | --- | --- | --- |
| Identity row list | `uint32` 0..N-1 in the workspace arena | N, `particleIndexGeneration()`, scratch pointer and extent | reorder, resize, compaction, migration, population replacement, foreign scratch re-preparation, test seam | one N-element checked fill per block, and the per-block 4N-byte arena growth (a fixed defect) |
| Local max mass | one `double` | N, `particleMassGeneration()`, mass-lane address | any dense-row change; any in-place mass write (notified); test seam | one O(N) mass-lane scan per block |

* Rollbacks keep both caches valid: they do not change row identity or mass.
* **Numerical equivalence.** A hit returns exactly the value a rescan would
  produce, because the key covers every mutation path. Cached and uncached
  runs must be bitwise identical; workflow test (11) checks this for
  homogeneous and mixed masses.
* **Not changed:** the global mass reduction is still done every block, and
  the core driver's per-block identity verification stays as a correctness
  guard.

## E. Owner validation manifest (do not run before compiling; NOT run here)

```bash
# 1. CPU-only compilation + unit/driver tests
cmake --preset cpu-only-debug && cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
# 2. HDF5/FFTW compilation
cmake --preset pm-hdf5-fftw-debug && cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
# 3. MPI compilation
cmake --preset mpi-hdf5-fftw-debug && cmake --build --preset build-mpi-hdf5-fftw-debug
# 4. Configuration and existing gravity tests
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'unit_config_parser|unit_tree_gravity|unit_tree_pm_split_kernel|integration_tree_gravity_vs_direct|integration_tree_pm_coupling_periodic|integration_tree_pm_fine_tick_caches' --output-on-failure
# 5. Rollback force-cache tests
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'unit_gravity_force_endpoint|integration_restart_checkpoint_roundtrip' --output-on-failure
# 6. MemoryGovernor failure tests
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'unit_memory_governor|unit_scratch_memory_governor|unit_retained_capacity_transaction|unit_memory_accounting' --output-on-failure
# 7. Hierarchical scheduler/cache tests
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'unit_particle_mass_generation|integration_hierarchical_time_bins|integration_hierarchical_timestep_regression|integration_hierarchical_multirate_convergence' --output-on-failure
# 8. Two- and four-rank MPI tests (workflow sections 3, 6-11 incl. refusal/failure injection)
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'integration_hierarchical_treepm_workflow|integration_pm_routing_admission_mpi|integration_tree_pm_remote_group_walk_mpi|integration_mpi_collective_safety|integration_dmo_empty_rank_smoke_mpi' --output-on-failure
# 9. FFTW/TreePM force accuracy
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'integration_pm_in_place_fft_layout|validation_gravity_invariant_suite|validation_tree_pm_ewald_accuracy|validation_tree_pm_recovery_ewald_accuracy' --output-on-failure
# 10. Multirate convergence and restart tests
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'integration_hierarchical_multirate_convergence|integration_restart_equivalence_harness|integration_parallel_two_rank_restart|integration_runtime_app_mpi_restart_smoke' --output-on-failure
# 11. Matched CHUI-AREPO DMO comparison (numpy/h5py + matched AREPO output)
ctest --test-dir build/mpi-hdf5-fftw-debug -R 'python_compare_dmo_arepo_identity|validation_dmo_zeldovich_workflow_single_rank|validation_hierarchical_dmo_growth|validation_hierarchical_dmo_rank_equivalence' --output-on-failure
python3 tools/compare_dmo_arepo.py --help
# 12. Memory and performance benchmarks (release)
cmake --preset mpi-hdf5-fftw-release && cmake --build --preset build-mpi-hdf5-fftw-release
# bench_tree_pm_coupling, bench_pm_solver, bench_tree_pm_sweep per docs/profiling.md
```

The deadlock bound for the workflow tests is
`COSMOSIM_CTEST_HIERARCHICAL_TIMEOUT_SECONDS` (default 1800 s); the other MPI
tests use the existing long timeout. The build directories assume the presets'
`build/<preset>` layout.

## F. Remaining risks (classified)

| Item | Classification |
| --- | --- |
| Every change in rounds 3–5 is uncompiled | Compilation pending |
| The new unit test replaces global `operator new`/`delete` (unaligned forms only); a toolchain that routes `std::vector<double>` through other forms would only weaken the count, not fail the build | Compilation pending |
| Order gates, bitwise cached/uncached equality, refusal-path behaviour under real MPI | Runtime validation pending |
| Multirate accuracy, force error, P(k), AREPO agreement | Scientific qualification pending |
| Allocator overhead beyond the 64 B/vector allowance; RSS vs accounted capacity | Runtime validation pending (documented model limit) |
| Interior-cycle restart; memory-dependent schedule under `single_block_fallback` | Deliberately unsupported (documented; `fail` policy for reproducibility) |
| LET delta exchange, asynchronous MPI, lossy precision, node-level prediction | Optional future optimization |
| Hydro/AMR/source physics in hierarchical KDK | Out-of-scope future physics |
| Unresolved source defects in the modified paths | **None known** |

The implementation is frozen. No further optimization work is planned in this
campaign.
