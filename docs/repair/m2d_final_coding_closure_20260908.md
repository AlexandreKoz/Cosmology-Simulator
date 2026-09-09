# CHUI / CosmoSim — M2D final coding follow-up

**Date:** 2026-09-08. **Mode:** focused repair.
**Base:** `Cosmology-Simulator-main(18)(3).zip`.
**Input SHA-256:** `6740570fd03cf4b1f4a2eb0e36eebc6cdef0909e27d89efd2610af1c06ff3fd2`.
**Verdict:** Concrete source-ID, I/O-copy and AMR active-index repairs are implemented and pass their focused available tests. **Complete M2D coding and production acceptance are not claimed.**

## Source authority and scope

The ZIP was integrity-checked and `AGENTS.md` read before source modification. The supplied repository is the base; no earlier patch was reapplied. Current v2 naming, typed `.param.txt` configuration, numerical order, and existing single-governor ownership are preserved. A temporary Git baseline provides exact changed-file fingerprints and a base-relative unified patch. This repair does not introduce a new executor, numerical method, precision, tolerance, source model, configuration key, or snapshot/restart schema.

## Implemented changes

### Source-ID ownership

The retained sharded ID set is replaced by a sorted vector of uint64 IDs. Lookup is binary search and sorted disjoint batches are merged backwards into already-admitted capacity, avoiding an additional full-shard merge buffer. Initial occupancy is validated in fixed local batches, and communication uses bounded peer rounds rather than retaining a population-sized MPI round plan. The same immutable birth-key validation, ID hash, collision ordinal, and deterministic winner ordering are retained.

The new `GovernedScratchArena` owns one physical allocation, a central-governor reservation, and a bounded PMR resource with no heap fallback. All new collision/validation temporary containers use it. Rank-local preparation failures are coordinated before subsequent collectives. Retained shard replacement is admitted before growth, and actual retained capacity is checked. Precommit result storage holds a lease until the source birth transaction finishes via an additive no-op-compatible interface hook. The internal registry factory permits direct tests without creating another production authority.

This does not make ID state free: the persistent cache requires approximately 8 bytes per retained ID, plus container metadata and allocator overhead. The scratch envelope is checked and depends on the admitted batch and global source-record count. It is not a claim of a fixed amount independent of all population sizes or of measured MPI peak memory. The existing source-population transaction is preserved.

### I/O metadata and payload copies

The snapshot writer uses a compact sorted index encoding particle index and sidecar row in one uint64 record. It removes three retained unordered maps and preserves arbitrary sidecar row ordering and duplicate/missing-index rejection. The exact record storage is `8 × (N_star + N_tracer + N_bh) + 256` bytes for the alignment allowance, with checked arithmetic and a physical arena lease supplied by the existing output owner. Standalone callers retain their previous optional-governor compatibility behavior. Restart module-sidecar payloads are written directly from canonical byte storage, eliminating the extra payload-sized vector without changing HDF5 data or schema.

### AMR active-set residency

The synchronized AMR hydro path replaces the sparse global active-row hash set with a sorted, governed uint32 index. Its declared storage is `4 × N_active + 256` bytes (zero when all cells are active). Patch-local membership uses the already ordered active-cell vector rather than another hash set. The solver receives the original active order, and no reconstruction, Riemann, ghost-fill, reflux, or refinement rule changes.

## Quantitative evidence and limitations

The source-ID representation is 8 bytes per stored ID rather than implementation-dependent hash nodes/buckets; no measured process saving is claimed. The snapshot index uses 8 bytes per optional-species row plus 256 bytes. A 1 GiB module payload previously required a second 1 GiB serialization vector; the direct writer removes that copy. These are source-derived allocation bounds, not measured RSS. The active lookup remains sparse in active count and has a fixed-size bounded arena. Existing allocator/external-library reserves remain separate.

### Completed validation

- `cmake --preset cpu-only-debug -B /mnt/data/chui_m2d_final_build_cpu` — PASS.
- Focused build of `cosmosim_workflows`, governor, source-runtime and registry tests — PASS.
- `OMP_NUM_THREADS=2 ctest --test-dir /mnt/data/chui_m2d_final_build_cpu -R '^(unit_memory_governor|integration_star_formation_source_runtime)$' --output-on-failure --timeout 90` — PASS 2/2.
- `cmake --preset hdf5-debug -B /mnt/data/chui_m2d_final_build_hdf5` — PASS; HDF5 1.14.5.
- HDF5 build of `test_unit_snapshot_hdf5_schema` and `test_unit_restart_checkpoint_schema` — PASS after incremental continuation.
- HDF5 CTest selection `^(unit_snapshot_hdf5_schema|unit_restart_checkpoint_schema)$` — PASS 2/2.
- `unit_particle_id_registry` — PASS, including deterministic collision rehash, duplicate rejection, retry and lease lifecycle.
- Focused CPU build of `test_unit_sidecar_row_lookup`, `test_unit_amr_hydro_geometry`, `test_unit_amr_ghost_fill`, `test_integration_reference_workflow`, `test_unit_particle_id_registry`, `test_unit_memory_governor`, and `test_integration_star_formation_source_runtime` — PASS.
- Corresponding aggregate CTest — six tests completed PASS; the command was interrupted before the reference workflow completed. The reference workflow was then rerun separately and PASS in 27.57 seconds, giving 7/7 completed passes across the invocations.
- `c++ -std=c++20 -fopenmp -I/mnt/data/chui_m2d_mpi_syntax -Iinclude -Isrc -fsyntax-only src/workflows/source_runtime.cpp` — PASS with a local declaration-only MPI stub and MPI-enabled generated feature header. This is only C++ syntax evidence, not real MPI compilation or execution.
- `bash scripts/ci/check_repo_hygiene.sh` — PASS.
- `git diff --check` — PASS.

### Environment-blocked and uncompleted validation

The authoritative `mpi-hdf5-fftw-debug` preset was attempted once: exit 1, `Could NOT find MPI (missing: MPI_CXX_FOUND CXX)`; `mpi-cxx` is unavailable. The `pm-hdf5-fftw-debug` preset was attempted once: exit 1, `COSMOSIM_ENABLE_FFTW=ON but FFTW3 serial double-precision library was not found`. An attempt to provision missing development dependencies did not complete. Production code was not altered to compensate.

Full CPU/HDF5 inventories, extracted-source completeness, real MPI np2/np3/np4 execution, distributed fault/restart equivalence, scientific batch-size comparisons, rank max/mean RSS/PSS, work imbalance and topology measurements were not completed. The earlier campaign reports are historical evidence and are not counted as fresh passes.

## Remaining coding boundary

The requested complete code closure is **not yet demonstrated**. Source-report metadata beyond the now-governed ID path, full AMR geometry/ghost/local-source/flux-register coexistence, and remaining snapshot/restart writer/readback metadata still lack complete owner-local physical bounds or complete peak certificates. The new sparse lookup does not solve those larger AMR allocations. The existing AMR regrid transaction, source population-growth transaction, output-local leases, and fail-closed task-contract completeness remain in force. This report does not relabel partial estimates as complete or enable speculative concurrency.

A future completion must replace or physically govern those remaining allocations, prove release/coexistence boundaries, and execute the dependency-complete distributed matrix. A controlled memory rejection is safety evidence but not proof of feasible migration progress. The whole-process workstation ceiling must be measured, not inferred from these local allocation estimates.

**Suggested branch:** `campaign-m2d-final-coding-closure`.
**Suggested PR title:** `Bound source-ID coordination and I/O staging under M2D memory contracts`.

## Packaging-only completion

**Generated UTC:** 2026-09-09T00:03:42.046858+00:00. This follow-up completes the artifact handoff only; it does not implement the remaining source/AMR/output contracts.

The persistent final worktree was compared with the preserved intermediate ZIP. The only previously unbundled differences were the three-line source-ID result-capacity guard in `src/workflows/source_runtime.cpp` and removal of trailing whitespace from the authority-date line in `CURRENT_STATUS.md`. The remaining 19 payloads were unchanged before this documentation update. The original ZIP was verified against all 666 files in the Git base, with no discrepancies.

The guard rejects a result-vector capacity exceeding the admitted birth-key count before filling or committing the occupied shard. It preserves the existing ID algorithm and collective preparation boundary. Its implementation had already been compiled and linked in the persistent CPU build. Fresh packaging-session validation:

```bash
cmake --build /mnt/data/chui_m2d_final_build_cpu --target test_unit_particle_id_registry test_integration_star_formation_source_runtime -j2
OMP_NUM_THREADS=2 ctest --test-dir /mnt/data/chui_m2d_final_build_cpu -R '^(unit_particle_id_registry|integration_star_formation_source_runtime)$' --output-on-failure --timeout 90 -j2
bash scripts/ci/check_repo_hygiene.sh
git diff --check
```

Build PASS (`ninja: no work to do`); CTest PASS 2/2, zero failures (7.73 seconds); hygiene PASS; whitespace check PASS. The existing report's other completed tests remain historical evidence from the preceding repair session and have not been rerun here. No new MPI/FFTW, full-inventory, topology, or whole-process memory result is claimed.

The final changed-files ZIP is regenerated against the exact supplied base. The unified patch is checked and applied to a clean extraction, every resulting changed file is compared byte-for-byte with the worktree, all bundled payload hashes are checked, and ZIP integrity is verified. The acceptance status remains **partial coding progress, not full M2D closure**.
