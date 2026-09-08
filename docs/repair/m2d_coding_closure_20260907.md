# CHUI / CosmoSim — M2D Coding Closure

**Mode:** Focused repair. **Base:** `Cosmology-Simulator-main(16)(5).zip`.
**Base SHA-256:** `7419c320ce7d97253ae3844195f8afa5ec17fdc951f6144bfd70d6d5c1e87aca`.
**Status:** Focused source implementation and CPU gates pass; full M2D acceptance remains provisional. No claim of a certified whole-process memory envelope.

## Authority and scope

The current `AGENTS.md`, runtime-truth, memory-governance, and M2D repair documents were read before editing. The supplied ZIP is byte-identical to the original base of the previously generated M2D acceptance-gates patch. That patch was applied only after its original-file fingerprints and unified patch were verified. Its 20 source payloads were then compared with the worktree. This preserves its owner-contract and task-admission repairs without recreating them.

## Implemented source closure

- A core retained-capacity transaction plans actual replacement allocations, reserves their coexistence before growth, and reconciles physical capacity into the existing governor baseline, including partial-growth exceptions. It uses O(number of lanes) metadata, checked uint64 arithmetic, and no population-scale duplication beyond each physical replacement.
- Source birth precommit now admits exact canonical, species-sidecar, species-index, and scheduler capacity before particle/cell mutation. Star and BH paths use the same existing ID-precommit authority. Empty-birth ranks participate in the growth admission gate. Standalone registries retain their compatibility behavior.
- The distributed SFC entry validation, local sorting/sampling and memory-refinement preparation now have collective error gates. The registered MPI fixture exercises the skewed hard-memory case, duplicate keys, infeasibility, disabled migration, rank-local configuration errors, and five preparation fault points. These cases are new code coverage, not claimed MPI runtime passes.
- No scientific equations, precision, time-bin policy, star/BH selection, ID construction, restart schema, or configuration key changes. Source iteration and source-stage ordering remain unchanged. Retained capacity is not charged again on later batches.
- The earlier M2D acceptance-gates patch adds state-dependent owner memory contracts, completeness/uncertainty, and owner-managed preflight without double charging. Incomplete contracts do not authorize overlap. Existing serial execution remains authoritative.

## Limits and reproducibility

The exact population-growth transaction is an owner-local physical admission mechanism, not proof that every source, AMR, or I/O allocation has a complete static peak. Existing AMR old/new transactions and output-local leases remain authoritative. Unmodeled source-ID coordination, AMR geometry/transaction scratch, and writer metadata must not be relabeled as fully certified. No speculative concurrent full-physics executor is enabled. Rejected source growth fails before accepted physical birth mutation; a partial capacity-growth failure preserves the old logical state and reconciles any already-retained capacity. This may change failure timing but not successful numerical results.

## Fresh validation

- `cmake --preset cpu-only-debug -S . -B /mnt/data/chui_m2d_coding_closure/build_cpu` — PASS, GCC 14.2.0 / CMake 3.31.6 / OpenMP 4.5.
- `cmake --build /mnt/data/chui_m2d_coding_closure/build_cpu --target cosmosim_workflows -j5` — PASS after incremental continuation; first command was interrupted after 85/99 objects without compiler errors.
- Focused build of `test_unit_retained_capacity_transaction`, `test_unit_runtime_module_registry`, and `test_integration_star_formation_source_runtime` — PASS.
- Focused CTest of those three targets — PASS 3/3, 0 failures.
- New injected partial-growth exception, exact-boundary, retained-reuse, consumed-plan rejection, and scheduler/sidecar tests — PASS 1/1.
- `cmake --build /mnt/data/chui_m2d_coding_closure/build_cpu --target test_unit_parallel_distributed_memory -j5` — PASS.
- `ctest --test-dir /mnt/data/chui_m2d_coding_closure/build_cpu -R '^unit_parallel_distributed_memory$' --output-on-failure --timeout 90` — PASS, including the existing exhaustive 86,016-case feasibility oracle.
- Focused five-target CPU build (retained-capacity, SFC, registry, reference workflow, source runtime) — PASS.
- `OMP_NUM_THREADS=2 ctest --test-dir /mnt/data/chui_m2d_coding_closure/build_cpu -R '^(unit_retained_capacity_transaction|unit_parallel_distributed_memory|unit_runtime_module_registry|integration_reference_workflow|integration_star_formation_source_runtime)$' --output-on-failure --timeout 90 -j2` — PASS 5/5, 0 failures, 26.79 s. Reference workflow 21.21 s; source runtime 5.57 s.
- Rebuilt the new transaction and source runtime after enforcing single-use plans/nonthrowing swap; their two focused tests PASS 2/2.

Additional validation and environment blockers are recorded below. The MPI-specific test code has not yet been compiled or executed with MPI enabled in this environment. The previous campaign's full CPU/HDF5 and topology claims are historical evidence, not fresh results.

## Remaining source work and acceptance boundary

This patch closes the concrete retained-capacity source-growth gap and adds the distributed preparation/fault regressions. It does **not** finish every source/AMR/output allocation model. In particular, the sharded ID registry and source metadata still have allocation paths without complete owner-wide admission; AMR already owns regrid transactions, but its total stage geometry/ghost bound is partial; snapshot/restart writer metadata and external-library allocations remain incompletely modeled. Do not label those contracts complete or enable concurrency based on this patch. A complete source/AMR/output certificate requires further owner-local bounds or verified physical leases, not a guessed aggregate constant. This is a remaining coding limitation, distinct from the unexecuted MPI/whole-process validation.

## Acceptance boundary and handoff

M2D's source-level safety foundation is usable for continued serial, locally governed development. Complete whole-task population/AMR/I/O peak certification, dependency-complete MPI fault/restart tests, rank max/mean RSS and work metrics, and representative full-physics envelope/topology measurements remain separate acceptance requirements. A controlled memory rejection is not proof of feasible migration progress. Do not promote this report to unconditional M2D acceptance or enable speculative full-physics overlap until those gates pass.

**Suggested branch:** `campaign-m2d-coding-closure`.
**Suggested PR title:** `Close retained population growth and M2D owner admission contracts`.

## Final additional validation and packaging

- `bash scripts/ci/check_repo_hygiene.sh` — PASS on the patched source tree; no build or runtime artifacts reside in that tree.
- `cmake --preset mpi-hdf5-fftw-debug -B /mnt/data/chui_m2d_coding_closure/build_mpi-hdf5-fftw-debug` — BLOCKED, exit 1: `Could NOT find MPI (missing: MPI_CXX_FOUND CXX)`; `mpi-cxx` unavailable. The correct preset was attempted once and not retried.
- `cmake --preset pm-hdf5-fftw-debug -B /mnt/data/chui_m2d_coding_closure/build_pm-hdf5-fftw-debug` — BLOCKED, exit 1: `COSMOSIM_ENABLE_FFTW=ON but FFTW3 serial double-precision library was not found.` Attempted once.
- `cmake --preset hdf5-debug -B /mnt/data/chui_m2d_coding_closure/build_hdf5-debug` — PASS; HDF5 1.14.5.
- `c++ -std=c++20 -fopenmp -DCOSMOSIM_ENABLE_MPI=0 -Iinclude -I/mnt/data/chui_m2d_coding_closure/build_cpu/generated/include -fsyntax-only tests/integration/test_distributed_sfc_rebalance_mpi.cpp` — PASS in non-MPI mode only. This does not compile the MPI-guarded branches or prove distributed runtime behavior.
- HDF5 focused build of retained-capacity, SFC, registry, analysis, reference-workflow and source-runtime targets — PASS after incremental continuation. The first command was interrupted after 70/111 steps without compiler errors; the remaining 41 steps compiled and linked successfully.
- HDF5 CTest selection of those six targets — four fast unit tests PASS. The aggregate invocation was externally interrupted while reference-workflow and source-runtime integrations were still running. Both are NOT COMPLETED in this configuration, not failed or passed.
- Full CPU/HDF5 inventories, extracted-source completeness, real MPI np2/np3/np4 execution, collective fault/restart equivalence, rank max/mean RSS, topology comparison, and whole-process 48 GiB envelope were not completed. No historical result is promoted to fresh evidence.

The authoritative changed-files ZIP is generated from the untouched input and final worktree. A throwaway Git index constructs the unified patch, which is applied to a clean extracted base with `git apply --check` followed by `git apply`. All resulting changed payloads are SHA-256 compared to the worktree and to the bundled `files/` entries. ZIP integrity is checked with `zipfile.testzip()`. The bundle contains source, tests, documentation, fingerprints, manifest, and acceptance report only; it excludes the original ZIP, build trees, binaries, caches and runtime outputs.

**Coding verdict:** The concrete retained-population-growth and distributed-preparation repairs are implemented and pass their available focused CPU tests. Complete M2D coding and production acceptance are **not** claimed: significant remaining source/AMR/output allocation contracts still require defensible owner-local bounds or governed physical leases. The current serial runtime and fail-closed incomplete-contract policy must remain in force. This is a narrower, honest handoff—not a recommendation to redesign the governor or introduce another scheduler.
