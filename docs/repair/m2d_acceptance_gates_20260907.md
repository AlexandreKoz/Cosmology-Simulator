# CHUÍ / CosmoSim — M2D Acceptance-Gates Source Repair

**Date:** 2026-09-07  
**Mode:** Focused repair  
**Base:** `Cosmology-Simulator-main(16)(4).zip`  
**Base SHA-256:** `7419c320ce7d97253ae3844195f8afa5ec17fdc951f6144bfd70d6d5c1e87aca`  
**Verdict:** Source repair and focused CPU gate pass; full M2D acceptance remains provisional.

## Scope and implemented changes

The current source was read after `AGENTS.md` and the repository's current
architecture, runtime-truth, memory, repair, and build contracts. No source
changes were imported from a previous worktree. The existing M2D SFC and
heavy-analysis fixes were preserved.

1. The existing runtime registry now accepts a state-dependent owner contract
   with a byte estimate, completeness, and uncertainty. A qualified read-only
   query exposes the model. Legacy callbacks and custom owner subclasses remain
   source compatible. Complete dispatcher-owned tasks hold RAII reservations;
   complete owner-managed tasks preflight without double charging. Incomplete
   models cannot authorize overlap or prematurely reject a shrinkable batch.
2. Gravity's existing conservative TreePM memory calculation is shared with
   its registry callback and physical admission. Retained gravity baseline is
   subtracted once. The existing collective gate and phase lease remain intact.
3. Source and hydro/AMR owners publish checked, capacity-aware known workspace
   estimates. Their complete-population/transaction peaks are not fabricated.
   The source no-op path is explicitly zero; mutating source paths remain
   incomplete. No source equations, batch order, or actual allocation contracts
   were changed.
4. Output/restart publishes a state- and cadence-dependent maximum of its
   snapshot and restart staging models, with the original planned-overlap
   subtraction shared with physical reservation. These phases have separate
   release boundaries. Writer/readback metadata and allocator growth remain
   explicitly unmodeled, so active output is not a complete overlap contract.
5. The registry regression exercises dynamic and partial contracts, checked
   arithmetic, RAII release, owner preflight, and the incomplete-owner safety
   boundary. The reference-workflow budget test recognizes the new earlier
   gravity admission while still requiring the actual one-byte hard ceiling.

## Resource and scientific evidence

The prior analytical 256³/512³ analysis-mesh bounds and 86,016-case SFC
feasibility regression remain historical context until independently rerun.
The current source preserves the accepted numerical operators, precision,
TreePM criteria, CFL/reconstruction, AMR refinement, source equations, and
physical restart schema. No configuration or provenance schema changes were
introduced. The serial stage order and optional-diagnostic cadence semantics
are unchanged. Dynamic admission can reject a required task earlier, but does
not alter the physical calculation or create a new numerical scheduling order.

## Fresh validation

The following commands were executed against this extracted base and its
patched worktree. Logs are outside the source bundle.

- `cmake --preset cpu-only-debug` — PASS (GCC 14.2.0, CMake 3.31.6).
- `cmake --build build/cpu-only-debug --target cosmosim_workflows -j5` — PASS.
- `cmake --build build/cpu-only-debug --target test_unit_runtime_module_registry -j5` — PASS.
- `ctest --test-dir build/cpu-only-debug -R '^unit_runtime_module_registry$' --output-on-failure --timeout 90` — PASS, 1/1.
- Focused six-target CPU build (analysis, governor, distributed-memory, registry,
  reference workflow, source runtime) — PASS.
- The first six-test run passed 5/6 and exposed a reference-test expectation
  limited to the old gravity admission label. The test was corrected to require
  the same hard_limit_bytes=1 and accept the new earlier gravity task-preflight
  label. No numerical or safety assertion was removed.
- `OMP_NUM_THREADS=2 ctest --test-dir build/cpu-only-debug -R '^integration_reference_workflow$' --output-on-failure --timeout 90` — PASS after that correction, 1/1.
- The source-runtime integration, analysis, governor, and SFC unit tests passed
  in the first focused invocation. A final combined focused run and extended
  dependency matrix are recorded in the bundle acceptance summary.

### Additional completed validation on the frozen source

- `cmake --preset hdf5-debug` — PASS; HDF5 1.14.5.
- `cmake --build build/hdf5-debug --target test_unit_analysis_diagnostics test_unit_memory_governor test_unit_parallel_distributed_memory test_unit_runtime_module_registry test_integration_reference_workflow test_integration_star_formation_source_runtime -j5` — PASS after incremental continuation. The earlier commands were interrupted externally without compiler errors.
- `OMP_NUM_THREADS=2 ctest --test-dir build/hdf5-debug -R '^(unit_analysis_diagnostics|unit_memory_governor|unit_parallel_distributed_memory|unit_runtime_module_registry|integration_reference_workflow|integration_star_formation_source_runtime)$' --output-on-failure --timeout 90 -j2` — PASS, 6/6, 0 failures, 34.19 seconds.
- `bash scripts/ci/check_repo_hygiene.sh` in a clean patched-source extraction — PASS.
- `cmake --preset mpi-hdf5-fftw-debug` — BLOCKED at configuration: `Could NOT find MPI (missing: MPI_CXX_FOUND CXX)`.
- `cmake --preset pm-hdf5-fftw-debug` — BLOCKED at configuration: required serial double-precision FFTW3 development library not found.

The unavailable dependency configurations were attempted once each. Full CPU/HDF5 inventories, extracted-source rebuild, real MPI/FFTW runtime, rank maximum/mean RSS, work-imbalance and topology measurements, and broad scientific equivalence remain unexecuted. No whole-process memory ceiling or production scaling result is claimed from the analytical owner models.

## Open acceptance gates and M2E handoff

A complete population-growth/AMR/output peak model was not established by this
patch. The existing local governor remains the admission authority for those
phases. Real MPI+FFTW execution, fault injection, migration/restart equivalence,
whole-process peak measurements, and topology tuning cannot be inferred from
CPU source inspection. The first M2E production acceptance must not rely on
unrestricted concurrency or an unproven 48 GiB envelope. Complete the remaining
owner models and run the authoritative distributed and whole-process matrix on
provisioned hardware. Do not restart the memory architecture or introduce a
second governor, predictive controller, mixed precision, or GPU scheduler.

**Suggested branch:** `campaign-m2d-acceptance-gates`  
**Suggested PR title:** `Complete major-task memory contracts and certify M2D runtime scheduling`
