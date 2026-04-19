# TreePM Phase 2 distributed integration closeout

Date: 2026-04-19 (UTC)

This closeout records the hard-gate evidence and explicit failure contracts for the
Phase 2 distributed TreePM implementation path.

## Hard-gate checklist (must pass as written)

- [ ] Configure MPI + HDF5 + FFTW:
  ```bash
  cmake --preset mpi-hdf5-fftw-debug
  ```
- [ ] Build MPI gravity and distributed-contract tests:
  ```bash
  cmake --build --preset build-mpi-hdf5-fftw-debug -j4 \
    --target test_validation_phase2_mpi_gravity test_unit_parallel_distributed_memory
  ```
- [ ] Run Phase 2 MPI gravity gate (single/two/three-rank CTest entries):
  ```bash
  ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure \
    -R "validation_phase2_mpi_gravity_single_rank|validation_phase2_mpi_gravity_two_rank|validation_phase2_mpi_gravity_three_rank"
  ```
- [ ] Run distributed contract/unit floor:
  ```bash
  ctest --test-dir build/mpi-hdf5-fftw-debug --output-on-failure \
    -R "unit_parallel_distributed_memory|integration_provenance_roundtrip|integration_docs_scaffold"
  ```
- [ ] Generate scaling artifacts (performance evidence, not correctness gate):
  ```bash
  cmake --build --preset build-mpi-hdf5-fftw-debug --target generate_mpi_gravity_scaling_artifacts
  ```

## Ownership/message contract summary

- **Particle owners** store authoritative particle state and own particle updates.
- **PM slab owners** own x-slab mesh cells for deposition, FFT/Poisson/gradient solve, and
  interpolation responses.
- **Tree short-range residual** uses active-target export/import with explicit per-peer request
  and response batches (`tree_exchange_batch_bytes` bounded).
- **Cadence contract** is global over gravity-kick opportunities; every rank must agree on
  refresh/reuse decisions and cadence metadata on each kick.

## Explicit failure contracts covered in code/tests

- Rank-count/config mismatch: runtime topology build throws when
  `parallel.mpi_ranks_expected` disagrees with actual communicator size.
- Unsupported decomposition mode: config parsing rejects non-`slab`
  `treepm_pm_decomposition_mode`.
- Communicator/layout mismatch: distributed PM solve throws when slab layout rank metadata
  does not match `MPI_COMM_WORLD`.
- Missing distributed restart metadata: restart decode throws if required slab ownership rows
  are absent.
- Inconsistent cadence state: restart decode throws when cadence metadata is contradictory
  (for example, non-zero refresh opportunity with zero long-range field version).

## Residual limitations (honest, non-silent)

- Decomposition mode is intentionally limited to `slab` in Phase 2.
- Restart long-range continuation policy remains `deterministic_rebuild`; cached PM long-range
  arrays are not serialized.
- Performance/load-balance sophistication (e.g., pencil decomposition, advanced migration
  heuristics, GPU overlap optimization) is intentionally deferred past Phase 2 closeout.
