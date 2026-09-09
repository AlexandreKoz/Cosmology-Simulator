# CHUÍ / CosmoSim — M2D final owner-lifetime closure

**Date:** 2026-09-08. **Mode:** focused repository repair. **Verdict:** PARTIAL CODING CLOSURE. The complete three-workstream gate is not satisfied. This is a code-bearing handoff, not a whole-task or workstation memory certificate.

## Authority and baseline

The supplied `Cosmology-Simulator-main(19)(2).zip` passed ZIP integrity and safe-entry checks. SHA-256: `17695a794eeb19ed36293511de66b1de96814f9199b5444a4fe27384b896da8b`. `AGENTS.md` was read first; the current architecture, runtime truth, memory-governance and M2D documents were inspected. The current v2 naming, `.param.txt` configuration, numerical methods, ID authority, population transaction, sparse AMR lookup, output leases, and serial dispatcher are preserved. An untouched temporary Git base supplies the exact patch fingerprints. Historical audits are not counted as fresh validation. No configuration, physics, precision, restart schema/version, or HDF5 field name changes were made.

## A. Source metadata implementation

Source parent membership is now a sorted unique uint64 index rather than transient hash nodes/buckets. The parent index has a null-upstream governed arena of `8 × N_patch + 256` bytes. Sparse active-star lookup uses sorted uint32 particle indices and preserves canonical star-row order; its arena is `4 × (N_active_particle + N_star) + 256` bytes. All-active execution uses canonical row order without materializing a full index. The old retained active-star vector and its reporting entry are removed.

Star-formation plans and immutable birth keys use a bounded PMR arena of `sizeof(StarBirthPlan) × N_batch + 8 × N_max_birth + 256` bytes, where `N_max_birth = N_batch × max_spawn_particles_per_cell_step`. The report's two output arrays are separately admitted at `12 × N_max_birth` bytes and reserved before ID precommit or birth mutation. Parent and birth-metadata preparation failures are coordinated before subsequent source collectives. The accepted ID construction, collision semantics, source selection, source order, and population-growth transaction are unchanged. The optional trailing PMR resource on `StarFormationModel::applyFromInputs` preserves existing standalone callers; a supplied resource is used without heap fallback. `starFormationBirthPlanBytes()` exposes only the private plan's exact width.

**Ownership limitation:** The report reservation is an admitted upper bound, not a physical allocator wrapper or an actual-capacity reconciliation. Other retained source staging, stellar-evolution and feedback reports/events, diffusion topology metadata, and applicable BH metadata remain outside a complete owner-wide contract. This patch does not establish collective-safe recovery for every source-side allocation failure.

## B. AMR geometry implementation and lifetime

`amrHydroGeometryCapacity` provides checked topology counts for a regular patch with dimensions `(nx,ny,nz)`. For `n=nx×ny×nz`, `g=2(nx×ny+nx×nz+ny×nz)` and `f=(nx−1)ny×nz+nx(ny−1)nz+nx×ny(nz−1)+g`, retained geometry storage is:

```
n × (sizeof(AmrHydroCellDescriptor) + 8 + 4)
+ f × (sizeof(HydroFace) + sizeof(HydroFluxRegisterFace) + sizeof(AmrHydroFaceDescriptor))
+ g × (sizeof(HydroGhostCell) + sizeof(AmrHydroGhostDescriptor))
```

The model also declares `8n` bytes of local construction-row scratch. Local and remote builders reserve the exact topology capacities; local row collection now scans the authoritative identity records once rather than materializing an additional global index. `AmrHydroPatchGeometry::ownedCapacityBytes()` sums actual capacities of the eight nested arrays. The synchronized orchestrator admits aggregate geometry and parent-object capacity before construction, separately admits the maximum construction-row scratch, checks actual capacity against the bound, and records `geometry_capacity_bytes` through existing AMR diagnostics and central memory telemetry. The geometry lease is released after geometry destruction; row scratch is released after construction. The public API is additive and existing callers require no migration. All geometry, ghost, source, and flux numerical semantics are unchanged.

**Ownership limitation:** This is a geometry-owner contract only. Prepared ghost snapshots, local source state, solver work arrays, flux accumulator/entry copies, coarse-fine synchronization, descriptor construction, and regrid old/new coexistence have not been combined into a complete simultaneous live-set contract. The existing regrid transaction and sparse-index admission remain authoritative. The geometry bound may retain conservative slack rather than reconciling all nested capacities into a new baseline. No complete AMR stage peak, collective failure guarantee, or overlap permission is claimed.

## C. Snapshot/restart and distributed metadata implementation

HDF5 one-dimensional readback now reads directly into its destination allocator, including aligned canonical lanes. It checks dataset dimension narrowing, byte multiplication and destination maximum size before resizing. String datasets write directly from the existing string's byte span and read directly into the final string; module sidecar payloads read directly into `vector<byte>`. This removes an avoidable full-sized intermediate for each affected read/write, preserving existing schema, full verification, integrity checks, and transactional `.part` publication.

Distributed restart decoding now consumes one line at a time instead of retaining a vector of all line strings. Declared item and slab counts are checked for narrowing, container maximum size, and consistency with the encoded input length before allocation. Existing key parsing, validation order, duplicate/missing entry checks, and serialized format are preserved. A new regression round-trips 4,096 ownership rows across three ranks and rejects oversized declared counts.

**Ownership limitation:** These changes remove copies but do not provide complete reader/writer physical admission. The distributed serializer still builds population-scale ownership conversion and text; parser result arrays, temporary validation maps, HDF5 metadata/descriptors, opaque library memory, and old/new canonical state coexistence still need exact bounds or physically governed transactions. No new full-readback or MPI-safe restart-commit certificate is claimed.

## Quantitative evidence

The source parent index has eight payload bytes per patch plus 256 bytes of arena allowance; sparse active metadata has four bytes per active particle plus four per star plus allowance. Birth reports are admitted at twelve bytes per maximum birth. Geometry counts are exact for the regular topology and actual vector capacities are queried. For a 2×2×2 patch the checked model gives 8 real cells, 24 ghosts, 36 faces, and 64 bytes of row scratch. The test verifies actual nested capacity equals the model on the tested allocator. Removing a full-sized readback/serialization copy eliminates one payload-sized intermediate: for a 1 GiB payload that is 1 GiB of logical additional residency. These are source-derived bounds and controlled local tests, **not measured RSS savings**. No 48 GiB whole-process envelope, large-scale throughput, or MPI rank-max/mean memory result is claimed.

## Executed validation

All CTest invocations use `OMP_NUM_THREADS=2` and `--output-on-failure`. Commands use the actual out-of-tree build paths. Initial builds were interrupted and resumed successfully; early signature/resource-fallback compilation errors were fixed before the passing results.

```
cmake --preset cpu-only-debug -B /mnt/data/chui_m2d_build_cpu
PASS
cmake --build /mnt/data/chui_m2d_build_cpu --target test_integration_star_formation_source_runtime test_unit_memory_governor -j4
PASS
ctest --test-dir /mnt/data/chui_m2d_build_cpu -R '^(unit_memory_governor|integration_star_formation_source_runtime)$' --output-on-failure --timeout 90 -j2
PASS 2/2, 5.66 s
cmake --build /mnt/data/chui_m2d_build_cpu --target test_unit_amr_hydro_geometry test_unit_amr_ghost_fill -j4
PASS
ctest --test-dir /mnt/data/chui_m2d_build_cpu -R '^(unit_amr_hydro_geometry|unit_amr_ghost_fill)$' --output-on-failure --timeout 90 -j2
PASS 2/2
cmake --build /mnt/data/chui_m2d_build_cpu --target test_unit_star_formation test_integration_star_formation_source_runtime test_unit_amr_hydro_geometry test_unit_amr_ghost_fill -j4
PASS
ctest --test-dir /mnt/data/chui_m2d_build_cpu -R '^(unit_star_formation|integration_star_formation_source_runtime|unit_amr_hydro_geometry|unit_amr_ghost_fill)$' --output-on-failure --timeout 90 -j2
PASS 4/4, 5.70 s
cmake --preset hdf5-debug -B /mnt/data/chui_m2d_build_hdf5
PASS
cmake --build /mnt/data/chui_m2d_build_hdf5 --target test_unit_restart_checkpoint_schema test_integration_restart_checkpoint_roundtrip test_unit_snapshot_hdf5_schema -j4
PASS
ctest --test-dir /mnt/data/chui_m2d_build_hdf5 -R '^(unit_restart_checkpoint_schema|integration_restart_checkpoint_roundtrip|unit_snapshot_hdf5_schema)$' --output-on-failure --timeout 120 -j2
PASS 3/3, 0.15 s
cmake --build /mnt/data/chui_m2d_build_cpu --target test_unit_parallel_distributed_memory test_unit_star_formation test_unit_amr_hydro_geometry -j4
PASS
ctest --test-dir /mnt/data/chui_m2d_build_cpu -R '^(unit_parallel_distributed_memory|unit_star_formation|unit_amr_hydro_geometry)$' --output-on-failure --timeout 90 -j2
PASS 3/3, 0.22 s
cmake --build /mnt/data/chui_m2d_build_cpu --target test_unit_parallel_distributed_memory -j4
PASS
ctest --test-dir /mnt/data/chui_m2d_build_cpu -R '^unit_parallel_distributed_memory$' --output-on-failure --timeout 90
PASS 1/1, 0.22 s (including new 4,096-row and malformed-count regression)
```

The new source test exercises tight-headroom rejection, bounded-arena failure before birth mutation, retry and deterministic IDs/birth keys/mass/metals over repeated runs. Geometry tests exercise exact counts, actual capacity, tight-headroom rejection and lease release, repeated stable capacities, degenerate and invalid dimensions. Existing ghost tests pass. HDF5 schema and roundtrip tests pass. These tests do not prove the remaining full-stage contracts.

Two attempted Ninja targets, `test_integration_amr_hydro` and `test_integration_amr_regrid`, did not exist; they were not counted as tests. The actual integration targets must be used for the broader matrix.

## Broader validation and environment

The following additional commands were executed after the first verified bundle. No missing dependencies were provisioned or worked around in production code.

```text
cmake --preset mpi-hdf5-fftw-debug -B /mnt/data/chui_m2d_build_mpi
BLOCKED, exit 1: Could NOT find MPI (missing: MPI_CXX_FOUND CXX); mpi-cxx and the MPI C++ development interface are unavailable.
cmake --preset pm-hdf5-fftw-debug -B /mnt/data/chui_m2d_build_pm
BLOCKED, exit 1: COSMOSIM_ENABLE_FFTW=ON but FFTW3 serial double-precision library was not found.
cmake --build /mnt/data/chui_m2d_build_hdf5 --target test_unit_restart_checkpoint_schema test_integration_restart_checkpoint_roundtrip test_unit_snapshot_hdf5_schema test_unit_parallel_distributed_memory -j4
PASS
OMP_NUM_THREADS=2 ctest --test-dir /mnt/data/chui_m2d_build_hdf5 -R '^(unit_restart_checkpoint_schema|integration_restart_checkpoint_roundtrip|unit_snapshot_hdf5_schema|unit_parallel_distributed_memory)$' --output-on-failure --timeout 120 -j2
PASS 4/4, 0.26 s
cmake --build /mnt/data/chui_m2d_build_cpu -j4
INCOMPLETE: interrupted after 88/333 build actions; no completed full-build claim.
cmake --build /mnt/data/chui_m2d_build_cpu --target test_integration_amr_conservative_refine test_integration_amr_conservative_derefine test_integration_amr_reflux_conservation test_integration_amr_production_hydro_integration test_integration_amr_hydro_subcycling test_integration_reference_workflow -j4
PASS
OMP_NUM_THREADS=2 ctest --test-dir /mnt/data/chui_m2d_build_cpu -R '^(integration_amr_conservative_refine|integration_amr_conservative_derefine|integration_amr_reflux_conservation|integration_amr_production_hydro_integration|integration_amr_hydro_subcycling|integration_reference_workflow)$' --output-on-failure --timeout 120 -j2
PASS 6/6, 21.17 s
```

Additional repository checks:

```text
bash scripts/ci/check_repo_hygiene.sh
PASS, exit 0 (all repository guardrails passed; terminal cleanup emitted a harmless TERM-not-set message).
git diff --check
PASS
bash scripts/ci/test_source_package_completeness.sh
INCOMPLETE: interrupted during the independently extracted source build after 81/101 actions. No source-package completeness pass or source compilation failure is claimed.
```

The remaining full CPU/HDF5 inventories, source-package-completeness completion, real MPI np2/np3/np4 runs, distributed failure/restart equivalence, legal-batch-size scientific comparisons, repeated whole-stage memory stability, rank max/mean RSS/PSS and work imbalance, representative topology comparisons, and full-physics workstation measurement are not completed. The full-build interruption is not a source compilation failure. Existing partial owner contracts remain fail-closed for overlap.

## Remaining coding boundaries and handoff

1. **Sources:** physically govern/reconcile retained report and contiguous staging capacity, stellar-evolution/feedback event and recipient metadata, diffusion topology, and any remaining significant source/BH allocations. Preserve the accepted ID and population transactions. Complete collective-safe preparation and partial-growth retry tests.
2. **AMR:** cover the simultaneously resident geometry, ghost snapshots, local source state, solver scratch, flux-register storage, synchronization and regrid old/new state; include actual capacities, release boundaries, failed-growth reconciliation, and distributed admission. The current geometry fix is not a substitute for this work.
3. **I/O:** bound or stream distributed serialization and remaining sidecar/descriptor/verification staging, reconcile retained capacities, account for opaque HDF5 runtime separately, and prove failure-safe restart loading and distributed readback under tight headroom.

The existing major-task registry remains partial and non-overlappable. Owner-held leases are not converted into dispatcher-owned duplicate reservations. No speculative task overlap, new governor, configuration system, or numerical relaxation was introduced. These unresolved paths are coding gaps, not merely unexecuted acceptance tests. Full source/AMR/I/O coding closure cannot honestly be declared from this patch.

**Acceptance handoff:** After these exact remaining boundaries are closed, execute the finite M2D CPU/HDF5/MPI/scientific/whole-process acceptance matrix and then proceed to M2E. Do not start another architectural redesign or infer workstation certification from analytical bounds. A hard-memory rejection demonstrates a safety boundary, not feasible workload completion. Suggested branch: `campaign-m2d-final-owner-lifetime-closure`. Suggested PR: `Close source, AMR, and I/O memory ownership gaps for M2D`.
