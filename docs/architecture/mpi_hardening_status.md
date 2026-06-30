# MPI hardening status

_MPI.1 executable lifecycle and MPI.2 PM routed-interpolation hardening status._

## Verdict

MPI support is still not accepted as complete production distributed execution.
MPI.1 gives the production `cosmosim_harness` executable ownership of the MPI
session boundary while preserving MPI-aware libraries as observers. MPI.2 now
repairs the PM reverse-routing ownership defect and hardens the request/response
protocol, but local runtime acceptance remains pending a successful real
multi-rank execution when an MPI launcher and MPI-enabled dependency stack are
unavailable. This does not claim broader distributed gravity, hydro, AMR,
migration, scaling, or rank-count-changing restart maturity.

## Supported now

- Serial and no-MPI builds remain supported by `COSMOSIM_ENABLE_MPI=OFF` and by
  serial fallbacks in `parallel::MpiContext` reductions and exchange helpers.
- `src/core/harness_main.cpp` owns executable-bound MPI lifecycle in
  MPI-enabled builds. It calls `MPI_Init_thread` only when MPI is not already
  initialized, requests `MPI_THREAD_FUNNELED`, records whether the executable
  owns finalization, and finalizes only when it initialized MPI and MPI has not
  already been finalized.
- `parallel::MpiContext` remains an observer/communicator wrapper. It does not
  initialize or finalize MPI.
- MPI-enabled tests initialize MPI explicitly in test binaries such as
  `tests/integration/test_reference_workflow_distributed_treepm_mpi.cpp`,
  `tests/integration/test_reference_workflow_distributed_isolated_pm_mpi.cpp`,
  `tests/integration/test_tree_pm_coupling_periodic.cpp`, and
  `tests/validation/test_validation_phase2_mpi_gravity.cpp`.
- `parallel::MpiContext` in
  `include/cosmosim/parallel/distributed_memory.hpp` and
  `src/parallel/distributed_memory.cpp` detects initialized
  `MPI_COMM_WORLD`, exposes `worldSize()`/`worldRank()`, validates
  `parallel.mpi_ranks_expected`, and provides sum/xor all-reduce wrappers.
- PM slab ownership is represented by `parallel::PmSlabLayout`,
  `pmOwnedXRangeForRank`, and `pmOwnerRankForGlobalCell`. Distributed periodic
  PM paths require `COSMOSIM_ENABLE_FFTW=ON` and `COSMOSIM_ENABLE_MPI=ON`.
- MPI.2 reverse interpolation is owner-local: particle-owner ranks emit opaque
  origin-local tokens (`origin_rank`, `origin_particle_index`,
  `request_sequence`, `exchange_epoch`); slab owners validate only sender,
  epoch, mesh routing, ownership, and finite mesh data, then echo those tokens
  in force/potential responses. Slab owners never use an origin-local particle
  index to address receiver-local particle arrays.
- Origin accumulation keeps an exchange-local request registry and accepts one
  and only one response per emitted request sequence. Duplicate, missing,
  stale/mismatched epoch, wrong-sender, wrong-slot, out-of-range, and non-finite
  responses are rejected before mutation of owner-local particle results.
- PM `MPI_Alltoallv(..., MPI_BYTE, ...)` paths for routed density records,
  interpolation requests, force responses, and potential responses now use a
  shared checked record-count/displacement-to-byte conversion. Negative values,
  size mismatches, non-representable record sizes, byte-count/displacement
  overflow, and cumulative layout overflow fail explicitly.
- `tests/integration/test_pm_periodic_mode.cpp` contains a deliberately uneven
  7-versus-1 owner-local particle split and compares distributed owner-order
  force/potential with a single-rank reference for CIC and TSC. The intended
  two-rank CTest is `integration_pm_periodic_mode_mpi_two_rank`.
- Runtime workflow code in `src/workflows/reference_workflow.cpp` builds a
  `DistributedExecutionTopology`, reduces global particle/cell identity
  summaries, rank-qualifies MPI run directories and snapshot/restart filenames
  when `world_size > 1`, and serializes `parallel::DistributedRestartState`.
- Restart schema in the inspected snapshot is `cosmosim_restart_v20`; restart
  payloads include particle scheduler state, gas-cell scheduler state keyed by
  `gas_cell_id`, gravity force cache state, output cadence state, stochastic
  state, and distributed gravity metadata.

## Known blockers

- **MPI.2 multi-rank runtime acceptance pending**: the receiver-local
  particle-index ownership bug is repaired and the reverse protocol now carries
  origin rank, origin particle index, request sequence, and a checked 64-bit
  exchange epoch with exact response coverage. A real successful launcher run
  of `integration_pm_periodic_mode_mpi_two_rank` is still required before
  claiming local multi-rank acceptance when MPI is unavailable in the current
  environment.
- **Distributed isolated/open PM maturity**: `PmSolver::solvePoissonIsolatedOpen`
  gathers/scatters through a root full-field solve for distributed slabs; this
  is not a scalable distributed open-boundary PM claim.
- **Scheduler identity migration**: local migration/rebuild paths use scheduler
  records, but multi-rank exact scheduler identity exchange remains a future
  contract and cannot be replaced by `ParticleMigrationRecord::time_bin`.
- **Distributed hydro and AMR acceptance**: hydro ghost refresh, conservative
  correction exchange, and AMR payload all-gathers exist, but no distributed
  hydro/AMR conservation or restart-continuation acceptance is claimed here.
- **Rank-count-changing restart**: `evaluateDistributedRestartCompatibility`
  requires runtime world size, PM grid, PM decomposition mode, local slab, PM
  cadence, kick state, and long-range field state compatibility.  Changing rank
  count on restart remains unsupported.

## State ownership observed

- Executable lifecycle authority belongs to the production harness. Tests may
  still own MPI initialization/finalization explicitly when compiled with MPI.
- `SimulationState` particle rows are authoritative only on
  `ParticleSidecar::owning_rank == world_rank`.  Imported ghost rows are
  non-authoritative and refreshed through `LocalGhostDescriptor`,
  `GhostLayerEpoch`, `GhostExchangeBufferSoA`, and
  `commitBlockingGhostRefreshResult`.
- `GasCellIdentityMap::gas_cell_id` is the gas-cell identity authority.  The
  workflow rebuilds particle-bound gas via gas particle IDs, but this does not
  promote legacy particle IDs into a general gas identity model.
- AMR patch ownership is carried by patch `owning_rank` lanes and exchanged as
  `AmrPatchPayloadRecord` / `AmrPatchCellPayloadRecord` for validation.
- Hydro ghost copies are read-only on consumer ranks.  Owner-rank conservative
  corrections are represented by `HydroConservativeFluxCorrectionRecord`.
- Restart payloads are continuation authority for scheduler, gas-cell
  scheduler, force-cache, output cadence, stochastic state, and distributed
  gravity metadata.

## Communicator and concurrency boundary

- MPI collectives and point-to-point exchanges use `MPI_COMM_WORLD` directly in
  the inspected production/library paths; there is no alternate communicator
  object or new global MPI singleton.
- The production executable requests `MPI_THREAD_FUNNELED`. The current runtime
  calls MPI from the main thread; configured OpenMP-style worker count is
  runtime configuration and does not authorize MPI calls from worker threads.
- Multi-rank collective paths include decomposition item all-gather, exact
  ownership identity reductions, PM density/interpolation all-to-alls, TreePM
  pseudo-particle hierarchy all-gather, TreePM short-range request/response
  all-to-all, particle migration all-to-all, AMR payload all-gathers, hydro
  conservative correction all-gather, and PM slab halo `MPI_Sendrecv`.
- Only owner ranks may mutate authoritative particle rows, PM slab cells, and
  AMR patch ownership.  Ghost and pseudo-particle records are derived or
  read-only on non-owner ranks.

## Numerical invariants and remaining evidence

MPI.2 now enforces its owner-local reverse-routing and one-response-per-request
protocol in source. A successful real multi-rank execution remains the required
command-backed acceptance evidence; the later MPI-stage invariants below remain
out of scope for this repair.

- Routed PM density contributions are accumulated only by the owning slab. For
  force and potential reverse interpolation, the source implementation now
  rejects missing or duplicate responses and mutates only origin-owned particle
  rows after origin-side validation; multi-rank runtime evidence remains
  required for acceptance.
- Particle identity, gas identity, sidecars, softening override authority, and
  scheduler identity records are conserved through migration.
- Finite-volume mass, momentum, and total energy are conserved across ghost
  refresh, ghost restore, and owner-side conservative corrections.
- Coarse-fine AMR reflux remains conservative across rank-owned patch
  boundaries.
- Restart/run continuation is equivalent to uninterrupted execution for the
  same rank count, topology, PM cadence, scheduler state, force cache, and
  output cadence state.

## Per-stage acceptance map

| Stage | Prerequisite | Evidence target |
| --- | --- | --- |
| MPI.1 executable lifecycle | `cosmosim_harness` owns MPI lifecycle; libraries observe initialized MPI only. | MPI build where the harness starts under `mpiexec`, initializes/finalizes exactly once, rank-qualifies runtime artifacts, and serial/no-MPI behavior is unchanged. |
| MPI.2 PM routing | Owner-local reverse routing with origin rank, origin particle index, request sequence, 64-bit exchange epoch, exact response coverage, and checked MPI_BYTE layouts; serial PM numerics unchanged. | Successful `integration_pm_periodic_mode_mpi_two_rank` run on the deliberately uneven 7-versus-1 partition, matching single-rank CIC/TSC force and potential references. |
| MPI.3 TreePM exchange | Preserve request/response `batch_token` and `request_id` coverage under stressed batching. | Multi-rank TreePM short-range and active-subset agreement with coverage diagnostics and no missing/duplicate response packets. |
| MPI.4 migration identity | Exchange full particle, gas, sidecar, softening, and scheduler identity payloads at a safe boundary. | Multi-rank migration test proving global ID partition, gas identity, sidecars, and scheduler records after migration/restart. |
| MPI.5 hydro ghosts | Accept conservative ghost/correction ownership. | Multi-rank finite-volume hydro conservation and ghost refresh equivalence with owner-side correction accounting. |
| MPI.6 AMR payloads | Define patch ownership, ghost fill, flux-register, and reflux rank boundaries. | Multi-rank AMR refine/derefine/reflux/restart equivalence with stable `gas_cell_id` ownership. |
| MPI.7 restart topology | Keep rank-count compatibility explicit. | Same-rank-count restart equivalence; rank-count-changing restart remains a documented rejection unless a new schema/versioned design is added. |
| MPI.8 CI closure | Dependency-enabled MPI runners available. | CI regexes include real two/three/four-rank validation, distributed workflow restart continuation, PM/TreePM routing, H2 gas migration, MPI.4 hydro, MPI.5/6 AMR boundary/reflux, not only unit smoke; failures report dependency or rank-capacity blockers. |

## Build and CI evidence map

- `CMakePresets.json` names `mpi-hdf5-fftw-debug` as the distributed gravity
  development preset and enables MPI through inheritance from
  `pm-hdf5-fftw-debug`.
- `CMakePresets.json` names `mpi-release` as an MPI-enabled release preset; it
  inherits CPU release settings and enables MPI, but not HDF5 or FFTW.
- `CMakeLists.txt` links `cosmosim_parallel` with `MPI::MPI_CXX` only when
  `COSMOSIM_ENABLE_MPI=ON`; distributed periodic PM additionally needs FFTW MPI
  linkage in `cosmosim_gravity`.
- CTest registers MPI-launched tests only when `COSMOSIM_ENABLE_MPI=ON`,
  including distributed workflow, PM periodic, TreePM periodic, distributed SFC
  rebalance, gas-cell migration, hydro interface conservation, and
  two/three/four-rank Phase 2 gravity validation tests. The production executable
  MPI TreePM smoke is registered only when MPI+HDF5+FFTW are all enabled, and
  invokes `cosmosim_harness` through `${MPIEXEC_EXECUTABLE}` with two ranks.
- `.github/workflows/ci.yml` includes an `mpi-hdf5-fftw-debug` matrix entry
  whose regex covers distributed PM/TreePM workflow, gas migration, hydro, AMR,
  and `validation_phase2_mpi_gravity_two_rank`,
  `validation_phase2_mpi_gravity_three_rank`, and
  `validation_phase2_mpi_gravity_four_rank`.  The optional MPI smoke job uses
  `mpi-release` and runs only `unit_parallel_distributed_memory` plus
  `integration_parallel_two_rank_restart`; that smoke path is not feature
  validation.

## MPI.8 restart topology contract

The executable restart contract is same-topology only. `ReferenceWorkflow` now validates a restart payload's
normalized config hash, MPI world size/rank, distributed restart schema, PM grid, decomposition mode, per-rank
slab table, owner table, and TreePM cadence/field metadata before resuming. Fresh-start decomposition is not
rerun for checkpoint payloads. Rank-count-changing restart, arbitrary topology remap, and treating rank-qualified
file stems as portability evidence remain prohibited claims. Per-rank HDF5 files are serial HDF5 artifacts; no
parallel HDF5/MPIO path is claimed here.

## Non-goals

- No solver refactor, parallel HDF5 claim, scale claim,
  rank-count-changing restart claim, or broad documentation rewrite.
- No rank-count-changing restart or arbitrary topology remap acceptance.
