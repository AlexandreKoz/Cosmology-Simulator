# MPI hardening status

_MPI.0--MPI.8 hardening status map, updated after directed AMR exchange repair._

## Verdict

MPI support is still not accepted as complete production distributed execution.
MPI.1 gives the production `cosmosim_harness` executable ownership of the MPI
session boundary while preserving MPI-aware libraries as observers. MPI.2 keeps
the PM reverse-routing ownership repair. This repair adds a directed AMR patch
payload exchange surface and routes AMR flux-register payloads to their owner
ranks instead of globally replicating live patch/cell/flux records. Local runtime
acceptance remains pending a successful real multi-rank execution when an MPI
launcher and MPI-enabled dependency stack are available; CPU-only builds cannot
be used as MPI acceptance evidence.


## Current hardening stage map

| Stage | Scope | Current status | Acceptance evidence still required |
| --- | --- | --- | --- |
| MPI.0 | Evidence/status truthfulness | Documentation and CI labels now distinguish CPU, serial AMR, and real MPI-launched coverage. | Real MPI lane execution on the dependency-complete preset. |
| MPI.1 | Executable lifecycle | Preserved: executable owns init/finalize; libraries observe. | Harness smoke under `mpiexec`. |
| MPI.2 | Routed PM | Preserved: origin-rank/request-sequence reverse routing and byte-safe all-to-all payload layouts. | Uneven periodic PM two-rank run. |
| MPI.3 | Gas-cell migration | Existing migration path remains present; this repair did not fully upgrade the MPI migration proof to all requested stable-ID/scheduler assertions. | Strengthened two-rank migration/restart test with named cross-rank patch/cell movement. |
| MPI.4 | Hydro ghosts | Existing distributed hydro ghost/correction test remains registered; this repair did not fully expand it to every requested conserved field assertion. | Two-rank workflow conservation/restart proof for mass, momentum x/y/z, total energy, and internal-energy behavior. |
| MPI.5 | AMR | Live production AMR patch/cell import now uses bounded rank envelopes plus candidate-peer payload exchange. Flux-register payloads route to authoritative owner ranks. Serial AMR boundary test is no longer counted as MPI coverage. | Real two-rank MPI execution of `integration_reference_workflow_distributed_amr_mpi_two_rank`; current test exercises the directed exchange protocol, not full AMR workflow/restart acceptance. |
| MPI.6 | Decomposition | Existing weighted SFC rebalance path preserved. | Two-/three-rank rebalance lane execution. |
| MPI.7 | Gravity scale-path policy | Periodic PM remains sparse-record `MPI_Alltoallv` transport and is accepted only as a small-cluster policy, not as neighbor-scalable proof. Isolated PM uses a root-workspace guard with explicit rejection when the workspace limit is exceeded. | Full MPI gravity lane and scaling artifacts. |
| MPI.8 | Restart/CI closure | CI manifest now names the new AMR MPI test and isolated-PM guard; serial AMR is not listed as MPI coverage. | Dependency-complete `mpi-hdf5-fftw-debug` lane. |

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
- **Distributed hydro and AMR acceptance**: hydro ghost refresh and conservative
  correction exchange remain existing pathways. Live AMR patch/cell exchange now
  uses directed candidate-peer transfer instead of full-record all-gather, and
  AMR flux-register payloads are delivered only to authoritative owner ranks.
  Full distributed hydro/AMR conservation and restart-continuation acceptance is
  still pending real MPI execution and stronger acceptance assertions.
- **Rank-count-changing restart**: `evaluateDistributedRestartCompatibility`
  requires runtime world size, PM grid, PM decomposition mode, local slab, PM
  cadence, kick state, and long-range field state compatibility.  Changing rank
  count on restart remains unsupported.
- Distributed AMR subcycling and remote temporal interpolation remain unsupported
  unless a later patch implements and validates them explicitly.

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
- AMR patch ownership is carried by patch `owning_rank` lanes. Live remote AMR
  patch ghosts are read-only consumer imports built from `AmrPatchPayloadRecord`
  and `AmrPatchCellPayloadRecord` received from candidate peers discovered by a
  fixed-size rank-envelope control plane; they are not duplicate authoritative
  state.
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
  all-to-all, particle migration all-to-all, directed AMR patch/cell exchange, hydro
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


## Directed AMR exchange repair note (2026-07-08)

The live AMR workflow no longer calls a full per-patch/per-cell AMR payload
`MPI_Allgather`/`MPI_Allgatherv` to construct remote patches. It first exchanges
one fixed-size rank envelope per rank as bounded control-plane metadata, then
uses candidate-peer point-to-point send/receive phases for patch descriptors and
cell payloads. The diagnostic surface reports candidate and neighbor peers,
patch descriptor records, cell records, flux records, control-plane bytes,
payload bytes, remote patch ghosts, remote interface candidates, and inbound/
outbound reflux counts.

AMR flux-register payloads are grouped by `owner_rank` and received only by the
authoritative rank that may apply reflux. Duplicate register keys are still
rejected before owner mutation. Remote AMR ghosts are read-only; no consumer rank
may mutate owner truth through imported cells or patch metadata.

Periodic PM is intentionally unchanged in this repair: it still uses sparse
record payloads over global all-to-all collectives. That is a documented
small-cluster transport policy, not a general neighbor-PM scalability claim.
