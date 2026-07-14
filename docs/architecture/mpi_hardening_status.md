# MPI hardening status

_MPI.0--MPI.8 hardening status map, updated after the 2026-07-13 gravity hardening run._

## Verdict

MPI support is still not accepted as complete production distributed execution.
MPI.1 gives the production `cosmosim_harness` executable ownership of the MPI
session boundary while preserving MPI-aware libraries as observers. MPI.2 keeps
the PM reverse-routing ownership repair. Periodic PM/TreePM now has real
MPI+HDF5+FFTW runtime evidence at np1, np2, np3, and np4, including strict
single-rank/distributed equivalence, empty/source-only/target-only ranks,
FFTW-compatible slabs, and periodic halo exchange. That is small-cluster gravity
acceptance, not general MPI production or scaling acceptance: hydro/AMR,
rank-changing restart, LET, sparse-peer PM, and large-rank evidence remain open.
The complete local MPI+HDF5+FFTW CTest inventory passes 150/150; this closes
the local feature-preset regression lane, not the stronger hydro/AMR oracle or
CI-runner requirements listed below.


## Current hardening stage map

| Stage | Scope | Current status | Acceptance evidence still required |
| --- | --- | --- | --- |
| MPI.0 | Evidence/status truthfulness | Documentation and CI labels distinguish CPU, serial AMR, and real MPI-launched coverage; the dependency-complete local lane passes 150/150. | Preserve execution on CI rather than relying only on local evidence. |
| MPI.1 | Executable lifecycle | Preserved: executable owns init/finalize; libraries observe; the registered harness MPI smoke passes locally. | Preserve the harness smoke on CI. |
| MPI.2 | Routed PM | Runtime-validated at np1--np4: origin-rank/request-sequence reverse routing, FFTW-compatible slab ownership, byte-safe all-to-all layouts, and exact response coverage. | Maintain the np1--np4 lane and add scale evidence beyond four ranks. |
| MPI.3 | Gas-cell migration | Existing migration path remains present; this repair did not fully upgrade the MPI migration proof to all requested stable-ID/scheduler assertions. | Strengthened two-rank migration/restart test with named cross-rank patch/cell movement. |
| MPI.4 | Hydro ghosts | Existing distributed hydro ghost/correction test remains registered; this repair did not fully expand it to every requested conserved field assertion. | Two-rank workflow conservation/restart proof for mass, momentum x/y/z, total energy, and internal-energy behavior. |
| MPI.5 | AMR | Live production AMR patch/cell import now uses bounded rank envelopes plus candidate-peer payload exchange. Flux-register payloads route to authoritative owner ranks. Serial AMR boundary test is no longer counted as MPI coverage. | Real two-rank MPI execution of `integration_reference_workflow_distributed_amr_mpi_two_rank`; current test exercises the directed exchange protocol, not full AMR workflow/restart acceptance. |
| MPI.6 | Decomposition | Distributed compact-cut SFC rebalance passes at np2/np3/np4, including a four-rank mixed-actionability case where only some ranks have outbound migrations. | Scaling and load-quality evidence beyond four ranks. |
| MPI.7 | Gravity scale-path policy | Periodic PM remains sparse-record `MPI_Alltoallv`; TreePM uses bounded hierarchy `MPI_Allgatherv` plus target/response `MPI_Alltoallv`. Accepted only as a small-cluster policy, not as LET/neighbor-scalable proof. Isolated PM retains a bounded root-workspace path. | LET or equivalent selective remote-node import, sparse-peer PM equivalence, and scaling artifacts beyond np4. |
| MPI.8 | Restart/CI closure | Focused MPI+HDF5+FFTW gravity gates pass at np1--np4 and the complete local feature inventory passes 150/150. CI manifest also names AMR and isolated-PM guards without counting serial AMR as MPI evidence. | Run the full feature lane on CI; the green local inventory does not strengthen the still-limited hydro/AMR scientific or restart oracles. |

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
  `src/parallel/distributed_memory.cpp` detects initialized, non-finalized
  `MPI_COMM_WORLD`, exposes `worldSize()`/`worldRank()`, validates
  `parallel.mpi_ranks_expected`, and provides sum/xor all-reduce wrappers.
  MPI-enabled libraries treat a process that has not called `MPI_Init` (or has
  already finalized) as serial; PM, TreePM, and PM-halo public entry points do
  not query communicator metadata in that state.
- PM slab ownership is represented by `parallel::PmSlabLayout`,
  `pmOwnedXRangeForRank`, and `pmOwnerRankForGlobalCell`. Distributed periodic
  PM paths require `COSMOSIM_ENABLE_FFTW=ON` and `COSMOSIM_ENABLE_MPI=ON`.
- PM x ownership matches FFTW-MPI's fixed input blocks exactly:
  `B=ceil(Nx/world_size)` and rank `r` owns the truncated interval
  `[min(rB,Nx), min((r+1)B,Nx))`. This admits legal zero-width high ranks and
  removes the former np3 ownership mismatch caused by a low-rank-remainder
  partition.
- A legal zero-width FFTW-MPI rank retains logical extent zero, uses a
  one-element backend-safe dummy allocation, and still enters plan creation,
  forward/inverse transforms, and all PM collectives. It is not reduced to a
  halo-only smoke participant.
- Periodic PM halo exchange derives peers from global plane ownership, skips
  empty slabs, handles a two-rank same-peer left/right topology, posts all
  `MPI_Irecv` operations before `MPI_Isend`, and completes with `MPI_Waitall`.
  Preparation/allocation failure is voted before point-to-point traffic, and
  all ranks fingerprint shape, requested depth, boundary mode, and sequence
  before proceeding. The np2/np3/np4 halo registrations passed in the
  2026-07-13 feature build.
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
- Distributed TreePM hierarchy records use an explicit version-1 little-endian
  wire format with decomposition/force epochs, exchange sequence, and periodic
  geometry frame. Every rank supplies exactly one root; a source-empty rank
  supplies an explicit zero-mass/zero-source sentinel rather than disappearing
  from coverage.
- Short-range target requests and responses use explicit versioned wire records
  with peer identity, batch/request identity, target identity, runtime epochs,
  finite-value checks, and exact expected/seen coverage. Explicit target
  coordinate lanes allow a source-empty rank to own targets without a dummy
  particle. Zero-record ranks still enter every globally coordinated batch.
- PM density contributions, force/potential requests, and force/potential
  responses also use explicit version-1 little-endian `PMW1` records rather
  than raw C++ object representations. Operation/layout, mesh shape, box,
  scale factor, gravity constant, assignment/deconvolution, execution,
  residency, decomposition, boundary, split, workspace, and entry-specific
  control consensus is voted from the actual initialized MPI world before
  assignment, periodic/open solve, or gather. A full-domain serial reference
  call inside MPI is legal only when all ranks enter the same call and layout
  mode.
- TreePM long-range reuse is guarded by a cache signature covering
  force epoch, force-evaluation scale factor, `G_code`, split and box geometry,
  assignment, boundary, decomposition mode, and deconvolution. Particle
  decomposition epoch is deliberately excluded because PM field ownership is
  the fixed FFT slab map. A cache mismatch makes an explicit reuse request fail
  coherently; divergent explicit refresh/reuse votes fail before PM
  collectives.
- Runtime gravity diagnostics expose local source/target/tree counts, global
  empty-rank counts, hierarchy packets, unique communicating peers, PM
  solve/reuse, halo values, local FFT slab dimensions, and tree
  build/multipole/opened-node work. These make empty-rank and communication
  behavior auditable but are not by themselves a scaling campaign.
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

- **Gravity scaling maturity**: owner-local PM routing and periodic TreePM
  correctness now have real np1--np4 runtime evidence. The transport still uses
  communicator-wide all-to-all count/data phases and a bounded hierarchy
  all-gather. Four ranks are not a strong/weak scaling campaign, and this is not
  LET or sparse-neighbor acceptance.
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
  conservative correction all-gather, and nonblocking PM slab halo
  `MPI_Irecv`/`MPI_Isend`/`MPI_Waitall`.
- Only owner ranks may mutate authoritative particle rows, PM slab cells, and
  AMR patch ownership.  Ghost and pseudo-particle records are derived or
  read-only on non-owner ranks.

## Numerical invariants and remaining evidence

MPI.2 now enforces its owner-local reverse-routing and one-response-per-request
protocol in source and has command-backed np1--np4 runtime evidence. The focused
gravity/MPI regression matrix passed 36/36 and the complete feature inventory
passed 150/150 in the 2026-07-13 MPI+HDF5+FFTW build. Later hydro/AMR/scale
invariants below remain separate and are not closed by gravity success.

- Routed PM density contributions are accumulated only by the owning slab. For
  force and potential reverse interpolation, the source implementation now
  rejects missing or duplicate responses and mutates only origin-owned particle
  rows after origin-side validation; uneven and empty-rank multi-rank fixtures
  now exercise this contract.
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
| MPI.2 PM routing | Owner-local reverse routing with origin rank, origin particle index, request sequence, 64-bit exchange epoch, exact response coverage, FFTW-compatible fixed-block slabs, and checked MPI_BYTE layouts; serial PM numerics unchanged. | Passed np1--np4 focused ownership/halo validation; retain as a regression gate and extend beyond four ranks for scaling evidence. |
| MPI.3 TreePM exchange | Versioned hierarchy and request/response records; explicit empty sentinel and target-only ranks; preserved `batch_token`/`request_id` coverage under stressed batching. Solve-entry metadata, PM halo-cache commit, compact active/PM/zoom target preparation, periodic unwrap/local-tree construction, worst-case traversal-stack reservation, reusable workspace growth, and each throwable protocol phase are failure-voted before peers advance. | Periodic TreePM and strict Phase-2 np1--np4 gates pass; LET/sparse-peer scaling remains open. |
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
  two/three/four-rank PM-halo, periodic TreePM, Phase 2 gravity, and DMO
  validation tests. The production executable
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
- Local feature-path evidence on 2026-07-13 used the
  `mpi-hdf5-fftw-debug` build. The final focused gravity/MPI matrix passed
  36/36, the production-workflow DMO matrix passed 5/5 including its rank
  comparator, and the complete feature inventory passed 150/150. This confirms
  local dependency availability and actual launcher execution for up to four
  ranks; it does not substitute for CI or scaling evidence.
- The DMO comparison also distinguishes a true non-MPI serial executable from
  MPI `np1`: 64 stable IDs compare with zero periodic position error, maximum
  velocity error `0`, and zero relative mass error. Their
  power-spectrum JSON is byte-identical (SHA-256
  `8acbf3d6826250ca4fbabb0761511ff1178cbb7c20472cb2d8c8073dd16d355c`).
  Across MPI np1--np4, maximum periodic position error is zero and maximum
  velocity error is `2.7755575615628914e-17`; this is decomposition-equivalence evidence,
  not a scaling result.

## MPI.8 restart topology contract

The executable restart contract is same-topology only. `ReferenceWorkflow` now validates a restart payload's
normalized config hash, MPI world size/rank, distributed restart schema, PM grid, decomposition mode, per-rank
slab table, owner table, and TreePM cadence/field metadata before resuming. Fresh-start decomposition is not
rerun for checkpoint payloads. Rank-count-changing restart, arbitrary topology remap, and treating rank-qualified
file stems as portability evidence remain prohibited claims. Per-rank HDF5 files are serial HDF5 artifacts; no
parallel HDF5/MPIO path is claimed here.

## Non-goals

- No parallel HDF5/MPIO claim, large-rank scaling claim, LET claim,
  rank-count-changing restart claim, or arbitrary topology-remap acceptance.
- No implication that periodic gravity np1--np4 closes distributed hydro or
  AMR production acceptance.


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

## Periodic gravity hardening note (2026-07-13)

This later gravity repair supersedes the earlier same-peer halo validation
pending status without rewriting that history. PM slab ownership now matches
FFTW-MPI's fixed ceil-block layout, and the halo protocol uses nonblocking
receive-first/send-second completion. The np2, np3, and np4 halo gates pass,
including zero-width slabs and the same-peer two-rank ring. A dedicated
`Nx=2,np=4` periodic solve proves that zero-width ranks participate in the
complete FFTW-MPI solve, not only in the halo protocol.

Distributed SFC rebalance also has np2/np3/np4 runtime coverage. Its actionable
migration vote is unconditional on every rank; the four-rank fixture includes
ranks with outbound work and ranks with no local payload, preventing a local
short circuit from skipping the collective while peers enter migration.

TreePM now builds a coherent largest-gap periodic tree frame, emits an explicit
zero-source hierarchy sentinel, carries version/geometry/decomposition/force/
exchange identity in explicit wire formats, and accepts independent active
target positions on source-empty ranks. Periodic TreePM and strict Phase-2
equivalence pass at np1--np4; the unchanged strict limits are
`relative_L2 <= 5e-6` and `max_relative <= 5e-5`.

The transport decision remains conservative: hierarchy summaries are globally
gathered and target/response records use communicator-wide all-to-all phases.
Correct zero-count participation and cutoff pruning improve correctness and
payload volume, but they do not turn this into a locally essential tree. LET or
selective remote-branch import remains the next major gravity scalability
stage.
