# Repair open issues (P01–P19 freeze ledger)

_Date captured: 2026-04-14 (UTC)_

## 2026-05-19 Stage 4 orchestrator dispatch status

| ID | Status | Area | Finding | Resolution / Remaining action |
| --- | --- | --- | --- | --- |
| P41-STAGE4-ORCHESTRATOR-DISPATCH-041 | Closed | Runtime orchestrator dispatch | The production orchestrator broadcast every stage to every callback, relying on callback-side self-filtering. | Repaired with typed per-stage handler registration and dispatch buckets, debug/contract failures for impossible off-stage direct calls, targeted unit coverage for registered-only delivery, unregistered stage absence, and deterministic per-stage order. Continue to reject future production broadcast/self-filter orchestration. |


## Stage 0 gate status update (2026-04-26 UTC)

- **Stage 0 status: NOT CLOSED.**
- Gate evidence: `docs/repair/stage0_runtime_truth_freeze.md`.
- New physics/performance work remains blocked pending closure of Stage 0 blockers below.

### Stage 0 closure blockers (this gate run)

| ID | Status | Area | Current limitation | Next-step evidence target |
|---|---|---|---|---|
| P0-10-STAGE0-GATE-028 | Closure candidate | Stage 0 closure gate | P0-09/P0-10 repair applied: missing test sources restored, schema drift repaired, runtime smoke normalized-config check aligned, and softening restart payload metadata fixed. Full closure still requires local completion of CPU/HDF5/PM preset suites. | Run the full command bundle from `docs/repair/stage0_runtime_truth_freeze.md`; mark Stage 0 closed after those suites pass locally. |
| P32-SOFTENING-PROVENANCE-029 | Closed | Softening + restart/provenance contract | Targeted CPU and HDF5 validation now pass for softening ownership, reorder/sidecar preservation, migration preservation, restart, and snapshot schema coverage; `integration_softening_ownership_invariants` builds a valid restart payload with matching normalized config hash/provenance hash and single-rank distributed restart state. | Continue to include `ctest --preset test-hdf5-debug --output-on-failure -R "restart|snapshot|softening"` in release-gate evidence. |
| P33-APP-CONFIG-ROUNDTRIP-030 | Closed pending local full-suite certification | Runtime app normalized config contract | `integration_runtime_app_smoke` now checks canonical axis-aware normalized config key `treepm_pm_grid_nx = 24`, matching the documented normalized-config contract. | Re-run runtime app smoke, config examples, and config parser tests locally. |
| P34-RELEASE-SCHEMA-DRIFT-031 | Closed pending local full-suite certification | Release/docs schema metadata drift | Release manifest, docs scaffold, and release readiness checks now agree on snapshot `gadget_arepo_v4`, restart `cosmosim_restart_v6`, and provenance `provenance_v4`. | Re-run docs scaffold, release readiness, and snapshot/restart schema tests locally. |
| P35-PM-VALIDATION-CONTRACT-032 | Closed | PM isolated/coupling/validation contract | PM+HDF5+FFTW validation now accepts either internal-node residual pruning or leaf-pair cutoff skips as cutoff traversal evidence, matching the TreePM diagnostic contract and closing the single-rank residual-cutoff expectation mismatch. | Full PM+HDF5+FFTW preset is green with `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure`; MPI+FFTW remains environment-blocked when `fftw3_mpi` is unavailable. |
| P36-STAGE2-TIMESTEP-RESTART-MIRROR-033 | Closed | Stage 2 timestep authority / restart mirror bypass | Adversarial Stage 2 review found that restart validation checked `CellSoa::time_bin` only when scheduler `bin_index` was cell-sized, allowing mixed particle/gas restarts with stale cell mirrors to pass if particle mirrors matched scheduler truth. | Patched restart validation/import to map cells through parent gas particle scheduler entries, added stale-cell-mirror rejection in `test_unit_restart_checkpoint_schema`, refreshed restart roundtrip fixtures to sync particle and cell mirrors from scheduler authority, documented compatibility behavior in `docs/output_schema.md`, and recorded the adversarial review in `docs/architecture/stage2_adversarial_timestep_review.md`. CPU and HDF5 debug gates passed locally. |


## Current blocker ledger after dependency-enabled validation

- TreePM Phase 3 implementation/closure work remains open by contract: see `docs/treepm_phase3_contract.md` hard-gate requirements and staged maturity map.
- Even with explicit PM kick-surface synchronization metadata, cadence-gated PM reuse is still an approximation and not full asynchronous multirate PM maturity.
- Zoom focused-PM no longer relies on configured spherical membership alone; file-driven high-resolution particle masks are now consumed directly. Remaining follow-up scope is richer irregular region metadata beyond particle-ID masks.
- Heavy gravity abnormality reference checks are intentionally policy-gated (`analysis.diagnostics_execution_policy = all_including_provisional`) and remain non-default for production cadence costs.


Environment blockers currently observed in this validation environment:


- MPI gravity gate now expects `validation_phase2_mpi_gravity_two_rank` and `validation_phase2_mpi_gravity_three_rank` in MPI-enabled presets; if MPI runners have fewer than 3 ranks available to `mpiexec`, this gate will fail by construction and should be treated as an environment-capacity blocker, not a silent skip.
- `cmake --preset pm-hdf5-fftw-debug` can fail when FFTW3 development files are unavailable.
- `cmake -S . -B build/mpi-fftw-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCOSMOSIM_ENABLE_TESTS=ON -DCOSMOSIM_ENABLE_MPI=ON -DCOSMOSIM_ENABLE_FFTW=ON -DCOSMOSIM_ENABLE_HDF5=OFF` fails in this container because MPI C++ tooling is missing (`Could NOT find MPI_CXX`), preventing two-rank PM distributed runtime closure in CI-like local validation.
- The `cpu-only-debug` preset used for quick repair validation sets `COSMOSIM_ENABLE_MPI=OFF`, so MPI two-rank tests such as `integration_reference_workflow_distributed_treepm_mpi_two_rank` are not registered in that preset.

Deferred to Phase 3 follow-up after axis-aware geometry repair:

- Full 2D pencil decomposition beyond FFTW transposed ownership remains open. Current `pencil` mode is now an end-to-end transposed PM path with explicit ownership contracts, but wider process-grid scaling and non-FFTW backend parity remain follow-up scope.
- Isolated/non-periodic PM is now implemented for the single-rank path; distributed isolated/open PM remains intentionally blocked until a multi-rank open-boundary strategy exists.
- CUDA PM CIC deposition/interpolation now accepts axis-aware periodic box lengths; isolated/open CUDA PM remains unsupported and fails fast.

No new in-tree CPU-path code blocker was observed on the repaired PM unit/integration single-rank checks.

Additional validation limitation for this pass:

- `test_integration_tree_pm_coupling_periodic` in `cpu-only-debug` remains expensive in this environment (naive DFT PM backend with MPI/FFTW disabled), so the performance-hardening patch evidence is currently anchored on `unit_pm_solver` plus benchmark runs (`bench_pm_solver`, `bench_tree_pm_coupling`, `bench_parallel_decomposition_exchange`) rather than full integration runtime closure in this container.
- `test_integration_tree_pm_coupling_periodic` also remains runtime-heavy after quadrupole + mature MAC tree upgrades; in this container it exceeded a 300s direct invocation timeout, so distributed short-range agreement evidence continues to rely on existing MPI-enabled validation gates (`validation_phase2_mpi_gravity_*`) outside the CPU-only preset.

## Recently closed items

| ID | Status | Area | Prior symptom | Closure evidence |
|---|---|---|---|---|
| P0-09 | Closed | Restart/reload runtime-truth round-trip invariants | Restart continuation safety needed explicit invariant checks for runtime-truth domains (particles/species/order, scheduler/bin/active set, softening overrides, gas identity, normalized config/provenance hashes) plus clear missing/legacy-field behavior. | Extended `test_integration_restart_checkpoint_roundtrip` with deterministic restart-truth checks, active-set reconstruction equivalence checks, required-field clear failure checks (`pending_bin_index`), and documented legacy compatibility coverage for missing per-particle softening dataset; updated `docs/restart_checkpointing.md` and `docs/architecture/adr_runtime_truth_ownership.md`; validated with `cmake --build --preset build-cpu-debug -j4 --target test_integration_restart_checkpoint_roundtrip` and `ctest --preset test-cpu-debug --output-on-failure -R "restart|checkpoint|roundtrip|provenance"`. |
| P0-08 | Closed | Active-set construction ownership, cache invalidation, and no-competing-authority tests | Active-set authority and solver-local derivation/invalidation rules were documented at high level but lacked explicit test-floor coverage for stale mirror divergence and mutable compact-view generation invalidation after reorder/resize. | Added deterministic active-set ownership tests in `test_unit_time_integration` (`testActiveSetAuthority`, `testActiveSetNoCompetingBuilders`) and mutable compact-view invalidation/derivation checks in `test_unit_simulation_state` (stale particle/cell view rejection after reorder/resize and solver-local view derivation from authoritative indices); updated `docs/architecture/adr_runtime_truth_ownership.md` with explicit owner/invalidation/generation/lifetime rules for active-set layers and forbidden module-local builders. |
| P0-07 | Closed | Config ownership: raw/normalized/derived/runtime/provenance lanes | Config-derived runtime values and continuation metadata ownership boundaries needed explicit classification + deterministic test floor for normalization/hash/restart metadata consistency and ambiguous time/scale semantics. | Added targeted ownership tests in `test_unit_config_parser` and `test_unit_units_cosmology_provenance`, plus restart normalized-text/hash mismatch rejection in `test_unit_restart_checkpoint_schema`; updated `docs/architecture/adr_runtime_truth_ownership.md` and `docs/configuration.md` with explicit ownership map and anti-ambiguity semantics for code-time vs scale-factor vs diagnostic redshift lanes. |
| P0-06 | Closed | Softening ownership/priority and override-preservation invariants | Softening ownership order (global/species/per-particle) and transform persistence rules were under-specified for resize/reorder/migration/restart/active extraction and diagnostics observer behavior. | Added `test_integration_softening_ownership_invariants` with `test_softening_priority_invariants`, `test_softening_override_resize_reorder_preservation`, and `test_softening_override_restart_roundtrip`; verified per-particle > species > global resolution, resize/reorder/migration persistence, restart roundtrip persistence, active-lane mirror behavior, and diagnostics non-authority; updated ADR softening preservation rules and mirror-authority notes. |
| P0-05 | Closed | Gas-cell identity and hydro-state ownership invariants | Gas identity semantics and temporary 1:1 gas-particle/gas-cell coupling were insufficiently explicit under reorder/resize/restart/active-hydro extraction paths. | Added `test_unit_gas_cell_identity_invariants` plus restart mapping checks in `test_integration_restart_checkpoint_roundtrip`; hardened the contract with `requireParticleBoundGasCellContract`, explicit row/ID helper seams, hydro active-view validation, migration-by-ID coverage in `test_integration_species_migration_invariants`, and reorder guardrails preventing silent gas-relative reorders; documented temporary contract + forbidden assumptions in `docs/architecture/adr_runtime_truth_ownership.md`. |
| P0-04 | Closed | Species migration + canonical grouping invariants | Species transitions were under-tested for canonical grouping continuity across sidecars/count ledgers/offset index maps, per-particle softening override survival, and scheduler-bin/active-set continuity after migration + reorder/resize. | Added deterministic migration invariant suite in `test_integration_species_migration_invariants` (single/batch migration, reorder/resize/active-extraction adversarial sequences, sidecar attach-detach checks, timestep-bin continuity, and softening-priority checks); repaired migration pack/commit lanes to carry optional per-particle softening overrides; updated runtime ownership docs for allowed migration path + sidecar rules + species-agnostic active-eligibility policy. |
| P0-03 | Closed | Particle ordering/resize/reorder/sidecar stale-index invariants | Index-based assumptions could silently desynchronize particle identity from species tags, sidecar payload rows, softening overrides, and cached active/kernel views across resize and reorder operations. | Hardened `reorderParticles` to move complete star/BH/tracer payload rows under `kMoveWithParent`, keep default parent-indirection rows stable while remapping `particle_index`, preserve optional softening sidecar values across particle resize, add index-generation invalidation and defensive species-sidecar ownership checks, and expand `test_integration_reorder_compaction_sidecars` coverage for reverse permutations, move-vs-indirection semantics, randomized reorder modes, compaction, and stale-view invalidation; validated with targeted build/ctest runs plus `rg "particle_index" src include tests` audit. |
| P32-TIMEBIN-AUTHORITY-027 | Closed | Timestep/bin ownership invariant coverage | Scheduler bin authority vs state mirror lanes lacked deterministic adversarial tests for unauthorized mutation, reassignment consistency, reorder pairing, and persistent-state restart equivalence. | Added targeted invariant coverage in `test_unit_time_integration` (`testTimestepBinAuthorityInvariant`, `testTimestepBinReassignmentAndRestartRoundTrip`, `testTimestepBinReorderIdentitySurvival`), plus new mirror-authority helper checks in `core/time_integration`; validated with `cmake --build --preset build-cpu-debug --target test_unit_time_integration` and `ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_time_integration|integration_hierarchical_time_bins|integration_restart_checkpoint_roundtrip"`. |
| P19-STATE-OWNERSHIP-005 | Closed | Core state ownership/invariants decomposition | `src/core/simulation_state.cpp` concentrated ownership rules, species indexing, transfer packing, metadata parsing, and active-view hot/cold logic in one implementation unit. | Split into focused implementation units; explicit gravity/hydro hot-field contract notes added; targeted invariant/contract checks added in `tests/unit/test_simulation_state.cpp`; verified by `cmake --build --preset build-cpu-debug -j4 --target test_unit_simulation_state` and `./build/cpu-only-debug/test_unit_simulation_state`. |
| P19-IO-CONTRACT-006 | Closed | Snapshot vs restart contract separation | Restart continuation-critical metadata/scheduler requirements and failure behavior were under-specified relative to snapshot interoperability boundaries. | Added shared continuation metadata contract helper, explicit restart completeness checklist, and negative tests for schema mismatch/missing required scheduler lane/finalize failure behavior; verified via `test_unit_snapshot_hdf5_schema`, `test_unit_restart_checkpoint_schema`, and `test_integration_restart_checkpoint_roundtrip`. |
| P19-GATE-FFTW-TEST-001 | Closed | PM FFTW analytic mode response | `unit_pm_solver` failed at `tests/unit/test_pm_solver.cpp:82` (`cosine_similarity > 0.98`). | `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R "unit_pm_solver|integration_tree_pm_coupling_periodic"` now passes. |
| P19-GATE-FFTW-TEST-002 | Closed | TreePM periodic PM long-range coupling (FFTW path) | `integration_tree_pm_coupling_periodic` failed with `rel_l2=18129.9` (required `<= 0.75`). | `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure` now passes (`36/36`). |
| P19-ARCH-CORE-BOUNDARY-003 | Closed | Core dependency direction | `include/cosmosim/core/reference_workflow.hpp` assembled analysis/I/O/physics workflow concerns inside `core/`. | Workflow assembly moved to `workflows/` and `integration_core_dependency_direction` guard added to fail on forbidden upward includes. |
| P19-CONFIG-CONTRACT-004 | Closed | Frozen configuration contract typing | Policy keys (`gravity_solver`, `hydro_solver`, `coordinate_frame`, boundaries, feedback mode/variant) remained string contracts after parse. | `unit_config_parser`, `unit_simulation_mode`, and `integration_simulation_mode_toy_runs` now validate enum-backed freeze behavior with deterministic normalized text/hash and centralized key+alias registry handling. |
| P19-OBS-EVENTS-007 | Closed | Operational observability | Runtime diagnostics lacked a structured run-health/event surface linked to provenance; failures were harder to audit from CI artifacts. | Added compact runtime event model in `core::ProfilerSession`, new `reference_operational_events.json` report, and test coverage in `test_unit_profiling` + `test_integration_reference_workflow`. |
| P19-DIAGNOSTICS-TIER-008 | Closed | Diagnostics maturity/scalability labeling | Diagnostics bundles mixed run-health outputs with heavy reference algorithms and did not encode maturity tier in code/output policy. | Added typed analysis execution policy, per-diagnostic tier/maturity/scalability metadata in bundle JSON, and tests proving provisional heavy diagnostics are non-default unless explicitly opted in (`all_including_provisional`). |
| P19-PARALLEL-CONTRACT-009 | Closed (R07/R07b/R07c hardened) | Distributed-memory ownership/determinism contracts | Owned-vs-ghost residency and reduction/config consensus expectations were partially implicit in ghost-owner vectors and ad-hoc checks. | Added explicit send/receive transfer-role descriptors, lifecycle intent + post-transfer residency typing, strict plan/descriptor consistency validation, unambiguous reduction baseline/measured naming, explicit reduction policy modes, and property-level config mismatch records with rank/value payloads; validated in `test_parallel_distributed_memory` and `test_parallel_two_rank_restart`. |
| P19-INFRA-GATE-010 | Closed | Infrastructure gate clarity + false-closure prevention | CPU-only pass signals could be misread as full closure when HDF5 or PM/HDF5/FFTW path validation was not surfaced as one explicit gate report. | Added `scripts/ci/enforce_infra_gates.sh` + CI `infrastructure_gates` job to run explicit CPU/HDF5/PM+HDF5+FFTW gate sets, write `infrastructure_gate_report.json`, and fail closure when any feature path gate fails; strengthened `integration_core_dependency_direction` with CMake target-link direction checks. |
| P20-CONFIG-ROUNDTRIP-011 | Closed | Frozen config canonicalization/provenance | Normalized config emitted non-parseable empty-value/self-hash combinations and could not round-trip through the same parser. | Removed self-hash embedding from normalized text, preserved empty-value parsing for `key =` canonical lines, hardened inline-comment stripping for URLs/path-like values, and added deterministic round-trip/hash checks in `test_unit_config_parser`. |
| P20-RESTART-STELLAR-012 | Closed | Restart completeness/integrity | Restart payloads omitted stellar-evolution bookkeeping lanes and the payload hash did not cover those lanes. | Bumped restart schema to `cosmosim_restart_v3`, serialized/read all stellar-evolution sidecar lanes, extended integrity hashing coverage, and verified round-trip in `test_integration_restart_checkpoint_roundtrip` + `test_unit_restart_checkpoint_schema`. |
| P20-SNAPSHOT-MASSTABLE-013 | Closed | Snapshot import safety | HDF5 snapshot import defaulted missing `Masses` datasets to zero instead of using `Header/MassTable`. | Snapshot reader now requires a positive `Header/MassTable` entry for fallback and reports `Masses=MassTable`; verified in `test_integration_snapshot_hdf5_roundtrip`. |
| P20-PARALLEL-HONESTY-014 | Closed | Distributed-memory scaffolding | Ghost descriptor planning treated local ghost rows as outbound payload rows, and distributed restart decode accepted missing ownership entries as rank 0. | Descriptor-only ghost planning now leaves outbound payload indices empty unless owned export rows are known, duplicate/missing restart ownership entries are rejected, and tests were updated in `test_unit_parallel_distributed_memory`. |
| P20-STATE-COVERAGE-015 | Closed | Species-sidecar invariants | Species-tagged star/BH/tracer particles could exist without exactly one corresponding sidecar row. | Ownership validation now enforces one-to-one species-tag/sidecar coverage and rejects duplicates; verified in `test_unit_simulation_state`. |
| P21-TREEPM-PHASE2-CONTRACT-016 | Closed | Distributed gravity contract freeze surface | Phase 2 distributed TreePM ownership and config/build/documentation surfaces were implicit and leaned on optional MPI release smoke. | Added typed numerics keys (`treepm_pm_decomposition_mode=slab`, `treepm_tree_exchange_batch_bytes`), provenance capture, `mpi-hdf5-fftw-debug` preset, CI matrix coverage for that preset, and explicit architecture freeze doc `docs/treepm_phase2_distributed_contract.md`; validated with `unit_config_parser`, `integration_docs_scaffold`, `integration_provenance_roundtrip`, and `unit_parallel_distributed_memory`. |
| P22-TREEPM-SLAB-LAYOUT-017 | Closed | PM slab ownership/storage contract | PM mesh ownership/indexing semantics for distributed slab decomposition were implicit and PM storage assumed globally replicated cube allocation. | Added typed `parallel::PmSlabLayout` ownership helpers (partitioning, owner lookup, global/local conversions), adapted `PmGridStorage` to allocate local slabs, kept one-rank full-domain behavior as degenerate slab, and added unit coverage in `test_unit_parallel_distributed_memory` + `test_unit_pm_solver` with explicit docs updates. |
| P23-RUNTIME-TOPOLOGY-018 | Closed | MPI/CUDA runtime infrastructure wiring | MPI slab ownership, `parallel.gpu_devices`, and CUDA PM execution policy existed as disconnected surfaces with no shared runtime topology/selection contract. | Added `core::CudaRuntimeInfo`, explicit `parallel::DistributedExecutionTopology` + rank/device assignment helpers, runtime world-size validation, CUDA device selection wiring in the reference workflow, new debug presets (`cuda-debug`, `mpi-cuda-hdf5-fftw-debug`), and targeted config/topology tests. |
| P24-TREEPM-DENSITY-DEPOSIT-019 | Closed | Phase 2 distributed PM density assignment | `PmSolver::assignDensity` previously rejected partial slabs and could not route owner contributions to slab owners. | Added slab-owner routed PM density contribution exchange with batched per-destination payloads and owner-side validation in `PmSolver::assignDensity`; added CIC/TSC wrap/boundary/mass-conservation and one-rank vs distributed agreement checks in PM tests; updated PM/distributed-memory contracts docs. |
| P25-TREEPM-INTERPOLATION-020 | Closed | Phase 2 distributed PM interpolation return path | PM interpolation previously required full-domain mesh residency on every rank and rejected partial slabs. | Added explicit reverse communication gather for force/potential where particle owners send stencil requests to slab owners and receive weighted contributions; validated CIC/TSC two-rank vs one-rank agreement and slab-boundary cases in `test_integration_pm_periodic_mode`. |
| P26-TREEPM-SHORT-RANGE-021 | Closed | Phase 2 distributed tree short-range residual | TreePM short-range residual previously evaluated only local sources and lacked peer-evaluated export/import tree participation. | Added explicit active-target export/import protocol with per-peer batching (`tree_exchange_batch_bytes`), remote peer tree evaluation, response coverage checks for duplicate/missing packets, workflow config wiring, and two-rank distributed-vs-single-rank coverage in `test_integration_tree_pm_coupling_periodic`. |
| P27-TREEPM-OWNERSHIP-022 | Closed | Gravity-aware decomposition + particle migration ownership commit | Decomposition cost used only generic work/memory terms, and migration intents existed without a committed owner-tag + sidecar rebuild boundary. | Added gravity-aware weighted decomposition terms (`owned_particle`, `active_target_count_recent`, `remote_tree_interactions_recent`, `memory_bytes`), per-rank decomposition component metrics, plus explicit `SimulationState` migration pack/commit APIs that rebuild ownership/species/sidecar mappings at a synchronization point; validated by `test_unit_parallel_distributed_memory` and `test_unit_simulation_state`. |
| P28-TREEPM-RESTART-CONTINUATION-023 | Closed | Distributed TreePM restart/provenance continuation safety | Restart payloads lacked explicit versioned distributed gravity continuation state, and rank/layout mismatches were hard to debug after resume. | Bumped restart schema to `cosmosim_restart_v4`; persisted `DistributedRestartState` (epoch, owning ranks, PM slab layout, cadence + field metadata, explicit restart policy); added restart compatibility diagnostics and mismatch coverage in `test_unit_parallel_distributed_memory`, plus restart round-trip validation in `test_integration_restart_checkpoint_roundtrip` and schema checks in `test_unit_restart_checkpoint_schema`. |
| P29-TREEPM-PHASE2-CLOSEOUT-024 | Closed | Phase 2 closeout coherence + explicit failure contracts | Final Phase 2 hard-gate checklist and some distributed failure conditions were not captured as one explicit closeout artifact with command-backed checks. | Added `docs/treepm_phase2_closeout.md`; hardened restart cadence-state validation and compatibility reporting (`pm_cadence_steps_match`); added explicit negative-contract coverage in `test_validation_phase2_mpi_gravity` and `test_unit_parallel_distributed_memory` for rank-count mismatch, unsupported decomposition mode, communicator/layout mismatch, missing distributed restart metadata, and cadence inconsistency. |


## Verified evidence (this run)

- CPU targeted repair floor passes:
  - `ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_config_parser|unit_simulation_state|unit_parallel_distributed_memory|unit_restart_checkpoint_schema"`
- HDF5 targeted repair floor passes:
  - `ctest --test-dir build/hdf5-debug --output-on-failure -R "unit_restart_checkpoint_schema|integration_snapshot_hdf5_roundtrip|integration_restart_checkpoint_roundtrip"`
- PM HDF5+FFTW configure remains blocked by missing FFTW3 development files in this environment:
  - `cmake --preset pm-hdf5-fftw-debug`

## Outcome tags

- **CPU-only preserved**
- **HDF5 path preserved**
- **PM HDF5+FFTW path environment-blocked**
- **P19 blocker ledger cleared**

## R08 hardening residual limitations (tracked, non-blocking)

- `integration_core_dependency_direction` target-link guard now consumes CMake Graphviz output from the configured build tree, which is materially stronger than top-level text scanning; however, it still relies on target-name pattern matching (`cosmosim_(analysis|io|physics|workflows|app)`) rather than full semantic layer tags. If naming conventions change, this check must be updated in lockstep.
- Infrastructure gate test scopes are now manifest-driven and converted into a ctest regex at runtime. This removes inline script fragility but still depends on stable test identifiers.

- Isolated/open PM remains single-rank in this stage; distributed isolated PM (MPI slabs) is intentionally blocked with runtime errors until a distributed open-boundary PM strategy is implemented.

- Tree softening maturity limitations after stage 2026-04-20:
  - Node-level accepted multipoles currently use `epsilon_node_max` as a conservative aggregate for mixed-source nodes; this favors safety/coherence over minimum-bias accuracy for strongly heterogeneous softening distributions.
  - Snapshot softening sidecar field is additive and optional; external consumers that assume only canonical GADGET datasets may ignore `GravitySofteningComoving` unless upgraded.

## New open follow-up items (post P29 closeout)

| ID | Status | Area | Current limitation | Next-step evidence target |
|---|---|---|---|---|
| P30-CLUSTERED-LOAD-MATURITY-025 | Open | Distributed gravity clustered-load maturity | Current deterministic contiguous SFC cuts reduce clustered overshoot but remain one-dimensional; extreme multi-cluster anisotropy can still create per-rank remote short-range pressure skew even when weighted load appears balanced. | Add clustered MPI scaling sweeps (2/3/4+ ranks) that jointly track weighted-load imbalance and remote request packet imbalance against acceptance thresholds in `bench_tree_only_scaling_mpi` + `bench_tree_pm_coupling`. |
| P31-STATE-TEST-DRIFT-026 | Closed | Snapshot schema contract drift | `test_unit_snapshot_hdf5_schema` had a stale `gadget_arepo_v2` expectation while runtime schema advertised `gadget_arepo_v4`; the unit contract and docs are now aligned to schema v4. | `test_unit_snapshot_hdf5_schema` rebuilt and passed in the 2026-04-26 CPU-only hardening pass. HDF5 roundtrip remains environment-dependent on HDF5-enabled preset availability. |
| P32-MIGRATION-IDENTITY-027 | Closed | Migration/compaction identity hardening | Migration, stale-ghost removal, ownership compaction, and rank-local gas rebuild needed one explicit identity-preserving contract with negative coverage for wrong-owner inbound records, stale ghosts, and stale sidecar indices. | Hardened `commitParticleMigration` stale-ghost validation; added targeted species migration negative/sidecar tests; documented `ParticleMigrationRecord` completeness and gas-cell rebuild-by-particle-ID implications. |
| P33-TRANSFORM-FUZZ-HARNESS-028 | Closed | Transform invariant hardening | Resize, reorder, compaction, migration, active-view invalidation, sparse optional sidecars, and HDF5 restart paths lacked one fixed-seed composed fuzz harness keyed by particle ID rather than row counts. | Added `integration_transform_fuzz_invariants` and `validation_transform_fuzz_invariants_long` with deterministic seed/trial failure context, all sidecar-family oracles, stale scatter negatives, gas rebuild checks, and HDF5 restart roundtrip when enabled. |
| P34-ACTIVE-VIEW-LIFETIME-029 | Closed | Active-view lifetime/generation hardening | Compact active/kernel views needed explicit transient lifetime documentation, source-generation stamps on read-only views, stale scatter coverage after migration/gas rebuild, and workspace reuse evidence. | Added generation fields to read-only `ParticleActiveView`/`CellActiveView`, kept mutable kernel scatter generation guards, documented non-ownership/hot-lane contracts in `TransientStepWorkspace`, and extended `test_unit_simulation_state` plus `test_integration_reorder_compaction_sidecars` for workspace reuse, migration invalidation, and gas rebuild invalidation. |


## 2026-04-20: Remaining Phase 3 campaign blockers

- Large-rank (beyond np2) strong/weak scaling certification sweep still required for Phase 3 closeout.
- External literature-target calibration for cosmological P(k) and halo profile systematics remains outside this stage; current campaign uses explicit in-repo reference contracts and envelopes.

## 2026-04-21: Phase 3 final integration closeout blockers

- `cmake --preset mpi-hdf5-fftw-debug` fails in this environment because `fftw3_mpi` is unavailable; distributed MPI+FFTW command bundle cannot be executed end-to-end in this cycle.
- `ctest --test-dir build/pm-hdf5-fftw-debug -R "integration_reference_workflow" --output-on-failure` fails with `runtime workflow schema compatibility validation failed`.
- `ctest --test-dir build/pm-hdf5-fftw-debug -R "integration_tree_pm_coupling_periodic" --output-on-failure` fails with `expected at least one low-res contaminant`.
- `validation_phase2_mpi_gravity_single_rank` residual-cutoff expectation mismatch is repaired; cutoff traversal evidence now accepts internal-node pruning or leaf-pair cutoff skips.
- `validation/artifacts/research_grade/phase3/correctness/power_spectrum_consistency.json` remains blocked (`missing_diagnostics_inputs`), so Phase 3 correctness/force-accuracy envelope is not fully evidenced.
- Multi-rank scaling evidence in the audited artifact set remains baseline-only and non-certifying (`phase2_baseline_scaling_summary.json`); checked-in CSV artifacts currently cover `np1` only.

Phase 3 status in this cycle: **incomplete** (see `docs/treepm_phase3_closeout.md`).

## 2026-05-07 restart sidecar exactness status

- Closed P0 blocker for v6 restart sidecar exactness in the HDF5 path: missing softening mask/value datasets and older restart schema versions are rejected, and targeted reorder/migration roundtrip tests cover identity-aligned sidecars.
- Remaining compatibility work, if needed, must be an explicit legacy importer with documented migration semantics and tests; no silent v5 compatibility is allowed.

## 2026-05-08 distributed ownership floor status

- P1 distributed ownership correctness floor is partially hardened in code by rank-local ID uniqueness summaries and reduced count/sum/xor partition checks in the reference workflow MPI path.
- Mature load balancing, pencil FFT, and production migration scheduling remain out of scope for this repair cycle.
- MPI+HDF5+FFTW closure remains dependent on the validation environment providing MPI FFTW (`fftw3_mpi`) and HDF5-enabled MPI presets; if unavailable, the MPI command bundle must be reported as environment-blocked rather than silently passed.


## 2026-05-08 notes

- P2 AMR/moving-mesh gas ownership readiness now has an RFC and isolated `GasCellIdentityMap` API guard. Remaining blocker: promotion to production requires explicit restart schema versioning, compatibility behavior, docs updates in `docs/output_schema.md`, and hydro remap tests keyed by stable `gas_cell_id`.

## 2026-05-08 Stage 1 runtime-truth gate status

- Closed CI-discipline gap for CPU-runnable P0 runtime-truth repairs: the P0 suite now has an exact CTest preset, labels, a label/feature-gating audit test, and one local/CI script under `scripts/ci/`.
- Feature-specific closure remains dependency-bound: HDF5, FFTW, MPI, CUDA, and Python presets are explicit and should fail at configure time when requested dependencies are absent; CPU-only Stage 1 success must not be reported as feature-path closure.


## Stage 2 timestep ownership follow-ups (opened 2026-05-10)

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| P35-STAGE2-TIMESTEP-OWNERSHIP-032 | Closed | Hierarchical timestep authority | Restart writer/hash/read now reject stale particle time-bin mirrors against scheduler `bin_index`, reader rebuilds mirrors from scheduler authority after validation, migration/transfer fields are documented as mirrors, and reorder/migration tests exercise stable-ID mirror preservation. | Remaining scheduler ownership work moved to P36 for new particle registration; no additional restart mirror guard is open for existing v6 lanes. |

## 2026-05-10 Scheduler ownership residuals

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| P36-SCHEDULER-SPAWN-OWNERSHIP-033 | Open | Particle creation + scheduler authority | Star-formation and black-hole creation no longer copy cell mirror bins into new particle mirror lanes, but production creation paths still need an explicit scheduler-owned element registration policy before spawned particles can participate in hierarchical activation without a rebuild. | Add scheduler element-registration/import APIs with tests covering spawned particles, mirror refresh, and restart continuity. |

## 2026-05-10 hierarchical timestep invariant status

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| P37-HIERARCHICAL-TIMESTEP-INVARIANTS-034 | Closed | Scheduler/debug invariants | Runtime checks now trap early activation, illegal sync-boundary bin transitions, skipped PM refresh commits, stale mirrors/descriptors, invalid restart timestep state, active-set mismatches, and non-monotonic tick overflow with deterministic context. | Feature-path closure still depends on HDF5/FFTW/MPI presets in environments with those dependencies available; CPU-only surrogate scheduler tests are in `test_unit_time_integration`. |


## 2026-05-10 scheduler identity exchange follow-up

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| P38-SCHEDULER-MPI-IDENTITY-035 | Open | MPI migration + scheduler authority | Local structural transforms now have stable-ID scheduler remap/rebuild helpers, but multi-rank migration does not yet have a fully specified scheduler identity-record exchange and rank-coordinated tick/max-bin contract. | Define and implement distributed scheduler identity payload exchange, compatibility checks, and MPI tests before claiming multi-rank exact scheduler continuation. |

## 2026-05-11 Stage 2 verification status

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| P39-STAGE2-GATE-INCLUSION-036 | Closed | CI/runtime-truth gates | Stage 2 scheduler-authority tests existed but the strongest CPU runtime-truth preset did not explicitly include `unit_time_integration`, and there was no Stage 2-named preset for scheduler authority, split-brain removal, PM sync legality, and restart equivalence. | `unit_time_integration` is now marked `runtime_truth;p0` with scheduler/restart/treepm labels, the label audit requires it, Stage 0/1 runtime-truth presets include it, and `test-stage2-runtime-truth-cpu-debug` names the Stage 2 critical gate directly. |
| P40-STAGE2-HDF5-RESTART-NEGATIVE-037 | Closed | Restart schema validation | Restart HDF5 roundtrip covered stale particle `time_bin` mirrors and payload-hash tamper, but did not corrupt a persisted scheduler hot lane directly. | `integration_restart_checkpoint_roundtrip` now mutates `/scheduler/active_flag` to an invalid value in-file and requires `readRestartCheckpointHdf5` to reject it before exact continuation is accepted. |


## 2026-05-11 Stage 2 documentation/release anti-drift status

| ID | Status | Area | Current blocker / ambiguity | Required follow-up |
| --- | --- | --- | --- | --- |
| P41-STAGE2-DOC-RELEASE-ANTI-DRIFT-038 | Closed | Docs/release claim discipline | Stage 2 timestep authority docs, restart/schema language, runtime truth map, ADR, repair-state docs, release notes, and docs scaffold now agree that scheduler objects own timestep truth and PM cadence legality while `time_bin` lanes are mirrors. | Keep `integration_docs_scaffold` green whenever changing Stage 2/Phase 3 language; do not claim Phase 3 multirate TreePM production closure without the hard-gate evidence in `docs/treepm_phase3_contract.md`. |

Command evidence for the documentation synchronization patch:

- `cmake --preset cpu-only-debug` → pass.
- `cmake --build --preset build-cpu-debug` → pass.
- `ctest --preset test-stage1-runtime-truth-cpu-debug --output-on-failure` → pass (17/17).
- `ctest --preset test-cpu-debug --output-on-failure` → pass (77/77).
- `cmake --preset hdf5-debug`, `cmake --build --preset build-hdf5-debug`, `ctest --preset test-hdf5-debug --output-on-failure` → pass (78/78).
- `cmake --preset pm-hdf5-fftw-debug`, `cmake --build --preset build-pm-hdf5-fftw-debug`, `ctest --preset test-pm-hdf5-fftw-debug --output-on-failure` → pass (78/78).
- `cmake --preset mpi-hdf5-fftw-debug` → blocked by missing `fftw3_mpi`; no MPI+HDF5+FFTW closure is claimed.

Remaining blockers are unchanged: P36 spawned-particle scheduler registration, P38 distributed scheduler identity exchange, and dependency-enabled MPI+HDF5+FFTW feature-path evidence are still required before broader closure claims.

## 2026-05-19 PM refresh stage-diagnostic update

No issue status changes in this patch. Open blockers remain P38 distributed scheduler identity exchange and dependency-enabled MPI+HDF5+FFTW feature-path evidence.
