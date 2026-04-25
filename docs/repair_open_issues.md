# Repair open issues (P01–P19 freeze ledger)

_Date captured: 2026-04-14 (UTC)_


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
| P0-05 | Closed | Gas-cell identity and hydro-state ownership invariants | Gas identity semantics and temporary 1:1 gas-particle/gas-cell coupling were insufficiently explicit under reorder/resize/restart/active-hydro extraction paths. | Added `test_unit_gas_cell_identity_invariants` plus restart mapping checks in `test_integration_restart_checkpoint_roundtrip`; added `debugAssertGasCellIdentityContract` and reorder guardrails preventing silent gas-relative reorders; documented temporary contract + forbidden assumptions in `docs/architecture/adr_runtime_truth_ownership.md`. |
| P0-04 | Closed | Species migration + canonical grouping invariants | Species transitions were under-tested for canonical grouping continuity across sidecars/count ledgers/offset index maps, per-particle softening override survival, and scheduler-bin/active-set continuity after migration + reorder/resize. | Added deterministic migration invariant suite in `test_integration_species_migration_invariants` (single/batch migration, reorder/resize/active-extraction adversarial sequences, sidecar attach-detach checks, timestep-bin continuity, and softening-priority checks); repaired migration pack/commit lanes to carry optional per-particle softening overrides; updated runtime ownership docs for allowed migration path + sidecar rules + species-agnostic active-eligibility policy. |
| P0-03 | Closed | Particle ordering/resize/reorder/sidecar stale-index invariants | Index-based assumptions could silently desynchronize particle identity from species tags, sidecar payload rows, softening overrides, and cached active/kernel views across resize and reorder operations. | Hardened `reorderParticles` to move full sidecar rows under `kMoveWithParent`, preserved optional softening sidecar values across particle resize, added index-generation invalidation checks for mutable particle/cell kernel views, and expanded integration coverage in `test_integration_reorder_compaction_sidecars` (`test_particle_resize_identity_invariants`, `test_particle_reorder_sidecar_invariants`, `test_particle_stale_view_invalidation`); validated with targeted and full CPU-debug ctest runs. |
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
| P31-STATE-TEST-DRIFT-026 | Open | Snapshot schema contract drift | `test_unit_snapshot_hdf5_schema` currently expects `gadget_arepo_v2` while runtime schema advertises `gadget_arepo_v4`, causing `ctest --preset test-cpu-debug --output-on-failure` failure during Stage 0 audit. | Align snapshot schema contract between implementation/tests/docs and re-run `test_unit_snapshot_hdf5_schema` plus `test_integration_snapshot_hdf5_roundtrip` with explicit migration/compatibility note if required. |


## 2026-04-20: Remaining Phase 3 campaign blockers

- Large-rank (beyond np2) strong/weak scaling certification sweep still required for Phase 3 closeout.
- External literature-target calibration for cosmological P(k) and halo profile systematics remains outside this stage; current campaign uses explicit in-repo reference contracts and envelopes.

## 2026-04-21: Phase 3 final integration closeout blockers

- `cmake --preset mpi-hdf5-fftw-debug` fails in this environment because `fftw3_mpi` is unavailable; distributed MPI+FFTW command bundle cannot be executed end-to-end in this cycle.
- `ctest --test-dir build/pm-hdf5-fftw-debug -R "integration_reference_workflow" --output-on-failure` fails with `runtime workflow schema compatibility validation failed`.
- `ctest --test-dir build/pm-hdf5-fftw-debug -R "integration_tree_pm_coupling_periodic" --output-on-failure` fails with `expected at least one low-res contaminant`.
- `ctest --test-dir build/pm-hdf5-fftw-debug -R "validation_phase2_mpi_gravity_single_rank" --output-on-failure` fails in communication-stress mode due residual-cutoff expectation mismatch.
- `validation/artifacts/research_grade/phase3/correctness/power_spectrum_consistency.json` remains blocked (`missing_diagnostics_inputs`), so Phase 3 correctness/force-accuracy envelope is not fully evidenced.
- Multi-rank scaling evidence in the audited artifact set remains baseline-only and non-certifying (`phase2_baseline_scaling_summary.json`); checked-in CSV artifacts currently cover `np1` only.

Phase 3 status in this cycle: **incomplete** (see `docs/treepm_phase3_closeout.md`).
