# Repair state recap (post-repair audit snapshot)

_Date captured: 2026-04-07 (UTC)_

This recap records **current command-backed audit evidence** for the emergency repair closeout pass.

## 0) Phase 2 distributed TreePM validation/gate suite wiring (2026-04-19 UTC)

Commands:

```bash
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug -j4 --target test_validation_phase2_mpi_gravity generate_mpi_gravity_scaling_artifacts
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -R "validation_phase2_mpi_gravity_two_rank|validation_phase2_mpi_gravity_three_rank"
```

Observed:

- Added an explicit MPI gravity validation gate (`validation_phase2_mpi_gravity_*`) covering:
  - distributed PM vs one-rank force reference (`rel_L2 <= 1e-10`),
  - distributed full TreePM vs one-rank force reference (`rel_L2 <= 5e-6`, `max_rel <= 5e-5`),
  - rank-count reproducibility sweep (`np=1,2,3`),
  - communication-stress path with tiny tree exchange batches,
  - MPI restart write/read continuation verification through reference workflow roundtrip.
- Added separate PM-only and tree-only MPI scaling artifact generators:
  - `validation/artifacts/pm_only_scaling_np{1,2}.csv`
  - `validation/artifacts/tree_only_scaling_np{1,2}.csv`
- CI now surfaces the MPI gravity gate in the MPI+HDF5+FFTW matrix row and uploads scaling artifacts.

Interpretation:

- Phase 2 evidence now has an auditable distributed gravity gate suite rather than relying on pseudo two-rank scaffolding alone.

## 0) TreePM Phase 2 distributed gravity performance hardening pass (2026-04-19 UTC)

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_pm_solver bench_pm_solver bench_tree_pm_coupling bench_parallel_decomposition_exchange
ctest --test-dir build/cpu-only-debug --output-on-failure -R unit_pm_solver
./build/cpu-only-debug/bench_pm_solver
./build/cpu-only-debug/bench_tree_pm_coupling
./build/cpu-only-debug/bench_parallel_decomposition_exchange
```

Observed:

- PM FFT plans and mesh-side scratch buffers remain cached by slab layout and are now paired with reused distributed interpolation send/receive buffers for force and potential reverse communication paths.
- TreePM short-range distributed exchange now reuses persistent payload/count buffers and overlaps communication with owner-local target evaluation by running local-local residual work while batched request exchange is in-flight.
- Active-set compact sidecar arrays in TreePM coupling now resize without per-step zero-fill churn.
- Bench output now includes repeated warmup/measured iteration fields and PM plan cache counters so reuse behavior is visible in profiler output instead of only one-shot timings.

Interpretation:

- Phase 2 distributed TreePM runtime overhead is reduced via explicit buffer/plan reuse and auditable overlap without changing force-split math or ownership semantics.

## 0) TreePM Phase 2 gravity-aware ownership decomposition + migration commit repair (2026-04-19 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory test_unit_simulation_state
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory|unit_simulation_state"
```

Observed:

- Morton SFC decomposition now supports explicit gravity-aware cost terms for:
  - owned particle count,
  - recent active target count,
  - recent remote tree interaction volume,
  - memory footprint (plus optional retained generic work term).
- Per-rank decomposition metrics now expose each component to make imbalance diagnostics auditable.
- `SimulationState` now provides explicit migration pack/commit boundaries:
  - packs hot lanes + metadata + species sidecars for migrating rows,
  - commits ownership only at one synchronization call,
  - rebuilds species counts and species index,
  - removes stale local ghost/import rows so post-commit ownership/sidecar state is unambiguous.

Interpretation:

- Distributed TreePM ownership has a real migration commit contract and gravity-aware decomposition input signal.
- One-rank baseline behavior is preserved when no migration is staged.

## 0) TreePM Phase 2 distributed workflow cadence-consensus + active-set wiring repair (2026-04-19 UTC)

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_integration_reference_workflow test_integration_reference_workflow_distributed_treepm_mpi test_integration_tree_pm_coupling_periodic
ctest --test-dir build/cpu-only-debug --output-on-failure -R "integration_reference_workflow"
./build/cpu-only-debug/test_integration_reference_workflow_distributed_treepm_mpi
```

Observed:

- The live gravity callback path now enforces rank-consensus cadence metadata for each gravity kick:
  - all ranks advance `gravity_kick_opportunity` every gravity kick stage,
  - refresh/reuse decisions are reduced and must agree on all ranks,
  - divergence fails loudly instead of silently drifting counters or field-version metadata.
- Early return on empty local active targets was removed so distributed PM cadence bookkeeping remains coherent even for ranks with no active targets on a given kick.
- Integration coverage now includes:
  - single-rank active-subset agreement against full solve in TreePM coupling tests,
  - distributed active-subset two-rank vs one-rank reference agreement in TreePM coupling tests,
  - distributed workflow MPI smoke test for cadence-record coherence and final digest agreement.

Interpretation:

- Distributed PM cadence and metadata contracts are now wired in the real reference workflow callback path.
- TreePM source/target split is preserved: PM refresh consumes all mass sources, while active-target export/import remains short-range-only work.

## 0) TreePM Phase 2 distributed PM interpolation return path repair (2026-04-19 UTC)

## 0) TreePM Phase 2 distributed short-range tree export/import repair (2026-04-19 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_integration_tree_pm_coupling_periodic
ctest --test-dir build/cpu-only-debug --output-on-failure -R "integration_tree_pm_coupling_periodic"
```

Observed:

- `TreePmCoordinator::evaluateShortRangeResidual` now runs an explicit active-target export/import protocol when `world_size>1`:
  - owner rank computes local-local residual,
  - owner exports target batches (bounded by `tree_exchange_batch_bytes`) to each peer,
  - peer evaluates requests against its local tree/source data and returns partial accelerations,
  - owner validates response coverage (`batch_token`, `request_id`) and accumulates remote partials.
- The workflow now wires `numerics.treepm_tree_exchange_batch_bytes` into `TreePmOptions`.
- MPI integration coverage was added (`integration_tree_pm_coupling_periodic_mpi_two_rank`) for two-rank distributed-vs-single-rank agreement including cutoff-boundary cross-rank peers.

Interpretation:

- Phase 2 short-range distributed TreePM now uses real peer participation instead of rank-local-only residual evaluation.
- The one-rank numerical contract remains the reference and is used as distributed comparison baseline.

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_integration_pm_periodic_mode
ctest --test-dir build/cpu-only-debug --output-on-failure -R "integration_pm_periodic_mode"
```

Observed:

- `PmSolver::interpolateForces` now supports slab-distributed reverse communication:
  particle-owner ranks send weighted stencil requests to slab owners, slab owners return weighted force contributions, and owners accumulate in local particle order.
- `PmSolver::interpolatePotential` uses the same reverse message contract for optional potential gather.
- Distributed PM integration coverage now includes CIC and TSC one-rank vs two-rank interpolation agreement with slab-boundary particle cases.

Interpretation:

- Phase 2 PM long-range path is now distributed through deposition, solve, and particle-force/potential return without replicated mesh gather assumptions.
- One-rank PM numerical contract remains the baseline reference and is used as the distributed comparison target.

## 0) TreePM Phase 2 distributed density assignment repair (2026-04-19 UTC)

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_pm_solver test_integration_pm_periodic_mode
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_pm_solver|integration_pm_periodic_mode"
```

Observed:

- `PmSolver::assignDensity` now supports slab-distributed owner-to-owner contribution routing with
  explicit slab ownership validation and no replicated full mesh requirement on each rank.
- One-rank path remains numerically consistent with prior contract and now also includes TSC wrap/mass checks.
- Integration tests include distributed density assignment agreement against one-rank reference for both CIC and TSC.

Interpretation:

- Phase 2 PM distributed deposition infrastructure is now implemented and tested at contract level.
- Distributed FFT solve support status is unchanged from prior repair state.

## 1) Preset and dependency gates

Commands:

```bash
cat CMakePresets.json
cmake --preset hdf5-debug
cmake --preset pm-hdf5-fftw-debug
```

Observed:

- Presets now include `hdf5-debug` and `pm-hdf5-fftw-debug`.
- Dependency failure behavior is explicit and actionable.
- In this environment, both dependency-enabled presets fail at configure with:
  - `COSMOSIM_ENABLE_HDF5=ON but HDF5 was not found.`

Interpretation:

- Preset-level closure for P05/P06 exists in-tree.
- Runtime closeout for HDF5/PM feature paths is blocked by missing dependency in this environment.

## 2) CPU-only baseline evidence

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4
ctest --preset test-cpu-debug --output-on-failure
```

Outcome:

- Configure: **PASS**
- Build: **PASS**
- Tests: **PASS** (`36/36`)

Interpretation:

- Baseline path remains healthy and no unrelated regressions were observed.

## 3) Guardrail scripts

Commands:

```bash
./scripts/ci/check_repo_hygiene.sh
./scripts/ci/guard_feature_paths.sh
```

Outcome:

- `check_repo_hygiene.sh`: **PASS**.
- `guard_feature_paths.sh`: CPU segment passes; feature segment fails at HDF5 configure gate, as expected when dependency is absent.

Interpretation:

- Hygiene and process checks are active and correctly prevent false success claims from CPU-only evidence.

## 4) IC reader gas import state evidence

Command:

```bash
nl -ba src/io/ic_reader.cpp | sed -n '562,592p'
```

Observed:

- For `PartType0`, IC import now writes:
  - `gas_cells.density_code[cell_i]`
  - `gas_cells.internal_energy_code[cell_i]`
- Unsupported reports are for unmapped metallicity/smoothing length, not thermodynamic bypass.

Related test artifact:

- `tests/unit/test_ic_reader.cpp` includes HDF5-gated checks for gas thermodynamic mapping and optional-density defaulting behavior.

Interpretation:

- P04 implementation appears closed in code/tests.
- Full runtime verification remains dependent on HDF5 availability.

## 5) Docs and naming/process hygiene

Commands:

```bash
test -f docs/build_instructions.md
test -f CONTRIBUTING.md
find . -maxdepth 1 -type f | sed 's#^./##' | sort
```

Outcome:

- `docs/build_instructions.md`: present.
- `CONTRIBUTING.md`: present.
- Root file naming is now policy-safe (no legacy violating spreadsheet at repo root).

Interpretation:

- P08/P09/P10 are closed at repository state level.

## 6) Audit conclusion pointer

See `docs/repair_closeout_report.md` for the stop/go decision.
Current state is **STOP** until HDF5 (and then PM/HDF5/FFTW) preset build+test commands pass on the intended feature paths.

## 7) Core boundary repair (reference workflow assembly)

Commands:

```bash
ctest --test-dir build/cpu-only-debug --output-on-failure -R "integration_core_dependency_direction|integration_reference_workflow"
```

Observed:

- Core dependency-direction guard passes and fails fast on forbidden upward includes in `include/cosmosim/core/**` and `src/core/**`.
- Reference workflow integration smoke test still passes after moving assembly to `workflows/`.

Interpretation:

- The architecture boundary leak from `core` into analysis/I/O/physics workflow assembly is repaired while keeping smoke behavior intact.

## 8) Configuration-contract hardening (typed freeze path)

Commands:

```bash
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_config_parser|unit_simulation_mode|integration_simulation_mode_toy_runs"
```

Observed:

- Policy-like config values are parsed from param strings, then frozen into typed enums for solver selection, coordinate frame, mode boundaries, and feedback mode/variant.
- Unknown key handling and deprecated alias mapping use a centralized key/alias registry in `src/core/config.cpp`.
- Unit/integration coverage validates unknown key rejection, alias behavior, invalid enum rejection, and deterministic normalized config/hash behavior.

Interpretation:

- String-driven runtime policy drift is reduced after freeze while preserving param-style UX and deterministic normalized-config provenance semantics.

## 9) SimulationState ownership/invariants decomposition repair

_Date captured: 2026-04-13 (UTC)_

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_simulation_state
./build/cpu-only-debug/test_unit_simulation_state
```

Observed:

- `src/core/simulation_state.cpp` responsibility concentration was split by invariant boundary:
  - storage-lane consistency (`simulation_state_structures.cpp`)
  - ownership checks (`simulation_state_ownership.cpp`)
  - species indexing and transfer packing (`simulation_state_species.cpp`)
  - metadata/module-sidecar serialization (`simulation_state_metadata.cpp`)
  - active-view assembly/scatter (`simulation_state_active_views.cpp`)
  - reorder/scratch logic remains in `simulation_state.cpp`
- Header-level hot-field contract notes now explicitly name allowed gravity/hydro active-view lanes.
- `tests/unit/test_simulation_state.cpp` now covers:
  - ownership invariant failure/recovery
  - unique-ID invariant failure/recovery
  - species-index rebuild correctness
  - transfer packet pack/unpack-equivalence checks
  - metadata serialize/deserialize round-trip
  - active-view hot-field writeback behavior and cold-lane protection
  - static assertion guardrails for gravity/hydro active-view compactness.

Interpretation:

- State ownership and active-view invariants are now reviewable in focused files with explicit hot/cold contracts and command-backed checks.

## 10) Snapshot/restart I/O contract boundary hardening

_Date captured: 2026-04-13 (UTC)_

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_snapshot_hdf5_schema test_unit_restart_checkpoint_schema test_integration_restart_checkpoint_roundtrip
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_snapshot_hdf5_schema|unit_restart_checkpoint_schema|integration_restart_checkpoint_roundtrip"
```

Observed:

- Shared continuation-metadata contract names/validation are centralized via `include/cosmosim/io/io_contract.hpp` + `src/io/internal/io_contract.cpp`.
- Restart docs and API now publish an exact-restart checklist (`exactRestartCompletenessChecklist()`), separating restart-completeness obligations from snapshot interoperability.
- New negative checks cover restart schema mismatch rejection, missing required scheduler dataset rejection, and finalize-failure behavior that leaves target path untouched while preserving the temporary artifact for diagnosis.

Interpretation:

- Snapshot and restart responsibilities are now explicitly separated and reviewable, with stronger error behavior on continuation-critical contract violations.

## 11) Operational observability repair (structured run events)

_Date captured: 2026-04-13 (UTC)_

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_profiling test_integration_reference_workflow
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_profiling|integration_reference_workflow"
```

Observed:

## 12) Phase-2 PM distributed FFT infrastructure repair state

_Date captured: 2026-04-19 (UTC)_

Commands:

```bash
cmake -S . -B build/mpi-fftw-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCOSMOSIM_ENABLE_TESTS=ON -DCOSMOSIM_ENABLE_MPI=ON -DCOSMOSIM_ENABLE_FFTW=ON -DCOSMOSIM_ENABLE_HDF5=OFF
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_pm_solver test_integration_pm_periodic_mode
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_pm_solver|integration_pm_periodic_mode"
```

Observed:

- `PmSolver::solvePoissonPeriodic` now accepts slab-owned local PM grids and uses communicator-aware FFTW MPI plans when MPI+FFTW are enabled.
- PM FFT plan/scratch buffers are cached by slab layout key and reused across calls; unit/integration checks now include plan cache reuse behavior.
- Single-rank PM periodic checks remain green on CPU-only debug path.
- This environment cannot complete MPI validation because CMake cannot find MPI (`Could NOT find MPI_CXX`).

Interpretation:

- Distributed PM FFT ownership/message contract is implemented in-tree with command-backed single-rank regression evidence.
- Multi-rank runtime closure remains blocked in this container by missing MPI runtime/development toolchain.

## 12) Infrastructure gate-bundle enforcement hardening

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
bash -n scripts/ci/enforce_infra_gates.sh
```

Observed:

- Added `scripts/ci/enforce_infra_gates.sh` to run three explicit infrastructure gates with command failure propagated per path:
  - `cpu_core_boundary_and_config_contract`
  - `hdf5_schema_and_exact_restart_contract`
  - `pm_hdf5_fftw_feature_path_validation`
- Each gate writes artifacts under `ci_artifacts/infrastructure_gates/<gate_id>/`.
- A machine-readable status file, `infrastructure_gate_report.json`, now records gate ids, presets, test regexes, and pass/fail status for CI interpretation.
- `.github/workflows/ci.yml` now includes a dedicated `infrastructure_gates` job and makes `reproducibility_gate` depend on both matrix coverage and the explicit gate bundle.
- `tests/integration/test_core_dependency_direction.cmake.in` now checks both source/header include direction and `CMakeLists.txt` target-link direction for `cosmosim_core`.

## 13) Phase 2 distributed TreePM contract freeze surface

_Date captured: 2026-04-19 (UTC)_

Commands:

```bash
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug -j4
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure -R "unit_config_parser|integration_docs_scaffold|integration_provenance_roundtrip|unit_parallel_distributed_memory"
```

Observed:

- Typed config now freezes Phase 2 distributed gravity controls:
  - `numerics.treepm_pm_decomposition_mode` (`slab` only),
  - `numerics.treepm_tree_exchange_batch_bytes` (`>0`).
- Provenance now records those controls for audit continuity.
- New architecture contract doc exists at `docs/treepm_phase2_distributed_contract.md`.
- CI matrix now includes `mpi-hdf5-fftw-debug`, so distributed gravity development is not only gated through optional `mpi-release` smoke.

Interpretation:

- Contract/documentation/build surfaces for Phase 2 are now explicit without claiming distributed TreePM algorithm completion.

Interpretation:

- CPU-only success is no longer treated as implicit feature-path closure because HDF5 and PM/HDF5/FFTW are evaluated in one explicit gate bundle with a consolidated status report.
- Architecture-boundary enforcement now covers both include-level and build-graph-level dependency direction.

## 12) Parallel contract hardening follow-up (R07)

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
cmake --build --preset build-cpu-debug -j4
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory|integration_parallel_two_rank_restart"
```

Observed:

- Distributed-memory transfer roles now expose typed `outbound_transfers` / `inbound_transfers` descriptors in addition to legacy send/recv index vectors.
- Plan/buffer invariants were tightened for invalid ownership combinations, zero-byte ghost payload shape, and decode shape checks.
- Reduction agreement fields now distinguish deterministic baseline vs measured sum, plus a policy helper for absolute/relative tolerance gating.
- Config consensus mismatch artifacts include property-level mismatch rows with baseline/offending rank-value pairs.

Interpretation:

- Pseudo-multi-rank contract clarity and diagnosability improved for review/CI artifacts without adding new MPI features or claiming production multi-rank closure.

## 13) Parallel contract hardening follow-up (R07c)

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory test_integration_parallel_two_rank_restart
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory|integration_parallel_two_rank_restart"
```

Observed:

- Ghost transfer descriptors now encode explicit lifecycle intent (`ghost refresh request`, `ghost refresh receive staging`) and post-transfer residency expectation, with reserved typed migration intents documented as out-of-scope scaffolding placeholders.
- `validateGhostExchangePlan()` now enforces neighbor-slot/peer-rank role correctness and exact descriptor-index equality against canonical send/recv vectors.
- Reduction agreement checks now use explicit mode selection (`absolute only`, `relative only`, `absolute and relative`, `absolute or relative`) rather than implicit boolean logic.
- Unit coverage now includes plan-drift rejection checks and per-policy-mode reduction agreement assertions.

Interpretation:

- Pseudo-multi-rank distributed-memory contracts are now more reviewable and diagnosable without claiming production MPI transfer correctness or migration commit semantics.

## 14) Diagnostics maturity-tier repair (analysis honesty/scalability guard)

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_config_parser test_unit_analysis_diagnostics test_integration_analysis_bundle
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_config_parser|unit_analysis_diagnostics|integration_analysis_bundle"
```

Observed:

- Analysis config now carries typed `diagnostics_execution_policy` with normalized-config persistence.
- Diagnostics bundles now emit per-diagnostic maturity metadata (`tier`, `maturity`, `scalability`, execution policy).
- Provisional heavy reference diagnostics (power spectrum direct summation) no longer run under default policy.
- Unit/integration tests now assert that heavy provisional diagnostics only run when explicitly opted in.

Interpretation:

## 13) Distributed-memory ownership/determinism contract hardening

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory test_integration_parallel_two_rank_restart
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory|integration_parallel_two_rank_restart"
```

Observed:

- Distributed-memory ghost planning now has an explicit typed ownership contract (`LocalGhostDescriptor` with `kOwned`/`kGhost`) and validation for invalid ownership combinations.
- Deterministic reduction agreement helpers now provide explicit deterministic reference sums and absolute/relative agreement reporting for reproducibility checks.
- Multi-rank config-freeze consensus checks now compare normalized config hash + rank-count expectation + deterministic-reduction mode across per-rank digests.
- Unit/integration tests cover typed ownership planning, deterministic reduction agreement checks, and config-consensus checks using pseudo-multi-rank vectors.

Interpretation:

- Reviewability and contract clarity improved for rank-owned vs ghost state and reduction determinism semantics without expanding scope into new decomposition/scaling features.
- Evidence remains CPU-only/pseudo-multi-rank in this environment; no claim is made that this alone closes full production MPI execution.

- Infrastructure run-health counters remain first-class and cheap.
- Validated light science diagnostics remain available in default runs.
- Reference/provisional heavy diagnostics are quarantined behind explicit non-default policy.

- `core::ProfilerSession` now carries a minimal structured runtime event model (kind, severity, subsystem, optional step/time/scale context, message, key/value payload).
- Reference workflow emits a machine-readable operational report (`reference_operational_events.json`) linked to deterministic config provenance via `provenance_config_hash_hex`.
- Key infrastructure lifecycle/failure surfaces are explicit in event records (config freeze validation, restart/snapshot write/read begin/complete/failure).

Interpretation:

- Operational troubleshooting and CI artifact review no longer depend on ad hoc text/exception surfaces alone.
- Reproducibility posture remains unchanged: operational events are additive and provenance-linked; no solver behavior or restart/snapshot schema semantics were changed.

## 13) Infrastructure gate hardening follow-up (R08 durability/diagnostics)

_Date captured: 2026-04-14 (UTC)_

Commands:

```bash
bash -n scripts/ci/run_preset_pipeline.sh
bash -n scripts/ci/enforce_infra_gates.sh
cmake --preset cpu-only-debug
ctest --test-dir build/cpu-only-debug --output-on-failure -R "integration_core_dependency_direction"
bash ./scripts/ci/enforce_infra_gates.sh ci_artifacts/local_infra_gates
```

Observed:

- Infrastructure gates are now defined in `scripts/ci/infrastructure_gates_manifest.tsv` with explicit fields for gate id, preset trio, intended test scope, metadata expectation, artifact subdir, and benchmark toggle.
- `run_preset_pipeline.sh` now emits `preset_pipeline_report-<preset>.json` with per-phase command metadata and failure phase (`configure`, `build`, `test`, `artifact_collection`, `benchmark`) to improve triage fidelity.
- `enforce_infra_gates.sh` now records richer per-gate entries in `infrastructure_gate_report.json`, including `build_preset`, `test_preset`, `test_scope`, `artifact_dir`, `failed_phase`, and explicit phase commands.
- `integration_core_dependency_direction` now validates target-link direction by inspecting CMake-generated Graphviz dependency metadata from the configured build tree instead of scanning only top-level `CMakeLists.txt` text.

Interpretation:

- Gate semantics are less rename-fragile and more reviewable due to a narrow declarative manifest.
- CI failure triage can identify failing phase and command per gate without log archaeology.
- Core target-link direction guard is stronger than single-file regex scanning while remaining narrow and infrastructure-focused.


## 11) 2026-04-14 targeted infrastructure repair pass

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j1 --target   test_unit_config_parser   test_unit_simulation_state   test_unit_parallel_distributed_memory   test_unit_restart_checkpoint_schema
ctest --test-dir build/cpu-only-debug --output-on-failure -R   "unit_config_parser|unit_simulation_state|unit_parallel_distributed_memory|unit_restart_checkpoint_schema"

cmake --preset hdf5-debug
ninja -C build/hdf5-debug -j8   test_unit_restart_checkpoint_schema   test_integration_snapshot_hdf5_roundtrip   test_integration_restart_checkpoint_roundtrip
ctest --test-dir build/hdf5-debug --output-on-failure -R   "unit_restart_checkpoint_schema|integration_snapshot_hdf5_roundtrip|integration_restart_checkpoint_roundtrip"

cmake --preset pm-hdf5-fftw-debug
```

Observed:

- CPU targeted repair tests pass after hardening config round-trip/hash semantics, species-sidecar coverage invariants, descriptor-only ghost-plan honesty, and distributed restart decode completeness checks.
- HDF5 targeted repair tests pass after fixing restart stellar-sidecar completeness, restart payload hash coverage, continuation metadata cross-checking, and `Header/MassTable` snapshot-mass fallback.
- `pm-hdf5-fftw-debug` configure is blocked in this environment because FFTW3 development files are unavailable.

Interpretation:

- The repaired invariants are command-backed on CPU and HDF5 feature paths.
- PM/FFTW path remains an environment blocker here, not a demonstrated code regression in this repair pass.

## 14) 2026-04-19 Phase 2 PM slab ownership/storage milestone

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory test_unit_pm_solver
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory|unit_pm_solver"
```

Observed:

- Added explicit PM slab ownership typing (`parallel::PmSlabLayout`, `PmSlabRange`) with deterministic uneven partitioning, owner-rank lookup, and validated global/local index conversions.
- `gravity::PmGridStorage` now accepts explicit slab layout and allocates only local slab storage (`local_nx * Ny * Nz`), while default construction remains the one-rank full-domain slab.
- PM solver entry points now fail fast on partial slabs, preventing pseudo-distributed use until distributed FFT and remote deposition/gather are implemented.
- Unit tests now cover uneven slab partitioning, ownership queries, index round-trips, one-rank slab equivalence, and explicit rejection of partial-slab use on the single-rank solver path.

Interpretation:

- This is an infrastructure-only Phase 2 milestone: ownership/storage contracts are now auditable without claiming distributed PM algorithm completion.
- Reproducibility posture is preserved for one-rank runs because default PM storage maps to the same full-domain indexing contract.


## 15) 2026-04-19 MPI/CUDA runtime topology milestone

Commands:

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug -j4 --target test_unit_parallel_distributed_memory test_unit_config_parser
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_parallel_distributed_memory|unit_config_parser"
```

Observed:

- Added explicit runtime topology assembly (`parallel::DistributedExecutionTopology`) that binds together MPI world size/rank, PM slab ownership, and rank-local CUDA device assignment.
- Added `core::CudaRuntimeInfo` and explicit CUDA device-selection helpers so `parallel.gpu_devices` becomes a real runtime contract rather than a dead config field.
- Reference workflow TreePM initialization now validates `parallel.mpi_ranks_expected` against the runtime world and enables the CUDA PM assignment/interpolation path only when the runtime request is valid.
- GPU requests now fail loudly when CUDA was requested but no visible devices are available; this avoids silent CPU fallback drift in distributed runs.

Interpretation:

- This remains infrastructure work, not a claim of distributed PM FFT or MPI+GPU overlap completion.
- The code path is now materially closer to an operational multi-rank/multi-device TreePM bring-up because rank-local execution intent is explicit and auditable.

## 16) 2026-04-19 distributed restart/provenance continuation-safety milestone

Commands:

```bash
cmake --build --preset build-cpu-debug -j4 --target test_unit_restart_checkpoint_schema test_unit_parallel_distributed_memory test_integration_restart_checkpoint_roundtrip test_integration_parallel_two_rank_restart
ctest --test-dir build/cpu-only-debug --output-on-failure -R "unit_restart_checkpoint_schema|unit_parallel_distributed_memory|integration_restart_checkpoint_roundtrip|integration_parallel_two_rank_restart"
```

Observed:

- Restart schema bumped to `cosmosim_restart_v4` with explicit `/distributed_gravity/state` payload (`parallel::DistributedRestartState`, schema_version=2).
- Distributed continuation payload now persists decomposition epoch, owning-rank table, PM slab layout metadata, gravity-kick cadence state, long-range field refresh/version metadata, and explicit `long_range_restart_policy`.
- Policy is explicit and enforced: `deterministic_rebuild` (cached long-range PM field arrays are not serialized; continuation rebuilds deterministically on cadence refresh).
- Provenance schema bumped to `provenance_v3` to include distributed restart diagnostics (epoch/world/grid/slab signature/kick opportunity/field version/restart policy).
- Added typed compatibility diagnostics (`evaluateDistributedRestartCompatibility`) plus negative mismatch coverage.

Interpretation:

- This is an infrastructure-repair continuation-safety improvement for Phase 2 distributed TreePM, not a new physics model.
- Reproducibility posture remains explicit: restart carries deterministic policy + auditable rank/layout metadata, and integrity hashing now includes distributed continuation state.
