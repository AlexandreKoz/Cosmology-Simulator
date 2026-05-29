# Stage 8 restart/output maturity map

This note records repository reconnaissance for Stage 8. It is an audit artifact only; it does not change production code, public APIs, schema versions, solvers, or numerics.

## 1. Repository inspection summary

Inspected restart, snapshot, scheduler, workflow, PM cadence, provenance/config, stochastic source-term, and existing test paths using `find`, `rg`, and direct `sed` reads. The current code already separates ordinary GADGET/AREPO-style snapshots from CosmoSim restart checkpoints: restart I/O uses `cosmosim_restart_v11` in `include/cosmosim/io/restart_checkpoint.hpp` and `src/io/restart_checkpoint.cpp`, while snapshot I/O uses `gadget_arepo_v4` in `include/cosmosim/io/snapshot_hdf5.hpp` and `src/io/snapshot_hdf5.cpp`.

Key architectural fact: `SimulationState::particles.time_bin` and `SimulationState::cells.time_bin` are diagnostic/compatibility mirrors, not scheduler authority. Exact continuation authority is `HierarchicalTimeBinScheduler` plus `IntegratorState`.

## 2. Restart writer/reader map

Restart write entry points:

- Public writer: `cosmosim::io::writeRestartCheckpointHdf5` in `src/io/restart_checkpoint.cpp`.
- Workflow call site: `maybeWriteOutputs` in `src/workflows/reference_workflow.cpp` builds `RestartWritePayload` and writes `report.restart_path`.
- Python bindings do not currently expose restart write/read; snapshot read/write is exposed separately.

Restart read entry points:

- Public reader: `cosmosim::io::readRestartCheckpointHdf5` in `src/io/restart_checkpoint.cpp`.
- Workflow readback check: `maybeWriteOutputs` immediately reads the just-written checkpoint and checks selected equality/compatibility fields.
- Integration/unit tests read checkpoints directly.

Restart HDF5 root attributes/datasets:

- Root attributes: `restart_schema_name`, `restart_schema_version`, `normalized_config_hash_hex`, `payload_integrity_hash_hex`, `payload_integrity_hash`.
- Root datasets through `sharedIoContractNames()`: normalized config text and serialized provenance record.
- `/state` group, including `particles`, `particle_sidecar`, `cells`, `gas_cells`, `patches`, `species_count_by_species`, `star_particles`, `black_holes`, `tracers`, `state_metadata`, `module_sidecars`, and `module_sidecar_names`.
- `/integrator` attributes: timeline, last KDK factors, `step_index`, `scheme`, boundary kind flags, time-bin context, `pm_long_range_field_valid`, and PM sync attributes.
- `/scheduler` attributes/datasets: `current_tick`, `max_bin`, `bin_index`, `next_activation_tick`, `active_flag`, `pending_bin_index`.
- `/distributed_gravity/state` serialized `parallel::DistributedRestartState`.

Persistent state serialized:

- Persistent particle hot lanes: positions, velocities, masses, time-bin mirror.
- Particle sidecar lanes: IDs, SFC keys, species, flags, owning rank, drift epoch, gravity softening values and override mask.
- Cell/gas/patch lanes: cell centers/mass/time-bin mirror/patch index, gas identity and thermodynamics, patch topology.
- Species ledgers and species-specific sidecars for stars, black holes, tracers.
- Opaque `ModuleSidecarRegistry` blocks with per-module schema versions.
- Integrator timeline/KDK boundary/PM cadence state.
- Scheduler persistent state.
- Distributed gravity decomposition and deterministic long-range restart policy.
- Normalized config and provenance.
- Payload integrity hash.

Restart validation already present:

- Writer requires non-null state, integrator, and scheduler; validates normalized config/provenance; validates `SimulationState::validateOwnershipInvariants()`; validates time-bin mirrors against scheduler state; and validates distributed gravity restart shape and `long_range_restart_policy == "deterministic_rebuild"` through the payload hash path.
- Reader requires exact schema version, required datasets/attributes through read helpers, validates mirrors against scheduler, rebuilds mirrors from scheduler, imports PM sync state through `PmSynchronizationState::importPersistentState`, validates ownership invariants, and verifies payload hash.

Notable gaps/risks:

- `writeRestartCheckpointHdf5` does not directly call `core::assertCanWriteCheckpointAtBoundary`; it trusts the workflow boundary discipline. Direct callers can write an unsafe half-step if they pass `IntegratorState::inside_kdk_step = true` unless hash/import validation happens to reject unrelated state.
- Workflow restart readback proves a file roundtrip and selected equality, not scientifically equivalent continuation over future steps.
- `restartRuntimeStateExactlyEquivalent` does not compare particle drift epoch lanes even though restart serializes them and the payload hash covers them.
- `RestartReadResult` returns `scheduler_state` but not an imported `HierarchicalTimeBinScheduler`; callers must reconstruct/import exactly and must not infer bins from `SimulationState::time_bin` mirrors.

## 3. Snapshot/output writer/reader map

Snapshot write entry points:

- Public writer: `cosmosim::io::writeGadgetArepoSnapshotHdf5` in `src/io/snapshot_hdf5.cpp`.
- Workflow call site: `maybeWriteOutputs` writes `report.snapshot_path` before optional restart checkpoint.
- Python bindings expose snapshot write through `src/python_bindings.cpp`.

Snapshot read entry points:

- Public reader: `cosmosim::io::readGadgetArepoSnapshotHdf5`.
- Workflow readback check only verifies particle count equality.
- Python bindings expose snapshot read.

Snapshot HDF5 groups/datasets/attributes:

- `/Header` with GADGET/AREPO-compatible header arrays plus `CosmoSimSchemaName`, `CosmoSimSchemaVersion`, and `CosmoSimBuild`.
- `/Config` attribute `normalized` with normalized config text.
- `/Provenance` attributes for schema/config/build/config hashes/raw/normalized/derived runtime state, TreePM geometry/cadence, softening, and PM FFT backend.
- `/PartType0` through `/PartType5` with canonical `Coordinates`, `Velocities`, `Masses`, `ParticleIDs`; optional `GravitySofteningComoving` and `GravitySofteningOverrideMask`; tracer datasets under `PartType3`; hard-link alias groups when enabled.

Snapshot semantics:

- Snapshot files are science/output artifacts and are not restart files. The snapshot reader initializes particle `time_bin` to zero and does not restore scheduler, integrator, PM sync, module sidecars, gas thermodynamics, output cadence, or restart metadata.
- Snapshot reader allows compatibility defaults for aliases, generated IDs when permitted, mass table fallback, optional softening, optional tracer fields, optional config/provenance groups.

Output boundary behavior:

- `OutputBoundaryCallback::onStage` returns silently when no output is pending or the boundary is not output/restart safe.
- `maybeWriteOutputs` also returns false for disabled/not-due output, step zero, or unsafe last boundary before calling `assertCanWriteSnapshotAtBoundary` / `assertCanWriteCheckpointAtBoundary`.
- This means workflow output can be silently skipped at unsafe boundaries rather than failing loudly; direct core boundary assertions do fail loudly.

## 4. Scheduler/KDK/PM/RNG/source-term ownership map

Scheduler ownership:

- `current_tick`: owned by `HierarchicalTimeBinScheduler::m_current_tick`; exported/imported via `TimeBinPersistentState::current_tick`.
- `bin_index`: owned by `HierarchicalTimeBinScheduler::m_hot.bin_index`; mirrored into `SimulationState::particles.time_bin` and `SimulationState::cells.time_bin` for diagnostics/I/O only.
- `next_activation_tick`: owned by `HierarchicalTimeBinScheduler::m_hot.next_activation_tick`.
- `active_flag`: owned by `HierarchicalTimeBinScheduler::m_hot.active_flag`; set during `beginSubstep`, cleared during `endSubstep`, and invalid outside open substeps.
- `pending_bin_index`: owned by `HierarchicalTimeBinScheduler::m_hot.pending_bin_index`; set by requested/candidate bin transitions and applied in `endSubstep`.
- Active sets: owned as derived, transient `m_active_elements` and `ActiveSetDescriptor`; not restart truth.
- Candidate transition scratch (`m_candidate_bin_index`, `m_candidate_source`, `m_candidate_label`), bin membership mirrors (`m_elements_by_bin`, `m_position_in_bin`), diagnostics, and open-substep flag are rebuilt or reset on scheduler import and are not serialized.

KDK/timeline ownership:

- `IntegratorState` owns current time/scale/redshift/Hubble rate, last drift/kick factors, `step_index`, `current_boundary_kind`, `last_completed_boundary_kind`, `inside_kdk_step`, and `last_completed_restart_safe`.
- `StepOrchestrator::executeSingleStep` sets `inside_kdk_step = true` before stage dispatch and resets it only after all KDK stages complete, then commits the timeline.
- `assertCanWriteSnapshotAtBoundary` and `assertCanWriteCheckpointAtBoundary` reject open KDK steps or unsafe last boundaries.

PM cadence and force validity:

- `IntegratorState::pm_sync_state` owns restartable PM cadence truth through `PmSynchronizationState`.
- Serialized PM fields are `cadence_steps`, `gravity_kick_opportunity`, `last_refresh_opportunity`, `field_version`, `last_refresh_step_index`, `last_refresh_scale_factor`, `refresh_commit_pending`, `pending_refresh_opportunity`, and `pending_refresh_field_version`.
- `IntegratorState::pm_long_range_field_valid` is serialized separately.
- TreePM solver buffers/long-range field arrays are not serialized; distributed restart policy is `deterministic_rebuild`.
- Direct restart after a pending PM refresh commit would be dangerous unless represented and reloadable; current serialization includes pending PM metadata, but restart write does not reject unsafe KDK phases globally.

RNG/stochastic source ownership:

- Star formation stochastic draws are stateless hash-based SplitMix64 values from `physics.sf_random_seed`, `step_index`, `cell_index`, and `rank_local_seed_offset`.
- Stellar feedback stochastic events are stateless hash-based SplitMix64 values from `physics.fb_random_seed`, star index, and step seed.
- No production `std::mt19937`/engine state was found outside tests. There is therefore no engine stream to serialize, but exact restart depends on preserving `step_index`, stable indices/IDs, rank-local offset semantics, and config seeds.
- Config seeds are `PhysicsConfig::sf_random_seed` and `PhysicsConfig::fb_random_seed`; deterministic reduction switch is `ParallelConfig::deterministic_reduction`.

Hydro/source-term persistent state:

- Gas thermodynamics live in `SimulationState::gas_cells` and are serialized in restart.
- Stellar evolution/feedback cumulative lanes live in `StarParticleSidecar` and module sidecars; restart serializes both.
- Black-hole and tracer persistent sidecars are serialized.
- Hydro scratch, PM scratch, compact active gravity views, workspace scratch, MPI staging buffers, and output scratch are not in `RestartPersistentStateView` and must remain excluded.

## 5. Existing tests and coverage gaps

Relevant tests:

- `tests/unit/test_restart_checkpoint_schema.cpp`: schema name/version, compatibility, payload root shape, checklist contents, integrity hash sensitivity to some state/integrator/distributed metadata, and missing payload/provenance validation.
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`: HDF5 restart roundtrip, scheduler persistent fields including `pending_bin_index`, PM sync fields, module/species/softening/gas identity, active-set reconstruction around restart, and required-field failures.
- `tests/integration/test_snapshot_hdf5_roundtrip.cpp`: GADGET/AREPO snapshot schema, aliases, config/provenance, softening/tracer fields, roundtrip count/state subset checks.
- `tests/integration/test_hierarchical_time_bins.cpp`: scheduler begin/end behavior, active bins, transitions, active flags, tick advance.
- `tests/integration/test_hierarchical_timestep_regression.cpp`: adaptive time-bin regression and scheduler/mirror invariants through workflow-like paths.
- `tests/integration/test_parallel_two_rank_restart.cpp`: distributed decomposition/restart metadata serialization and deterministic reduction checks.

Coverage gaps for Stage 8:

- No continuation-equivalence test that runs N steps uninterrupted versus checkpoint/restart/resume and compares final state with strict tolerances.
- No direct test that `writeRestartCheckpointHdf5` rejects `inside_kdk_step` or unsafe boundary states.
- No test that output cadence state survives a real resume; cadence currently appears as pending-output workflow state, not a restartable typed contract.
- Snapshot roundtrip is intentionally not restart-equivalence coverage.
- PM deterministic rebuild/resume coverage is metadata-heavy; it should be extended to future-step equivalence when HDF5/FFTW are available.
- Stochastic modules rely on stateless seeds; tests should pin that restarts preserve exact stochastic outcomes when star formation/feedback is enabled.

## 6. Stage 8 task-by-task changed-file target map

Prompt 8.1 — Restart boundary hardening:

- Targets: `include/cosmosim/io/restart_checkpoint.hpp`, `src/io/restart_checkpoint.cpp`, `src/core/time_integration.cpp`, `tests/unit/test_restart_checkpoint_schema.cpp`, `tests/integration/test_restart_checkpoint_roundtrip.cpp`, `docs/restart_checkpointing.md`, `docs/output_schema.md`.
- Goal: make unsafe `IntegratorState` a hard restart writer error and document boundary contract without changing solver behavior.

Prompt 8.2 — Scheduler resume contract:

- Targets: `include/cosmosim/core/time_integration.hpp`, `src/core/time_integration.cpp`, `src/io/restart_checkpoint.cpp`, restart tests.
- Goal: provide explicit helper(s) to import `RestartReadResult::scheduler_state` into `HierarchicalTimeBinScheduler`, prove no heuristic reconstruction from mirrors, and validate `current_tick/bin_index/next_activation_tick/active_flag/pending_bin_index` exactly.

Prompt 8.3 — Integrator/KDK continuation validation:

- Targets: `src/io/restart_checkpoint.cpp`, `src/core/time_integration.cpp`, continuation tests.
- Goal: validate `current_boundary_kind`, `last_completed_boundary_kind`, `inside_kdk_step`, `last_completed_restart_safe`, and timeline fields on load/write with actionable messages.

Prompt 8.4 — PM cadence/force-validity restart equivalence:

- Targets: `include/cosmosim/core/time_integration.hpp`, `src/core/time_integration.cpp`, `src/io/restart_checkpoint.cpp`, `src/gravity/tree_pm_coupling.cpp`, PM restart tests.
- Goal: prove `pm_sync_state` and `pm_long_range_field_valid` are sufficient for deterministic rebuild and cadence continuation; reject unreloadable pending refresh states if they cannot be resumed safely.

Prompt 8.5 — Output cadence persistent contract:

- Targets: `include/cosmosim/workflows/reference_workflow.hpp`, `src/workflows/reference_workflow.cpp`, restart schema/doc/tests if cadence becomes restart truth.
- Goal: replace transient `PendingOutputBoundary` assumptions with typed output cadence state if resumed runs must preserve next output/checkpoint decisions exactly.

Prompt 8.6 — Stochastic source restart equivalence:

- Targets: `src/physics/star_formation.cpp`, `src/physics/stellar_feedback.cpp`, source-term tests, restart continuation tests.
- Goal: document and test stateless hash RNG contract; ensure stable identifiers/step/rank inputs are sufficient, or add explicit persistent stochastic module sidecar if not.

Prompt 8.7 — Snapshot/restart boundary separation:

- Targets: `src/io/snapshot_hdf5.cpp`, `src/io/restart_checkpoint.cpp`, `docs/output_schema.md`, `docs/restart_checkpointing.md`, snapshot tests.
- Goal: strengthen docs/tests that snapshots are not restarts and optional snapshot defaults must not leak into restart continuation.

Prompt 8.8 — End-to-end restart continuation test harness:

- Targets: new/updated integration tests under `tests/integration`, `CMakeLists.txt`, maybe `src/workflows/reference_workflow.cpp` only if a narrow resume entry point is needed.
- Goal: strict uninterrupted-vs-restarted comparison over scheduler, KDK, PM cadence, stochastic source, and output cadence state.

Prompt 8.9 — Parallel restart compatibility and rank output boundaries:

- Targets: `src/parallel/distributed_memory.cpp`, `include/cosmosim/parallel/distributed_memory.hpp`, `src/workflows/reference_workflow.cpp`, `tests/integration/test_parallel_two_rank_restart.cpp`, MPI validation tests.
- Goal: prove rank-qualified restart names, decomposition metadata, deterministic reductions, and PM slab metadata resume safely.

Prompt 8.10 — Documentation and release-gate alignment:

- Targets: `docs/restart_checkpointing.md`, `docs/output_schema.md`, `docs/time_integration.md`, `docs/architecture/runtime_truth_map.md`, `docs/code_review.md`, `docs/repair_state_recap.md`, `docs/repair_open_issues.md` if assumptions/issues change.
- Goal: align docs with code-backed invariants and name exact tests/presets required for closure.

## 7. Risk register

1. Unsafe direct restart write: `writeRestartCheckpointHdf5` can be called outside workflow boundary checks.
2. Snapshot/restart confusion: snapshot reader defaults missing restart-critical state and must never be used as restart continuation input.
3. Scheduler mirror confusion: `SimulationState::time_bin` mirrors are serialized but are not authority; future code must import scheduler state exactly.
4. Silent skipped output: workflow output callback returns instead of failing on unsafe boundary; useful operationally, but it can hide missing checkpoint cadence unless logged/tested.
5. PM pending refresh state: serialization exists, but safe reload semantics need explicit continuation proof or writer rejection at unsafe phases.
6. Output cadence state: `PendingOutputBoundary` is transient and not serialized; resumed runs may duplicate/skip outputs unless cadence is reconstructed exactly from step index and config or persisted.
7. Stateless RNG index inputs: star formation/feedback reproducibility depends on stable step indices, cell/star indices, and rank-local seed offsets; migration/reordering can change stochastic outcomes if not identity-based.
8. Roundtrip-only validation: existing workflow and tests can prove read/write integrity without proving future scientific equivalence.
9. Optional snapshot fields: reader fallback/default behavior is correct for compatibility but unsafe as restart truth.
10. HDF5/FFTW/MPI feature paths: CPU-only evidence is insufficient for Stage 8 closure when PM/restart paths are in scope.
