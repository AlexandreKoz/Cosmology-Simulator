# Architecture decision log

## 2026-04-19 — ADR-FEATURE-GRAVITY-010: Freeze Phase 2 distributed TreePM contract and evidence floor

### Status
Accepted (feature-branch distributed gravity contract freeze)

### Context
Phase 1 TreePM contract is implemented and documented, but Phase 2 distributed ownership/communication semantics were not frozen in typed config and architecture docs. CI surface also centered on CPU/HDF5/PM paths plus an optional MPI release smoke job, which is insufficient for distributed gravity development cadence.

### Decision
- Freeze Phase 2 config-facing contract with two typed keys:
  - `numerics.treepm_pm_decomposition_mode` (currently only `slab`)
  - `numerics.treepm_tree_exchange_batch_bytes` (batch cap for tree export/import payloads)
- Add a dedicated distributed-gravity developer preset:
  - `mpi-hdf5-fftw-debug`
- Publish repository-specific Phase 2 architecture contract in `docs/treepm_phase2_distributed_contract.md`, including:
  - PM slab ownership semantics,
  - particle-owner vs slab-owner roles,
  - distributed FFT stage ordering,
  - tree export/import batch contract,
  - active-set and migration timing constraints,
  - restart continuity and deterministic limits.
- Explicitly record that pseudo-multi-rank tests are not completion evidence for distributed TreePM implementation.

### Consequences
- Positive: distributed gravity work now has an auditable typed/config/docs baseline without claiming algorithm completion.
- Positive: CI/build surface can exercise MPI+HDF5+FFTW debug configuration directly rather than only optional MPI release smoke.
- Tradeoff: only slab decomposition is valid in this phase; future modes require explicit ADR + typed/config/docs/test expansion.

### Evidence references
- `include/cosmosim/core/config.hpp`
- `src/core/config.cpp`
- `CMakePresets.json`
- `docs/treepm_phase2_distributed_contract.md`
- `.github/workflows/ci.yml`

## 2026-04-12 — ADR-REPAIR-BOUNDARY-005: Keep reference workflow assembly outside `core/`

### Status
Accepted (infrastructure boundary repair)

### Context
`include/cosmosim/core/reference_workflow.hpp` assembled diagnostics, physics callbacks, and I/O roundtrips directly from the `core/` include tree. This inverted the dependency direction guardrail (`core` foundational; higher layers depend on it).

### Decision
- Move concrete reference workflow assembly into a dedicated integration layer:
  - `include/cosmosim/workflows/reference_workflow.hpp`
  - `src/workflows/reference_workflow.cpp`
- Keep `core/` focused on abstract stage contracts and typed simulation state.
- Add a dependency-direction guard test that fails when `include/cosmosim/core/` or `src/core/` includes `analysis/`, `io/`, `physics/`, `workflows/`, or `app/` headers.

### Consequences
- Positive: Restores one-way dependency flow from higher orchestration layers into `core`.
- Positive: Makes boundary regressions auditable in CI via an explicit guard test.
- Tradeoff: Callers should migrate toward `cosmosim::workflows` naming; transitional aliases remain for now.

## 2026-04-07 — ADR-REPAIR-FREEZE-001: Freeze repair evidence before touching production code

### Status
Accepted (stabilization phase)

### Context
The repository already contains significant architecture work (SoA storage, hot/cold split, sidecars, typed config, provenance, staged integration contracts). Entering emergency stabilization without a frozen baseline risks accidental churn and invalid claims of closure.

### Decision
Before any further repair implementation:

1. Capture a command-backed repair-state recap.
2. Capture a command-backed open-issues ledger.
3. Freeze the current public-interface baseline from `include/cosmosim/**`.
4. Record both default-path green status and feature-path red status explicitly.

### Consequences
- Positive: Later prompts can cite a single factual baseline instead of re-discovering repo state.
- Positive: Prevents CPU-only green path from being mistaken for HDF5 feature-path closure.
- Positive: Supports minimal, auditable fixes by tying each open issue to file+command+symptom.
- Tradeoff: Adds documentation overhead but no production-code churn.

### Evidence references
- `docs/repair_state_recap.md`
- `docs/repair_open_issues.md`

## 2026-04-07 — ADR-REPAIR-AUDIT-002: Hold progression to P20 until dependency-enabled presets are proven

### Status
Accepted (emergency stabilization gate)

### Context
Final repair audit was run against required scope. CPU-only path is green, but this environment lacks HDF5 development libraries, preventing configure/build/test execution on required HDF5 and PM/HDF5/FFTW paths.

### Decision
Do **not** authorize progression to P20+ until all required dependency-enabled preset commands pass in an environment with HDF5 (and FFTW for PM path).

### Consequences
- Positive: Prevents invalid closure claims based on CPU-only evidence.
- Positive: Keeps dependency failure behavior explicit and actionable.
- Tradeoff: Progression is delayed until proper feature-path environment is available.

### Evidence references
- `docs/repair_state_recap.md`
- `docs/repair_closeout_report.md`

## 2026-04-07 — ADR-REPAIR-BLOCK-003: Keep P20 progression blocked until PM HDF5+FFTW validation passes

### Status
Accepted (stabilization gate enforcement)

### Context
Dependency-enabled stabilization validation was re-run with presets in this environment.

- CPU-only configure/build/test passed.
- HDF5 configure/build/test passed.
- PM HDF5+FFTW configure/build passed, but test preset failed (`2/36` failures):
  - `unit_pm_solver` assertion failure in `tests/unit/test_pm_solver.cpp:82`.
  - `integration_tree_pm_coupling_periodic` runtime failure with `rel_l2=18129.9` against required `<= 0.75`.

### Decision
Progression beyond P19 remains blocked. Do not authorize P20 until `test-pm-hdf5-fftw-debug` passes completely.

### Consequences
- Positive: Enforces dependency-enabled acceptance criteria and avoids CPU-only false closure.
- Positive: Keeps blocker scope narrow (PM/FFTW test path) without speculative refactoring.
- Tradeoff: Delivery remains gated pending PM HDF5+FFTW test stabilization.

### Evidence references
- `docs/repair_closeout_report.md`
- `docs/repair_open_issues.md`

## 2026-04-07 — ADR-REPAIR-GATE-004: Clear P20 gate only when all three preset validation paths pass

### Status
Accepted (post-repair gate policy)

### Context
P19 stabilization was previously blocked by PM FFTW-path numerical convention defects affecting:
- `unit_pm_solver` analytic mode directionality check.
- `integration_tree_pm_coupling_periodic` long-range amplitude consistency (`rel_l2` blow-up).

After targeted PM FFTW-path correction, required validation commands were re-run and all three required presets passed:
- `test-pm-hdf5-fftw-debug`
- `test-hdf5-debug`
- `test-cpu-debug`

### Decision
P20 gate is considered **CLEARED only if all three preset paths pass in the same validation cycle**:
1. CPU-only preset path
2. HDF5 preset path
3. PM HDF5+FFTW preset path

Any future regression in one path re-blocks the gate.

### Consequences
- Positive: Prevents partial-path closure and preserves dependency-enabled confidence.
- Positive: Keeps PM FFTW path first-class in stabilization acceptance.
- Tradeoff: Slightly longer validation cycle per gate decision.

### Evidence references
- `docs/repair_closeout_report.md`
- `docs/repair_open_issues.md`


## 2026-04-14 — ADR-REPAIR-CONTINUATION-005: Remove normalized-config self-hash and reject incomplete restart continuations

### Status
Accepted (targeted infrastructure repair)

### Context
The normalized config artifact embedded its own hash, which made the canonical text internally inconsistent and non-round-trippable. Separately, restart payloads/schema omitted stellar-evolution continuation lanes while the integrity hash and compatibility contract still claimed exact continuation.

### Decision
- Define the normalized config hash over the normalized config text alone; do not embed a self-hash line inside the canonical text.
- Require continuation metadata agreement: normalized text hash == stored normalized-config hash == provenance config hash.
- Bump restart schema to `cosmosim_restart_v3` and reject older incomplete restart artifacts.
- Persist and hash the full stellar-evolution star sidecar payload.

### Consequences
- Positive: normalized config dumps are deterministic, parseable, and hash-consistent.
- Positive: restart artifacts once again match their exact-continuation claim for currently serialized state.
- Tradeoff: `cosmosim_restart_v2` artifacts are intentionally incompatible because they omit required continuation state.

### Evidence references
- `tests/unit/test_config_parser.cpp`
- `tests/unit/test_restart_checkpoint_schema.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/integration/test_snapshot_hdf5_roundtrip.cpp`

## 2026-04-18 — ADR-FEATURE-GRAVITY-006: Feature-branch exception for PM periodic solver contract hardening

### Status
Accepted (feature-branch gravity upgrade)

### Context
Repository repair guardrails include a repair-only prohibition on solver-behavior edits. This branch is explicitly scoped as a gravity-upgrade feature branch requiring PM solver contract hardening, including potential output formalization and stricter periodic Poisson convention tests.

### Decision
- Authorize PM solver behavior/documentation/test changes in this branch only, while preserving existing architecture direction:
  - no dependency-direction inversion,
  - explicit typed/configured runtime inputs,
  - auditable solver sign/normalization conventions,
  - evidence-backed closure via preset configure/build/test commands.

### Consequences
- Positive: Enables scientifically auditable PM potential+force contract without violating boundary discipline.
- Positive: Removes ambiguity around `δρ`, zero mode, Fourier `k` mapping, and force sign/normalization.
- Tradeoff: Branch intentionally changes numerical-solver contract surface and therefore requires same-patch docs/tests updates.

### Evidence references
- `include/cosmosim/gravity/pm_solver.hpp`
- `src/gravity/pm_solver.cpp`
- `docs/pm_gravity_solver.md`
- `tests/unit/test_pm_solver.cpp`
- `tests/integration/test_pm_periodic_mode.cpp`
- `tests/validation/test_validation_integration.cpp`


## 2026-04-18 — ADR-FEATURE-GRAVITY-007: Enable PM assignment/gather upgrade (CIC+TSC) on gravity feature branch

### Status
Accepted (feature-branch gravity upgrade)

### Context
The gravity-upgrade branch requires end-to-end PM assignment/gather configurability (`cic`, `tsc`) with auditable transfer-kernel behavior. This intentionally changes solver behavior versus prior CIC-only runtime gating.

### Decision
- Keep PM assignment and gather strictly matched by selected runtime scheme (no mixed silent modes).
- Implement TSC alongside CIC in deposition and interpolation paths.
- Make k-space deconvolution scheme-aware and apply it to the combined transfer operator with a protective floor.
- Keep conservative defaults (`cic`, deconvolution off) and record behavior in same-patch docs/tests.

### Consequences
- Positive: `kTsc` is no longer a placeholder; config/runtime/docs/tests align.
- Positive: PM transfer behavior is explicit and auditable across assignment choices.
- Tradeoff: TSC increases stencil work (27-point 3D) compared to CIC (8-point 3D).
- Tradeoff: Explicit CUDA PM path remains CIC-only in this build and rejects `tsc` rather than silently diverging.

### Evidence references
- `include/cosmosim/gravity/pm_solver.hpp`
- `src/gravity/pm_solver.cpp`
- `src/workflows/reference_workflow.cpp`
- `docs/pm_gravity_solver.md`
- `docs/configuration.md`
- `tests/unit/test_pm_solver.cpp`
- `tests/integration/test_pm_periodic_mode.cpp`
- `tests/validation/test_validation_integration.cpp`

## 2026-04-18 — ADR-FEATURE-GRAVITY-008: TreePM split/cutoff contract upgrade on gravity feature branch

### Status
Accepted (feature-branch gravity upgrade)

### Context
This branch stage requires solver-behavior edits in the TreePM coupling layer to make the AREPO/GADGET-style split contract explicit and runtime-configurable from normalized config controls (`asmth_cells`, `rcut_cells`). Existing residual traversal lacked explicit RCUT truncation behavior.

### Decision
- Treat mesh-cell controls as authoritative and derive:
  - `r_s = asmth_cells * Δmesh`
  - `r_cut = rcut_cells * Δmesh`
  - `Δmesh = box_size / PMGRID`
- Keep Gaussian split family:
  - PM: `exp(-k^2 r_s^2)`
  - residual tree: Gaussian complementary short-range factor.
- Enforce explicit residual cutoff in traversal:
  - node-level AABB pruning,
  - node-acceptance guard requiring full in-cutoff enclosure,
  - leaf pair skipping beyond `r_cut`.
- Expand diagnostics/docs/tests so split parameters, composition checks, and cutoff pruning are auditable.

### Consequences
- Positive: Removes hidden split heuristics and dual-source ambiguity.
- Positive: `rcut_cells` now changes actual residual traversal work and force contribution domain.
- Positive: Improves reproducibility/audibility via runtime diagnostics and same-patch docs/tests.
- Tradeoff: Tree residual traversal now performs additional cheap AABB distance bounds to guarantee safe pruning.

### Evidence references
- `include/cosmosim/gravity/tree_pm_split_kernel.hpp`
- `include/cosmosim/gravity/tree_pm_coupling.hpp`
- `src/gravity/tree_pm_coupling.cpp`
- `src/workflows/reference_workflow.cpp`
- `docs/tree_pm_coupling.md`
- `tests/unit/test_tree_pm_split_kernel.cpp`
- `tests/integration/test_tree_pm_coupling_periodic.cpp`

## 2026-04-18 — ADR-FEATURE-GRAVITY-009: Enable explicit PM long-range cadence/reuse in reference runtime

### Status
Accepted (feature-branch gravity upgrade)

### Context
`numerics.treepm_update_cadence_steps` existed in typed/normalized config but runtime still hard-rejected values above `1`, leaving no auditable PM refresh/reuse behavior despite hierarchical stepping and active-set gravity callbacks.

### Decision
- Treat `treepm_update_cadence_steps` as an authoritative runtime control in the reference workflow.
- Apply cadence per gravity kick opportunity (`gravity_kick_pre`, `gravity_kick_post`):
  - refresh PM long-range mesh solve every `N` opportunities,
  - reuse cached long-range field otherwise,
  - always recompute short-range tree residual each kick.
- Keep default conservative (`N=1`) to preserve immediate-refresh baseline behavior.
- Record refresh/reuse metadata (field version, build step, build scale factor) in workflow report and operational events.

### Consequences
- Positive: cadence key now changes real runtime behavior and is diagnosable.
- Positive: PM reuse is explicit, deterministic, and traceable instead of hidden in callback internals.
- Tradeoff: Phase 1 cadence is single-rank and intentionally simple; it is not yet a full multirate integrator.

### Evidence references
- `src/workflows/reference_workflow.cpp`
- `include/cosmosim/workflows/reference_workflow.hpp`
- `include/cosmosim/gravity/tree_pm_coupling.hpp`
- `src/gravity/tree_pm_coupling.cpp`
- `tests/integration/test_reference_workflow_end_to_end.cpp`
- `docs/time_integration.md`
- `docs/reference_workflow.md`
- `docs/pm_gravity_solver.md`

## 2026-04-18 — ADR-FEATURE-GRAVITY-010: Persist TreePM provenance contract in snapshot/restart/log outputs

### Status
Accepted (feature-branch gravity upgrade)

### Context
TreePM controls and derived scales materially affect force accuracy/reproducibility, but output artifacts did not carry a complete gravity contract across snapshots, restarts, and runtime operational events.

### Decision
- Extend `core::ProvenanceRecord` with gravity provenance fields covering:
  - controls (`pm_grid`, assignment, deconvolution, `asmth_cells`, `rcut_cells`, cadence),
  - derived scales (`Δmesh`, `r_s`, `r_cut`),
  - softening policy/kernel/epsilon and PM FFT backend.
- Bump snapshot schema identity to `gadget_arepo_v2` / `schema_version=2` and persist the new provenance attributes in `/Provenance`.
- Keep restart schema `cosmosim_restart_v3`; serialize gravity provenance through existing provenance dataset and integrity hash path.
- Emit structured runtime events (`gravity.treepm_setup`, enriched `gravity.pm_long_range_field`) that expose the same gravity contract at setup/refresh points.

### Consequences
- Positive: Snapshot/restart/runtime logs can reconstruct the exact TreePM runtime contract used for a run.
- Positive: Derived physical/comoving scales are auditable without re-deriving from partial controls.
- Tradeoff: Snapshot schema version increased; downstream readers that hardcode schema version `1` must accept `2`.

### Evidence references
- `include/cosmosim/core/provenance.hpp`
- `src/core/provenance.cpp`
- `include/cosmosim/io/snapshot_hdf5.hpp`
- `src/io/snapshot_hdf5.cpp`
- `src/workflows/reference_workflow.cpp`
- `tests/integration/test_snapshot_hdf5_roundtrip.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/integration/test_provenance_roundtrip.cpp`
- `tests/integration/test_reference_workflow_end_to_end.cpp`

## 2026-04-19 — ADR-FEATURE-GRAVITY-011: TreePM Phase 1 integration hard-gate closure contract

### Status
Accepted (feature-branch gravity upgrade)

### Context
Final TreePM Phase 1 integration requires coherence across runtime wiring, config normalization, provenance/output audibility, validation docs, and command-backed evidence. Earlier stage changes introduced PMGRID/ASMTH/RCUT/cadence controls and TreePM runtime behavior, but Phase 1 closure requires a single explicit contract and hard-gate evidence checklist.

### Decision
- Adopt an explicit Phase 1 integration hard gate that must be evidenced together:
  1. no hard-coded PM mesh in workflow (PMGRID is config-driven),
  2. no hidden split assumptions (ASMTH/RCUT drive runtime `r_s`/`r_cut`),
  3. deterministic repeated-run reproducibility evidence on identical config,
  4. force-error map artifact over PMGRID/ASMTH/RCUT sweep.
- Record the final contract and evidence procedure in `docs/treepm_phase1_closeout.md` and align cross-doc language (`configuration`, `tree_pm_coupling`, `validation_ladder`, `validation_plan`).

### Consequences
- Positive: Phase 1 completion claims are tied to auditable artifacts and commands rather than narrative-only status.
- Positive: Runtime/config/docs/output semantics are synchronized for TreePM controls.
- Tradeoff: Phase 1 closeout remains scoped to single-rank periodic TreePM and non-Ewald references.

### Evidence references
- `docs/treepm_phase1_closeout.md`
- `tests/integration/test_reference_workflow_end_to_end.cpp`
- `bench/bench_tree_pm_force_error_map.cpp`
- `validation/artifacts/tree_pm_force_error_map.csv`

## 2026-04-19 — ADR-INFRA-RESTART-012: Version distributed TreePM continuation state in restart/provenance

### Status
Accepted (infrastructure repair)

### Context
Phase 2 distributed TreePM restart continuation lacked an explicit, versioned contract for rank ownership, PM slab layout, cadence state, and long-range field refresh metadata. This made rank/layout mismatch debugging ambiguous and allowed implicit continuation assumptions.

### Decision
- Bump restart schema to `cosmosim_restart_v4` and require a new restart group `/distributed_gravity/state` carrying serialized `parallel::DistributedRestartState` (schema_version `2`).
- Persist explicit distributed continuation fields:
  decomposition epoch, owning-rank table, PM grid/layout metadata, kick-opportunity cadence state, long-range refresh/version metadata, and explicit long-range restart policy.
- Adopt and enforce restart policy `deterministic_rebuild` for PM long-range field continuation.
  Cached long-range arrays are not serialized; restart resumes with deterministic rebuild at the next cadence-triggered refresh.
- Extend provenance (`provenance_v3`) with distributed restart diagnostics (epoch/world/grid/slab signature/kick+field version/policy).
- Add compatibility diagnostics API `evaluateDistributedRestartCompatibility(...)` to emit explicit mismatch reasons.

### Consequences
- Positive: distributed restart semantics are explicit, versioned, and integrity-hashed.
- Positive: rank/layout/config mismatch debugging has typed reporting instead of opaque failures.
- Tradeoff: `cosmosim_restart_v3` artifacts are intentionally incompatible with `v4`.

### Evidence references
- `include/cosmosim/io/restart_checkpoint.hpp`
- `src/io/restart_checkpoint.cpp`
- `include/cosmosim/parallel/distributed_memory.hpp`
- `src/parallel/distributed_memory.cpp`
- `include/cosmosim/core/provenance.hpp`
- `src/core/provenance.cpp`
- `tests/unit/test_restart_checkpoint_schema.cpp`
- `tests/integration/test_restart_checkpoint_roundtrip.cpp`
- `tests/unit/test_parallel_distributed_memory.cpp`
