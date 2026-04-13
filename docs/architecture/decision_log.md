# Architecture decision log

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
