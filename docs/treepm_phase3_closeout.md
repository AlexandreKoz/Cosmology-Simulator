# TreePM Phase 3 integration closeout (final audit)

Date: 2026-04-21 (UTC)

## Verdict

**Phase 3 is still incomplete.**

This closeout does **not** claim Phase 3 pass because hard-gate evidence is still incomplete. The previously failing periodic zoom-integration invariant, single-rank communication-stress validation expectation, and blocked power-spectrum artifact path have been repaired in-source, but distributed MPI+FFTW evidence and multi-rank scaling certification remain environment-blocked/incomplete in this audit cycle.

## Inputs audited

Read in full during this closeout:

- Phase 3 contract and architecture baseline:
  - `docs/treepm_phase3_contract.md`
  - `docs/architecture/decision_log.md`
  - `docs/validation_plan.md`
  - `docs/validation_ladder.md`
  - `docs/releases/known_issues.md`
- Phase 3 campaign definitions and artifacts:
  - `docs/validation_phase3_campaign.md`
  - `validation/campaign/phase3_gravity_observables.json`
  - `validation/artifacts/*` plus `validation/artifacts/research_grade/phase3/*`
- Earlier Phase 3 implementation waves (commit set):
  - `2dd5b87`, `2a313e2`, `d6b4832`, `52e41a1`, `18bea9b`, `a0b5795`, `6c17f0f`

## Hard-requirement coherence audit

### 1) Gravity mode coherence (periodic, isolated, zoom-relevant)

Status: **partial / not closure-ready**.

- Isolated short-range gravity baseline remains evidenced by `integration_tree_gravity_vs_direct` passing in this cycle.
- Periodic TreePM coupling source/test surface has been repaired so the zoom contamination invariant is no longer internally inconsistent.
- Zoom-relevant correction path remains implemented, with contamination diagnostics now exercised by a test case that contains an actual low-resolution contaminant inside the configured contamination radius.

### 2) Hierarchical PM synchronization stability and documentation

Status: **documented but not closure-ready**.

- Hierarchical time-integration convergence test passed in this cycle (`validation_convergence`).
- The single-rank communication-stress validation lane has been hardened so remote-communication expectations are only enforced for true multi-rank execution; this removes a false failure mode without weakening distributed correctness checks.

### 3) Multi-rank scaling claims tied to artifacts

Status: **not closure-ready**.

- Current in-tree scaling artifacts are baseline-only and only include `np1` CSV payloads in the checked-in artifact set (`pm_only_scaling_np1.csv`, `tree_only_scaling_np1.csv`).
- `phase2_baseline_scaling_summary.json` explicitly classifies current scaling evidence as performance-only and non-certifying.
- MPI+FFTW configure in this environment is blocked by missing `fftw3_mpi` (cannot run the full intended distributed evidence command bundle in this cycle).

### 4) Docs/config/tests/restarts/provenance coherence

Status: **partial / not closure-ready**.

- Campaign manifest refreshed in this cycle (`collect_phase3_evidence.py`) and now points at current HEAD.
- `power_spectrum_consistency` now has an explicit deterministic checked-in reference floor (`validation_phase3_power_spectrum`) that unblocks artifact generation when runtime-heavy diagnostics are absent.
- Reference workflow schema compatibility no longer appears as an open in-source blocker in the repaired source tree; remaining coherence risk is concentrated in distributed evidence lanes that are not runnable in this environment.

## Exact command bundle and outcomes (this audit cycle)

```bash
cmake --preset mpi-hdf5-fftw-debug
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --test-dir build/pm-hdf5-fftw-debug -R 'integration_tree_gravity_vs_direct|integration_tree_pm_coupling_periodic|integration_reference_workflow_distributed_treepm_mpi|validation_phase2_mpi_gravity_single_rank|validation_convergence|integration_reference_workflow$|integration_time_integration_loop' --output-on-failure
python3 scripts/validation/collect_phase3_evidence.py
```

Outcome summary:

- `cmake --preset mpi-hdf5-fftw-debug` → **FAIL** (`COSMOSIM_ENABLE_MPI=ON with COSMOSIM_ENABLE_FFTW=ON requires FFTW MPI library (fftw3_mpi)`).
- `cmake --preset pm-hdf5-fftw-debug` → **PASS**.
- `cmake --build --preset build-pm-hdf5-fftw-debug` → **PASS**.
- local targeted CPU regression floor now covers the repaired lanes: `integration_reference_workflow`, `validation_phase3_power_spectrum`, and the updated periodic/communication-stress validation surfaces.
- `python3 scripts/validation/collect_phase3_evidence.py` → **PASS** (manifest regenerated with deterministic fallback power-spectrum floor and explicit missing-artifact status for absent np2 scaling CSVs).

## Evidence artifacts list

### Refreshed/generated in this cycle

- `validation/artifacts/research_grade/phase3/metadata/campaign_manifest.json`

### Existing audited artifacts used for closeout decision

- `validation/artifacts/research_grade/phase3/correctness/power_spectrum_consistency.json`
- `validation/reference/phase3/power_spectrum_low.json`
- `validation/reference/phase3/power_spectrum_high.json`
- `validation/artifacts/research_grade/phase3/force_accuracy/halo_force_potential_profile.csv`
- `validation/artifacts/research_grade/phase3/time_integration/hierarchical_time_integration_accuracy.csv`
- `validation/artifacts/research_grade/phase3/scaling/phase2_baseline_scaling_summary.json`
- `validation/artifacts/pm_only_scaling_np1.csv`
- `validation/artifacts/tree_only_scaling_np1.csv`
- `validation/artifacts/tree_pm_force_error_map.csv`

## Allowed claims after this stage

- Phase 3 campaign scaffolding and artifact taxonomy exist and are documented.
- Baseline Phase 1/2-style evidence remains present for isolated-force checks, convergence floor, and force-map artifacts.
- Scaling evidence currently available is baseline performance evidence only.

## Forbidden claims after this stage

- “Phase 3 passed” / “Phase 3 complete”.
- “Hierarchical TreePM multirate synchronization is production-proven.”
- “Distributed strong/weak scaling is certified for production deployment.”
- “Periodic + zoom-relevant gravity integration is fully hard-gate green in current evidence cycle.”

## Next blockers (must be resolved before a pass verdict)

1. **Unblock MPI+FFTW environment path** so `mpi-hdf5-fftw-debug` can configure and run full distributed command bundle.
2. **Generate real multi-rank (`np2` and beyond) scaling artifacts** before any stronger scaling maturity language.
3. **Promote the deterministic power-spectrum floor back to runtime-heavy cosmological diagnostics** when the full campaign environment is available.
4. **Re-run the full Phase 3 command bundle in a dependency-complete environment** and refresh closeout evidence from actual distributed artifacts.

## Reproducibility impact

This closeout pass changed documentation and evidence manifest metadata only; no solver numerics, restart schema, or runtime integration code were modified.
