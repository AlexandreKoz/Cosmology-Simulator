# Known issues (initial release)

## GPU PM acceleration is experimental

- Scope: CUDA PM path.
- Impact: not release-validated for production science runs.
- Mitigation: treat CUDA runs as exploratory; use CPU/MPI production path for published results.

## Black-hole AGN model calibration is incomplete

- Scope: black-hole AGN toy model.
- Impact: physically plausible exploratory behavior but not yet a validated production calibration.
- Mitigation: disable AGN for baseline release reference runs unless explicitly studying the model's exploratory behavior.

## HDF5 interoperability requires feature-enabled build

- Scope: snapshot/restart I/O workflows.
- Impact: users may assume HDF5 works in all presets, but it is gated by `COSMOSIM_ENABLE_HDF5`.
- Mitigation: use HDF5-enabled presets and follow `docs/build_instructions.md` dependency requirements.

## Large-cluster scaling evidence not included in this initial release

- Scope: high-rank production deployment.
- Impact: architecture supports scale-up, but this release package does not claim full strong/weak scaling certification on large HPC systems.
- Mitigation: use desktop/small-cluster scale for this release; add dedicated scaling campaign in a follow-up release.


## TreePM Phase 3 maturity is not yet implemented

- Scope: integrator-grade hierarchical cadence/PM synchronization maturity and Phase 3 closure evidence.
- Impact: current release evidence supports Phase 1/2 baseline correctness contracts, but not Phase 3 completion claims.
- Mitigation: follow `docs/treepm_phase3_contract.md` hard-gate requirements before using Phase 3-complete language.


## Phase 3 scaling evidence remains baseline-only in this stage

- Scope: strong/weak scaling certification.
- Impact: campaign now exports explicit baseline scaling artifacts and provenance, but certification-grade large-rank sweeps are still pending.
- Mitigation: treat `validation/artifacts/research_grade/phase3/scaling/phase2_baseline_scaling_summary.json` as baseline evidence only and complete the larger-rank follow-on campaign before certification claims.
