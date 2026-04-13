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
