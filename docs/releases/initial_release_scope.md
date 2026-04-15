# Initial serious release scope (v0.1.0-initial)

## Supported operating modes (release-ready)

- `zoom_in`
- `cosmo_cube`
- `isolated_galaxy`
- `isolated_cluster`

These modes are release-ready for CPU-only workstations and small CPU clusters when run through preset-based build and test workflows.

## Numerical architecture in scope

- TreePM gravity.
- Finite-volume Godunov hydrodynamics in comoving variables.
- Patch-based AMR.
- Hierarchical timestepping.
- Modular baryonic physics (cooling/heating, star formation, stellar feedback, stellar evolution, black-hole AGN as documented below).
- HDF5-centered interoperability path where enabled.

## Dependency and build support matrix

### Production-ready paths

- CPU-only build via `cpu-only-debug` and `cpu-only-release` presets.
- MPI-enabled decomposition/restart path via `mpi-debug` preset.
- HDF5 snapshot/restart workflow via `hdf5-debug` preset.

### Conditionally supported path

- FFTW-enabled PM path via `fftw-debug` preset when FFTW is installed.

### Experimental path (not validated for production science)

- CUDA/GPU PM acceleration remains **experimental** and limited to smoke-level verification.

## Physics module status

### Production-ready in this release

- Cooling/heating source-term integration.
- Star-formation spawning and bookkeeping.
- Stellar evolution + stellar feedback deposition.

### Experimental or limited

- Black-hole AGN model is available for toy and exploratory studies only; calibration and production validation are not yet complete.
- Tracer support remains an opt-in feature path and is not a release gate for scientific production runs.

## Reproducibility and schema commitments

- Normalized, typed configuration is required before run execution.
- Provenance output is mandatory for release reference runs (`provenance_v1`).
- Snapshot schema compatibility: `gadget_arepo_v1`.
- Restart schema compatibility: `cosmosim_restart_v3`.
- Release metadata schema and packaging manifest: `release_manifest_v1`.
