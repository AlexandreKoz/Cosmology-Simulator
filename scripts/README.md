# Scripts

Place developer and CI helper scripts here using `lower_snake` file names.

## ci/

- `check_repo_hygiene.sh`: validates repository-root naming safety and required preset presence.
- `guard_feature_paths.sh`: runs baseline and PM/HDF5/FFTW guarded configure/build/test flows.
- `run_preset_pipeline.sh`: configure/build/test pipeline runner with artifact collection, optional benchmark sentinel execution, and per-preset phase report JSON.
- `check_build_metadata.py`: validates generated build metadata against versioned CI regression expectations.
- `run_reproducibility_gate.sh`: executes reproducibility-focused tests and records baseline hashes/report artifacts.
- `infrastructure_gates_manifest.tsv`: tabular manifest for explicit infrastructure-gate definitions (preset trio, test scope, metadata contract, artifact subdir, benchmark toggle).
