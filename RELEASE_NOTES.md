# CosmoSim Release Notes

## v0.1.0-initial (2026-04-08)

This initial serious release is scoped to desktop-first and small-cluster-first science workflows, with a clear scale-up path to larger HPC systems.

### Release scope and evidence bundle

- Release scope and support matrix: [`docs/releases/initial_release_scope.md`](docs/releases/initial_release_scope.md)
- Validation evidence summary: [`docs/releases/initial_release_validation_summary.md`](docs/releases/initial_release_validation_summary.md)
- Benchmark summary and environment disclosure: [`docs/releases/initial_release_benchmark_summary.md`](docs/releases/initial_release_benchmark_summary.md)
- Known issues and mitigations: [`docs/releases/known_issues.md`](docs/releases/known_issues.md)
- Release checklist and gate criteria: [`docs/releases/release_checklist.md`](docs/releases/release_checklist.md)
- Release artifact metadata manifest: [`release/release_manifest_v1.json`](release/release_manifest_v1.json)

### Included example configurations

- Reference release smoke workflow: [`configs/release/release_smoke_zoom_in.param.txt`](configs/release/release_smoke_zoom_in.param.txt)
- Isolated galaxy release smoke workflow: [`configs/release/release_smoke_isolated_galaxy.param.txt`](configs/release/release_smoke_isolated_galaxy.param.txt)
- Cosmological cube release smoke workflow: [`configs/release/release_smoke_cosmo_cube.param.txt`](configs/release/release_smoke_cosmo_cube.param.txt)

### Runtime entry point in this release

The shipped executable `cosmosim_harness` is a real config-driven runtime application in this release. A valid `param.txt` now produces a concrete run directory outcome, including normalized config and operational reporting artifacts, rather than only printing a banner.


### Stage 2 / Phase 3 claim boundary

Stage 2 timestep-authority repair evidence is limited to infrastructure invariants: scheduler-owned timestep bins, active-set provenance, restart mirror validation, and PM cadence legality. These release notes do **not** claim production-proven Phase 3 hierarchical TreePM multirate synchronization, distributed scaling certification, or Phase 3 closure.

### Schema and provenance commitments in this release

- Snapshot schema family remains `gadget_arepo_v4` for external interoperability.
- Restart schema family remains `cosmosim_restart_v6`.
- Provenance schema remains `provenance_v4` and includes deterministic normalized config hash metadata.
- Naming-rule policy follows `lower_snake_v1` for files/directories and explicit unit/frame suffixing in ambiguous internal variables.
