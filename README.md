# CosmoSim

CosmoSim is a desktop-first and small-cluster-first cosmological simulation framework with a scale-up path to larger HPC systems.

The flagship science target is physically credible zoom-in galaxy formation from cosmological initial conditions, with cosmological cubes, isolated galaxies, and isolated clusters maintained as supported operating modes.

## Documentation map

- **User workflow**
  - Build and dependency setup: [`docs/build_instructions.md`](docs/build_instructions.md)
  - Runtime parameter reference (GADGET/AREPO-style `param.txt`): [`docs/configuration.md`](docs/configuration.md)
  - Snapshot/restart/provenance schema guide: [`docs/output_schema.md`](docs/output_schema.md)
  - Validation ladder and reference decks: [`docs/validation_plan.md`](docs/validation_plan.md)
  - Profiling and benchmark workflow: [`docs/profiling.md`](docs/profiling.md)
  - Initial release package and checklist: [`RELEASE_NOTES.md`](RELEASE_NOTES.md)

- **Developer architecture and workflow**
  - Architecture index and ownership map: [`docs/architecture/overview.md`](docs/architecture/overview.md)
  - AI-agent task interface and repo-specific guardrails: [`AGENTS.md`](AGENTS.md), [`docs/architecture/agent_task_interface.md`](docs/architecture/agent_task_interface.md)
  - Design decisions and compatibility notes: [`docs/architecture/decision_log.md`](docs/architecture/decision_log.md)
  - Prompt/review contract for Codex + human contributors: [`docs/architecture/developer_workflow_contract.md`](docs/architecture/developer_workflow_contract.md)
  - Contribution rules: [`CONTRIBUTING.md`](CONTRIBUTING.md)


## Current maturity boundary

Stage 2 timestep-authority repairs document scheduler-owned timestep bins, scheduler-built active sets, restart mirror validation, and PM cadence legality as infrastructure invariants. This is not a claim that Phase 3 hierarchical TreePM multirate synchronization is production-proven; Phase 3 closure remains gated by the evidence contract in [`docs/treepm_phase3_contract.md`](docs/treepm_phase3_contract.md).

## Repository layout

- `include/cosmosim/` — public API headers by module.
- `src/` — implementation and ownership boundaries.
- `docs/` — user and developer documentation.
- `configs/` — normalized example parameter files.
- `validation/` — validation input decks and tolerances.
- `tests/unit`, `tests/integration`, `tests/validation` — test ladder.
- `bench/` — benchmark and profiling hooks.

## Quickstart (CPU-only)

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

For HDF5/FFTW/MPI/GPU paths, use the preset matrix documented in [`docs/build_instructions.md`](docs/build_instructions.md).

## Reproducibility discipline

CosmoSim treats reproducibility as a first-order requirement:

- typed + validated + normalized config state,
- deterministic config hash and provenance record,
- explicit schema versioning for snapshot/restart payloads,
- stable naming constraints for output stems.

These rules are implemented in the core and I/O interfaces and documented in [`docs/configuration.md`](docs/configuration.md) and [`docs/output_schema.md`](docs/output_schema.md).
