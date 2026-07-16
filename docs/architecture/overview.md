# Architecture overview

This file is the architecture index for contributors.

## Module ownership map

1. `core`: configuration, provenance, units, time integration, shared state layout.
2. `gravity`: TreePM, PM, tree ordering/gravity, coupling contracts.
3. `hydro`: Godunov finite-volume hydro kernels and interfaces.
4. `amr`: patch-based AMR lifecycle and synchronization boundaries.
5. `physics`: cooling, SF, feedback, stellar evolution, AGN, tracers.
6. `io`: IC ingest, snapshot I/O, restart checkpointing.
7. `analysis`: diagnostics and halo/merger workflow scaffolding.
8. `parallel`: distributed-memory exchange and backend seams.
9. `utils`: narrow cross-module support only.

## Workflow runtime services

The workflow composition root creates one `workflows::RuntimeServices` bundle
per process. It borrows the process MPI context and profiler and carries the
deterministic-execution policy. Drift, gravity, and hydro runtime callbacks
receive this bundle (or an explicit borrowed dependency from it); they do not
construct replacement MPI contexts. `workflows::FailureCoordinator` provides
the phase gate used to reduce rank-local failures before any rank enters the
following collective phase. The first integrated gate covers restart-topology
validation; additional collective protocols should use the same authority.

`core::moduleDescriptors()` is the compile-time composition catalog. Each
descriptor names its typed config fragment, state/sidecar requirements,
stage/task and timestep contributions, restart and migration participation,
diagnostics, capability prerequisites, and explicit incompatibilities. The
legacy `moduleNames()` surface remains a compatibility view over the same
ordered catalog; new composition checks should consume descriptors.

## Architectural invariants

- SoA-first and hot/cold separated state organization.
- Explicit ownership boundaries; avoid duplicate systems.
- Typed, normalized, validated config before execution.
- Stable output naming and explicit schema versioning.
- Deterministic behavior by default unless controlled stochasticity is explicitly configured.

## Cross-links

- Decision log: [`decision_log.md`](decision_log.md)
- Developer workflow contract: [`developer_workflow_contract.md`](developer_workflow_contract.md)
- Agent task interface: [`agent_task_interface.md`](agent_task_interface.md)
- Config reference: [`../configuration.md`](../configuration.md)
- Output schema guide: [`../output_schema.md`](../output_schema.md)
- Validation ladder: [`../validation_plan.md`](../validation_plan.md)
- Profiling workflow: [`../profiling.md`](../profiling.md)
