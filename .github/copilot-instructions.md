# Copilot instructions for CHUI / CosmoSim

Read `AGENTS.md` before proposing edits. For non-trivial repository work, also read `docs/architecture/agent_task_interface.md` and `docs/architecture/developer_workflow_contract.md`.

Key constraints:

- The current user-facing project name is CHUI, but existing `CosmoSim`/`cosmosim` code and schema names are compatibility identifiers. Do not rename them opportunistically.
- Do not add pseudocode, placeholder core logic, or broad unrelated rewrites.
- Preserve single-source-of-truth runtime ownership, restart/schema compatibility, typed config validation, and deterministic evidence discipline.
- New or changed public interfaces under `include/cosmosim/**` require same-patch docs and migration notes.
- HDF5 GADGET/AREPO-style group and dataset names are canonical; do not restyle them.
- For clean handoff artifacts, exclude build trees, caches, generated outputs, binaries, and `CMakeUserPresets.json`.
