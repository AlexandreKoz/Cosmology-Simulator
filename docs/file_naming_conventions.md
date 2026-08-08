# CHUI / CosmoSim file and artifact naming conventions

**Version 2.0 — 2026-08-08**

This repository policy supersedes older examples that used YAML configuration
or `.cc` as the default implementation suffix.

- Source and directory names are ASCII `lower_snake` unless an externally
  mandated filename says otherwise (`CMakeLists.txt`, `AGENTS.md`, etc.).
- Public C++ headers live under `include/cosmosim/<module>/` and use `.hpp`.
- C++ implementations live under `src/<module>/` and use the repository's
  established `.cpp` suffix (CUDA translation units use `.cu`).
- Simulation configuration files use `.param.txt`; do not add YAML/JSON input
  aliases.
- Snapshots use `snap_###.hdf5`; restarts use `restart_###.hdf5` unless a
  workflow's documented rank suffix is required.
- Files being written atomically use the `.part` suffix and are renamed to the
  final target only after a successful close/flush sequence. `.tmp` is not the
  simulation-output partial-file convention.
- Scratch/runtime artifacts remain outside version-controlled source content.
- Canonical external HDF5 group/dataset names are schema, not filesystem style,
  and are preserved exactly.

Changes to these rules require a deliberate governance edit; do not drift them
through isolated implementation patches.
