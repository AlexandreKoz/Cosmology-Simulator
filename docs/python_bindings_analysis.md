# Python bindings and analysis interface

## Scope and stability

The Python module exposes a **stable, high-value analysis surface** for post-processing snapshots while keeping solver kernels in C++.

Initial public bindings:

- `load_frozen_config(path)`
- `read_snapshot(path, config)`
- `write_snapshot(path, state, config, normalized_config_text, git_sha="unknown")`
- `make_uniform_dark_matter_state(particle_count, mass_code, box_size_mpc_comov)`
- `DiagnosticsEngine.compute_run_health(...)`

## Ownership and copy semantics

- Snapshot and state array accessors return **read-only NumPy views** over C++ SoA storage.
- Views are zero-copy and lifetime-tied to the owning C++ object held by Python.
- Mutable views are intentionally not exposed in this first interface to preserve invariants.

### Zero-copy paths

- `SimulationState.particle_id()`
- `SimulationState.position_x_comoving()`
- `SimulationState.position_y_comoving()`
- `SimulationState.position_z_comoving()`
- `SimulationState.velocity_x_peculiar()`
- `SimulationState.mass_code()`

These methods produce NumPy arrays without duplicating backing data.

### Copying paths

- `DiagnosticsEngine.compute_gas_xy_projection_density(...)` returns a copied Python list (from C++ `std::vector<double>`).
- Bundle summary helper fields are copied into Python dictionaries.

## Build gating

Python bindings are optional and disabled by default.

```bash
cmake -S . -B build -DCOSMOSIM_ENABLE_PYTHON=ON -DCOSMOSIM_ENABLE_HDF5=ON
```

If `COSMOSIM_ENABLE_PYTHON=ON`, the configure stage requires:

- `Python3` interpreter + development module support
- `pybind11` CMake package

Core engine libraries remain buildable with `COSMOSIM_ENABLE_PYTHON=OFF`.

## Current limitations (conservative)

- The Python path is currently aimed at **single-rank post-processing workflows** and notebook analysis.
- Live multi-rank control, scheduler mutation, and in-situ distributed orchestration are not yet exposed.
- Snapshot loading still follows the existing HDF5 backend contract and requires `COSMOSIM_ENABLE_HDF5=ON`.

## Reproducibility and schema notes

- No snapshot schema changes were introduced.
- No config key changes were introduced.
- No restart schema changes were introduced.
- `write_snapshot(...)` requires explicit normalized config text to keep provenance auditable.
