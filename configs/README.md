# Simulation configuration profiles

CosmoSim simulation inputs use the authoritative lower-snake-case
`*.param.txt` format. YAML and JSON are not accepted simulation-input formats.
JSON remains appropriate for generated machine-readable artifacts such as
provenance, manifests, capability summaries, and counters.

## Canonical profiles and compatibility entry points

The scenario subdirectories contain the canonical scientific profiles. Three
root-level files remain as compatibility entry points because existing scripts
and tests refer to them and repository symlinks are not portable to all Windows
workflows:

| Compatibility entry point | Canonical profile |
|---|---|
| `adaptive_bound_jeans_isolated_galaxy.param.txt` | `isolated_galaxy/disk_adaptive_bound_jeans_v01.param.txt` |
| `zoom_in.param.txt` | `zoom_in/zoom_adaptive_bound_jeans_v01.param.txt` |
| `cosmo_cube.param.txt` | `cosmo_cube/cube_effective_multiphase_tng_like_v01.param.txt` |

`integration_config_alias_consistency` uses CMake's cross-platform
`compare_files` command to reject byte-level drift between each compatibility
entry point and its canonical profile. Update the canonical profile first, then
copy it to the compatibility path.

Files under `configs/release/` are independent release-smoke fixtures. They may
currently match a scientific profile, but they are intentionally not aliases:
their purpose, run names, and future acceptance constraints may diverge.
