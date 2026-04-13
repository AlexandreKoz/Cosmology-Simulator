# Halo/Subhalo and Merger-Tree Workflow (v1 planning scaffold)

## Scope and ownership boundaries

- **On-the-fly mode** (`analysis.halo_on_the_fly=true`) is intended to run analysis products during the simulation loop, but in v1 the planner is a standalone callable integrated at snapshot boundaries.
- **Post-processing mode** (`analysis.halo_on_the_fly=false`) is the conservative default and assumes a snapshot-aligned invocation.
- Simulation state ownership remains in `core::SimulationState`; halo workflow owns only derived catalog artifacts and schema-governed metadata.

## v1 schema decisions (explicit and centralized)

- Halo catalog format: `cosmosim_halo_catalog_v1`.
- Merger-tree plan format: `cosmosim_merger_tree_plan_v1`.
- Shared schema version: `1`.
- Required provenance fields:
  - `run_name`
  - `normalized_config_hash`
  - `snapshot_step_index`
  - `snapshot_scale_factor`
- Stable output naming:
  - `${analysis.halo_catalog_stem}_step_${step}.json`
  - `${analysis.merger_tree_stem}_step_${step}.json`

## Baseline finder and limitations

- Finder: Friends-of-Friends with periodic minimum-image metric in comoving coordinates.
- Candidate species: dark matter always; gas, stars, BH toggled by config; tracers excluded.
- Group IDs are deterministic per snapshot by sorting accepted groups with minimum particle ID.
- Complexity is O(N^2) pair checks in v1; this is intentionally explicit as a planning scaffold and not a distributed production halo finder.
- Subhalo entries are placeholders derived from host halos; no bound-particle subhalo finder is claimed.
- Merger-tree entries are a **plan scaffold** and do not yet provide validated progenitor/descendant matching.

## Config keys

- `analysis.enable_halo_workflow` (bool)
- `analysis.halo_on_the_fly` (bool)
- `analysis.halo_catalog_stem` (stable stem)
- `analysis.merger_tree_stem` (stable stem)
- `analysis.halo_fof_linking_length_factor` in `(0,1]`
- `analysis.halo_fof_min_group_size >= 2`
- `analysis.halo_include_gas` (bool)
- `analysis.halo_include_stars` (bool)
- `analysis.halo_include_black_holes` (bool)

## Restart and interoperability notes

- No persistent solver state is added to restart payloads in v1.
- Catalog outputs are external analysis products and include snapshot linkage and config provenance for reproducibility.
- HDF5 catalog output is intentionally deferred; JSON schema is explicit for fast iteration and testability.
