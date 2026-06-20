# H2 gas-cell identity repair acceptance

Date: 2026-06-19  
Suggested PR: `fix-h2-gas-cell-identity-workflow-scheduler-migration-acceptance`

## Scope and authority

This addendum supersedes the earlier optimistic H2 closeout language where it conflicts
with the current implementation. `SimulationState::gas_cell_identity` is the authoritative
identity-to-row contract. `gas_cell_id` is the stable key; `local_cell_row` is transient;
`parent_particle_id` is optional lineage/compatibility metadata only.

Ordinary production updates use `SimulationState::replaceGasCellIdentityRecords(...)` or
`restoreGasCellIdentityRecords(...)`. Those methods update the map, derive compatibility
lanes (`gas_cells.gas_cell_id`, `gas_cells.parent_particle_id`, and `cells.patch_index` when
PatchSoa exists), and reject map/mirror divergence at the next guarded boundary.

Importing a legacy state from mirror lanes remains deliberately explicit:
`refreshGasCellIdentityMapFromSidecarLanes()` is a bootstrap/import bridge, not a normal
runtime mutation path. Root-patch materialization updates canonical identity records and
then derives mirrors; it no longer reconstructs the map from mutable sidecars.

## Production workflow result

The non-AMR Cartesian compatibility path builds an explicit Cartesian geometry-to-dense-row
mapping from cell centers. It accepts 1D, 2D, and 3D uniform Cartesian grids independent of
dense storage order, rejects duplicate/missing/off-lattice/non-rectilinear inputs, and
scatters through stable gas-cell identity. It does not sort `SimulationState` into a preferred
physical order.

The reference workflow owns two schedules:

- a particle scheduler indexed by particle rows; and
- a gas-cell scheduler indexed by dense gas-cell rows but persisted/remapped by stable
  `gas_cell_id`.

`CellSoa::time_bin` is a derived mirror of the gas-cell scheduler. It is not scheduler
authority. The workflow independently submits adaptive candidates, advances both schedules
at the same tick, and synchronizes their respective mirrors only after the scheduler boundary.

## Restart contract

Restart schema v19 (`cosmosim_restart_v19`) adds `/gas_cell_scheduler` with:

- `identity_key = gas_cell_id`;
- `current_tick` and `max_bin` attributes;
- `gas_cell_id`, `bin_index`, `next_activation_tick`, `active_flag`, and
  `pending_bin_index` datasets.

The payload hash and restart diagnostics include the dedicated gas-cell scheduler. Modern
write/read validation requires its rows, stable IDs, and `CellSoa::time_bin` mirrors to agree
with the canonical map. A documented legacy reconstruction path remains for v14--v18 files;
it derives cells only from persisted cell mirrors and never from parent particle rows.

Patchless legacy/non-AMR cell states are legal restart inputs when their `patch_index` values
are the zero sentinel and their identity records use `owning_patch_id = 0`. AMR-backed states
remain required to provide complete, non-overlapping PatchSoa coverage.

## Evidence executed in this repair environment

HDF5 debug tests passed:

- `unit_gas_cell_identity_invariants`
- `integration_gas_cell_identity_migration`
- `integration_gas_cell_split_merge_remap`
- `integration_hydro_decoupled_gas_cells`
- `integration_amr_production_hydro_integration`
- `integration_amr_hydro_subcycling`
- `integration_restart_checkpoint_roundtrip`
- `integration_reference_workflow`
- `integration_reorder_compaction_sidecars`
- `integration_species_migration_invariants`
- `integration_transform_fuzz_invariants`
- `unit_time_integration`
- restart-equivalence harness, multirate-bin, AMR-hydro, pending-flux-register, and
  temporal-ghost scenarios.

`integration_reference_workflow` now runs two full workflow branches with an injected but
otherwise production-path state: canonical rows and shuffled rows. Each branch contains one
parentless gas cell, two gas cells sharing the same optional parent particle, explicit patch
geometry, independent gas scheduling, HDF5 output, and restart readback. The final restart
states compare equal by `gas_cell_id`; the parentless relation and both sibling records survive.

## Remaining acceptance boundary

The local migration and gas-cell commit tests are meaningful single-process evidence, but they
are not evidence of real rank-to-rank migration. This environment has no MPI compiler/launcher,
so MPI configuration and two-/three-rank H2 migration execution were not performed. H2 must
not be described as MPI-accepted until an MPI-enabled build runs a test that transfers
parentless and shared-parent gas-cell records across ranks, rebuilds the gas-cell scheduler by
`gas_cell_id`, and verifies restart continuation after that transfer.
