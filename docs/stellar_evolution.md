# Stellar Evolution and Enrichment Bookkeeping

This module tracks table-driven stellar mass return, metal yield, and feedback-energy budgets without performing deposition.

## Configuration keys

Under `[physics]` in param.txt:

- `enable_stellar_evolution` (`true|false`): toggles bookkeeping updates.
- `stellar_evolution_table_path` (string path): optional resource path to a table file.
  - Empty path uses the built-in reference table.
- `stellar_evolution_hubble_time_years` (float): age proxy conversion for `dt_code` and `a/a_birth` to years.

## Table format

Whitespace-delimited rows with 13 numeric columns:

1. `age_yr`
2. `return_fraction_total`
3. `metal_yield_fraction_total`
4. `energy_erg_per_initial_mass_code`
5. `return_winds`
6. `return_ccsn`
7. `return_snia`
8. `metal_winds`
9. `metal_ccsn`
10. `metal_snia`
11. `energy_winds`
12. `energy_ccsn`
13. `energy_snia`

Optional metadata comments in file header:

- `# table_id = <id>`
- `# table_version = <version>`

A sample is provided at `resources/stellar_evolution/default_v1.txt`.

## Provenance and outputs

Each step updates the `stellar_evolution` module sidecar with:

- `schema_version`
- `table_id`
- `table_version`
- `table_source_path`
- aggregate returned mass/metals/energy counters

Per-star cumulative bookkeeping is stored in `StarParticleSidecar` lanes, including channel-specific cumulative mass/metals/energy.
