# Stellar evolution and enrichment bookkeeping

The production implementation and scientific rationale are documented in
[`metals_enrichment_and_mixing.md`](metals_enrichment_and_mixing.md).

## Configuration

Under `[physics]`:

- `enable_stellar_evolution` enables delayed SSP bookkeeping and source-runtime deposition;
- `stellar_evolution_table_path` selects a v2 age-by-birth-metallicity table;
- `stellar_evolution_require_production_table=true` rejects the built-in zero-yield compatibility table;
- `stellar_evolution_hubble_time_years` is retained only for configuration compatibility and legacy isolated API calls. Production cosmological ages use the FLRW integral.

## V2 table rows

Each whitespace-delimited row contains 22 values:

1. birth metallicity;
2. SSP age in years;
3. cumulative returned mass fraction;
4. cumulative total ejected metal fraction;
5. cumulative newly synthesized metal fraction;
6. cumulative event count per initial code mass;
7. cumulative energy in erg per initial code mass;
8-10. returned mass fractions for winds/AGB, CCSN, SNIa;
11-13. total ejected metal fractions by channel;
14-16. newly synthesized metal fractions by channel;
17-19. event counts by channel;
20-22. energies by channel.

Metadata comments use `# key = value`; supported keys include table ID/version,
source papers/repository, redistribution license, SHA-256, IMF, stellar mass
range, solar abundance reference, and `production_calibrated`.

The deterministic test table is `resources/stellar_evolution/test_synthetic_v2.txt`.
It is synthetic and non-production. `resources/stellar_evolution/default_v1.txt`
is retained only as historical source data and is not accepted by the v2 production loader.
