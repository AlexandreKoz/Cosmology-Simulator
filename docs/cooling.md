# Cooling and effective-ISM interaction

The authoritative cooling implementation is `CoolingRateProvider` plus `CoolingSourceIntegrator`. The effective multiphase EOS uses that provider to construct a deterministic reference equilibrium table once, with its cooling assumptions and hash recorded in provenance.

Below the effective density threshold, gas follows ordinary cooling/heating. Above threshold, the thermodynamic closure establishes a lower effective branch. Dense gas below that branch receives an explicit positive EOS source adjustment; dense gas far above the branch remains hot and is excluded from star formation until it approaches equilibrium. Finite-timescale relaxation is available; instantaneous closure remains explicit in the source ledger.

The following quantities are intentionally distinct:

- radiative cooling/heating energy change;
- effective-ISM energy added/removed;
- star-birth internal-energy sink;
- stellar-feedback energy injection.

No local metallicity-dependent multiphase table is rebuilt in the hot loop. Local metallicity continues to affect the ordinary cooling and enrichment pathways.
