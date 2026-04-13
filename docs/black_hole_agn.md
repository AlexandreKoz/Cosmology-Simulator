# Black-Hole and AGN baseline module

This document defines the baseline `black_hole_agn` model used by CosmoSim for the first full-physics BH/AGN path.

## Scope and model choices

- This is a **single explicit baseline model** (no alternative accretion or feedback closures mixed into this path).
- The module is controlled by `physics.enable_black_hole_agn` and defaults to **off** for lightweight modes.
- All calibration-sensitive parameters are explicit in config normalization and module sidecar metadata.

## Equations and conventions

The module uses:

- Bondi-like rate:
  \[
  \dot M_\mathrm{Bondi} = \alpha\,4\pi G^2 M_\mathrm{BH}^2\rho\left(c_s^2 + v_\mathrm{rel}^2\right)^{-3/2}
  \]
- Eddington rate:
  \[
  \dot M_\mathrm{Edd} = \frac{4\pi G M_\mathrm{BH} m_p}{\epsilon_r\sigma_T c}
  \]
- Accretion cap:
  \[
  \dot M_\mathrm{acc} = \min(\dot M_\mathrm{Bondi},\dot M_\mathrm{Edd})
  \]
- Feedback power:
  \[
  \dot E_\mathrm{fb} = \epsilon_f\epsilon_r\dot M_\mathrm{acc}c^2
  \]

### Conserved and non-conserved quantities

- Conserved by explicit bookkeeping in this baseline:
  - BH mass growth added to `black_holes.subgrid_mass_code` and synchronized into `particles.mass_code`.
  - Integrated feedback budget tracked in BH sidecar counters.
- Not strictly globally conserved in this baseline:
  - Gas mass removal by BH accretion is not yet modeled as a conservative sink from host gas cells.
  - Feedback deposition is local host-cell internal-energy increment via coupling efficiency; hydrodynamic redistribution is handled by hydro later.

## Restart/schema/provenance implications

- Restart `state/black_holes` now persists:
  - `host_cell_index`
  - `eddington_ratio`
  - `cumulative_accreted_mass_code`
  - `cumulative_feedback_energy_code`
  - `duty_cycle_active_time_code`
  - `duty_cycle_total_time_code`
- `physics.*` BH config keys are included in normalized config text and hash.
- Per-step module metadata sidecar key: `black_hole_agn` with counters and parameter echo.

## Future extension boundary

- MPI ownership/migration policy remains explicit via `particle_sidecar.owning_rank`; domain migration hooks are future work.
- Relative velocity closure currently uses `v_rel = 0` unless external estimators are provided.
