# Cooling and heating module

The cooling/heating module is split into three auditable layers:

1. **Table loading** (`MetalLineCoolingTable`): loads optional metal-line cooling tables and records provenance (`source_path`, `table_label`, `row_count`).
2. **Rate lookup** (`CoolingRateProvider`): combines primordial rates, optional metal-line rates, UV background choice, and explicit self-shielding approximation.
3. **Source integration** (`CoolingSourceIntegrator`): integrates specific internal energy with bounded fractional substeps and exposes stiffness diagnostics.

## Config and provenance hooks

Physics parameters now include:
- `physics.cooling_model`
- `physics.uv_background_model`
- `physics.self_shielding_model`
- `physics.metal_line_table_path`
- `physics.temperature_floor_k`

The normalized config dump preserves these keys in stable order for restart/snapshot reproducibility.

## Assumptions

- Primordial baseline uses a compact analytic cooling/heating fit intended as a stable reference implementation, not a full collisional network.
- Metal-line contribution is optional and table-driven; missing table assets fail at load time if table cooling is requested.
- Self-shielding currently exposes `none` and `rahmati13_like` as explicit, parameterized approximations.
