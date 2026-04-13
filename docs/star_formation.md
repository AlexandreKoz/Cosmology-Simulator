# Star formation module

The baseline star-formation module implements a Schmidt-Kennicutt-inspired volumetric law for eligible gas:

\[
\dot{\rho}_{*} = \epsilon_{\mathrm{ff}} \frac{\rho_{\mathrm{gas}}}{t_{\mathrm{ff}}},\qquad
 t_{\mathrm{ff}} = \sqrt{\frac{3\pi}{32 G \rho_{\mathrm{gas}}}}
\]

## Eligibility criteria

A gas cell is eligible when all checks pass:

- `density_code >= physics.sf_density_threshold_code`
- `temperature_k <= physics.sf_temperature_threshold_k`
- `velocity_divergence_code <= physics.sf_min_converging_flow_rate_code`
- positive gas mass

These checks are separated from rate evaluation and spawn logic in `StarFormationModel`.

## Time-stepping integration

`StarFormationCallback` hooks the module into the `kSourceTerms` stage of
`core::StepOrchestrator`, so star formation can run inside the standard kick-drift-kick loop
instead of only in standalone tests/tools.

## Spawn modes

- **Deterministic** (`physics.sf_stochastic_spawning=false`):
  - formed mass over `dt` is transferred directly from gas to one spawned star particle.
- **Stochastic** (`physics.sf_stochastic_spawning=true`):
  - expected mass is converted to a Poisson-like draw around
    `lambda = expected_mass / physics.sf_min_star_particle_mass_code`.
  - spawn count uses `floor(lambda)` plus one Bernoulli trial from a reproducible hash RNG.

## Conservation and bookkeeping

For every spawn event:

- gas mass is reduced in `cells.mass_code`.
- gas density is scaled by the same mass fraction in `gas_cells.density_code`.
- one star particle is appended with species tag `kStar`.
- star sidecar metadata is written:
  - `formation_scale_factor`
  - `birth_mass_code`
  - `metallicity_mass_fraction`

A module sidecar block named `star_formation` stores run counters and normalized parameters for restart/provenance workflows.

## Config and provenance implications

New normalized config fields in `[physics]`:

- `sf_density_threshold_code`
- `sf_temperature_threshold_k`
- `sf_min_converging_flow_rate_code`
- `sf_epsilon_ff`
- `sf_min_star_particle_mass_code`
- `sf_stochastic_spawning`
- `sf_random_seed`

Assumptions (conservative baseline):

1. `density_code` and `mass_code` are internally consistent code units.
2. newborn stars inherit cell center position and are initialized with zero peculiar velocity.
3. one star sidecar row is emitted per spawn event to keep downstream stellar-evolution hooks explicit.
