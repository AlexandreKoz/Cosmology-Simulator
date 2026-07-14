# Power-spectrum diagnostics

## Scope

`DiagnosticsEngine::computePowerSpectrumEstimate` is the explicit deterministic
three-dimensional estimator used by controlled validation cases. It is an
additive interface; the existing compact `computePowerSpectrum` overloads keep
their historical NGP/no-deconvolution behavior and continue to omit empty bins.

The detailed estimator accepts explicit mesh, bin-count, mass-assignment,
assignment-window, and Poisson shot-noise policies. It returns every requested
bin, including empty bins, together with its policy and normalization metadata.
It remains a direct-DFT reference implementation with `O(N_mesh^6)` work and is
not a high-dynamic-range or distributed production-analysis claim.

## Estimator contract

For a periodic cube of code volume `V`:

1. Particle mass is deposited deterministically with NGP or CIC assignment.
2. The density contrast is `delta = rho / mean(rho) - 1`.
3. The discrete transform is normalized as
   `delta_k = sum(delta_grid exp(-i k x_grid)) / N_mesh`.
4. The reported spectrum is `P(k) = V |delta_k|^2`.
5. CIC deconvolution divides power by the square of
   `W_CIC = product_axis sinc(pi n_axis / mesh_n)^2`. NGP uses the corresponding
   first-power assignment window.
6. Shells are linearly spaced from zero to the three-dimensional corner mode,
   and the DC mode is excluded. `mode_count` counts every represented discrete
   mesh mode. Empty requested shells remain present with `empty=true`,
   `mode_count=0`, and zero power.
7. The weighted Poisson level is `V sum(m_i^2) / sum(m_i)^2`. The selected policy
   either reports this value without subtraction or subtracts it from each
   non-empty bin.

`k` is reported in `code_length^-1`; `P(k)` and the Poisson level are reported
in `code_length^3`. For configurations whose code length is Mpc comoving, these
map to `Mpc_comoving^-1` and `Mpc_comoving^3`.

## Validation use and limits

The controlled DMO Zel'dovich gate uses a `12^3` CIC mesh, 32 linear shells,
assignment-window deconvolution, and reported-but-not-subtracted Poisson noise.
The narrow first-mode shell measures linear growth from
`sqrt(P_evolved / P_initial)`. Unit coverage checks policy labels, units,
normalization, total mode count, explicit empty bins, and subtraction behavior.
The DMO validation additionally checks initial amplitude, evolved growth, and
restart-equivalent spectra and writes a deterministic JSON artifact. That JSON
retains empty bins and maps their internal zero-power sentinel to JSON `null` so
absence of measured power cannot be confused with a measured zero.

The growth acceptance rule is expressed relative to the expected increment
`Delta D=a_final/a_initial-1`, not the order-unity growth factor. Both the
direct displacement-mode and `sqrt(P_final/P_initial)` increments must be positive
and agree with `Delta D` within 7.5 percent; the two measured increments must
agree with each other within relative `1e-3`. This prevents a zero-growth
spectrum from passing a loose order-unity comparison.

The 2026-07-13 controlled gate measured fundamental powers
`2.1315894140899176e-6` and `2.1827523417862502e-6`, giving
`sqrt(P_evolved/P_initial)=1.011929959672577` versus expected scale-factor
growth `1.0119399775244009`. After normalizing only `world_size`, MPI np1--np4
JSON artifacts had SHA-256
`53ee518b63d21cd28790415b8f076614497536bb6a62c95276b8ed6c14984bf3`.
A true non-MPI serial run and MPI-world-one run emitted byte-identical JSON
(SHA-256
`8acbf3d6826250ca4fbabb0761511ff1178cbb7c20472cb2d8c8073dd16d355c`).
These are reproducibility checks for the controlled fixture, not a performance
or high-dynamic-range claim.

No config key, snapshot schema, restart schema, or provenance schema changes are
introduced by this additive analysis interface.
