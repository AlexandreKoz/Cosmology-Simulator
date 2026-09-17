# DMO first-light correctness closure — 2026-09-16

Scope: first cosmological dark-matter-only run only. This repair does not reopen
hydro, AMR, source physics, TreePM opening calibration, or the M1/M2 memory
architecture.

## Closed defects

- **COS-2:** `.param.txt` production cosmology is explicitly flat and
  radiation-free; `omega_matter + omega_lambda` must equal one within `1e-10`.
  The lower-level FLRW background also enforces closure across its explicit
  matter, Lambda, radiation, and curvature terms, so accepted backgrounds have
  `H(a=1)=H0`.
- **IC-2:** added `gadget_arepo_stored_peculiar` for standard GADGET/AREPO
  cosmological stored velocities (`v_peculiar = v_stored * sqrt(a)`). The
  pre-existing `sqrt_a_scaled_peculiar` semantic is preserved unchanged for
  compatibility. Canonical cosmological bridge configs now select the explicit
  standard convention.
- **RST-2:** exact runtime restart comparison now includes
  `last_drift_time_code` and `last_drift_scale_factor`. Serialization already
  persisted both lanes, so no schema bump is required.
- **PERF-3:** removed `scale_factor` from the invariant PM spectral-operator
  cache key only. MPI collective-entry consensus still includes scale factor.
  Profiling now exposes `spectral_operator_rebuilds` for regression coverage.

## Reproducibility / compatibility

- Invalid non-closed flat cosmologies now fail validation instead of silently
  changing the present-day Hubble rate.
- Existing custom IC decks using `sqrt_a_scaled_peculiar` retain their previous
  interpretation. Standard GADGET/AREPO cosmological decks should use the new
  explicit convention.
- Restart schema and HDF5 dataset names are unchanged.
- PM force equations, FFT normalization, Nyquist treatment, TreePM split, and
  rank epoch consensus are unchanged.
