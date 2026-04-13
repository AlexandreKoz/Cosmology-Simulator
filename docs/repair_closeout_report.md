# Repair Closeout Report (P01–P19)

_Date: 2026-04-07 (UTC)_

## Gate decision

**GO — P19 stabilization criteria satisfied; progression to P20 is permitted.**

All three required preset paths now pass in this environment, including the previously blocked PM HDF5+FFTW path.

## Validation outcomes

- **CPU-only preset path:** PASS (`36/36` tests).
- **HDF5 preset path:** PASS (`36/36` tests).
- **PM HDF5+FFTW preset path:** PASS (`36/36` tests).

## Exact commands executed

```bash
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure -R "unit_pm_solver|integration_tree_pm_coupling_periodic"
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
ctest --preset test-hdf5-debug --output-on-failure
ctest --preset test-cpu-debug --output-on-failure
```

## PM/FFTW blocker resolution evidence

### Previously failing test: `unit_pm_solver`
- **Current status:** PASS
- **Prior symptom:** assertion failure at `tests/unit/test_pm_solver.cpp:82` (`cosine_similarity > 0.98`).
- **Resolution summary:** PM Fourier-space gradient sign was corrected to return acceleration (`-∇φ`) instead of `+∇φ`.

### Previously failing test: `integration_tree_pm_coupling_periodic`
- **Current status:** PASS
- **Prior symptom:** `rel_l2=18129.9`, required `<= 0.75`.
- **Resolution summary:** FFTW inverse-transform output is now normalized by total cell count (`Nx*Ny*Nz`) before writing PM potential/forces, restoring consistent PM long-range amplitude for TreePM composition.

## Stabilization gate conclusion

With CPU-only, HDF5, and PM HDF5+FFTW preset paths all green, the P19 stabilization blocker is resolved.
