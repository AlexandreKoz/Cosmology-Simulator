# GPU enablement strategy (phase 1)

This document records the first auditable GPU path for PM gravity in CosmoSim.

## Scope and boundaries

- **Targeted kernels in this phase:** PM CIC assignment and CIC force interpolation.
- **Solver stage intentionally kept on host in this phase:** FFT + Poisson + spectral gradient.
- **Execution policy is explicit, never implicit:** use `PmSolveOptions::execution_policy` and `PmSolveOptions::data_residency`.

## Ownership and synchronization policy

- `execution_policy = host_serial` keeps all state in host memory.
- `execution_policy = cuda` requires `data_residency = kPreferDevice`; this avoids hidden offload side effects.
- Runtime workflow integration now maps `parallel.gpu_devices` onto explicit rank-local CUDA device assignment; device selection is round-robin across the configured visible-device pool.
- Synchronization points are explicit and auditable:
  1. host → device particle upload,
  2. device assignment kernel completion,
  3. device → host density download,
  4. host Poisson solve completion,
  5. host → device force upload,
  6. device interpolation kernel completion,
  7. device → host acceleration download.

## Numerical conventions and invariants

- Variable conventions remain unchanged from CPU PM path:
  - positions in `box_size_mpc_comoving`, periodic wrap,
  - CIC assignment/interpolation,
  - comoving Poisson prefactor `-4π G_code a^2`.
- The CUDA path is designed to be numerically comparable to host within tolerance and keeps the same assignment convention.
- Mean-density subtraction, Fourier-space Poisson solve, and gradient construction remain the host reference implementation.

## Profiling fields

`PmProfileEvent` now separates:

- `transfer_h2d_ms`
- `transfer_d2h_ms`
- `device_kernel_ms`

These complement the existing PM timing breakdown and are intended for transfer-vs-kernel cost attribution.

## Build and fallback behavior

- CUDA support uses modern CMake CUDA language support and `CUDA::cudart`/`CUDA::cufft` imported targets.
- When `COSMOSIM_ENABLE_CUDA=OFF`, the CPU reference path remains fully available.
- Requesting `execution_policy = cuda` in a non-CUDA build throws a clear runtime error.

## Current limitations (explicit)

- FFT/Poisson is host-only in this phase; no cuFFT Poisson path yet.
- No persistent long-lived PM device residency across timesteps yet; this phase focuses on explicit ownership and correctness.
- Multi-GPU and MPI+GPU overlap are out of scope in phase 1 and intentionally not implied.
