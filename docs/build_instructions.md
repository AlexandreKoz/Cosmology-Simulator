# Build instructions

This page is the authoritative build and dependency workflow for CosmoSim.

## Supported toolchain baseline

- CMake >= 3.24
- C++20 compiler (Clang/GCC/MSVC with C++20 mode)
- Ninja (recommended generator)

Optional dependencies are feature-gated and preset-driven:

- MPI (`COSMOSIM_ENABLE_MPI`)
- HDF5 (`COSMOSIM_ENABLE_HDF5`)
- FFTW (`COSMOSIM_ENABLE_FFTW`)
- CUDA/cuFFT (`COSMOSIM_ENABLE_CUDA`)
- Python + pybind11 (`COSMOSIM_ENABLE_PYTHON`)

## Preset matrix (recommended)

| Preset | Intended use |
|---|---|
| `cpu-only-debug` | Default local development and most CI-facing checks |
| `cpu-debug` | Compatibility configure alias for `cpu-only-debug`, matching the `build-cpu-debug`/`test-cpu-debug` shorthand |
| `cpu-only-release` | Performance sanity on CPU-only builds |
| `hdf5-debug` | Snapshot/restart/provenance schema and I/O work |
| `pm-hdf5-fftw-debug` | PM/TreePM validation with HDF5+FFTW |
| `mpi-hdf5-fftw-debug` | Required distributed feature-complete validation surface (MPI+HDF5+FFTW/FFTW-MPI) |
| `cuda-debug` | Single-rank CUDA PM development surface |
| `mpi-cuda-hdf5-fftw-debug` | Distributed gravity + CUDA infrastructure surface |
| `asan-debug` | Address-sanitizer safety checks |
| `mpi-release` | MPI-only release smoke; not HDF5/FFTW/TreePM restart coverage |
| `cuda-release` | GPU-capable build path |

## Standard developer path (CPU-only)

```bash
cmake --preset cpu-only-debug
cmake --build --preset build-cpu-debug
ctest --preset test-cpu-debug --output-on-failure
```

## HDF5 path

```bash
cmake --preset hdf5-debug
cmake --build --preset build-hdf5-debug
ctest --preset test-hdf5-debug --output-on-failure
```

## PM/TreePM + HDF5 + FFTW path

```bash
cmake --preset pm-hdf5-fftw-debug
cmake --build --preset build-pm-hdf5-fftw-debug
ctest --preset test-pm-hdf5-fftw-debug --output-on-failure
```

## Distributed gravity development path (MPI + HDF5 + FFTW)

```bash
cmake --preset mpi-hdf5-fftw-debug
cmake --build --preset build-mpi-hdf5-fftw-debug
ctest --preset test-mpi-hdf5-fftw-debug --output-on-failure
```

This preset is the required CI lane for accepted distributed PM/TreePM, workflow restart, gas migration,
hydro interface, AMR boundary/reflux, and multi-rank gravity validation. It requires MPI C++ tooling, HDF5,
serial FFTW, and FFTW-MPI (`fftw3_mpi`; Debian/Ubuntu package `libfftw3-mpi-dev`). The `mpi-release` preset
is intentionally MPI-only smoke coverage and must not be cited as FFT-PM, HDF5, restart-topology, or
parallel-HDF5 evidence.

## MPI and GPU paths

`parallel.gpu_devices > 0` now acts as an explicit runtime request for the CUDA PM assignment/interpolation path. This request is **not** silently downgraded: the runtime validates visible devices and raises a clear error if CUDA was requested but no compatible runtime devices are present.

```bash
cmake --preset mpi-release
cmake --build --preset build-mpi-release
ctest --preset test-mpi-release --output-on-failure

cmake --preset cuda-release
cmake --build --preset build-cuda-release
ctest --preset test-cuda-release --output-on-failure
```

> GPU preset availability depends on CUDA toolkit/compiler/runtime environment.

## Python bindings path

```bash
cmake -S . -B build/py \
  -DCOSMOSIM_ENABLE_PYTHON=ON \
  -DCOSMOSIM_ENABLE_HDF5=ON \
  -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"
cmake --build build/py --target cosmosim_python_package
ctest --test-dir build/py --output-on-failure
```

## Generated build metadata

Every configure run emits:

- `cosmosim_feature_summary.txt`
- `cosmosim_build_metadata.json`

in the build directory. Keep these artifacts for reproducibility and incident debugging.

## Common troubleshooting

- Missing HDF5: set `-DHDF5_ROOT=/path/to/hdf5`.
- Missing FFTW in feature preset: set `-DPKG_CONFIG_PATH=/path/to/fftw/lib/pkgconfig`.
- Custom toolchain paths: copy `CMakeUserPresets.json.example` to `CMakeUserPresets.json` and edit locally (do not commit).
