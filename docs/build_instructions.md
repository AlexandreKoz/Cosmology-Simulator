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

## Running CHUÍ

The repository-root `chui` launcher is the normal POSIX/Linux/WSL operator interface. It uses only the Python 3 standard library and remains a thin orchestration layer over the existing `cosmosim_harness`; scientific configuration stays exclusively in the typed `.param.txt` path.

After building one runnable preset, use:

```bash
./chui run configs/chui_firstlight_32_smoke.param.txt
```

If multiple preset build trees contain `cosmosim_harness`, the launcher fails with the available choices rather than executing an arbitrary binary. Select one deterministically with `--preset`, or bypass preset discovery with `--exe`. For non-quiet launches, a single `[CHUI][LAUNCH]` line records the selected executable, build directory, serial/MPI mode, and requested rank count before `exec` replaces the wrapper process:

```bash
./chui run CONFIG --preset pm-hdf5-fftw-debug
./chui run CONFIG --exe /absolute/path/to/cosmosim_harness
```

For MPI runs, `--mpi N` prepends the MPI launcher recorded in the selected build's `CMakeCache.txt` when available, otherwise it uses `mpiexec`/`mpirun` from `PATH`. If the selected build metadata explicitly says `COSMOSIM_ENABLE_MPI=OFF`, the launcher fails instead of starting multiple independent serial processes. It does not inspect or rewrite `parallel.mpi_ranks_expected`; the authoritative runtime still validates communicator size before expensive simulation work. Arguments after `--` are passed to the MPI launcher:

```bash
./chui run CONFIG --preset mpi-hdf5-fftw-release --mpi 8
./chui run CONFIG --preset mpi-hdf5-fftw-debug --mpi 2 -- --bind-to core
```

Native progress reporting is owned by the C++ runtime, not the launcher. The default direct harness and launcher paths emit bounded rank-0 status. Under MPI, `run_directory` is the shared/logical run location while `rank_directory` identifies rank 0's rank-local artifacts. A committed snapshot line reports the logical member count and `.complete` manifest; the rank-0 member is labeled explicitly rather than being presented as the complete snapshot. Presentation controls are:

```text
--quiet
--status-every N
--status-seconds SEC
```

For example:

```bash
./chui run CONFIG --status-every 5 --status-seconds 20
./chui run CONFIG --quiet
```

The direct executable remains fully supported:

```bash
./build/<preset>/cosmosim_harness CONFIG
```

`./chui --help` and `./chui run --help` document the complete launcher surface. The root launcher is currently POSIX/Linux/WSL-first; platforms that do not support it can invoke `cosmosim_harness` directly. External `chui_telemetry.py` process telemetry remains optional diagnostic tooling and is never launched implicitly.

## HDF5 path

The supported HDF5 source/API range is **1.10.x through 1.14.x**. Object inspection uses stable handle/type queries rather than version-sensitive unversioned `H5Oget_info_by_name` signatures. CMake fails closed below 1.10 and on unqualified HDF5 2.x.

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


## OpenMP shared-memory execution

OpenMP is enabled by default when CMake can find `OpenMP::OpenMP_CXX` and may be
disabled explicitly with `-DCOSMOSIM_ENABLE_OPENMP=OFF`. In an OpenMP build,
`parallel.omp_threads = N` controls the actual runtime team size; `0` requests
the runtime default. A non-OpenMP build remains functional but rejects thread
requests above one. Runtime capability output reports the compiled, requested,
and active thread state. Current real OpenMP work includes the active Tree
gravity traversal and production FFT analysis.

## Install/export consumer

CosmoSim installs public headers, the generated `build_config.hpp`, libraries,
and CMake package metadata. A downstream CMake project can use:

```cmake
find_package(CosmoSim CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE CosmoSim::core CosmoSim::analysis)
```

A typical local validation is:

```bash
cmake --install build/cpu-only-debug --prefix /tmp/cosmosim-install
cmake -S /path/to/consumer -B /tmp/cosmosim-consumer \
  -DCMAKE_PREFIX_PATH=/tmp/cosmosim-install
cmake --build /tmp/cosmosim-consumer
```
