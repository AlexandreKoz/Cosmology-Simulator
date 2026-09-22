# Build instructions

This page is the authoritative build and dependency workflow for CosmoSim.

## Supported toolchain baseline

- CMake >= 3.24
- C++20 compiler **and standard library** with the facilities used by CHUÍ. Configure probes `#include <span>` / `std::span` instead of trusting `CMAKE_CXX_STANDARD` alone. GCC 9 / libstdc++ 9 is therefore unsupported for the current tree; newer toolchains are accepted by capability rather than by a hard-coded compiler brand/version.
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

For a known CMake build, `./chui run CONFIG` asks the C++ harness preflight to parse the authoritative `.param.txt` and return `parallel.mpi_ranks_expected`. If the value is greater than one, the launcher composes the MPI command automatically; Python does not implement a second scientific parser and never rewrites the config. An explicit `--mpi N` remains supported but must agree with the typed config. The MPI launcher is taken from the selected build's `CMakeCache.txt` when available, otherwise `mpiexec`/`mpirun` is resolved from `PATH`. If the selected build explicitly has `COSMOSIM_ENABLE_MPI=OFF`, automatic or explicit multi-rank launch fails before starting multiple independent serial processes. Arguments after `--` remain expert MPI-launcher passthrough:

```bash
./chui run CONFIG --preset mpi-hdf5-fftw-release --mpi 8
./chui run CONFIG --preset mpi-hdf5-fftw-debug --mpi 2 -- --bind-to core
```

Native progress reporting is owned by the C++ runtime, not the launcher. The default direct harness and launcher paths emit bounded rank-0 status. Under MPI, `run_directory` is the shared/logical run location while `rank_directory` identifies rank 0's rank-local artifacts. With the normal Parallel-HDF5 science layout a committed snapshot record names the one ordinary `snap_###.hdf5` analysis product plus its transactional `.complete` marker. Explicit legacy `sharded` mode instead reports the logical member count and labels rank 0's member as a member, never as the complete science state. Presentation controls are:

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

The supported HDF5 source/API range is **1.10.x through 1.14.x** and CHUÍ consumes the **HDF5 C API**; the unused HDF5 C++ bindings are not a source or installed-package dependency. Object inspection uses stable handle/type queries rather than version-sensitive unversioned `H5Oget_info_by_name` signatures. CMake fails closed below 1.10 and on unqualified HDF5 2.x.

`HDF5_IS_PARALLEL` is recorded only as provider metadata. Production capability truth comes from a compile/link probe against the exact selected headers/libraries and MPI stack using `H5Pset_fapl_mpio` and `H5Pset_dxpl_mpio`. On distro systems that install serial and OpenMPI HDF5 together, CHUÍ prefers the parallel C wrapper (`h5pcc`, `h5pcc.openmpi`, or `h5pcc.mpich`) unless the user explicitly selected `HDF5_ROOT`, `HDF5_DIR`, or `HDF5_C_COMPILER_EXECUTABLE`. `feature_hdf5_parallel=true` therefore means the MPIO calls actually linked, not merely that aggregate CMake metadata claimed parallel support. Serial HDF5 is valid for serial output and the explicit legacy `output.snapshot_layout=sharded` MPI compatibility backend.

The MPIO `CheckCXXSourceCompiles` result is deliberately invalidated immediately before each configure-time capability check. This prevents a cached `TRUE`/`FALSE` from a previously selected HDF5 prefix or MPI toolchain from suppressing the new proof. The configure regression poisons the cached result to `FALSE`, requires the `Performing Test COSMOSIM_HDF5_MPIO_CAPABILITY_PROBE` test to execute and succeed on a known Parallel-HDF5 preset, and—where a serial HDF5 C wrapper is installed—verifies that `COSMOSIM_REQUIRE_PARALLEL_HDF5=ON` rejects it.

`COSMOSIM_REQUIRE_PARALLEL_HDF5=ON` turns that capability into a configure-time hard requirement. The production `mpi-hdf5-fftw-debug` and `mpi-hdf5-fftw-release` presets enable it, so they cannot configure successfully against serial HDF5 and later fail during a normal `snapshot_layout=auto` run. `mpi-serial-hdf5-fftw-debug` is the explicitly named compatibility configure surface when sharded output is intentionally being tested.

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
serial FFTW, and FFTW-MPI (`fftw3_mpi`; Debian/Ubuntu package `libfftw3-mpi-dev`). To qualify the default
analysis-ready one-file MPI science snapshot backend, the HDF5 selected by this preset must additionally be
MPI-enabled Parallel HDF5; a serial-HDF5 MPI build is reported explicitly and cannot silently masquerade as
that capability. The `mpi-release` preset is intentionally MPI-only smoke coverage and must not be cited as
FFT-PM, HDF5, restart-topology, or Parallel-HDF5 evidence.


## Toolchain and dependency provenance

Every configure reports the selected C++ compiler path/version and the resolved MPI, HDF5, FFTW, and FFTW-MPI providers that are relevant to the chosen preset. `cosmosim_build_metadata.json` also records `feature_hdf5_parallel`, so launcher/runtime diagnostics can distinguish a normal HDF5 build from a build capable of the one-file MPI science backend. If a Conda MPI wrapper is combined with a system compiler/HDF5 stack, CMake emits a mixed-prefix warning instead of silently overriding the user's selection. This is a risk diagnostic, not a blanket Conda prohibition.

For a normal native Linux stack, prefer one coherent compiler/MPI/HDF5/FFTW provider family (distribution packages or one HPC module stack). If an explicit custom stack is required, use local `CMakeUserPresets.json` overrides and verify the configure provenance before building.

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
