# CHUÍ — Consolidated source patch (2026-10-28)

**Base:** Uploaded source archive `Cosmology-Simulator-main(24)(5).zip` (repository-relative `src/`, `tests/` paths). **This is an overlay patch, not a full repository.** Every file from all three previous patches in the conversation is included, alongside the new CTest-related corrections. The Round-5 previous patch files are included even where byte-identical to the uploaded source ZIP, as requested.

## Apply (Linux)

Save and optionally back up local changes first. Run from the **CHUÍ repository root**:

```sh
cd ~/CHUI
unzip -l ~/Downloads/CHUI_Consolidated_Round5_and_CTest_Fixes_2026-10-28.zip
unzip -o ~/Downloads/CHUI_Consolidated_Round5_and_CTest_Fixes_2026-10-28.zip
```

The overlay uses `src/...` and `tests/...` paths; there is no `Cosmology-Simulator-main/` nesting. It overwrites those paths in place. Adjust the ZIP path if necessary.

If `CMAKE_SUPPRESS_REGENERATION=ON` is set to avoid the observed CMake regeneration loop, **reconfigure explicitly** after applying these source changes:

```sh
cmake -S . -B build/mpi-hdf5-fftw-debug \
  -DCMAKE_C_COMPILER=/usr/bin/gcc-13 \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++-13 \
  -DMPI_CXX_COMPILER=/usr/bin/mpicxx \
  -DCMAKE_SUPPRESS_REGENERATION=ON
cmake --build build/mpi-hdf5-fftw-debug --parallel 4
OMP_NUM_THREADS=1 ctest --test-dir build/mpi-hdf5-fftw-debug \
  --parallel 2 --output-on-failure -E '(_8_rank|_eight_rank|_np8)' \
  2>&1 | tee chui_ctest_after_consolidated_patch.log
```

Run to completion; if MPI launch becomes too busy run `--parallel 1`. Do not judge correctness only from a green CTest summary: physical convergence remains required. Eight-rank tests are omitted on the reported local machine because OpenMPI reported insufficient slots, not because those tests are inherently invalid.

## Previously issued patches (all included verbatim except the tree-gravity patch, which is *corrected further*)

**Round 5 / four fixes:**
- `src/physics/black_hole_agn.cpp`
- `src/physics/stellar_evolution.cpp`
- `src/physics/tracer_support.cpp`
- `src/workflows/time_coordinator.cpp`
- `tests/unit/test_black_hole_agn.cpp`
- `tests/unit/test_gravity_force_endpoint.cpp`
- `tests/unit/test_stellar_evolution.cpp`
- `tests/unit/test_time_integration.cpp`
- `tests/unit/test_tracer_support.cpp`

**Separate mass-generation test fixture:**
- `tests/unit/test_particle_mass_generation.cpp` (species counts in fixture).

**Separate in-place FFT overflow test fixture:**
- `tests/unit/test_tree_gravity.cpp` (updated origin patch + additionally corrected memory estimate and overflow message expectations, based on reproduced local failure).

## New changes and reasons

| File | Change |
|---|---|
| `tests/validation/test_dmo_zeldovich_workflow.cpp` | Gather a concrete time-bin vector populated through `particleTimeBin()` for compact DMO; prevents null MPI send buffer. |
| `src/gravity/pm_solver.cpp` | Correct isolated PM origin Green function: analytic split-kernel limit, midpoint-cell average for unsplit `1/r`; fixes non-attractive nearest-cell finite-difference behavior. **Needs numerical cross-validation with MPI/FFTW.** |
| `tests/integration/test_tree_pm_fine_tick_caches.cpp` | Enable group walks before asserting grouped topology cache reuse. |
| `tests/unit/test_output_cadence.cpp` | Establish completed snapshot step 4 before checking checkpoint step 6. |
| `tests/integration/test_restart_equivalence_dm_only.cpp` | Respect compact DMO drift sidecar storage; materialize for explicit lane corruption test. |
| `tests/integration/test_distributed_ic_reader_mpi.cpp` | Permit zero routing batches for serial one-rank fallback; preserve strict multi-rank collective budget and print actual counters on failure. |
| `tests/integration/test_configure_presets.cmake.in` | Forward working C and C++ compiler paths into nested CMake preset checks. |
| `tests/integration/test_installed_package_consumer.cmake.in` | Forward compiler/MPI toolchain hints and prefer parallel HDF5 in installed consumer probe. |
| `tests/integration/test_runtime_truth_ctest_labels.cmake.in` | Accept current CTest registration name quoting / bracket syntax; preserve all label checks. |
| `tests/integration/test_star_formation_source_runtime_mpi.cpp` | Synchronize gas identity through authoritative replacement API; lower fixture softening to satisfy certified epsilon/split-scale envelope. |
| `tests/integration/test_restart_checkpoint_roundtrip.cpp` | Supply an authoritative gas-cell time-bin scheduler when writing checkpoints that contain gas cells. |
| `tests/integration/test_distributed_gas_cell_migration_mpi.cpp` | Give ranks the same collective checkpoint output roots. |
| `tests/unit/test_config_parser.cpp` | Use source-free cosmo-cube DMO fixture for hierarchical policy test; fix stale message substring. |
| `tests/integration/test_hierarchical_multirate_convergence.cpp` | Hold the cosmological midpoint fixed over the whole PM cycle rather than recomputing it each fine block (previously reversed drift interval). |

## Verification carried out here

- GCC 14, CPU-only Debug configuration: **8/8 focused CTest cases passed**: `unit_output_cadence`, `unit_config_parser`, `unit_tree_gravity`, `unit_particle_mass_generation`, `unit_gravity_force_endpoint`, `integration_tree_pm_fine_tick_caches`, `integration_hierarchical_multirate_convergence`, `integration_runtime_truth_ctest_labels`.
- Isolated-PM CPU-only linked smoke using the modified implementation: left/right nearest-cell forces were attractive and antisymmetric (`+0.156582` and `-0.156582`, relative antisymmetry error ~`2.5e-16`). **This is not a complete FFTW+MPI regression test and does not establish convergence across meshes.**
- Original full MPI+HDF5+FFTW user build and CTest results refer to an **earlier** source state, before the new changes. I cannot reproduce full MPI/HDF5/FFTW tests in this container (dependencies unavailable).

## Still unresolved — intentionally not papered over

- MPI `MPI_ERR_TRUNCATE` in distributed reference TreePM and isolated PM workflow (`MPI_Allreduce` mismatch). No verified code-level cause yet.
- 24 multi-rank distributed-IC collective count failures; added diagnostic values, not relaxed count.
- Hierarchical TreePM trajectory RMS scientific gate; PM refinement nonmonotonic two-body accuracy; distributed grouped-walk direct-pair accounting; distributed SFC fault-injection collective handling; AMR patch migration scheduler mismatch; potential routing request telemetry; distributed TreePM request batch budget.
- Eight-rank launch needs adequate MPI slots or explicit (careful) oversubscription; legacy multi-rank Cartesian hydro requires a supported production AMR fixture.
- Any new MPI/HDF5 behavior changed in this overlay requires real local compilation and tests before production scientific runs.

**Do not mark the full suite green or CHUÍ production-qualified based only on this patch.**
