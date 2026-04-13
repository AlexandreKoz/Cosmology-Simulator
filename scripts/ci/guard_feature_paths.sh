#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

run_path() {
  local configure_preset="$1"
  local build_preset="$2"
  local test_preset="$3"
  local test_regex="$4"

  echo "[guard] configure: $configure_preset"
  cmake --preset "$configure_preset"

  echo "[guard] build: $build_preset"
  cmake --build --preset "$build_preset"

  if [[ -n "$test_regex" ]]; then
    echo "[guard] targeted tests: $test_preset (regex: $test_regex)"
    ctest --preset "$test_preset" -R "$test_regex"
  else
    echo "[guard] full tests: $test_preset"
    ctest --preset "$test_preset"
  fi
}

# Baseline CPU path.
run_path "cpu-only-debug" "build-cpu-debug" "test-cpu-debug" ""

# Real feature path: PM + HDF5 + FFTW.
run_path "pm-hdf5-fftw-debug" "build-pm-hdf5-fftw-debug" "test-pm-hdf5-fftw-debug" "integration_pm_periodic_mode|unit_snapshot_hdf5_schema"
