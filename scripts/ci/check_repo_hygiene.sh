#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

echo "[hygiene] checking repository-root filename safety"

bad_root_entries=()
while IFS= read -r entry; do
  name="${entry#./}"
  if [[ ! "$name" =~ ^[A-Za-z0-9._-]+$ ]]; then
    bad_root_entries+=("$name")
  fi
done < <(find . -maxdepth 1 -mindepth 1 \( -type f -o -type d \) -print)

if (( ${#bad_root_entries[@]} > 0 )); then
  echo "[hygiene] ERROR: root entries with non-automation-safe names detected:" >&2
  printf '  - %s\n' "${bad_root_entries[@]}" >&2
  echo "[hygiene] Move/rename these entries or quarantine with explicit rationale." >&2
  exit 1
fi

echo "[hygiene] checking recursive filename and artifact safety"

bad_nested_entries=()
while IFS= read -r entry; do
  name="${entry##*/}"
  if [[ ! "$name" =~ ^[A-Za-z0-9._-]+$ ]] ||
     [[ "$name" == *Zone.Identifier* ]] ||
     [[ "$name" == *#Uf03aZone.Identifier* ]]; then
    bad_nested_entries+=("${entry#./}")
  fi
done < <(find . -mindepth 1 -path './.git' -prune -o \
  \( -type f -o -type d \) -print)

if (( ${#bad_nested_entries[@]} > 0 )); then
  echo "[hygiene] ERROR: nested unsafe/ADS filenames detected:" >&2
  printf '  - %s\n' "${bad_nested_entries[@]}" >&2
  exit 1
fi

echo "[hygiene] checking recursively for generated build/runtime artifacts"

bad_artifacts=()
while IFS= read -r entry; do
  bad_artifacts+=("${entry#./}")
done < <(find . -path './.git' -prune -o \
  \( -type d \( \
       -name build -o -name 'build-*' -o -name 'cmake-build-*' -o \
       -name CMakeFiles -o -name Testing -o -name _deps -o \
       -name integration_outputs -o -name validation_outputs -o \
       -name test_outputs -o -name output -o -name .pytest_cache -o \
       -name __pycache__ -o -name .mypy_cache -o -name .ruff_cache -o \
       -name .ipynb_checkpoints -o -name vcpkg_installed \
     \) -print -prune \) -o \
  \( -type f \( \
       -name CMakeCache.txt -o -name CTestTestfile.cmake -o \
       -name CMakeUserPresets.json -o -name '*.o' -o -name '*.obj' -o \
       -name '*.a' -o -name '*.so' -o -name '*.so.*' -o -name '*.dll' -o \
       -name '*.dylib' -o -name '*.exe' -o -name '*.pyc' -o -name '*.pyo' -o \
       -name '*.gcda' -o -name '*.gcno' -o -name '*.profraw' -o \
       -name '*.profdata' -o -name '*.part' -o -name core -o -name 'core.*' \
     \) -print \))

if (( ${#bad_artifacts[@]} > 0 )); then
  echo "[hygiene] ERROR: generated build/runtime artifacts are present:" >&2
  printf '  - %s\n' "${bad_artifacts[@]}" >&2
  exit 1
fi

echo "[hygiene] checking IC subsystem structural guardrails"
if find src/io -maxdepth 2 -type f \
    \( -name 'ic_utils.cpp' -o -name 'ic_utils.hpp' -o -name 'ic_utils.hh' \) \
    | grep -q .; then
  echo "[hygiene] ERROR: generic IC utility dumping-ground file detected" >&2
  exit 1
fi
ic_reader_lines="$(wc -l < src/io/ic_reader_file_set.cpp)"
if (( ic_reader_lines > 4200 )); then
  echo "[hygiene] ERROR: src/io/ic_reader_file_set.cpp has $ic_reader_lines lines; split responsibilities before exceeding 4200" >&2
  exit 1
fi
for ic_unit in \
  src/io/ic_canonical_bundle.cpp \
  src/io/ic_conversion_catalog.cpp \
  src/io/ic_reader_session.cpp \
  src/io/ic_record_codec.cpp \
  src/io/ic_byte_codec.cpp \
  src/io/ic_sha256.cpp; do
  unit_lines="$(wc -l < "$ic_unit")"
  if (( unit_lines > 1200 )); then
    echo "[hygiene] ERROR: $ic_unit has $unit_lines lines and is becoming a replacement monolith" >&2
    exit 1
  fi
done

echo "[hygiene] checking required preset names"
required_presets=(
  "cpu-only-debug"
  "hdf5-debug"
  "pm-hdf5-fftw-debug"
  "mpi-hdf5-fftw-debug"
  "mpi-release"
  "build-cpu-debug"
  "build-hdf5-debug"
  "build-pm-hdf5-fftw-debug"
  "build-mpi-hdf5-fftw-debug"
  "build-mpi-release"
  "test-cpu-debug"
  "test-hdf5-debug"
  "test-pm-hdf5-fftw-debug"
  "test-mpi-hdf5-fftw-debug"
  "test-mpi-release"
)

preset_list="$(cmake --list-presets=all 2>&1)"
for preset in "${required_presets[@]}"; do
  if ! grep -Fq "\"$preset\"" <<<"$preset_list"; then
    echo "[hygiene] ERROR: required preset '$preset' is missing" >&2
    exit 1
  fi
done



echo "[hygiene] checking FFTW MPI CI dependency wiring"
if ! grep -Fq 'configure_preset: mpi-hdf5-fftw-debug' .github/workflows/ci.yml; then
  echo "[hygiene] ERROR: mpi-hdf5-fftw-debug CI row is missing" >&2
  exit 1
fi
if ! grep -Fq 'libfftw3-mpi-dev' .github/workflows/ci.yml; then
  echo "[hygiene] ERROR: CI FFTW MPI paths must install libfftw3-mpi-dev, not only libfftw3-dev" >&2
  exit 1
fi
if ! grep -Fq 'REQUIRE_MPI' cmake/cosmosim_find_fftw.cmake; then
  echo "[hygiene] ERROR: FFTW discovery must expose a REQUIRE_MPI path for distributed PM presets" >&2
  exit 1
fi

echo "[hygiene] repository hygiene checks passed"

echo "[hygiene] checking CI shell script syntax"
bash -n ./scripts/ci/run_preset_pipeline.sh
bash -n ./scripts/ci/enforce_infra_gates.sh
bash -n ./scripts/ci/run_reproducibility_gate.sh
bash -n ./scripts/ci/run_stage1_runtime_truth_gate.sh

exit 0
