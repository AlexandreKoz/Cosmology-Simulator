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
       -name '*.a' -o -name '*.lib' -o -name '*.so' -o -name '*.so.*' -o \
       -name '*.dll' -o -name '*.dylib' -o -name '*.exe' -o -name '*.pdb' -o \
       -name '*.pyc' -o -name '*.pyo' -o \
       -name '*.gcda' -o -name '*.gcno' -o -name '*.profraw' -o \
       -name '*.profdata' -o -name '*.part' -o -name core -o -name 'core.*' \
     \) -print \))

if (( ${#bad_artifacts[@]} > 0 )); then
  echo "[hygiene] ERROR: generated build/runtime artifacts are present:" >&2
  printf '  - %s\n' "${bad_artifacts[@]}" >&2
  exit 1
fi

echo "[hygiene] checking authoritative simulation configuration format"
while IFS= read -r config_file; do
  config_name="${config_file##*/}"
  if [[ ! "$config_name" =~ ^[a-z0-9]+(_[a-z0-9]+)*\.param\.txt$ ]]; then
    echo "[hygiene] ERROR: simulation config must use lower_snake_case .param.txt: ${config_file#./}" >&2
    exit 1
  fi
done < <(find configs -type f ! -name '*.md' -print | sort)

if find configs -type f \( -iname '*.yaml' -o -iname '*.yml' -o -iname '*.json' \) | grep -q .; then
  echo "[hygiene] ERROR: YAML/JSON simulation input found under configs/; use .param.txt" >&2
  exit 1
fi
cmake -DCOSMOSIM_SOURCE_DIR="$repo_root" \
  -P "$repo_root/cmake/check_config_aliases.cmake"

echo "[hygiene] checking IC subsystem structural guardrails"
if find src/io -maxdepth 2 -type f \
    \( -name 'ic_utils.cpp' -o -name 'ic_utils.hpp' -o -name 'ic_utils.hh' \) \
    | grep -q .; then
  echo "[hygiene] ERROR: generic IC utility dumping-ground file detected" >&2
  exit 1
fi
while IFS= read -r ic_unit; do
  unit_lines="$(wc -l < "$ic_unit")"
  if (( unit_lines > 2000 )); then
    echo "[hygiene] ERROR: $ic_unit has $unit_lines lines; IC source units must remain at or below 2000 lines" >&2
    exit 1
  fi
done < <(find src/io -maxdepth 1 -type f -name 'ic_*.cpp' | sort)

common_header_users="$(
  grep -RIl --exclude='ic_file_set_common.hpp' \
    'io/internal/ic_file_set_common.hpp' include src tools 2>/dev/null || true
)"
while IFS= read -r user; do
  [[ -z "$user" ]] && continue
  case "$user" in
    src/io/ic_file_set_manifest.cpp|src/io/ic_stream_ingestion.cpp|\
    src/io/ic_serial_ingestion.cpp|src/io/ic_distributed_audit.cpp|\
    src/io/ic_distributed_ingestion.cpp) ;;
    *)
      echo "[hygiene] ERROR: private IC file-set header leaked into $user" >&2
      exit 1
      ;;
  esac
done <<<"$common_header_users"

if [[ ! -f src/io/ic_mpi_collectives.cpp ]] ||
   [[ ! -f src/io/internal/ic_mpi_collectives.hpp ]]; then
  echo "[hygiene] ERROR: centralized IC MPI collective instrumentation is missing" >&2
  exit 1
fi
for distributed_source in   src/io/ic_distributed_ingestion.cpp   src/io/ic_distributed_audit.cpp   src/io/ic_failure_protocol.cpp; do
  if grep -Eq 'MPI_(Allreduce|Bcast|Gather|Gatherv|Alltoall|Alltoallv)[[:space:]]*\('       "$distributed_source"; then
    echo "[hygiene] ERROR: $distributed_source bypasses IC MPI collective instrumentation" >&2
    exit 1
  fi
done

echo "[hygiene] checking source-package guardrails"
for package_script in \
  scripts/ci/source_package_common.sh \
  scripts/ci/package_source_zip.sh \
  scripts/ci/test_source_package_completeness.sh; do
  if [[ ! -f "$package_script" ]]; then
    echo "[hygiene] ERROR: required source-package script is missing: $package_script" >&2
    exit 1
  fi
  bash -n "$package_script"
done
if grep -Eq -- "--exclude=(['\"])?core(\.\*)?(['\"])?([[:space:]]|$)" \
    scripts/ci/source_package_common.sh scripts/ci/package_source_zip.sh; then
  echo "[hygiene] ERROR: source packaging must not exclude generic core/core.* paths" >&2
  exit 1
fi
if ! grep -Fq -- "-type f" scripts/ci/source_package_common.sh ||
   ! grep -Fq -- "-name core" scripts/ci/source_package_common.sh ||
   ! grep -Fq -- "-name 'core.*'" scripts/ci/source_package_common.sh; then
  echo "[hygiene] ERROR: crash-dump-name checks must be restricted to regular files" >&2
  exit 1
fi

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
