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

echo "[hygiene] repository hygiene checks passed"

echo "[hygiene] checking CI shell script syntax"
bash -n ./scripts/ci/run_preset_pipeline.sh
bash -n ./scripts/ci/enforce_infra_gates.sh
bash -n ./scripts/ci/run_reproducibility_gate.sh

exit 0
