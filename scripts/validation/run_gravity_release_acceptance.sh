#!/usr/bin/env bash
set -euo pipefail

build_dir="${1:-build/cpu-only-debug}"
poster_ranks="${COSMOSIM_POSTER_RANKS:-4}"
manifest_out="${COSMOSIM_GRAVITY_ACCEPTANCE_MANIFEST:-gravity_acceptance_manifest.json}"
evidence_dir="${COSMOSIM_GRAVITY_EVIDENCE_DIR:-gravity_acceptance_evidence}"

if [[ ! -d "$build_dir" ]]; then
  echo "gravity acceptance: build directory not found: $build_dir" >&2
  exit 2
fi
mkdir -p "$evidence_dir"

test_list="$(ctest --test-dir "$build_dir" -N)"
has_test() {
  local name="$1"
  grep -Eq "Test +#[0-9]+: ${name}$" <<<"$test_list"
}

run_ctest_regex() {
  local regex="$1"
  ctest --test-dir "$build_dir" --output-on-failure -R "$regex"
}

evidence_args=()
rank_args=(--rank-count "1")

echo "== serial gravity contract gates =="
run_ctest_regex '^(unit_tree_gravity|unit_pm_solver|integration_tree_gravity_vs_direct|integration_tree_pm_coupling_periodic|integration_star_formation_source_runtime)$'

fftw_enabled=false
pm_backend="unknown"
pm_backend_capability="unavailable"
fftw_version="unavailable"
if grep -q '^COSMOSIM_ENABLE_FFTW:BOOL=ON$' "$build_dir/CMakeCache.txt" 2>/dev/null; then
  fftw_enabled=true
  pm_backend="fftw"
  pm_backend_capability="production_fftw"
  fftw_version="enabled_version_not_extracted"
elif grep -q '^COSMOSIM_ENABLE_FFTW:BOOL=OFF$' "$build_dir/CMakeCache.txt" 2>/dev/null; then
  pm_backend="naive_dft"
  pm_backend_capability="diagnostic_naive_dft"
fi

if $fftw_enabled && has_test validation_tree_pm_ewald_accuracy; then
  echo "== TreePM Ewald/softening envelope evidence =="
  for ratio in 0 0.05 0.10 0.20; do
    evidence_file="$evidence_dir/ewald_eps_over_rs_${ratio//./p}.json"
    python3 tools/gravity_ctest_evidence.py \
      --build-dir "$build_dir" \
      --test validation_tree_pm_ewald_accuracy \
      --output "$evidence_file" \
      --profile "${COSMOSIM_ACCEPTANCE_PROFILE_ID:-unknown_not_extracted}" \
      --rank-count 1 \
      --fft-backend "$pm_backend" \
      --softening-ratio "$ratio" \
      --env "COSMOSIM_EWALD_SOFTENING_RATIO=$ratio"
    evidence_args+=(--evidence "$evidence_file")
  done
else
  echo "Production FFTW Ewald lane unavailable in this build; no FFTW/Ewald certification is claimed."
fi

mpi_version="unavailable"
if command -v mpiexec >/dev/null 2>&1 && grep -q '^COSMOSIM_ENABLE_MPI:BOOL=ON$' "$build_dir/CMakeCache.txt" 2>/dev/null; then
  mpi_version="$(mpiexec --version 2>/dev/null | head -n1 || true)"
  [[ -n "$mpi_version" ]] || mpi_version="enabled_version_not_extracted"
  echo "== registered MPI gravity/PM/restart gates =="
  mapfile -t mpi_tests < <(
    awk '/Test +#[0-9]+:/ {print $NF}' <<<"$test_list" | \
      grep -E '(tree_pm|treepm|pm.*slab|zeldovich|restart|gas).*(mpi_np[234]|two[_-]rank|three[_-]rank|four[_-]rank)|(mpi_np[234]|two[_-]rank|three[_-]rank|four[_-]rank).*(tree_pm|treepm|pm|zeldovich|restart|gas)' || true
  )
  if [[ ${#mpi_tests[@]} -eq 0 ]]; then
    echo "MPI is enabled but no matching gravity/PM/restart CTest lanes were discovered." >&2
    exit 3
  fi
  declare -A rank_seen=()
  for test_name in "${mpi_tests[@]}"; do
    rank=0
    if [[ "$test_name" =~ _mpi_np([0-9]+) ]]; then
      rank="${BASH_REMATCH[1]}"
    elif [[ "$test_name" =~ two[_-]rank ]]; then
      rank=2
    elif [[ "$test_name" =~ three[_-]rank ]]; then
      rank=3
    elif [[ "$test_name" =~ four[_-]rank ]]; then
      rank=4
    fi
    evidence_file="$evidence_dir/${test_name}.json"
    python3 tools/gravity_ctest_evidence.py \
      --build-dir "$build_dir" \
      --test "$test_name" \
      --output "$evidence_file" \
      --profile "${COSMOSIM_ACCEPTANCE_PROFILE_ID:-unknown_not_extracted}" \
      --rank-count "$rank" \
      --fft-backend "$pm_backend"
    evidence_args+=(--evidence "$evidence_file")
    if [[ "$rank" -gt 0 && -z "${rank_seen[$rank]:-}" ]]; then
      rank_args+=(--rank-count "$rank")
      rank_seen[$rank]=1
    fi
  done
  if [[ "$poster_ranks" =~ ^[0-9]+$ && "$poster_ranks" -gt 1 && -z "${rank_seen[$poster_ranks]:-}" ]]; then
    echo "Poster rank count ${poster_ranks} was not covered by a discovered MPI CTest lane; production acceptance remains incomplete." >&2
  fi
else
  echo "MPI runtime/tests unavailable in this build; no MPI acceptance is claimed."
fi

hdf5_version="unavailable"
if grep -q '^COSMOSIM_ENABLE_HDF5:BOOL=ON$' "$build_dir/CMakeCache.txt" 2>/dev/null; then
  hdf5_version="enabled_version_not_extracted"
fi
cuda_version="unavailable"
if grep -q '^COSMOSIM_ENABLE_CUDA:BOOL=ON$' "$build_dir/CMakeCache.txt" 2>/dev/null; then
  cuda_version="enabled_version_not_extracted"
fi
source_revision="unavailable"
if command -v git >/dev/null 2>&1 && git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  source_revision="$(git rev-parse HEAD)"
fi

python3 tools/gravity_acceptance_manifest.py create \
  --output "$manifest_out" \
  --source-revision "$source_revision" \
  --build-profile "$build_dir" \
  --mpi-version "$mpi_version" \
  --fftw-version "$fftw_version" \
  --hdf5-version "$hdf5_version" \
  --cuda-version "$cuda_version" \
  --profile-id "${COSMOSIM_ACCEPTANCE_PROFILE_ID:-unknown_not_extracted}" \
  --pm-backend "$pm_backend" \
  --pm-backend-capability "$pm_backend_capability" \
  --pm-grid "${COSMOSIM_ACCEPTANCE_PMGRID:-unknown_not_extracted}" \
  --assignment "${COSMOSIM_ACCEPTANCE_ASSIGNMENT:-unknown_not_extracted}" \
  --deconvolution "${COSMOSIM_ACCEPTANCE_DECONVOLUTION:-unknown_not_extracted}" \
  --pm-decomposition "${COSMOSIM_ACCEPTANCE_PM_DECOMPOSITION:-unknown_not_extracted}" \
  --theta-mac "${COSMOSIM_ACCEPTANCE_THETA_MAC:-unknown_not_extracted}" \
  --asmth-cells "${COSMOSIM_ACCEPTANCE_ASMTH_CELLS:-unknown_not_extracted}" \
  --rcut-cells "${COSMOSIM_ACCEPTANCE_RCUT_CELLS:-unknown_not_extracted}" \
  --softening-profile "${COSMOSIM_ACCEPTANCE_SOFTENING_PROFILE:-unknown_not_extracted}" \
  "${rank_args[@]}" \
  "${evidence_args[@]}" \
  --note "Only attached evidence is certified by this manifest; unavailable feature lanes remain unexecuted."
python3 tools/gravity_acceptance_manifest.py verify "$manifest_out"
echo "Acceptance manifest written to $manifest_out and bound to the current gravity-contract source fingerprint."
