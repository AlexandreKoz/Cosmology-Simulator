#!/usr/bin/env bash
set -euo pipefail

build_dir="${1:-build/cpu-only-debug}"
poster_ranks="${COSMOSIM_POSTER_RANKS:-4}"
manifest_out="${COSMOSIM_GRAVITY_ACCEPTANCE_MANIFEST:-gravity_acceptance_manifest.json}"

if [[ ! -d "$build_dir" ]]; then
  echo "gravity acceptance: build directory not found: $build_dir" >&2
  exit 2
fi

run_ctest_regex() {
  local regex="$1"
  ctest --test-dir "$build_dir" --output-on-failure -R "$regex"
}

echo "== serial gravity contract gates =="
run_ctest_regex '^(unit_tree_gravity|unit_pm_solver|integration_tree_gravity_vs_direct|integration_tree_pm_coupling_periodic|integration_star_formation_source_runtime)$'

if ctest --test-dir "$build_dir" -N | grep -q 'validation_tree_pm_ewald_accuracy'; then
  echo "== Ewald executable present =="
  echo "Run scripts/validation/run_treepm_softening_envelope.sh $build_dir/test_validation_tree_pm_ewald_accuracy on a production FFTW build."
fi

if command -v mpiexec >/dev/null 2>&1 && ctest --test-dir "$build_dir" -N | grep -qi 'mpi'; then
  echo "== MPI gravity/PM/restart gates registered in this build =="
  # CTest owns the exact mpiexec command encoded by CMake. A registered lane
  # failure must fail this script; do not mask it with a best-effort fallback.
  run_ctest_regex '(tree_pm|treepm|pm.*slab|zeldovich|restart).*(two|three|four)[_-]rank|(two|three|four)[_-]rank.*(tree_pm|treepm|pm|zeldovich|restart)'
  case "$poster_ranks" in
    1|2|3|4)
      echo "Poster-rank request (${poster_ranks}) is covered only if a matching registered CTest lane appears above."
      ;;
    *)
      echo "Poster-rank count ${poster_ranks} has no generic synthetic CTest mapping; run the production DMO acceptance executable/config explicitly at that rank count."
      ;;
  esac
else
  echo "MPI runtime/tests unavailable in this build; no MPI acceptance is claimed."
fi

if command -v python3 >/dev/null 2>&1; then
  pm_backend="unknown"
  pm_backend_capability="unavailable"
  if grep -q '^COSMOSIM_ENABLE_FFTW:BOOL=ON$' "$build_dir/CMakeCache.txt" 2>/dev/null; then
    pm_backend="fftw"
    pm_backend_capability="production_fftw"
  elif grep -q '^COSMOSIM_ENABLE_FFTW:BOOL=OFF$' "$build_dir/CMakeCache.txt" 2>/dev/null; then
    pm_backend="naive_dft"
    pm_backend_capability="diagnostic_naive_dft"
  fi

  python3 tools/gravity_acceptance_manifest.py create \
    --output "$manifest_out" \
    --build-profile "$build_dir" \
    --profile-id "unverified-release-candidate" \
    --pm-backend "$pm_backend" \
    --pm-backend-capability "$pm_backend_capability" \
    --pm-grid "from-run-config" \
    --assignment "from-run-config" \
    --deconvolution "from-run-config" \
    --pm-decomposition "from-run-config" \
    --theta-mac "from-run-config" \
    --asmth-cells "from-run-config" \
    --rcut-cells "from-run-config" \
    --softening-profile "from-run-config" \
    --rank-count "1" \
    --note "MPI/FFTW/HDF5/CUDA evidence must be attached separately for the exact source fingerprint"
  python3 tools/gravity_acceptance_manifest.py verify "$manifest_out"
  echo "Acceptance manifest written to $manifest_out. It records source identity; it is not scientific certification by itself."
fi
