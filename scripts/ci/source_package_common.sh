#!/usr/bin/env bash

# Shared source-package invariants. Callers must enable their own strict mode.

source_package_sentinels=(
  "AGENTS.md"
  "CMakeLists.txt"
  "CMakePresets.json"
  "src/core/version.cpp"
  "src/core/harness_main.cpp"
  "include/cosmosim/core/config.hpp"
  "src/io/ic_distributed_ingestion.cpp"
  "src/io/ic_distributed_audit.cpp"
  "src/io/ic_mpi_collectives.cpp"
  "src/io/internal/ic_mpi_collectives.hpp"
  "tests/integration/test_distributed_ic_reader_mpi.cpp"
  "scripts/ci/check_repo_hygiene.sh"
  "scripts/ci/source_package_common.sh"
  "scripts/ci/package_source_zip.sh"
  "scripts/ci/test_source_package_completeness.sh"
)

source_package_rsync_excludes=(
  ".git/"
  "build/"
  "build-*/"
  "cmake-build-*/"
  "CMakeFiles/"
  "Testing/"
  "_deps/"
  "integration_outputs/"
  "validation_outputs/"
  "test_outputs/"
  "output/"
  ".pytest_cache/"
  "__pycache__/"
  ".mypy_cache/"
  ".ruff_cache/"
  ".ipynb_checkpoints/"
  "vcpkg_installed/"
  "CMakeCache.txt"
  "CTestTestfile.cmake"
  "CMakeUserPresets.json"
  "*.zip"
  "*.pyc"
  "*.pyo"
  "*.part"
  "*.o"
  "*.obj"
  "*.a"
  "*.lib"
  "*.so"
  "*.so.*"
  "*.dll"
  "*.dylib"
  "*.exe"
  "*.pdb"
  "*.gcda"
  "*.gcno"
  "*.profraw"
  "*.profdata"
)

source_package_require_commands() {
  local command_name
  for command_name in bash cmake find rsync sha256sum sort unzip zip; do
    if ! command -v "$command_name" >/dev/null 2>&1; then
      echo "[package] ERROR: required command '$command_name' is unavailable" >&2
      return 1
    fi
  done
}

source_package_verify_sentinels() {
  local tree_root="$1"
  local tree_label="$2"
  local sentinel
  local missing=0

  for sentinel in "${source_package_sentinels[@]}"; do
    if [[ ! -f "$tree_root/$sentinel" ]]; then
      echo "[package] ERROR: sentinel missing from $tree_label: $sentinel" >&2
      missing=1
    fi
  done

  if (( missing != 0 )); then
    return 1
  fi
}

source_package_assert_no_ads() {
  local tree_root="$1"
  local tree_label="$2"
  local matches

  matches="$(find "$tree_root" -path "$tree_root/.git" -prune -o \
    \( -name '*Zone.Identifier*' -o -name '*#Uf03aZone.Identifier*' \) \
    -print)"
  if [[ -n "$matches" ]]; then
    echo "[package] ERROR: ADS content found in $tree_label:" >&2
    printf '%s\n' "$matches" >&2
    return 1
  fi
}

source_package_stage_tree() {
  local source_root="$1"
  local staging_root="$2"
  local rsync_args=(-a --delete)
  local exclude_pattern

  for exclude_pattern in "${source_package_rsync_excludes[@]}"; do
    rsync_args+=("--exclude=$exclude_pattern")
  done

  mkdir -p "$staging_root"
  rsync "${rsync_args[@]}" "$source_root/" "$staging_root/"
}

source_package_write_regular_hash_manifest() {
  local tree_root="$1"
  local output_manifest="$2"
  local relative_path
  local digest

  : > "$output_manifest"
  while IFS= read -r -d '' relative_path; do
    relative_path="${relative_path#./}"
    digest="$(sha256sum "$tree_root/$relative_path" | awk '{print $1}')"
    printf '%s\t%s\n' "$relative_path" "$digest" >> "$output_manifest"
  done < <(cd "$tree_root" && find . -type f -print0 | LC_ALL=C sort -z)
}

source_package_write_symlink_manifest() {
  local tree_root="$1"
  local output_manifest="$2"
  local relative_path
  local target

  : > "$output_manifest"
  while IFS= read -r -d '' relative_path; do
    relative_path="${relative_path#./}"
    target="$(readlink "$tree_root/$relative_path")"
    printf '%s\t%s\n' "$relative_path" "$target" >> "$output_manifest"
  done < <(cd "$tree_root" && find . -type l -print0 | LC_ALL=C sort -z)
}

source_package_verify_no_forbidden_artifacts() {
  local tree_root="$1"
  local tree_label="$2"
  local matches

  matches="$(find "$tree_root" -path "$tree_root/.git" -prune -o \
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
       -name '*.pyc' -o -name '*.pyo' -o -name '*.gcda' -o -name '*.gcno' -o \
       -name '*.profraw' -o -name '*.profdata' -o -name '*.part' -o \
       -name '*.zip' -o -name core -o -name 'core.*' -o \
       -name '*Zone.Identifier*' -o -name '*#Uf03aZone.Identifier*' \
     \) -print \))"

  if [[ -n "$matches" ]]; then
    echo "[package] ERROR: generated, crash-dump, archive, or ADS content found in $tree_label:" >&2
    printf '%s\n' "$matches" >&2
    return 1
  fi
}
