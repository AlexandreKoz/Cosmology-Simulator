#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=source_package_common.sh
source "$repo_root/scripts/ci/source_package_common.sh"

output_zip="${1:-$repo_root/../cosmosim_source.zip}"
archive_root="${2:-$(basename "$repo_root")}"

if [[ "$output_zip" != /* ]]; then
  output_parent="$(dirname "$output_zip")"
  mkdir -p "$output_parent"
  output_parent="$(cd "$output_parent" && pwd)"
  output_zip="$output_parent/$(basename "$output_zip")"
fi

staging_parent="$(mktemp -d)"
validation_parent="$(mktemp -d)"
package_success=0
cleanup() {
  rm -rf "$staging_parent" "$validation_parent"
  if (( package_success == 0 )); then
    rm -f "$output_zip"
  fi
}
trap cleanup EXIT

staging_root="$staging_parent/$archive_root"
extracted_root="$validation_parent/$archive_root"
staged_file_manifest="$validation_parent/staged_files.tsv"
extracted_file_manifest="$validation_parent/extracted_files.tsv"
staged_symlink_manifest="$validation_parent/staged_symlinks.tsv"
extracted_symlink_manifest="$validation_parent/extracted_symlinks.tsv"

source_package_require_commands

if [[ ! "$archive_root" =~ ^[A-Za-z0-9._-]+$ ]] ||
   [[ "$archive_root" == "." || "$archive_root" == ".." ]]; then
  echo "[package] ERROR: archive root must be one automation-safe path component" >&2
  exit 1
fi

mkdir -p "$(dirname "$output_zip")"
rm -f "$output_zip"

(
  cd "$repo_root"
  bash scripts/ci/check_repo_hygiene.sh
)
source_package_assert_no_ads "$repo_root" "source repository"
source_package_verify_no_forbidden_artifacts "$repo_root" "source repository"
source_hygiene_result="pass"

source_package_stage_tree "$repo_root" "$staging_root"
(
  cd "$staging_root"
  bash scripts/ci/check_repo_hygiene.sh
)
source_package_assert_no_ads "$staging_root" "staged tree"
source_package_verify_no_forbidden_artifacts "$staging_root" "staged tree"
source_package_verify_sentinels "$staging_root" "staged tree"
staged_sentinel_result="pass"

source_package_write_regular_hash_manifest "$staging_root" "$staged_file_manifest"
source_package_write_symlink_manifest "$staging_root" "$staged_symlink_manifest"

(
  cd "$staging_parent"
  zip -q -y -r "$output_zip" "$archive_root"
)
unzip -t "$output_zip" >/dev/null
unzip_test_result="pass"

archive_roots="$(unzip -Z1 "$output_zip" | awk -F/ 'NF {print $1}' | LC_ALL=C sort -u)"
if [[ "$archive_roots" != "$archive_root" ]]; then
  echo "[package] ERROR: archive must contain exactly one repository root '$archive_root'" >&2
  printf '[package] observed roots:\n%s\n' "$archive_roots" >&2
  exit 1
fi

if unzip -Z1 "$output_zip" | grep -E \
  '(^|/)(build|build-[^/]*|cmake-build-[^/]*|CMakeFiles|Testing|_deps|integration_outputs|validation_outputs|test_outputs|output|\.pytest_cache|__pycache__|\.mypy_cache|\.ruff_cache|\.ipynb_checkpoints|vcpkg_installed)(/|$)|(^|/)(CMakeCache\.txt|CTestTestfile\.cmake|CMakeUserPresets\.json)$|Zone\.Identifier|#Uf03aZone\.Identifier|\.(o|obj|a|lib|so|dll|dylib|exe|pdb|pyc|pyo|gcda|gcno|profraw|profdata|part|zip)$' \
  >/dev/null; then
  echo "[package] ERROR: archive name denylist matched generated or ADS content" >&2
  exit 1
fi

unzip -q "$output_zip" -d "$validation_parent"
if [[ ! -d "$extracted_root" ]]; then
  echo "[package] ERROR: extracted repository root is missing: $extracted_root" >&2
  exit 1
fi

source_package_verify_sentinels "$extracted_root" "extracted archive"
extracted_sentinel_result="pass"
source_package_assert_no_ads "$extracted_root" "extracted archive"
source_package_verify_no_forbidden_artifacts "$extracted_root" "extracted archive"
archive_denylist_result="pass"

source_package_write_regular_hash_manifest "$extracted_root" "$extracted_file_manifest"
source_package_write_symlink_manifest "$extracted_root" "$extracted_symlink_manifest"

if ! diff -u \
    <(cut -f1 "$staged_file_manifest") \
    <(cut -f1 "$extracted_file_manifest"); then
  echo "[package] ERROR: staged/extracted regular-file inventory differs" >&2
  exit 1
fi
if ! diff -u "$staged_symlink_manifest" "$extracted_symlink_manifest"; then
  echo "[package] ERROR: staged/extracted symlink inventory or targets differ" >&2
  exit 1
fi
inventory_result="pass"

if ! diff -u "$staged_file_manifest" "$extracted_file_manifest"; then
  echo "[package] ERROR: staged/extracted regular-file SHA-256 values differ" >&2
  exit 1
fi
hash_result="pass"

(
  cd "$extracted_root"
  bash scripts/ci/check_repo_hygiene.sh
)
extracted_hygiene_result="pass"

(
  cd "$extracted_root"
  cmake --preset cpu-only-debug
)
configure_result="pass"
(
  cd "$extracted_root"
  cmake --build --preset build-cpu-debug --parallel "${COSMOSIM_PACKAGE_BUILD_JOBS:-4}" --target cosmosim_core cosmosim_harness
)
build_result="pass"

entry_count="$(unzip -Z1 "$output_zip" | wc -l | tr -d '[:space:]')"
file_count="$(wc -l < "$extracted_file_manifest" | tr -d '[:space:]')"
archive_size="$(wc -c < "$output_zip" | tr -d '[:space:]')"
archive_sha256="$(sha256sum "$output_zip" | awk '{print $1}')"

package_success=1

echo "[package] source ZIP ready: $output_zip"
echo "[package] archive path: $output_zip"
echo "[package] archive root: $archive_root"
echo "[package] archive entry count: $entry_count"
echo "[package] archive file count: $file_count"
echo "[package] archive size bytes: $archive_size"
echo "[package] SHA-256: $archive_sha256"
echo "[package] source hygiene result: $source_hygiene_result"
echo "[package] staged sentinel result: $staged_sentinel_result"
echo "[package] extracted sentinel result: $extracted_sentinel_result"
echo "[package] staged/extracted inventory result: $inventory_result"
echo "[package] staged/extracted hash result: $hash_result"
echo "[package] extracted hygiene result: $extracted_hygiene_result"
echo "[package] extracted CMake configure result: $configure_result"
echo "[package] extracted build result: $build_result (cosmosim_core, cosmosim_harness)"
echo "[package] archive denylist result: $archive_denylist_result"
echo "[package] unzip test result: $unzip_test_result"
