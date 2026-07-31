#!/usr/bin/env bash
set -euo pipefail

source_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=source_package_common.sh
source "$source_root/scripts/ci/source_package_common.sh"

source_package_require_commands
source_package_verify_sentinels "$source_root" "regression-test source tree"

work_root="$(mktemp -d)"
trap 'rm -rf "$work_root"' EXIT
candidate_root="$work_root/source_candidate"
archive_path="$work_root/source_package.zip"
package_log="$work_root/package.log"
extract_parent="$work_root/independent_extract"
extracted_root="$extract_parent/Cosmology-Simulator-main"
candidate_manifest="$work_root/candidate_files.tsv"
extracted_manifest="$work_root/extracted_files.tsv"
candidate_symlinks="$work_root/candidate_symlinks.tsv"
extracted_symlinks="$work_root/extracted_symlinks.tsv"

source_package_stage_tree "$source_root" "$candidate_root"
(
  cd "$candidate_root"
  bash scripts/ci/check_repo_hygiene.sh
)
source_package_verify_sentinels "$candidate_root" "regression-test staged candidate"
source_package_assert_no_ads "$candidate_root" "regression-test staged candidate"
source_package_verify_no_forbidden_artifacts "$candidate_root" "regression-test staged candidate"

# A regular crash-dump-style file must fail closed, while source directories
# named core remain legitimate and are checked after successful extraction.
crash_probe="$candidate_root/core.1234"
rejected_archive="$work_root/rejected_core_dump.zip"
touch "$crash_probe"
if bash "$candidate_root/scripts/ci/package_source_zip.sh" \
    "$rejected_archive" Cosmology-Simulator-main >/dev/null 2>&1; then
  echo "[package-test] ERROR: packaging accepted a regular core.<suffix> crash-dump candidate" >&2
  exit 1
fi
if [[ -e "$rejected_archive" ]]; then
  echo "[package-test] ERROR: failed packaging retained a user-facing archive" >&2
  exit 1
fi
rm -f "$crash_probe"

bash "$candidate_root/scripts/ci/package_source_zip.sh" \
  "$archive_path" Cosmology-Simulator-main | tee "$package_log"
grep -Fq "[package] extracted CMake configure result: pass" "$package_log"
grep -Fq "[package] extracted build result: pass (cosmosim_core, cosmosim_harness)" "$package_log"

unzip -t "$archive_path" >/dev/null
mkdir -p "$extract_parent"
unzip -q "$archive_path" -d "$extract_parent"

archive_roots="$(unzip -Z1 "$archive_path" | awk -F/ 'NF {print $1}' | LC_ALL=C sort -u)"
if [[ "$archive_roots" != "Cosmology-Simulator-main" ]]; then
  echo "[package-test] ERROR: expected one archive root, observed:" >&2
  printf '%s\n' "$archive_roots" >&2
  exit 1
fi

source_package_verify_sentinels "$extracted_root" "regression-test extracted archive"
[[ -d "$extracted_root/src/core" ]]
[[ -d "$extracted_root/include/cosmosim/core" ]]
source_package_assert_no_ads "$extracted_root" "regression-test extracted archive"
source_package_verify_no_forbidden_artifacts "$extracted_root" "regression-test extracted archive"

source_package_write_regular_hash_manifest "$candidate_root" "$candidate_manifest"
source_package_write_regular_hash_manifest "$extracted_root" "$extracted_manifest"
source_package_write_symlink_manifest "$candidate_root" "$candidate_symlinks"
source_package_write_symlink_manifest "$extracted_root" "$extracted_symlinks"

diff -u <(cut -f1 "$candidate_manifest") <(cut -f1 "$extracted_manifest")
diff -u "$candidate_manifest" "$extracted_manifest"
diff -u "$candidate_symlinks" "$extracted_symlinks"

(
  cd "$extracted_root"
  bash scripts/ci/check_repo_hygiene.sh
)

echo "[package-test] PASS: source package is clean, source-complete, inventory-identical, hash-identical, and the real package path builds required core targets"
