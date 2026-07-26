#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
output_zip="${1:-$repo_root/../cosmosim_source.zip}"
archive_root="${2:-$(basename "$repo_root")}" 
staging_parent="$(mktemp -d)"
trap 'rm -rf "$staging_parent"' EXIT
staging_root="$staging_parent/$archive_root"
mkdir -p "$staging_root"

rsync -a --delete \
  --exclude='.git/' \
  --exclude='build/' --exclude='build-*/' --exclude='cmake-build-*/' \
  --exclude='integration_outputs/' --exclude='validation_outputs/' \
  --exclude='test_outputs/' --exclude='output/' \
  --exclude='.pytest_cache/' --exclude='__pycache__/' \
  --exclude='.mypy_cache/' --exclude='.ruff_cache/' \
  --exclude='*.pyc' --exclude='*.pyo' --exclude='*.part' \
  --exclude='*.o' --exclude='*.obj' --exclude='*.a' --exclude='*.lib' \
  --exclude='*.so' --exclude='*.so.*' --exclude='*.dll' --exclude='*.dylib' \
  --exclude='*.exe' --exclude='*.pdb' --exclude='core' --exclude='core.*' \
  --exclude='*Zone.Identifier*' --exclude='*#Uf03aZone.Identifier*' \
  "$repo_root/" "$staging_root/"

(
  cd "$staging_root"
  bash scripts/ci/check_repo_hygiene.sh
)

rm -f "$output_zip"
mkdir -p "$(dirname "$output_zip")"
(
  cd "$staging_parent"
  zip -q -r "$output_zip" "$archive_root"
)
unzip -t "$output_zip" >/dev/null

if unzip -Z1 "$output_zip" | grep -E \
  '(^|/)(build|build-[^/]*|cmake-build-[^/]*|CMakeFiles|Testing|_deps|integration_outputs|validation_outputs|test_outputs|output|\.pytest_cache|__pycache__|\.mypy_cache|\.ruff_cache|\.ipynb_checkpoints|vcpkg_installed)(/|$)|(^|/)(CMakeCache\.txt|CTestTestfile\.cmake|CMakeUserPresets\.json|core|core\.[^/]+)$|Zone\.Identifier|#Uf03aZone\.Identifier|\.(o|obj|a|lib|so|dll|dylib|exe|pdb|pyc|pyo|gcda|gcno|profraw|profdata|part)$' \
  >/dev/null; then
  echo "[package] ERROR: archive denylist matched generated or ADS content" >&2
  exit 1
fi

echo "[package] source ZIP ready: $output_zip"
sha256sum "$output_zip"
