#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <artifact_dir>" >&2
  exit 2
fi

artifact_dir="$1"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"
mkdir -p "$artifact_dir"

./scripts/ci/run_preset_pipeline.sh cpu-only-debug build-cpu-debug test-cpu-debug "unit_config_parser|integration_provenance_roundtrip|integration_feature_summary|validation_regression" "$artifact_dir" 1

./scripts/ci/check_build_metadata.py \
  "$artifact_dir/cosmosim_build_metadata-cpu-only-debug.json" \
  "validation/reference/ci_build_metadata_expectations_v1.json" \
  "cpu-only-debug"

sha256sum validation/reference/validation_tolerances_v1.txt > "$artifact_dir/validation_tolerances_v1.sha256"
sha256sum bench/baselines/benchmark_sizes_v1.txt > "$artifact_dir/benchmark_sizes_v1.sha256"

cat > "$artifact_dir/reproducibility_report.txt" <<REPORT
reproducibility_gate=pass
preset=cpu-only-debug
checked_tests=unit_config_parser,integration_feature_summary,validation_regression
metadata_contract=validation/reference/ci_build_metadata_expectations_v1.json
validation_tolerance_hash_file=validation_tolerances_v1.sha256
benchmark_size_hash_file=benchmark_sizes_v1.sha256
REPORT
