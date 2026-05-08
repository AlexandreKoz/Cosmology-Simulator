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

stage1_tests=(
  unit_simulation_state
  unit_gas_cell_identity_invariants
  unit_hot_cold_sidecar_layout
  integration_reorder_compaction_sidecars
  integration_species_migration_invariants
  integration_transform_fuzz_invariants
  integration_softening_ownership_invariants
  integration_hierarchical_time_bins
  integration_hierarchical_timestep_regression
  integration_parallel_two_rank_restart
  unit_snapshot_hdf5_schema
  unit_restart_checkpoint_schema
  integration_snapshot_hdf5_roundtrip
  integration_restart_checkpoint_roundtrip
  integration_provenance_roundtrip
  integration_runtime_truth_ctest_labels
)

test_regex="$(IFS='|'; echo "${stage1_tests[*]}")"

bash ./scripts/ci/run_preset_pipeline.sh \
  cpu-only-debug \
  build-cpu-debug \
  test-stage1-runtime-truth-cpu-debug \
  "$test_regex" \
  "$artifact_dir" \
  0

python3 ./scripts/ci/check_build_metadata.py \
  "$artifact_dir/cosmosim_build_metadata-cpu-only-debug.json" \
  validation/reference/ci_build_metadata_expectations_v1.json \
  cpu-only-debug

cat > "$artifact_dir/stage1_runtime_truth_gate_report.txt" <<REPORT
stage1_runtime_truth_gate=pass
preset=cpu-only-debug
test_preset=test-stage1-runtime-truth-cpu-debug
checked_tests=$(IFS=,; echo "${stage1_tests[*]}")
metadata_contract=validation/reference/ci_build_metadata_expectations_v1.json
optional_dependency_negative_check=integration_runtime_truth_ctest_labels
REPORT
