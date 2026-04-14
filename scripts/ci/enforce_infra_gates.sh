#!/usr/bin/env bash
set -u -o pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <artifact_dir>" >&2
  exit 2
fi

artifact_dir="$1"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

mkdir -p "$artifact_dir"

report_json="$artifact_dir/infrastructure_gate_report.json"
report_md="$artifact_dir/infrastructure_gate_report.md"

declare -a gate_rows=()
overall_status="pass"

run_gate() {
  local gate_id="$1"
  local configure_preset="$2"
  local build_preset="$3"
  local test_preset="$4"
  local test_regex="$5"
  local run_bench="$6"
  local gate_artifact_dir="$artifact_dir/$gate_id"
  local status="pass"
  local message="ok"

  mkdir -p "$gate_artifact_dir"
  echo "[infra-gate] running ${gate_id} (${configure_preset})"

  if ! bash ./scripts/ci/run_preset_pipeline.sh \
      "$configure_preset" \
      "$build_preset" \
      "$test_preset" \
      "$test_regex" \
      "$gate_artifact_dir" \
      "$run_bench"; then
    status="fail"
    message="preset pipeline failed"
    overall_status="fail"
  fi

  if [[ "$status" == "pass" ]]; then
    if ! python3 ./scripts/ci/check_build_metadata.py \
      "$gate_artifact_dir/cosmosim_build_metadata-${configure_preset}.json" \
      "validation/reference/ci_build_metadata_expectations_v1.json" \
      "$configure_preset"; then
      status="fail"
      message="metadata contract validation failed"
      overall_status="fail"
    fi
  fi

  gate_rows+=("{\"gate_id\":\"${gate_id}\",\"configure_preset\":\"${configure_preset}\",\"status\":\"${status}\",\"message\":\"${message}\",\"test_regex\":\"${test_regex}\"}")
}

run_gate \
  "cpu_core_boundary_and_config_contract" \
  "cpu-only-debug" \
  "build-cpu-debug" \
  "test-cpu-debug" \
  "integration_core_dependency_direction|unit_config_parser|integration_release_readiness_artifacts" \
  "0"

run_gate \
  "hdf5_schema_and_exact_restart_contract" \
  "hdf5-debug" \
  "build-hdf5-debug" \
  "test-hdf5-debug" \
  "unit_snapshot_hdf5_schema|unit_restart_checkpoint_schema|integration_snapshot_hdf5_roundtrip|integration_restart_checkpoint_roundtrip" \
  "0"

run_gate \
  "pm_hdf5_fftw_feature_path_validation" \
  "pm-hdf5-fftw-debug" \
  "build-pm-hdf5-fftw-debug" \
  "test-pm-hdf5-fftw-debug" \
  "unit_pm_solver|integration_pm_periodic_mode|integration_tree_pm_coupling_periodic|validation_regression" \
  "1"

{
  echo "{"
  echo "  \"overall_status\": \"${overall_status}\","
  echo "  \"generated_by\": \"scripts/ci/enforce_infra_gates.sh\","
  echo "  \"gates\": ["
  for i in "${!gate_rows[@]}"; do
    if [[ "$i" -gt 0 ]]; then
      echo ","
    fi
    echo -n "    ${gate_rows[$i]}"
  done
  echo
  echo "  ]"
  echo "}"
} > "$report_json"

{
  echo "# Infrastructure gate report"
  echo
  echo "- overall_status: ${overall_status}"
  for row in "${gate_rows[@]}"; do
    echo "- ${row}"
  done
} > "$report_md"

if [[ "$overall_status" != "pass" ]]; then
  echo "[infra-gate] one or more infrastructure gates failed" >&2
  exit 1
fi

echo "[infra-gate] all infrastructure gates passed"
