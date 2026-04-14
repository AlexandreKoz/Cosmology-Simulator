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
manifest_path="scripts/ci/infrastructure_gates_manifest.tsv"

if [[ ! -f "$manifest_path" ]]; then
  echo "[infra-gate] missing manifest: ${manifest_path}" >&2
  exit 2
fi

declare -a gate_rows=()
overall_status="pass"

json_escape() {
  python3 -c 'import json,sys; print(json.dumps(sys.stdin.read()))'
}

run_gate() {
  local gate_id="$1"
  local configure_preset="$2"
  local build_preset="$3"
  local test_preset="$4"
  local test_scope_names_csv="$5"
  local metadata_expectation="$6"
  local artifact_subdir="$7"
  local run_bench="$8"

  local test_scope_regex
  test_scope_regex="${test_scope_names_csv//,/|}"

  local gate_artifact_dir="$artifact_dir/$artifact_subdir"
  local status="pass"
  local failed_phase=""
  local message="ok"

  mkdir -p "$gate_artifact_dir"
  echo "[infra-gate] running ${gate_id} (${configure_preset})"

  local configure_cmd="cmake --preset ${configure_preset}"
  local build_cmd="cmake --build --preset ${build_preset}"
  local test_cmd
  if [[ "$test_scope_regex" == "-" ]]; then
    test_cmd="ctest --preset ${test_preset} --output-junit ${gate_artifact_dir}/ctest-${test_preset}.xml"
  else
    test_cmd="ctest --preset ${test_preset} -R ${test_scope_regex} --output-junit ${gate_artifact_dir}/ctest-${test_preset}.xml"
  fi
  local metadata_cmd="python3 ./scripts/ci/check_build_metadata.py ${gate_artifact_dir}/cosmosim_build_metadata-${configure_preset}.json ${metadata_expectation} ${configure_preset}"

  local pipeline_status_path="$gate_artifact_dir/preset_pipeline_report-${configure_preset}.json"
  if ! bash ./scripts/ci/run_preset_pipeline.sh \
      "$configure_preset" \
      "$build_preset" \
      "$test_preset" \
      "$test_scope_regex" \
      "$gate_artifact_dir" \
      "$run_bench"; then
    status="fail"
    if [[ -f "$pipeline_status_path" ]]; then
      failed_phase="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("failed_phase","unknown"))' "$pipeline_status_path")"
      message="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("message","preset pipeline failed"))' "$pipeline_status_path")"
    else
      failed_phase="preset_pipeline"
      message="preset pipeline failed before report generation"
    fi
    overall_status="fail"
  fi

  if [[ "$status" == "pass" ]]; then
    if ! python3 ./scripts/ci/check_build_metadata.py \
      "$gate_artifact_dir/cosmosim_build_metadata-${configure_preset}.json" \
      "$metadata_expectation" \
      "$configure_preset"; then
      status="fail"
      failed_phase="metadata_validation"
      message="metadata contract validation failed"
      overall_status="fail"
    fi
  fi

  local message_json
  message_json="$(printf '%s' "$message" | json_escape)"

  gate_rows+=("{\"gate_id\":\"${gate_id}\",\"configure_preset\":\"${configure_preset}\",\"build_preset\":\"${build_preset}\",\"test_preset\":\"${test_preset}\",\"test_scope\":\"${test_scope_names_csv}\",\"test_scope_regex\":\"${test_scope_regex}\",\"metadata_expectation\":\"${metadata_expectation}\",\"artifact_dir\":\"${gate_artifact_dir}\",\"status\":\"${status}\",\"failed_phase\":\"${failed_phase}\",\"message\":${message_json},\"commands\":{\"configure\":\"${configure_cmd}\",\"build\":\"${build_cmd}\",\"test\":\"${test_cmd}\",\"metadata_validation\":\"${metadata_cmd}\"}}")
}

while IFS=$'\t' read -r gate_id configure_preset build_preset test_preset test_scope_names_csv metadata_expectation artifact_subdir run_bench; do
  [[ -z "$gate_id" || "$gate_id" == \#* ]] && continue
  run_gate "$gate_id" "$configure_preset" "$build_preset" "$test_preset" "$test_scope_names_csv" "$metadata_expectation" "$artifact_subdir" "$run_bench"
done < "$manifest_path"

{
  echo "{"
  echo "  \"overall_status\": \"${overall_status}\"," 
  echo "  \"generated_by\": \"scripts/ci/enforce_infra_gates.sh\"," 
  echo "  \"manifest\": \"${manifest_path}\"," 
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
  echo "- manifest: ${manifest_path}"
  for row in "${gate_rows[@]}"; do
    echo "- ${row}"
  done
} > "$report_md"

if [[ "$overall_status" != "pass" ]]; then
  echo "[infra-gate] one or more infrastructure gates failed" >&2
  exit 1
fi

echo "[infra-gate] all infrastructure gates passed"
