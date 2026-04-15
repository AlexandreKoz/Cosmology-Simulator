#!/usr/bin/env bash
set -u -o pipefail

if [[ $# -lt 6 ]]; then
  echo "usage: $0 <configure_preset> <build_preset> <test_preset> <ctest_regex_or_dash> <artifact_dir> <run_bench_0_or_1>" >&2
  exit 2
fi

configure_preset="$1"
build_preset="$2"
test_preset="$3"
test_regex="$4"
artifact_dir="$5"
run_bench="$6"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

mkdir -p "$artifact_dir"

status="pass"
failed_phase=""
message="ok"

ctest_junit="$artifact_dir/ctest-${test_preset}.xml"
configure_cmd="cmake --preset ${configure_preset}"
build_cmd="cmake --build --preset ${build_preset}"
if [[ "$test_regex" == "-" ]]; then
  test_cmd="ctest --preset ${test_preset} --output-junit ${ctest_junit}"
else
  test_cmd="ctest --preset ${test_preset} -R '${test_regex}' --output-junit ${ctest_junit}"
fi
feature_copy_cmd="cp build/${configure_preset}/cosmosim_feature_summary.txt ${artifact_dir}/cosmosim_feature_summary-${configure_preset}.txt"
metadata_copy_cmd="cp build/${configure_preset}/cosmosim_build_metadata.json ${artifact_dir}/cosmosim_build_metadata-${configure_preset}.json"

run_cmd() {
  local phase="$1"
  local display_cmd="$2"
  shift 2
  echo "[preset-pipeline] phase=${phase} cmd=${display_cmd}"
  if ! "$@"; then
    status="fail"
    failed_phase="$phase"
    message="${phase} phase failed"
    return 1
  fi
  return 0
}

require_expected_tests_registered() {
  local test_listing_file="$artifact_dir/ctest-registered-${test_preset}.txt"
  local parsed_listing_file="$artifact_dir/ctest-registered-names-${test_preset}.txt"

  echo "[preset-pipeline] phase=test_registration cmd=ctest --preset ${test_preset} -N"
  if ! ctest --preset "$test_preset" -N > "$test_listing_file"; then
    status="fail"
    failed_phase="test_registration"
    message="failed to enumerate registered tests"
    return 1
  fi

  sed -n 's/^  Test[[:space:]]*#[0-9][0-9]*: //p' "$test_listing_file" > "$parsed_listing_file"

  if [[ "$test_regex" == "-" ]]; then
    return 0
  fi

  IFS='|' read -r -a expected_tests <<< "$test_regex"
  local expected_name
  for expected_name in "${expected_tests[@]}"; do
    [[ -z "$expected_name" ]] && continue
    if ! grep -Fxq "$expected_name" "$parsed_listing_file"; then
      status="fail"
      failed_phase="test_selection"
      message="expected test not registered: ${expected_name}"
      return 1
    fi
  done

  return 0
}

if run_cmd "configure" "$configure_cmd" cmake --preset "$configure_preset" \
  && run_cmd "build" "$build_cmd" cmake --build --preset "$build_preset" \
  && require_expected_tests_registered; then
  if [[ "$test_regex" == "-" ]]; then
    run_cmd "test" "$test_cmd" ctest --preset "$test_preset" --output-junit "$ctest_junit"
  else
    run_cmd "test" "$test_cmd" ctest --preset "$test_preset" -R "$test_regex" --output-junit "$ctest_junit"
  fi
fi

if [[ "$status" == "pass" ]]; then
  build_dir="build/${configure_preset}"
  if [[ ! -f "$build_dir/cosmosim_feature_summary.txt" || ! -f "$build_dir/cosmosim_build_metadata.json" ]]; then
    status="fail"
    failed_phase="artifact_collection"
    message="required build artifacts missing from ${build_dir}"
  else
    run_cmd "artifact_collection_feature_summary" "$feature_copy_cmd" cp "$build_dir/cosmosim_feature_summary.txt" "$artifact_dir/cosmosim_feature_summary-${configure_preset}.txt" || true
    if [[ "$status" == "pass" ]]; then
      run_cmd "artifact_collection_metadata" "$metadata_copy_cmd" cp "$build_dir/cosmosim_build_metadata.json" "$artifact_dir/cosmosim_build_metadata-${configure_preset}.json" || true
    fi
  fi
fi

if [[ "$status" == "pass" && "$run_bench" == "1" ]]; then
  bench_exe="build/${configure_preset}/bench_layout_smoke"
  bench_output="$artifact_dir/bench_layout_smoke-${configure_preset}.txt"
  if [[ -x "$bench_exe" ]]; then
    echo "[preset-pipeline] phase=benchmark cmd=${bench_exe} | tee ${bench_output}"
    if ! "$bench_exe" | tee "$bench_output"; then
      status="fail"
      failed_phase="benchmark"
      message="benchmark phase failed"
    fi
  else
    echo "bench_layout_smoke unavailable for preset ${configure_preset}" | tee "$bench_output"
  fi
fi

python3 - <<PY > "$artifact_dir/preset_pipeline_report-${configure_preset}.json"
import json
print(json.dumps({
  "configure_preset": "${configure_preset}",
  "build_preset": "${build_preset}",
  "test_preset": "${test_preset}",
  "test_scope": "${test_regex}",
  "artifact_dir": "${artifact_dir}",
  "status": "${status}",
  "failed_phase": "${failed_phase}",
  "message": "${message}",
  "commands": {
    "configure": "${configure_cmd}",
    "build": "${build_cmd}",
    "test": "${test_cmd}",
    "artifact_collection_feature_summary": "${feature_copy_cmd}",
    "artifact_collection_metadata": "${metadata_copy_cmd}"
  }
}, separators=(",", ":")))
PY

if [[ "$status" != "pass" ]]; then
  echo "[preset-pipeline] failed in phase: ${failed_phase}" >&2
  exit 1
fi

exit 0
