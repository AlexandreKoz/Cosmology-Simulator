#!/usr/bin/env bash
set -euo pipefail

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

cmake --preset "$configure_preset"
cmake --build --preset "$build_preset"

ctest_junit="$artifact_dir/ctest-${test_preset}.xml"
if [[ "$test_regex" == "-" ]]; then
  ctest --preset "$test_preset" --output-junit "$ctest_junit"
else
  ctest --preset "$test_preset" -R "$test_regex" --output-junit "$ctest_junit"
fi

build_dir="build/${configure_preset}"
cp "$build_dir/cosmosim_feature_summary.txt" "$artifact_dir/cosmosim_feature_summary-${configure_preset}.txt"
cp "$build_dir/cosmosim_build_metadata.json" "$artifact_dir/cosmosim_build_metadata-${configure_preset}.json"

if [[ "$run_bench" == "1" ]]; then
  bench_exe="$build_dir/bench_layout_smoke"
  if [[ -x "$bench_exe" ]]; then
    "$bench_exe" | tee "$artifact_dir/bench_layout_smoke-${configure_preset}.txt"
  else
    echo "bench_layout_smoke unavailable for preset ${configure_preset}" | tee "$artifact_dir/bench_layout_smoke-${configure_preset}.txt"
  fi
fi
