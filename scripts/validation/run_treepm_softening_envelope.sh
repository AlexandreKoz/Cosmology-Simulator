#!/usr/bin/env bash
set -euo pipefail

exe="${1:-./build/mpi-fftw-release/test_validation_tree_pm_ewald_accuracy}"
for ratio in 0 0.05 0.10 0.20; do
  echo "=== TreePM Ewald epsilon/r_s=${ratio} ==="
  COSMOSIM_EWALD_SOFTENING_RATIO="${ratio}" "${exe}"
done
