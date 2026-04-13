#!/usr/bin/env python3
import json
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 4:
        print(
            "usage: check_build_metadata.py <metadata_json> <expectations_json> <preset>",
            file=sys.stderr,
        )
        return 2

    metadata_path = Path(sys.argv[1])
    expectations_path = Path(sys.argv[2])
    preset = sys.argv[3]

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    expectations = json.loads(expectations_path.read_text(encoding="utf-8"))

    if preset not in expectations:
      raise RuntimeError(f"missing preset expectations for {preset}")

    expected = expectations[preset]

    if metadata.get("project") != "cosmosim":
      raise RuntimeError("metadata project must be cosmosim")

    if metadata.get("preset") != preset:
      raise RuntimeError(f"metadata preset mismatch: expected {preset}, got {metadata.get('preset')}")

    features = metadata.get("features", {})
    for key, value in expected.get("features", {}).items():
      if features.get(key) != value:
        raise RuntimeError(f"feature {key} mismatch for {preset}: expected {value}, got {features.get(key)}")

    targets = set(metadata.get("targets", []))
    required_targets = set(expected.get("required_targets", []))
    missing = sorted(required_targets - targets)
    if missing:
      raise RuntimeError(f"missing required targets for {preset}: {', '.join(missing)}")

    dependencies = metadata.get("dependencies", {})
    for key, expected_state in expected.get("dependencies", {}).items():
      got = dependencies.get(key)
      if got != expected_state:
        raise RuntimeError(f"dependency {key} mismatch for {preset}: expected {expected_state}, got {got}")

    print(f"metadata validation passed for preset {preset}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
