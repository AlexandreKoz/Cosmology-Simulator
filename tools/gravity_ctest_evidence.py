#!/usr/bin/env python3
"""Run one registered CTest lane and capture machine-readable gravity evidence."""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import pathlib
import re
import subprocess
import sys

METRIC = re.compile(r"([A-Za-z][A-Za-z0-9_]*)=([-+0-9.eE]+|not_measured)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--build-dir", type=pathlib.Path, required=True)
    ap.add_argument("--test", required=True)
    ap.add_argument("--output", type=pathlib.Path, required=True)
    ap.add_argument("--profile", default="unknown")
    ap.add_argument("--rank-count", type=int, default=1)
    ap.add_argument("--fft-backend", default="unknown")
    ap.add_argument("--softening-ratio", default="unknown")
    ap.add_argument("--env", action="append", default=[])
    args = ap.parse_args()

    env = os.environ.copy()
    for item in args.env:
        if "=" not in item:
            ap.error(f"--env must be KEY=VALUE, got {item!r}")
        key, value = item.split("=", 1)
        env[key] = value

    command = ["ctest", "--test-dir", str(args.build_dir), "-V", "--output-on-failure", "-R", f"^{re.escape(args.test)}$"]
    proc = subprocess.run(command, text=True, capture_output=True, env=env, check=False)
    combined = (proc.stdout or "") + (proc.stderr or "")
    metric_lines: list[dict[str, object]] = []
    for line in combined.splitlines():
        values: dict[str, object] = {}
        for match in METRIC.finditer(line):
            raw = match.group(2)
            if raw == "not_measured":
                values[match.group(1)] = raw
            else:
                try:
                    values[match.group(1)] = float(raw)
                except ValueError:
                    values[match.group(1)] = raw
        if values:
            metric_lines.append({"line": line.strip(), "metrics": values})

    payload = {
        "schema": "chui.gravity.ctest.evidence.v1",
        "generated_utc": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
        "test_identifier": args.test,
        "command": command,
        "profile": args.profile,
        "rank_count": args.rank_count,
        "fft_backend": args.fft_backend,
        "softening_ratio_epsilon_over_rs": args.softening_ratio,
        "result": "PASS" if proc.returncode == 0 else "FAIL",
        "return_code": proc.returncode,
        "metrics": metric_lines,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    sys.stdout.write(combined)
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
