#!/usr/bin/env python3
"""Create/verify a source-bound gravity acceptance manifest.

The tool never claims certification by itself. It binds externally produced test
results to the exact gravity-relevant source tree and records the build/runtime
profile used to obtain them.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import pathlib
import platform
import subprocess
import sys
from typing import Any

ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCE_GLOBS = (
    # Gravity numerical implementation and public contracts.
    "src/gravity/**/*",
    "include/cosmosim/gravity/**/*",
    # Distributed ownership/transport used by PM and sparse short-range work.
    "src/parallel/**/*",
    "include/cosmosim/parallel/**/*",
    # Authoritative gravity source assembly and source-state mutation owners.
    "src/workflows/gravity_runtime.cpp",
    "src/workflows/hydro_amr_runtime.cpp",
    "src/workflows/source_runtime.cpp",
    "src/workflows/time_coordinator.cpp",
    "src/workflows/reference_workflow.cpp",
    "src/workflows/output_restart_runtime.cpp",
    "src/workflows/runtime_capabilities.cpp",
    "src/workflows/internal/gas_cell_ownership.cpp",
    "include/cosmosim/workflows/gravity_source_ownership.hpp",
    # Force/KDK/epoch/configuration semantics that can change scientific meaning.
    "src/core/config.cpp",
    "src/core/time_integration.cpp",
    "src/core/simulation_state_ownership.cpp",
    "include/cosmosim/core/config.hpp",
    "include/cosmosim/core/time_integration.hpp",
    "include/cosmosim/core/simulation_state.hpp",
    # Build-feature definitions can change available numerical backends.
    "CMakeLists.txt",
    "cmake/templates/build_config.hpp.in",
)


def files_for_fingerprint() -> list[pathlib.Path]:
    files: set[pathlib.Path] = set()
    for pattern in SOURCE_GLOBS:
        for path in ROOT.glob(pattern):
            if path.is_file():
                files.add(path)
    return sorted(files, key=lambda p: p.relative_to(ROOT).as_posix())


def source_fingerprint() -> tuple[str, list[str]]:
    digest = hashlib.sha256()
    rels: list[str] = []
    for path in files_for_fingerprint():
        rel = path.relative_to(ROOT).as_posix()
        rels.append(rel)
        digest.update(rel.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest(), rels


def command_output(argv: list[str]) -> str:
    try:
        proc = subprocess.run(argv, cwd=ROOT, check=False, capture_output=True, text=True, timeout=10)
        text = (proc.stdout or proc.stderr).strip().splitlines()
        return text[0] if text else "unavailable"
    except (OSError, subprocess.SubprocessError):
        return "unavailable"


def load_evidence(paths: list[pathlib.Path]) -> list[Any]:
    out: list[Any] = []
    for path in paths:
        value = json.loads(path.read_text(encoding="utf-8"))
        out.append({"source": str(path), "payload": value})
    return out


def build_manifest(args: argparse.Namespace) -> dict[str, Any]:
    fingerprint, source_files = source_fingerprint()
    return {
        "schema": "chui.gravity.acceptance.v1",
        "status": "evidence_only_not_self_certifying",
        "generated_utc": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
        "source": {
            "gravity_contract_sha256": fingerprint,
            "source_revision": args.source_revision,
            "source_file_count": len(source_files),
            "invalidation_scope": [
                "src/gravity/**",
                "include/cosmosim/gravity/**",
                "gravity workflow source ownership/integration semantics",
                "src/parallel/** and include/cosmosim/parallel/**",
                "PM/MPI decomposition and gravity configuration contracts",
                "core force/KDK/source-generation semantics",
                "gravity-relevant build feature definitions",
            ],
        },
        "build": {
            "profile": args.build_profile,
            "compiler": args.compiler or command_output([os.environ.get("CXX", "c++"), "--version"]),
            "cmake": command_output(["cmake", "--version"]),
            "mpi": args.mpi_version,
            "fftw": args.fftw_version,
            "hdf5": args.hdf5_version,
            "cuda": args.cuda_version,
            "platform": platform.platform(),
        },
        "gravity_profile": {
            "profile_id": args.profile_id,
            "pm_backend": args.pm_backend,
            "pm_backend_capability": args.pm_backend_capability,
            "pm_grid": args.pm_grid,
            "assignment": args.assignment,
            "deconvolution": args.deconvolution,
            "pm_decomposition": args.pm_decomposition,
            "theta_mac": args.theta_mac,
            "asmth_cells": args.asmth_cells,
            "rcut_cells": args.rcut_cells,
            "softening_profile": args.softening_profile,
            "rank_counts_tested": args.rank_count,
        },
        "evidence": load_evidence(args.evidence),
        "notes": args.note,
    }


def cmd_create(args: argparse.Namespace) -> int:
    manifest = build_manifest(args)
    text = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(args.output)
    return 0


def cmd_verify(args: argparse.Namespace) -> int:
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    current, _ = source_fingerprint()
    recorded = manifest.get("source", {}).get("gravity_contract_sha256", "")
    if current != recorded:
        print(
            f"STALE gravity acceptance manifest: recorded={recorded} current={current}",
            file=sys.stderr,
        )
        return 2
    print(f"CURRENT gravity acceptance source fingerprint: {current}")
    return 0


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    create = sub.add_parser("create", help="write a source-bound evidence manifest")
    create.add_argument("--output", type=pathlib.Path, required=True)
    create.add_argument("--source-revision", default="unavailable")
    create.add_argument("--build-profile", required=True)
    create.add_argument("--compiler", default="")
    create.add_argument("--mpi-version", default="unavailable")
    create.add_argument("--fftw-version", default="unavailable")
    create.add_argument("--hdf5-version", default="unavailable")
    create.add_argument("--cuda-version", default="unavailable")
    create.add_argument("--profile-id", required=True)
    create.add_argument("--pm-backend", required=True)
    create.add_argument("--pm-backend-capability", required=True)
    create.add_argument("--pm-grid", required=True)
    create.add_argument("--assignment", required=True)
    create.add_argument("--deconvolution", required=True)
    create.add_argument("--pm-decomposition", required=True)
    create.add_argument("--theta-mac", required=True)
    create.add_argument("--asmth-cells", required=True)
    create.add_argument("--rcut-cells", required=True)
    create.add_argument("--softening-profile", required=True)
    create.add_argument("--rank-count", action="append", default=[])
    create.add_argument("--evidence", type=pathlib.Path, action="append", default=[])
    create.add_argument("--note", action="append", default=[])
    create.set_defaults(func=cmd_create)

    verify = sub.add_parser("verify", help="fail if a manifest no longer matches gravity-relevant source")
    verify.add_argument("manifest", type=pathlib.Path)
    verify.set_defaults(func=cmd_verify)
    return p


def main() -> int:
    args = parser().parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
