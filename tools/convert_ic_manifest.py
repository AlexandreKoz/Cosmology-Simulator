#!/usr/bin/env python3
"""Create a small JSON manifest documenting expected IC aliases.

This helper is intentionally lightweight so workflows can validate schema mapping
without requiring the full C++ binary toolchain.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ALIAS_MAP = {
    "coordinates": ["Coordinates", "Position", "POS"],
    "velocities": ["Velocities", "Velocity", "VEL"],
    "particle_ids": ["ParticleIDs", "ParticleID", "ID"],
    "masses": ["Masses", "Mass"],
}


def build_manifest(ic_path: str) -> dict:
    return {
        "ic_path": ic_path,
        "supported_container": "hdf5",
        "header_requirements": ["NumPart_ThisFile", "Time"],
        "header_optional": ["MassTable"],
        "alias_map": ALIAS_MAP,
        "assumed_units": {"length": "kpc", "mass": "msun", "velocity": "km_s"},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Emit a converter alias manifest for CosmoSim IC ingestion")
    parser.add_argument("ic_path", help="Path to IC file (for metadata only)")
    parser.add_argument("--output", "-o", default="ic_manifest.json", help="Output JSON path")
    args = parser.parse_args()

    manifest = build_manifest(args.ic_path)
    output_path = Path(args.output)
    output_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote manifest: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
