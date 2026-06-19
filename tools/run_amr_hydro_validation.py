#!/usr/bin/env python3
"""Analyze reproducible CHUI AMR hydro profile campaigns and external references.

This tool never runs another code implicitly.  It consumes explicit profile CSVs and
versioned JSON manifests, validates setup compatibility, resamples to a common
coordinate grid, and writes one structured JSON result.  It is usable for CHUI
resolution ladders and for actual RAMSES/ENZO/AREPO/Athena++ outputs once a user
provides a converted canonical profile and provenance manifest.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

SCHEMA = "cosmosim_amr_profile_v1"

REQUIRED_METADATA = {
    "schema_version", "case_name", "geometry", "source_code_name", "source_code_version",
    "gamma", "final_time_code", "units", "boundary_conditions", "gravity_enabled",
    "cooling_enabled", "cosmological_expansion_enabled", "provenance",
}

@dataclass(frozen=True)
class Point:
    coordinate: float
    rho: float
    velocity: float
    pressure: float


def fail(message: str) -> None:
    raise ValueError(message)


def load_manifest(path: Path) -> dict:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    missing = REQUIRED_METADATA - manifest.keys()
    if missing:
        fail(f"{path}: manifest missing required keys: {sorted(missing)}")
    if manifest["schema_version"] != SCHEMA:
        fail(f"{path}: unsupported schema_version {manifest['schema_version']!r}")
    if manifest["geometry"] not in {"line", "radial"}:
        fail(f"{path}: geometry must be 'line' or 'radial'")
    if not math.isfinite(float(manifest["gamma"])) or float(manifest["gamma"]) <= 1.0:
        fail(f"{path}: gamma must be finite and > 1")
    if not math.isfinite(float(manifest["final_time_code"])) or float(manifest["final_time_code"]) < 0.0:
        fail(f"{path}: final_time_code must be finite and non-negative")
    return manifest


def profile_fields(geometry: str) -> tuple[str, str, str, str]:
    return (("x_code", "rho_code", "vel_x_code", "pressure_code") if geometry == "line"
            else ("radius_code", "rho_code", "radial_velocity_code", "pressure_code"))


def load_profile(path: Path, geometry: str) -> list[Point]:
    coordinate_key, rho_key, velocity_key, pressure_key = profile_fields(geometry)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            fail(f"{path}: CSV header is required")
        missing = {coordinate_key, rho_key, velocity_key, pressure_key} - set(reader.fieldnames)
        if missing:
            fail(f"{path}: CSV missing fields {sorted(missing)}")
        points = [Point(float(row[coordinate_key]), float(row[rho_key]), float(row[velocity_key]), float(row[pressure_key]))
                  for row in reader]
    if len(points) < 2:
        fail(f"{path}: profile must contain at least two rows")
    previous = -math.inf
    for point in points:
        if not all(math.isfinite(value) for value in (point.coordinate, point.rho, point.velocity, point.pressure)):
            fail(f"{path}: profile contains non-finite values")
        if point.rho <= 0.0 or point.pressure <= 0.0:
            fail(f"{path}: profile contains non-positive density/pressure")
        if point.coordinate <= previous:
            fail(f"{path}: coordinate must be strictly increasing")
        previous = point.coordinate
    return points


def require_compatible(lhs: dict, rhs: dict) -> None:
    for key in ("case_name", "geometry", "units", "boundary_conditions", "gravity_enabled", "cooling_enabled",
                "cosmological_expansion_enabled"):
        if lhs[key] != rhs[key]:
            fail(f"incompatible comparison metadata: {key}: {lhs[key]!r} != {rhs[key]!r}")
    if abs(float(lhs["gamma"]) - float(rhs["gamma"])) > 1.0e-14:
        fail("incompatible comparison metadata: gamma")
    if abs(float(lhs["final_time_code"]) - float(rhs["final_time_code"])) > 1.0e-12:
        fail("incompatible comparison metadata: final_time_code")


def interpolate(points: list[Point], coordinate: float) -> Point:
    if coordinate < points[0].coordinate or coordinate > points[-1].coordinate:
        fail("common coordinate is outside profile support")
    if coordinate == points[0].coordinate:
        return points[0]
    if coordinate == points[-1].coordinate:
        return points[-1]
    for upper_index in range(1, len(points)):
        upper = points[upper_index]
        if coordinate <= upper.coordinate:
            lower = points[upper_index - 1]
            alpha = (coordinate - lower.coordinate) / (upper.coordinate - lower.coordinate)
            return Point(coordinate,
                         (1.0-alpha)*lower.rho + alpha*upper.rho,
                         (1.0-alpha)*lower.velocity + alpha*upper.velocity,
                         (1.0-alpha)*lower.pressure + alpha*upper.pressure)
    raise AssertionError("unreachable")


def common_grid(lhs: list[Point], rhs: list[Point], samples: int) -> list[float]:
    start = max(lhs[0].coordinate, rhs[0].coordinate)
    end = min(lhs[-1].coordinate, rhs[-1].coordinate)
    if not end > start:
        fail("profiles do not overlap in physical coordinate")
    if samples < 2:
        fail("common-grid sample count must be >= 2")
    return [start + (end-start)*index/(samples-1) for index in range(samples)]


def norms(lhs: list[Point], rhs: list[Point], grid: Iterable[float]) -> dict[str, float]:
    grid = list(grid)
    dl_rho = dl_vel = dl_press = l2_rho = 0.0
    for coordinate in grid:
        a = interpolate(lhs, coordinate)
        b = interpolate(rhs, coordinate)
        drho = a.rho - b.rho
        dl_rho += abs(drho)
        dl_vel += abs(a.velocity - b.velocity)
        dl_press += abs(a.pressure - b.pressure)
        l2_rho += drho*drho
    count = float(len(grid))
    return {"l1_density": dl_rho/count, "l1_velocity": dl_vel/count,
            "l1_pressure": dl_press/count, "l2_density": math.sqrt(l2_rho/count)}


def strongest_gradient_coordinate(points: list[Point]) -> float:
    return max((0.5*(points[i-1].coordinate+points[i].coordinate),
                abs((points[i].rho-points[i-1].rho)/(points[i].coordinate-points[i-1].coordinate)))
               for i in range(1, len(points)))[0]


def analyze_pair(numerical_manifest: Path, numerical_csv: Path, reference_manifest: Path,
                 reference_csv: Path, samples: int) -> dict:
    numerical_meta = load_manifest(numerical_manifest)
    reference_meta = load_manifest(reference_manifest)
    require_compatible(numerical_meta, reference_meta)
    numerical = load_profile(numerical_csv, numerical_meta["geometry"])
    reference = load_profile(reference_csv, reference_meta["geometry"])
    grid = common_grid(numerical, reference, samples)
    result = {
        "schema_version": "cosmosim_amr_comparison_result_v1",
        "case_name": numerical_meta["case_name"],
        "geometry": numerical_meta["geometry"],
        "numerical_source": {"name": numerical_meta["source_code_name"], "version": numerical_meta["source_code_version"]},
        "reference_source": {"name": reference_meta["source_code_name"], "version": reference_meta["source_code_version"]},
        "final_time_code": numerical_meta["final_time_code"],
        "gamma": numerical_meta["gamma"],
        "common_grid_samples": len(grid),
        "metrics": norms(numerical, reference, grid),
        "numerical_feature_coordinate": strongest_gradient_coordinate(numerical),
        "reference_feature_coordinate": strongest_gradient_coordinate(reference),
        "provenance": {"numerical": numerical_meta["provenance"], "reference": reference_meta["provenance"]},
    }
    result["feature_coordinate_absolute_difference"] = abs(
        result["numerical_feature_coordinate"] - result["reference_feature_coordinate"])
    return result


def analyze_campaign(campaign_path: Path, samples: int) -> dict:
    campaign = json.loads(campaign_path.read_text(encoding="utf-8"))
    if campaign.get("schema_version") != "cosmosim_amr_hydro_campaign_v1":
        fail("campaign schema_version must be cosmosim_amr_hydro_campaign_v1")
    entries = campaign.get("resolution_runs")
    if not isinstance(entries, list) or len(entries) < 3:
        fail("campaign needs at least three resolution_runs")
    results = []
    for entry in entries:
        result = analyze_pair(Path(entry["numerical_manifest"]), Path(entry["numerical_profile"]),
                              Path(campaign["reference_manifest"]), Path(campaign["reference_profile"]), samples)
        result["effective_resolution"] = int(entry["effective_resolution"])
        results.append(result)
    results.sort(key=lambda item: item["effective_resolution"])
    for previous, current in zip(results, results[1:]):
        ratio = current["effective_resolution"] / previous["effective_resolution"]
        if ratio <= 1:
            fail("campaign resolutions must be strictly increasing")
        current["observed_density_order_from_previous"] = math.log(
            previous["metrics"]["l1_density"] / current["metrics"]["l1_density"]) / math.log(ratio)
    return {"schema_version": "cosmosim_amr_hydro_campaign_result_v1", "campaign": campaign.get("campaign"),
            "results": results, "note": "Metrics only; interpret shock-limited order physically. This tool does not certify science readiness."}


def main() -> int:
    parser = argparse.ArgumentParser()
    subcommands = parser.add_subparsers(dest="command", required=True)
    for command in ("compare", "shock", "sedov"):
        sub = subcommands.add_parser(command)
        sub.add_argument("--numerical-manifest", type=Path, required=True)
        sub.add_argument("--numerical-profile", type=Path, required=True)
        sub.add_argument("--reference-manifest", type=Path, required=True)
        sub.add_argument("--reference-profile", type=Path, required=True)
        sub.add_argument("--output", type=Path, required=True)
        sub.add_argument("--samples", type=int, default=512)
    campaign = subcommands.add_parser("campaign")
    campaign.add_argument("--campaign", type=Path, required=True)
    campaign.add_argument("--output", type=Path, required=True)
    campaign.add_argument("--samples", type=int, default=512)
    args = parser.parse_args()
    if args.command == "campaign":
        result = analyze_campaign(args.campaign, args.samples)
    else:
        result = analyze_pair(args.numerical_manifest, args.numerical_profile,
                              args.reference_manifest, args.reference_profile, args.samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        raise SystemExit(f"validation error: {error}")
