#!/usr/bin/env python3
"""Stream GADGET/AREPO HDF5 ICs into canonical CHUÍ IC v1 plus audit manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sqlite3
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable

import h5py
import numpy as np

KPC_TO_SI = 3.0856775814913673e19
MSUN_TO_SI = 1.98847e30
KM_S_TO_SI = 1.0e3
TYPE_COUNT = 6
ALIASES: dict[str, tuple[str, ...]] = {
    "Coordinates": ("Coordinates", "Position", "POS"),
    "Velocities": ("Velocities", "Velocity", "VEL"),
    "ParticleIDs": ("ParticleIDs", "ParticleID", "ID"),
    "Masses": ("Masses", "Mass"),
    "InternalEnergy": ("InternalEnergy", "U", "Internal_Energy"),
    "Density": ("Density", "Rho"),
    "Metallicity": ("Metallicity", "GFM_Metallicity"),
    "SmoothingLength": ("SmoothingLength", "Hsml", "Smoothing_Length"),
    "StellarFormationTime": ("GFM_StellarFormationTime", "StellarFormationTime", "BirthTime"),
    "InitialMass": ("GFM_InitialMass", "InitialMass", "BirthMass"),
    "BH_Mass": ("BH_Mass",),
    "BH_Mdot": ("BH_Mdot",),
    "ParentParticleIDs": ("ParentParticleIDs", "TracerParentIDs"),
    "InjectionStep": ("InjectionStep",),
    "HostCellIndex": ("HostCellIndex",),
    "MassFractionOfHost": ("MassFractionOfHost",),
    "LastHostMass": ("LastHostMass",),
    "CumulativeExchangedMass": ("CumulativeExchangedMass",),
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, help="first member or single source HDF5")
    p.add_argument("--output", required=True, type=Path, help="canonical CHUÍ HDF5 output")
    p.add_argument("--manifest", required=True, type=Path, help="audit manifest JSON output")
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--source-convention", choices=["gadget_arepo_bridge_v1"])
    source.add_argument("--source-manifest", type=Path, help="strict CHUÍ audit manifest describing the source file set and conversions")
    p.add_argument("--chunk-particles", type=int, default=65536)
    p.add_argument("--source-length-unit-to-si", type=float, default=KPC_TO_SI)
    p.add_argument("--source-mass-unit-to-si", type=float, default=MSUN_TO_SI)
    p.add_argument("--source-velocity-unit-to-si", type=float, default=KM_S_TO_SI)
    p.add_argument("--target-length-unit-to-si", type=float, default=KPC_TO_SI)
    p.add_argument("--target-mass-unit-to-si", type=float, default=MSUN_TO_SI)
    p.add_argument("--target-velocity-unit-to-si", type=float, default=KM_S_TO_SI)
    p.add_argument("--coordinate-frame", choices=["comoving", "physical"])
    p.add_argument(
        "--velocity-convention",
        choices=["physical_peculiar", "sqrt_a_scaled_peculiar", "comoving_coordinate_rate"],
    )
    p.add_argument("--length-h-exponent", type=float, default=0.0)
    p.add_argument("--length-a-exponent", type=float, default=0.0)
    p.add_argument("--mass-h-exponent", type=float, default=0.0)
    p.add_argument("--mass-a-exponent", type=float, default=0.0)
    p.add_argument("--velocity-h-exponent", type=float, default=0.0)
    p.add_argument("--velocity-a-exponent", type=float, default=0.0)
    p.add_argument("--part-type2-policy", choices=["reject", "dark_matter", "star", "black_hole", "tracer"], default="reject")
    p.add_argument("--part-type3-policy", choices=["reject", "dark_matter", "star", "black_hole", "tracer"], default="reject")
    return p.parse_args()


def read_header(path: Path) -> dict[str, Any]:
    with h5py.File(path, "r") as f:
        if "Header" not in f:
            raise ValueError(f"{path}: missing /Header")
        h = f["Header"].attrs
        required = ["NumPart_ThisFile", "NumPart_Total", "MassTable", "Time", "Redshift", "BoxSize", "Omega0", "OmegaLambda", "HubbleParam", "NumFilesPerSnapshot"]
        for key in required:
            if key not in h:
                raise ValueError(f"{path}: missing Header/{key}")
        low = np.asarray(h["NumPart_Total"], dtype=np.uint64)
        high = np.asarray(h.get("NumPart_Total_HighWord", np.zeros(TYPE_COUNT, dtype=np.uint32)), dtype=np.uint64)
        if low.shape != (TYPE_COUNT,) or high.shape != (TYPE_COUNT,):
            raise ValueError(f"{path}: particle total arrays must have length six")
        total = low | (high << np.uint64(32))
        return {
            "this": np.asarray(h["NumPart_ThisFile"], dtype=np.uint64),
            "low": low,
            "high": high,
            "total": total,
            "mass_table": np.asarray(h["MassTable"], dtype=np.float64),
            "time": float(h["Time"]),
            "redshift": float(h["Redshift"]),
            "box": float(h["BoxSize"]),
            "omega0": float(h["Omega0"]),
            "omega_lambda": float(h["OmegaLambda"]),
            "hubble": float(h["HubbleParam"]),
            "num_files": int(h["NumFilesPerSnapshot"]),
        }


def discover_files(first: Path, count: int) -> list[Path]:
    first = first.resolve()
    if count == 1:
        return [first]
    match = re.match(r"^(.*?)(?:\.\d+)?(\.hdf5|\.h5)$", first.name)
    if not match:
        raise ValueError("multifile input must end in .hdf5 or .h5")
    prefix, suffix = match.groups()
    files = [first.parent / f"{prefix}.{i}{suffix}" for i in range(count)]
    missing = [str(path) for path in files if not path.is_file()]
    if missing:
        raise ValueError("missing multifile members: " + ", ".join(missing))
    return files


def close(a: float, b: float) -> bool:
    return abs(a - b) <= 1.0e-10 * max(1.0, abs(a), abs(b))


def validate_headers(files: list[Path]) -> list[dict[str, Any]]:
    headers = [read_header(path) for path in files]
    first = headers[0]
    summed = np.zeros(TYPE_COUNT, dtype=np.uint64)
    for path, h in zip(files, headers, strict=True):
        if h["num_files"] != len(files):
            raise ValueError(f"{path}: NumFilesPerSnapshot mismatch")
        for key in ("total", "high", "mass_table"):
            if not np.array_equal(h[key], first[key]):
                raise ValueError(f"{path}: inconsistent Header/{key}")
        for key in ("time", "redshift", "box", "omega0", "omega_lambda", "hubble"):
            if not close(h[key], first[key]):
                raise ValueError(f"{path}: inconsistent Header/{key}")
        summed += h["this"]
    if not np.array_equal(summed, first["total"]):
        raise ValueError("sum of NumPart_ThisFile does not match 64-bit NumPart_Total")
    if not close(first["redshift"], 1.0 / first["time"] - 1.0):
        raise ValueError("Header/Time and Header/Redshift are inconsistent")
    return headers


def select(group: h5py.Group, canonical: str, required: bool) -> str | None:
    for alias in ALIASES[canonical]:
        if alias in group:
            return alias
    if required:
        raise ValueError(f"{group.name}/{canonical}: required dataset absent")
    return None


def dtype_info(dtype: np.dtype[Any]) -> tuple[str, str, int, bool, str]:
    dtype = np.dtype(dtype)
    if dtype.kind == "f":
        cls, signed, name = "floating_point", True, f"float{dtype.itemsize * 8}"
    elif dtype.kind in "iu":
        signed = dtype.kind == "i"
        cls = "integer"
        name = f"{'int' if signed else 'uint'}{dtype.itemsize * 8}"
    else:
        raise ValueError(f"unsupported HDF5 scalar dtype {dtype}")
    order = {"<": "little_endian", ">": "big_endian", "=": "native", "|": "not_applicable"}.get(dtype.byteorder, "native")
    return name, cls, dtype.itemsize, signed, order


def field_manifest(
    file_index: int,
    canonical_path: str,
    selected_alias: str,
    dataset: h5py.Dataset,
    semantics: str,
    args: argparse.Namespace,
    source_unit: str,
    target_unit: str,
    base_unit: float,
    h_exp: float = 0.0,
    a_exp: float = 0.0,
    velocity_convention: str = "not_velocity",
) -> dict[str, Any]:
    scalar_type, scalar_class, width, signed, order = dtype_info(dataset.dtype)
    return {
        "source_file_index": file_index,
        "dataset_path": canonical_path,
        "selected_alias": selected_alias,
        "scalar_type": scalar_type,
        "scalar_class": scalar_class,
        "byte_width": width,
        "is_signed": signed,
        "byte_order": order,
        "rank": len(dataset.shape),
        "dimensions": list(map(int, dataset.shape)),
        "record_count": int(dataset.shape[0]),
        "base_unit_to_si": base_unit,
        "hubble_exponent": h_exp,
        "scale_factor_exponent": a_exp,
        "coordinate_frame": args.coordinate_frame,
        "velocity_convention": velocity_convention,
        "semantics": semantics,
        "disposition": "converted",
        "source_unit": source_unit,
        "target_unit": target_unit,
        "conversion_equation": "target = stored * base_unit_to_si * h^hubble_exponent * a^scale_factor_exponent / target_si_per_code",
    }


def attribute_manifest(file_index: int, name: str, value: np.ndarray[Any, Any], semantics: str, args: argparse.Namespace, base: float, source: str, target: str, h_exp: float = 0.0, a_exp: float = 0.0) -> dict[str, Any]:
    value = np.asarray(value)
    scalar_type, scalar_class, width, signed, order = dtype_info(value.dtype)
    return {
        "source_file_index": file_index,
        "dataset_path": f"/Header/{name}",
        "selected_alias": name,
        "scalar_type": scalar_type,
        "scalar_class": scalar_class,
        "byte_width": width,
        "is_signed": signed,
        "byte_order": order,
        "rank": value.ndim,
        "dimensions": list(map(int, value.shape)),
        "record_count": 1 if value.ndim == 0 else int(value.shape[0]),
        "base_unit_to_si": base,
        "hubble_exponent": h_exp,
        "scale_factor_exponent": a_exp,
        "coordinate_frame": args.coordinate_frame,
        "velocity_convention": "not_velocity",
        "semantics": semantics,
        "disposition": "converted",
        "source_unit": source,
        "target_unit": target,
        "conversion_equation": "target = stored * base_unit_to_si * h^hubble_exponent * a^scale_factor_exponent / target_si_per_code",
    }


def species_policies(args: argparse.Namespace) -> list[str]:
    def map_family(value: str, family: int) -> str:
        if value == "dark_matter":
            return f"family{family}_as_dark_matter"
        return {"reject": "reject", "star": "star", "black_hole": "black_hole", "tracer": "tracer"}[value]
    return ["gas", "dark_matter", map_family(args.part_type2_policy, 2), map_family(args.part_type3_policy, 3), "star", "black_hole"]


def inspect_sources(files: list[Path], headers: list[dict[str, Any]], args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[str], list[str], list[str]]:
    fields: list[dict[str, Any]] = []
    defaults: list[str] = []
    dropped: list[str] = []
    warnings: list[str] = []
    policies = species_policies(args)
    for file_index, (path, header) in enumerate(zip(files, headers, strict=True)):
        with h5py.File(path, "r") as f:
            fields.append(attribute_manifest(file_index, "MassTable", np.asarray(f["Header"].attrs["MassTable"]), "extensive", args, args.source_mass_unit_to_si, "source_mass", "runtime_mass", args.mass_h_exponent, args.mass_a_exponent))
            fields.append(attribute_manifest(file_index, "BoxSize", np.asarray(f["Header"].attrs["BoxSize"]), "coordinate", args, args.source_length_unit_to_si, "source_length", "runtime_length", args.length_h_exponent, args.length_a_exponent))
            for ptype in range(TYPE_COUNT):
                count = int(header["this"][ptype])
                if not count:
                    continue
                if policies[ptype] == "reject":
                    raise ValueError(f"{path}: populated PartType{ptype} has reject policy")
                group_name = f"PartType{ptype}"
                if group_name not in f:
                    raise ValueError(f"{path}: missing /{group_name}")
                group = f[group_name]
                required = {"Coordinates": True, "Velocities": True, "ParticleIDs": True}
                selected: dict[str, str | None] = {}
                for canonical, is_required in required.items():
                    selected[canonical] = select(group, canonical, is_required)
                selected["Masses"] = select(group, "Masses", False)
                if selected["Masses"] is None and header["mass_table"][ptype] <= 0:
                    raise ValueError(f"{group.name}: Masses absent and MassTable is zero")
                field_specs = [
                    ("Coordinates", "coordinate", args.source_length_unit_to_si, "source_length", "runtime_length", args.length_h_exponent, args.length_a_exponent, "not_velocity"),
                    ("Velocities", "velocity", args.source_velocity_unit_to_si, "source_velocity", "runtime_velocity", args.velocity_h_exponent, args.velocity_a_exponent, args.velocity_convention),
                    ("ParticleIDs", "identifier", 1.0, "identifier", "identifier", 0.0, 0.0, "not_velocity"),
                    ("Masses", "extensive", args.source_mass_unit_to_si, "source_mass", "runtime_mass", args.mass_h_exponent, args.mass_a_exponent, "not_velocity"),
                ]
                for canonical, semantics, base, src, dst, h_exp, a_exp, vel in field_specs:
                    alias = selected.get(canonical)
                    if alias is None:
                        continue
                    ds = group[alias]
                    if ds.shape[0] != count:
                        raise ValueError(f"{ds.name}: first dimension disagrees with NumPart_ThisFile")
                    if canonical in ("Coordinates", "Velocities") and ds.shape != (count, 3):
                        raise ValueError(f"{ds.name}: expected [N,3]")
                    if canonical not in ("Coordinates", "Velocities") and ds.shape != (count,):
                        raise ValueError(f"{ds.name}: expected [N]")
                    fields.append(field_manifest(file_index, f"/{group_name}/{canonical}", alias, ds, semantics, args, src, dst, base, h_exp, a_exp, vel))
                optional: list[tuple[str, str, float, str, str]] = []
                if ptype == 0:
                    optional += [
                        ("InternalEnergy", "specific", args.source_velocity_unit_to_si**2, "source_velocity_squared", "runtime_specific_energy"),
                        ("Density", "intensive", args.source_mass_unit_to_si / args.source_length_unit_to_si**3, "source_density", "runtime_density"),
                        ("Metallicity", "intensive", 1.0, "mass_fraction", "mass_fraction"),
                        ("SmoothingLength", "coordinate", args.source_length_unit_to_si, "source_length", "runtime_length"),
                    ]
                if policies[ptype] == "star":
                    optional += [
                        ("StellarFormationTime", "intensive", 1.0, "scale_factor", "scale_factor"),
                        ("InitialMass", "extensive", args.source_mass_unit_to_si, "source_mass", "runtime_mass"),
                        ("Metallicity", "intensive", 1.0, "mass_fraction", "mass_fraction"),
                    ]
                if policies[ptype] == "black_hole":
                    optional += [
                        ("BH_Mass", "extensive", args.source_mass_unit_to_si, "source_mass", "runtime_mass"),
                        ("BH_Mdot", "intensive", args.source_mass_unit_to_si / (args.source_length_unit_to_si / args.source_velocity_unit_to_si), "source_mass/source_time", "runtime_mass/runtime_time"),
                    ]
                if policies[ptype] == "tracer":
                    optional += [
                        ("ParentParticleIDs", "identifier", 1.0, "identifier", "identifier"),
                        ("InjectionStep", "identifier", 1.0, "step", "step"),
                        ("HostCellIndex", "identifier", 1.0, "index", "index"),
                        ("MassFractionOfHost", "intensive", 1.0, "fraction", "fraction"),
                        ("LastHostMass", "extensive", args.source_mass_unit_to_si, "source_mass", "runtime_mass"),
                        ("CumulativeExchangedMass", "extensive", args.source_mass_unit_to_si, "source_mass", "runtime_mass"),
                    ]
                seen: set[str] = {alias for alias in selected.values() if alias}
                for canonical, semantics, base, src, dst in optional:
                    alias = select(group, canonical, False)
                    required_for_sidecar = policies[ptype] == "tracer" or (policies[ptype] == "black_hole" and canonical == "BH_Mass")
                    if alias is None:
                        if required_for_sidecar:
                            raise ValueError(f"{group.name}/{canonical}: required by selected sidecar policy")
                        default = {
                            "InternalEnergy": "zero", "Density": "zero", "StellarFormationTime": "Header/Time",
                            "InitialMass": "particle_mass", "Metallicity": "zero", "BH_Mdot": "zero",
                        }.get(canonical)
                        if default:
                            defaults.append(f"/{group_name}/{canonical}={default}")
                        continue
                    seen.add(alias)
                    ds = group[alias]
                    if ds.shape != (count,):
                        raise ValueError(f"{ds.name}: expected [N]")
                    h_exp = args.mass_h_exponent if canonical in ("InitialMass", "BH_Mass", "LastHostMass", "CumulativeExchangedMass") else 0.0
                    a_exp = args.mass_a_exponent if h_exp else 0.0
                    fields.append(field_manifest(file_index, f"/{group_name}/{canonical}", alias, ds, semantics, args, src, dst, base, h_exp, a_exp))
                    if ptype == 0 and canonical in ("Metallicity", "SmoothingLength"):
                        dropped.append(f"/{group_name}/{canonical}: canonical source retained only in audit; runtime has no target lane")
                for name in group.keys():
                    if name not in seen:
                        dropped.append(f"/{group_name}/{name}: unrecognized source field")
    if len(files) > 1:
        baseline = {field["dataset_path"]: field for field in fields if field["source_file_index"] == 0}
        for file_index in range(1, len(files)):
            current = {field["dataset_path"]: field for field in fields if field["source_file_index"] == file_index}
            if set(current) != set(baseline):
                raise ValueError(f"source file {file_index}: inconsistent dataset schema")
            for path, expected in baseline.items():
                actual = current[path]
                expected_tail = expected["dimensions"] if path.startswith("/Header/") else expected["dimensions"][1:]
                actual_tail = actual["dimensions"] if path.startswith("/Header/") else actual["dimensions"][1:]
                signature = ("selected_alias", "scalar_type", "scalar_class", "byte_width", "is_signed", "byte_order", "rank")
                if any(actual[key] != expected[key] for key in signature) or actual_tail != expected_tail:
                    raise ValueError(f"source file {file_index}: inconsistent dataset schema for {path}")
    if dropped:
        warnings.append("Source contains fields not represented in canonical CHUÍ runtime state; see dropped_fields.")
    return fields, defaults, dropped, warnings


def factor(base: float, h_exp: float, a_exp: float, h: float, a: float, target: float) -> float:
    return base * h**h_exp * a**a_exp / target



def hash_file(path: Path, chunk_bytes: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(chunk_bytes):
            digest.update(block)
    return digest.hexdigest()



def policy_argument(value: str, family: int) -> str:
    mapping = {
        "reject": "reject",
        f"family{family}_as_dark_matter": "dark_matter",
        "star": "star",
        "black_hole": "black_hole",
        "tracer": "tracer",
    }
    if value not in mapping:
        raise ValueError(f"source manifest has unsupported PartType{family} policy: {value}")
    return mapping[value]


def resolve_manifest_source(path: Path, value: str) -> Path:
    source = Path(value)
    return source.resolve() if source.is_absolute() else (path.parent / source).resolve()


def load_source_manifest(path: Path, args: argparse.Namespace) -> dict[str, Any]:
    path = path.resolve()
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read source manifest {path}: {exc}") from exc
    if not isinstance(manifest, dict):
        raise ValueError("source manifest root must be a JSON object")
    if manifest.get("schema_name") != "chui_ic_audit_manifest" or manifest.get("schema_version") != 2:
        raise ValueError("source manifest must use chui_ic_audit_manifest schema version 2")
    if manifest.get("dialect") != "gadget_arepo_bridge_v1" or manifest.get("dialect_version") != "1":
        raise ValueError("source manifest dialect must be gadget_arepo_bridge_v1 version 1")
    source_files = manifest.get("source_files")
    source_hashes = manifest.get("source_sha256")
    source_sizes = manifest.get("source_file_sizes_bytes")
    fields = manifest.get("fields")
    policies = manifest.get("species_policy")
    parameters = manifest.get("conversion_parameters")
    if not isinstance(source_files, list) or not source_files or not all(isinstance(value, str) and value for value in source_files):
        raise ValueError("source manifest source_files must be a non-empty string array")
    if not isinstance(source_hashes, list) or len(source_hashes) != len(source_files):
        raise ValueError("source manifest source_sha256 must match source_files")
    if not isinstance(source_sizes, list) or len(source_sizes) != len(source_files):
        raise ValueError("source manifest source_file_sizes_bytes must match source_files")
    if not isinstance(fields, list) or not fields:
        raise ValueError("source manifest fields must be a non-empty array")
    if not isinstance(policies, list) or len(policies) != TYPE_COUNT:
        raise ValueError("source manifest species_policy must contain six entries")
    if policies[0] != "gas" or policies[1] != "dark_matter" or policies[4] != "star" or policies[5] != "black_hole":
        raise ValueError("source manifest has invalid canonical species policies")
    if not isinstance(parameters, dict):
        raise ValueError("source manifest is missing conversion_parameters")
    required_parameters = (
        "coordinate_frame", "velocity_convention",
        "source_length_unit_to_si", "source_mass_unit_to_si", "source_velocity_unit_to_si",
        "target_length_unit_to_si", "target_mass_unit_to_si", "target_velocity_unit_to_si",
        "length_h_exponent", "length_a_exponent", "mass_h_exponent", "mass_a_exponent",
        "velocity_h_exponent", "velocity_a_exponent",
    )
    missing = [key for key in required_parameters if key not in parameters]
    if missing:
        raise ValueError("source manifest conversion_parameters missing: " + ", ".join(missing))
    args.source_convention = "gadget_arepo_bridge_v1"
    args.input = resolve_manifest_source(path, source_files[0])
    args.coordinate_frame = parameters["coordinate_frame"]
    args.velocity_convention = parameters["velocity_convention"]
    if args.coordinate_frame not in ("comoving", "physical"):
        raise ValueError("source manifest has unsupported coordinate_frame")
    if args.velocity_convention not in ("physical_peculiar", "sqrt_a_scaled_peculiar", "comoving_coordinate_rate"):
        raise ValueError("source manifest has unsupported velocity_convention")
    for key in required_parameters[2:]:
        value = parameters[key]
        if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            raise ValueError(f"source manifest conversion parameter {key} must be finite")
        setattr(args, key, float(value))
    args.part_type2_policy = policy_argument(str(policies[2]), 2)
    args.part_type3_policy = policy_argument(str(policies[3]), 3)
    manifest["_resolved_source_files"] = [resolve_manifest_source(path, value) for value in source_files]
    manifest["_source_manifest_path"] = path
    return manifest


def field_signature(field: dict[str, Any]) -> tuple[Any, ...]:
    required = (
        "source_file_index", "dataset_path", "selected_alias", "scalar_type", "scalar_class",
        "byte_width", "is_signed", "byte_order", "rank", "dimensions", "record_count",
    )
    missing = [key for key in required if key not in field]
    if missing:
        raise ValueError("source manifest field missing: " + ", ".join(missing))
    dimensions = field["dimensions"]
    if not isinstance(dimensions, list):
        raise ValueError("source manifest field dimensions must be an array")
    return tuple(field[key] if key != "dimensions" else tuple(dimensions) for key in required)


def validate_supplied_manifest(
    supplied: dict[str, Any],
    files: list[Path],
    headers: list[dict[str, Any]],
    fields: list[dict[str, Any]],
    hashes: list[str],
    policies: list[str],
) -> None:
    resolved_files = supplied["_resolved_source_files"]
    if [path.resolve() for path in files] != resolved_files:
        raise ValueError("discovered source-file order does not match the supplied manifest")
    if int(supplied.get("num_files_per_snapshot", -1)) != len(files):
        raise ValueError("source manifest NumFilesPerSnapshot does not match the source set")
    expected_sizes = [int(value) for value in supplied["source_file_sizes_bytes"]]
    actual_sizes = [path.stat().st_size for path in files]
    if actual_sizes != expected_sizes:
        raise ValueError("source file size differs from the supplied manifest")
    if hashes != supplied["source_sha256"]:
        raise ValueError("source file SHA-256 differs from the supplied manifest")
    if policies != supplied["species_policy"]:
        raise ValueError("resolved species policies differ from the supplied manifest")
    if [list(map(int, header["this"])) for header in headers] != supplied.get("num_part_this_file"):
        raise ValueError("source NumPart_ThisFile differs from the supplied manifest")
    if list(map(int, headers[0]["total"])) != supplied.get("num_part_total"):
        raise ValueError("source NumPart_Total differs from the supplied manifest")
    if list(map(int, headers[0]["high"])) != supplied.get("num_part_total_high_word"):
        raise ValueError("source NumPart_Total_HighWord differs from the supplied manifest")
    expected_fields = {field_signature(field): field for field in supplied["fields"]}
    actual_fields = {field_signature(field): field for field in fields}
    if set(actual_fields) != set(expected_fields):
        missing = len(set(expected_fields) - set(actual_fields))
        unexpected = len(set(actual_fields) - set(expected_fields))
        raise ValueError(f"actual HDF5 schema differs from the supplied manifest (missing={missing}, unexpected={unexpected})")

def id_blob(value: int) -> sqlite3.Binary:
    if value < 0 or value > 0xFFFFFFFFFFFFFFFF:
        raise ValueError(f"particle ID is outside uint64 range: {value}")
    return sqlite3.Binary(value.to_bytes(8, byteorder="little", signed=False))

def velocity_extra(convention: str, a: float) -> float:
    return {"physical_peculiar": 1.0, "sqrt_a_scaled_peculiar": 1.0 / math.sqrt(a), "comoving_coordinate_rate": a}[convention]


def chunks(count: int, size: int) -> Iterable[tuple[int, int]]:
    for start in range(0, count, size):
        yield start, min(size, count - start)


def canonical_dataset_plan(ptype: int, policy: str) -> list[tuple[str, np.dtype[Any]]]:
    plan: list[tuple[str, np.dtype[Any]]] = [
        ("Coordinates", np.dtype("<f8")), ("Velocities", np.dtype("<f8")),
        ("Masses", np.dtype("<f8")), ("ParticleIDs", np.dtype("<u8")),
    ]
    if ptype == 0:
        plan += [("InternalEnergy", np.dtype("<f8")), ("Density", np.dtype("<f8"))]
    if policy == "star":
        plan += [("StellarFormationTime", np.dtype("<f8")), ("InitialMass", np.dtype("<f8")), ("Metallicity", np.dtype("<f8"))]
    if policy == "black_hole":
        plan += [("BH_Mass", np.dtype("<f8")), ("BH_Mdot", np.dtype("<f8"))]
    if policy == "tracer":
        plan += [("ParentParticleIDs", np.dtype("<u8")), ("InjectionStep", np.dtype("<u8")), ("HostCellIndex", np.dtype("<u8")), ("MassFractionOfHost", np.dtype("<f8")), ("LastHostMass", np.dtype("<f8")), ("CumulativeExchangedMass", np.dtype("<f8"))]
    return plan


def convert(files: list[Path], headers: list[dict[str, Any]], args: argparse.Namespace, manifest: dict[str, Any]) -> None:
    total = headers[0]["total"].astype(np.uint64)
    policies = species_policies(args)
    tmp_output = args.output.with_suffix(args.output.suffix + ".part")
    tmp_manifest = args.manifest.with_suffix(args.manifest.suffix + ".part")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    for path in (tmp_output, tmp_manifest):
        path.unlink(missing_ok=True)
    db_fd, db_name = tempfile.mkstemp(prefix="chui_ic_ids_", suffix=".sqlite", dir=args.output.parent)
    os.close(db_fd)
    db = sqlite3.connect(db_name)
    try:
        db.execute("PRAGMA journal_mode=OFF")
        db.execute("PRAGMA synchronous=OFF")
        db.execute("CREATE TABLE ids(id BLOB PRIMARY KEY) WITHOUT ROWID")
        with h5py.File(tmp_output, "w") as out:
            header = out.create_group("Header")
            low = (total & np.uint64(0xFFFFFFFF)).astype(np.uint32)
            high = (total >> np.uint64(32)).astype(np.uint32)
            header.attrs["NumPart_ThisFile"] = low
            header.attrs["NumPart_Total"] = low
            header.attrs["NumPart_Total_HighWord"] = high
            header.attrs["MassTable"] = np.zeros(TYPE_COUNT, dtype=np.float64)
            header.attrs["Time"] = headers[0]["time"]
            header.attrs["Redshift"] = headers[0]["redshift"]
            box_factor = factor(args.source_length_unit_to_si, args.length_h_exponent, args.length_a_exponent, headers[0]["hubble"], headers[0]["time"], args.target_length_unit_to_si)
            if args.coordinate_frame == "physical":
                box_factor /= headers[0]["time"]
            header.attrs["BoxSize"] = headers[0]["box"] * box_factor
            header.attrs["Omega0"] = headers[0]["omega0"]
            header.attrs["OmegaLambda"] = headers[0]["omega_lambda"]
            header.attrs["HubbleParam"] = headers[0]["hubble"]
            header.attrs["NumFilesPerSnapshot"] = np.uint32(1)
            header.attrs["ChuiIcSchemaName"] = "chui_canonical_v1"
            header.attrs["ChuiIcSchemaVersion"] = np.uint32(1)
            header.attrs["ChuiLengthUnitToSI"] = args.target_length_unit_to_si
            header.attrs["ChuiMassUnitToSI"] = args.target_mass_unit_to_si
            header.attrs["ChuiVelocityUnitToSI"] = args.target_velocity_unit_to_si
            header.attrs["ChuiCoordinateFrame"] = "comoving"
            header.attrs["ChuiVelocityConvention"] = "physical_peculiar"
            header.attrs["ConversionManifestSha256"] = "0" * 64
            outputs: dict[int, dict[str, h5py.Dataset]] = {}
            for ptype in range(TYPE_COUNT):
                n = int(total[ptype])
                if not n:
                    continue
                group = out.create_group(f"PartType{ptype}")
                outputs[ptype] = {}
                for name, dtype in canonical_dataset_plan(ptype, policies[ptype]):
                    shape = (n, 3) if name in ("Coordinates", "Velocities") else (n,)
                    chunk_shape = (min(args.chunk_particles, n), 3) if len(shape) == 2 else (min(args.chunk_particles, n),)
                    outputs[ptype][name] = group.create_dataset(name, shape=shape, dtype=dtype, chunks=chunk_shape)
            offsets = np.zeros(TYPE_COUNT, dtype=np.uint64)
            for path, source_header in zip(files, headers, strict=True):
                with h5py.File(path, "r") as src:
                    for ptype in range(TYPE_COUNT):
                        count = int(source_header["this"][ptype])
                        if not count:
                            continue
                        group = src[f"PartType{ptype}"]
                        aliases = {name: select(group, name, False) for name, _ in canonical_dataset_plan(ptype, policies[ptype])}
                        aliases["Coordinates"] = select(group, "Coordinates", True)
                        aliases["Velocities"] = select(group, "Velocities", True)
                        aliases["ParticleIDs"] = select(group, "ParticleIDs", True)
                        aliases["Masses"] = select(group, "Masses", False)
                        for start, amount in chunks(count, args.chunk_particles):
                            dst = int(offsets[ptype]) + start
                            ids = np.asarray(group[aliases["ParticleIDs"]][start:start+amount], dtype=np.uint64)
                            if np.any(ids == 0):
                                raise ValueError(f"{path}: particle IDs must be nonzero")
                            try:
                                db.executemany("INSERT INTO ids(id) VALUES (?)", ((id_blob(int(x)),) for x in ids))
                            except sqlite3.IntegrityError as exc:
                                raise ValueError(f"duplicate particle ID detected across source files: {exc}") from exc
                            pos = np.asarray(group[aliases["Coordinates"]][start:start+amount], dtype=np.float64)
                            vel = np.asarray(group[aliases["Velocities"]][start:start+amount], dtype=np.float64)
                            lf = factor(args.source_length_unit_to_si, args.length_h_exponent, args.length_a_exponent, source_header["hubble"], source_header["time"], args.target_length_unit_to_si)
                            if args.coordinate_frame == "physical":
                                lf /= source_header["time"]
                            vf = factor(args.source_velocity_unit_to_si, args.velocity_h_exponent, args.velocity_a_exponent, source_header["hubble"], source_header["time"], args.target_velocity_unit_to_si) * velocity_extra(args.velocity_convention, source_header["time"])
                            pos *= lf
                            vel *= vf
                            if not np.all(np.isfinite(pos)) or not np.all(np.isfinite(vel)):
                                raise ValueError(f"{path}: non-finite converted coordinates/velocities")
                            mass_alias = aliases["Masses"]
                            if mass_alias:
                                masses = np.asarray(group[mass_alias][start:start+amount], dtype=np.float64)
                            else:
                                masses = np.full(amount, source_header["mass_table"][ptype], dtype=np.float64)
                            masses *= factor(args.source_mass_unit_to_si, args.mass_h_exponent, args.mass_a_exponent, source_header["hubble"], source_header["time"], args.target_mass_unit_to_si)
                            if np.any(~np.isfinite(masses)) or np.any(masses <= 0):
                                raise ValueError(f"{path}: masses must be finite and positive")
                            outputs[ptype]["Coordinates"][dst:dst+amount] = pos
                            outputs[ptype]["Velocities"][dst:dst+amount] = vel
                            outputs[ptype]["ParticleIDs"][dst:dst+amount] = ids
                            outputs[ptype]["Masses"][dst:dst+amount] = masses
                            def scalar(name: str, default: float | np.ndarray[Any, Any], multiplier: float = 1.0, dtype: Any = np.float64) -> np.ndarray[Any, Any]:
                                alias = aliases.get(name)
                                values = np.asarray(group[alias][start:start+amount], dtype=dtype) if alias else np.asarray(default, dtype=dtype)
                                if values.ndim == 0:
                                    values = np.full(amount, values, dtype=dtype)
                                return values * multiplier
                            if ptype == 0:
                                outputs[ptype]["InternalEnergy"][dst:dst+amount] = scalar("InternalEnergy", 0.0, vf*vf)
                                density_factor = factor(args.source_mass_unit_to_si, args.mass_h_exponent, args.mass_a_exponent, source_header["hubble"], source_header["time"], args.target_mass_unit_to_si) / lf**3
                                outputs[ptype]["Density"][dst:dst+amount] = scalar("Density", 0.0, density_factor)
                            if policies[ptype] == "star":
                                outputs[ptype]["StellarFormationTime"][dst:dst+amount] = scalar("StellarFormationTime", source_header["time"])
                                outputs[ptype]["InitialMass"][dst:dst+amount] = scalar("InitialMass", masses, factor(args.source_mass_unit_to_si, args.mass_h_exponent, args.mass_a_exponent, source_header["hubble"], source_header["time"], args.target_mass_unit_to_si)) if aliases.get("InitialMass") else masses
                                outputs[ptype]["Metallicity"][dst:dst+amount] = scalar("Metallicity", 0.0)
                            if policies[ptype] == "black_hole":
                                outputs[ptype]["BH_Mass"][dst:dst+amount] = scalar("BH_Mass", 0.0, factor(args.source_mass_unit_to_si, args.mass_h_exponent, args.mass_a_exponent, source_header["hubble"], source_header["time"], args.target_mass_unit_to_si))
                                mdot_factor = factor(args.source_mass_unit_to_si, args.mass_h_exponent, args.mass_a_exponent, source_header["hubble"], source_header["time"], args.target_mass_unit_to_si) / (args.source_length_unit_to_si / args.source_velocity_unit_to_si) * (args.target_length_unit_to_si / args.target_velocity_unit_to_si)
                                outputs[ptype]["BH_Mdot"][dst:dst+amount] = scalar("BH_Mdot", 0.0, mdot_factor)
                            if policies[ptype] == "tracer":
                                for name in ("ParentParticleIDs", "InjectionStep", "HostCellIndex"):
                                    outputs[ptype][name][dst:dst+amount] = scalar(name, 0, dtype=np.uint64)
                                outputs[ptype]["MassFractionOfHost"][dst:dst+amount] = scalar("MassFractionOfHost", 0.0)
                                mf = factor(args.source_mass_unit_to_si, args.mass_h_exponent, args.mass_a_exponent, source_header["hubble"], source_header["time"], args.target_mass_unit_to_si)
                                outputs[ptype]["LastHostMass"][dst:dst+amount] = scalar("LastHostMass", 0.0, mf)
                                outputs[ptype]["CumulativeExchangedMass"][dst:dst+amount] = scalar("CumulativeExchangedMass", 0.0, mf)
                        offsets[ptype] += count
            if not np.array_equal(offsets, total):
                raise ValueError("canonical output write coverage does not match global counts")
        manifest["canonical_output_file"] = str(args.output.resolve())
        manifest_text = json.dumps(manifest, indent=2, sort_keys=False) + "\n"
        manifest_hash = hashlib.sha256(manifest_text.encode()).hexdigest()
        with h5py.File(tmp_output, "r+") as finalized_output:
            finalized_output["Header"].attrs.modify(
                "ConversionManifestSha256", manifest_hash)
        tmp_manifest.write_text(manifest_text, encoding="utf-8")
        os.replace(tmp_output, args.output)
        os.replace(tmp_manifest, args.manifest)
    finally:
        db.close()
        Path(db_name).unlink(missing_ok=True)
        tmp_output.unlink(missing_ok=True)
        tmp_manifest.unlink(missing_ok=True)


def main() -> int:
    args = parse_args()
    supplied_manifest: dict[str, Any] | None = None
    if args.source_manifest is not None:
        if args.source_manifest.resolve() == args.manifest.resolve():
            raise ValueError("--source-manifest and --manifest must be different files")
        supplied_manifest = load_source_manifest(args.source_manifest, args)
    else:
        if args.input is None:
            raise ValueError("--input is required with --source-convention")
        if args.coordinate_frame is None:
            raise ValueError("--coordinate-frame is required with --source-convention")
        if args.velocity_convention is None:
            raise ValueError("--velocity-convention is required with --source-convention")
    if args.chunk_particles <= 0:
        raise ValueError("--chunk-particles must be positive")
    for name in ("source_length_unit_to_si", "source_mass_unit_to_si", "source_velocity_unit_to_si", "target_length_unit_to_si", "target_mass_unit_to_si", "target_velocity_unit_to_si"):
        if not math.isfinite(getattr(args, name)) or getattr(args, name) <= 0:
            raise ValueError(f"--{name.replace('_', '-')} must be finite and positive")
    first_header = read_header(args.input)
    files = discover_files(args.input, first_header["num_files"])
    if args.output.resolve() in {path.resolve() for path in files}:
        raise ValueError("--output must not replace a source IC member")
    headers = validate_headers(files)
    fields, defaults, dropped, warnings = inspect_sources(files, headers, args)
    hashes = [hash_file(path) for path in files]
    policies = species_policies(args)
    if supplied_manifest is not None:
        validate_supplied_manifest(supplied_manifest, files, headers, fields, hashes, policies)
    original_headers = []
    for h in headers:
        original_headers.append(
            ";".join([
                f"NumFilesPerSnapshot={h['num_files']}",
                "NumPart_ThisFile=" + ",".join(map(str, h["this"])),
                "NumPart_Total=" + ",".join(map(str, h["total"])),
                "NumPart_Total_HighWord=" + ",".join(map(str, h["high"])),
                "MassTable=" + ",".join(map(repr, h["mass_table"])),
                f"BoxSize={h['box']}", f"Time={h['time']}", f"Redshift={h['redshift']}",
                f"Omega0={h['omega0']}", f"OmegaLambda={h['omega_lambda']}", f"HubbleParam={h['hubble']}",
            ])
        )
    equations = sorted({field["conversion_equation"] for field in fields})
    manifest: dict[str, Any] = {
        "schema_name": "chui_ic_audit_manifest",
        "schema_version": 2,
        "converter_version": "tools/convert_ic.py_v1",
        "dialect": args.source_convention,
        "dialect_version": "1",
        "num_files_per_snapshot": len(files),
        "source_files": [str(path.resolve()) for path in files],
        "source_provenance_ids": [f"sha256:{value}" for value in hashes],
        "source_file_sizes_bytes": [path.stat().st_size for path in files],
        "source_sha256": hashes,
        "original_header_attributes": original_headers,
        "num_part_this_file": [list(map(int, h["this"])) for h in headers],
        "num_part_total": list(map(int, headers[0]["total"])),
        "num_part_total_high_word": list(map(int, headers[0]["high"])),
        "mass_table": list(map(float, headers[0]["mass_table"])),
        "box_size": headers[0]["box"], "time": headers[0]["time"], "redshift": headers[0]["redshift"],
        "omega_matter": headers[0]["omega0"], "omega_lambda": headers[0]["omega_lambda"], "hubble_param": headers[0]["hubble"],
        "species_policy": policies,
        "fields": fields,
        "defaulted_fields": defaults,
        "converted_fields": [field["dataset_path"] for field in fields],
        "dropped_fields": dropped,
        "rejected_fields": [],
        "preserved_auxiliary_fields": [],
        "conversion_equations": equations,
        "warnings": warnings,
        "conversion_parameters": {
            "coordinate_frame": args.coordinate_frame,
            "velocity_convention": args.velocity_convention,
            "source_length_unit_to_si": args.source_length_unit_to_si,
            "source_mass_unit_to_si": args.source_mass_unit_to_si,
            "source_velocity_unit_to_si": args.source_velocity_unit_to_si,
            "target_length_unit_to_si": args.target_length_unit_to_si,
            "target_mass_unit_to_si": args.target_mass_unit_to_si,
            "target_velocity_unit_to_si": args.target_velocity_unit_to_si,
            "length_h_exponent": args.length_h_exponent,
            "length_a_exponent": args.length_a_exponent,
            "mass_h_exponent": args.mass_h_exponent,
            "mass_a_exponent": args.mass_a_exponent,
            "velocity_h_exponent": args.velocity_h_exponent,
            "velocity_a_exponent": args.velocity_a_exponent,
        },
    }
    if supplied_manifest is not None:
        source_manifest_path = supplied_manifest["_source_manifest_path"]
        manifest["source_manifest_file"] = str(source_manifest_path)
        manifest["source_manifest_sha256"] = hash_file(source_manifest_path)
    convert(files, headers, args, manifest)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # fail closed with one actionable diagnostic
        print(f"convert_ic.py: error: {exc}", file=sys.stderr)
        raise SystemExit(2)
