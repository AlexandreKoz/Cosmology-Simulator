#!/usr/bin/env python3
"""End-to-end regression for the streaming external-to-canonical IC converter."""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError:
    raise SystemExit(77)


def write_member(path: Path, member: int, duplicate: bool = False) -> None:
    local = np.zeros(6, dtype=np.uint32)
    total = np.zeros(6, dtype=np.uint32)
    local[1] = 2
    total[1] = 4
    with h5py.File(path, "w") as f:
        h = f.create_group("Header")
        h.attrs["NumPart_ThisFile"] = local
        h.attrs["NumPart_Total"] = total
        h.attrs["NumPart_Total_HighWord"] = np.zeros(6, dtype=np.uint32)
        h.attrs["MassTable"] = np.zeros(6, dtype=np.float64)
        h.attrs["Time"] = 0.5
        h.attrs["Redshift"] = 1.0
        h.attrs["BoxSize"] = 50.0
        h.attrs["Omega0"] = 0.315
        h.attrs["OmegaLambda"] = 0.685
        h.attrs["HubbleParam"] = 0.5
        h.attrs["NumFilesPerSnapshot"] = np.uint32(2)
        g = f.create_group("PartType1")
        first = member * 2
        g.create_dataset(
            "Coordinates",
            data=np.array(
                [[1.0 + first, 2.0, 3.0], [2.0 + first, 2.0, 3.0]],
                dtype=np.float32,
            ),
        )
        g.create_dataset(
            "Velocities",
            data=np.array([[4.0, 0.0, 0.0], [8.0, 0.0, 0.0]], dtype=np.float32),
        )
        g.create_dataset("Masses", data=np.array([2.0, 3.0], dtype=np.float32))
        ids = np.array(
            [0x8000000000000001 + first, 0x8000000000000002 + first],
            dtype=np.uint64,
        )
        if duplicate and member == 1:
            ids[1] = np.uint64(0x8000000000000001)
        g.create_dataset("ParticleIDs", data=ids)


def run_converter(converter: Path, source: Path, output: Path, manifest: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(converter),
            "--input",
            str(source),
            "--output",
            str(output),
            "--manifest",
            str(manifest),
            "--source-convention",
            "gadget_arepo_bridge_v1",
            "--coordinate-frame",
            "physical",
            "--velocity-convention",
            "sqrt_a_scaled_peculiar",
            "--length-h-exponent",
            "-1",
            "--mass-h-exponent",
            "-1",
            "--chunk-particles",
            "1",
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def run_converter_from_manifest(
    converter: Path, source_manifest: Path, output: Path, manifest: Path
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(converter),
            "--source-manifest",
            str(source_manifest),
            "--output",
            str(output),
            "--manifest",
            str(manifest),
            "--chunk-particles",
            "2",
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_convert_ic.py <converter>")
    converter = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="chui_convert_ic_test_") as temp:
        root = Path(temp)
        first = root / "source.0.hdf5"
        second = root / "source.1.hdf5"
        write_member(first, 0)
        write_member(second, 1)
        output = root / "canonical.hdf5"
        manifest_path = root / "canonical.audit.json"
        completed = run_converter(converter, first, output, manifest_path)
        assert completed.returncode == 0, completed.stderr
        manifest_text = manifest_path.read_text(encoding="utf-8")
        manifest = json.loads(manifest_text)
        assert manifest["schema_name"] == "chui_ic_audit_manifest"
        assert manifest["schema_version"] == 2
        assert len(manifest["source_files"]) == 2
        assert all(len(value) == 64 for value in manifest["source_sha256"])
        assert any(
            field["dataset_path"] == "/PartType1/Coordinates"
            and field["scalar_type"] == "float32"
            and field["dimensions"] == [2, 3]
            for field in manifest["fields"]
        )
        assert any(
            field["dataset_path"] == "/Header/BoxSize"
            and field["rank"] == 0
            and field["dimensions"] == []
            for field in manifest["fields"]
        )
        with h5py.File(output, "r") as f:
            header = f["Header"].attrs
            assert header["ChuiIcSchemaName"] == "chui_canonical_v1"
            assert header["NumFilesPerSnapshot"] == 1
            assert header["ConversionManifestSha256"] == hashlib.sha256(
                manifest_text.encode("utf-8")
            ).hexdigest()
            ids = f["PartType1/ParticleIDs"][:]
            assert ids.dtype == np.dtype("uint64")
            assert int(ids[0]) > 2**63
            coordinates = f["PartType1/Coordinates"][:]
            velocities = f["PartType1/Velocities"][:]
            masses = f["PartType1/Masses"][:]
            # physical -> comoving divides by a, while h^-1 doubles again.
            assert math.isclose(coordinates[0, 0], 4.0)
            # sqrt(a)-scaled peculiar -> physical peculiar divides by sqrt(a).
            assert math.isclose(velocities[0, 0], 4.0 / math.sqrt(0.5))
            assert math.isclose(masses[0], 4.0)

        manifest_output = root / "canonical_from_manifest.hdf5"
        manifest_audit = root / "canonical_from_manifest.audit.json"
        manifest_completed = run_converter_from_manifest(
            converter, manifest_path, manifest_output, manifest_audit
        )
        assert manifest_completed.returncode == 0, manifest_completed.stderr
        second_manifest = json.loads(manifest_audit.read_text(encoding="utf-8"))
        assert second_manifest["source_manifest_file"] == str(manifest_path.resolve())
        assert len(second_manifest["source_manifest_sha256"]) == 64
        with h5py.File(output, "r") as expected, h5py.File(manifest_output, "r") as actual:
            for name in ("Coordinates", "Velocities", "Masses", "ParticleIDs"):
                np.testing.assert_array_equal(
                    expected[f"PartType1/{name}"][:], actual[f"PartType1/{name}"][:]
                )
        tampered_manifest = root / "tampered_source.audit.json"
        tampered = json.loads(manifest_text)
        tampered["source_sha256"][0] = "0" * 64
        tampered_manifest.write_text(json.dumps(tampered), encoding="utf-8")
        tampered_result = run_converter_from_manifest(
            converter, tampered_manifest, root / "tampered.hdf5", root / "tampered.audit.json"
        )
        assert tampered_result.returncode != 0
        assert "SHA-256" in tampered_result.stderr

        duplicate_first = root / "duplicate.0.hdf5"
        duplicate_second = root / "duplicate.1.hdf5"
        write_member(duplicate_first, 0, duplicate=True)
        write_member(duplicate_second, 1, duplicate=True)
        failed = run_converter(
            converter,
            duplicate_first,
            root / "duplicate.canonical.hdf5",
            root / "duplicate.audit.json",
        )
        assert failed.returncode != 0
        assert "duplicate particle ID" in failed.stderr
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
