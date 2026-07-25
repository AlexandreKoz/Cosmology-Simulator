#!/usr/bin/env python3
"""Compatibility launcher for the authoritative C++ CHUÍ IC converter."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence


def locate_converter() -> str:
    override = os.environ.get("COSMOSIM_CONVERT_IC_EXECUTABLE")
    if override:
        return override
    found = shutil.which("cosmosim_convert_ic")
    if found:
        return found
    root = Path(__file__).resolve().parents[1]
    candidates = [
        root / "build" / "hdf5-debug" / "cosmosim_convert_ic",
        root / "build" / "mpi-hdf5-fftw-debug" / "cosmosim_convert_ic",
    ]
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    raise FileNotFoundError(
        "cosmosim_convert_ic executable not found; build the HDF5 preset or "
        "set COSMOSIM_CONVERT_IC_EXECUTABLE"
    )


def print_launcher_help() -> None:
    print(
        "CHUÍ canonical IC converter launcher\n\n"
        "Build the authoritative converter with:\n"
        "  cmake --preset hdf5-debug\n"
        "  cmake --build --preset build-hdf5-debug --target cosmosim_convert_ic\n\n"
        "Then run:\n"
        "  tools/convert_ic.py [cosmosim_convert_ic options]\n\n"
        "Set COSMOSIM_CONVERT_IC_EXECUTABLE to use an executable from another "
        "build tree. Scientific conversion is implemented only in the linked "
        "C++ converter and runtime I/O library."
    )


def main(argv: Sequence[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    try:
        executable = locate_converter()
    except FileNotFoundError as error:
        if any(argument in {"-h", "--help"} for argument in arguments):
            print_launcher_help()
            return 0
        print(str(error), file=sys.stderr)
        return 127
    return subprocess.call([executable, *arguments])


if __name__ == "__main__":
    raise SystemExit(main())
