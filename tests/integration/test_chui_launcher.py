#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tempfile
import time


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def make_executable(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_presets(root: Path, names: list[str]) -> None:
    data = {
        "version": 6,
        "configurePresets": [
            {"name": "base", "hidden": True, "binaryDir": "${sourceDir}/build/${presetName}"},
            *[{"name": name, "inherits": "base"} for name in names],
        ],
        "buildPresets": [
            {"name": f"build-{name}", "configurePreset": name}
            for name in names
        ],
    }
    (root / "CMakePresets.json").write_text(json.dumps(data), encoding="utf-8")


def run(*args: str, cwd: Path, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(args), cwd=cwd, env=env, text=True, capture_output=True, check=False
    )


def main() -> int:
    source_launcher = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="chui_launcher_test_") as temp:
        root = Path(temp)
        launcher = root / "chui"
        shutil.copy2(source_launcher, launcher)
        launcher.chmod(launcher.stat().st_mode | stat.S_IXUSR)
        write_presets(root, ["alpha"])

        help_result = run(str(launcher), "--help", cwd=root)
        require(help_result.returncode == 0, "./chui --help should succeed")
        require("{run}" in help_result.stdout, "top-level help should expose run")
        run_help = run(str(launcher), "run", "--help", cwd=root)
        require(run_help.returncode == 0, "./chui run --help should succeed")
        require("--status-every" in run_help.stdout, "run help should document status controls")

        fake = root / "fake harness"
        make_executable(fake, "#!/bin/sh\nprintf '%s\\n' \"$@\"\nexit \"${FAKE_EXIT:-0}\"\n")
        config = Path("configs/config with spaces.param.txt")
        (root / config.parent).mkdir(parents=True, exist_ok=True)
        (root / config).write_text("schema_version = 1\n", encoding="utf-8")
        explicit = run(
            str(launcher), "run", str(config), "--exe", str(fake),
            "--quiet", "--status-every", "5", "--status-seconds", "2.5",
            cwd=root,
        )
        require(explicit.returncode == 0, "explicit --exe launch should propagate success")
        forwarded = explicit.stdout.splitlines()
        require(forwarded[0] == str(config), "config path must be forwarded without relocation")
        require("--quiet" in forwarded, "quiet flag must be forwarded to harness")
        require("--status-every" in forwarded and "5" in forwarded,
                "step cadence must be forwarded to harness")
        require("--status-seconds" in forwarded and "2.5" in forwarded,
                "wall-clock cadence must be forwarded to harness")
        require("[CHUI][LAUNCH]" not in explicit.stdout,
                "quiet mode must suppress launcher provenance")

        absolute = run(
            str(launcher), "run", str((root / config).resolve()), "--exe", str(fake), "--quiet",
            cwd=root,
        )
        require(absolute.returncode == 0, "absolute config path should be accepted")
        require(absolute.stdout.splitlines()[0] == str((root / config).resolve()),
                "absolute config path must be forwarded unchanged")

        exit_env = os.environ.copy()
        exit_env["FAKE_EXIT"] = "7"
        exit_result = run(str(launcher), "run", str(config), "--exe", str(fake), cwd=root, env=exit_env)
        require(exit_result.returncode == 7, "launcher must preserve child exit status through exec")

        signal_fake = root / "signal harness"
        make_executable(
            signal_fake,
            "#!/usr/bin/env python3\nimport time\ntime.sleep(30)\n",
        )
        signal_proc = subprocess.Popen(
            [str(launcher), "run", str(config), "--exe", str(signal_fake), "--quiet"],
            cwd=root,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        time.sleep(0.2)
        signal_proc.terminate()
        signal_returncode = signal_proc.wait(timeout=5)
        require(signal_returncode == -15,
                "process replacement should deliver SIGTERM directly to the launched runtime")

        preset_exe = root / "build" / "alpha" / "cosmosim_harness"
        make_executable(preset_exe, "#!/bin/sh\nprintf 'preset:%s\\n' \"$1\"\n")
        preset_result = run(str(launcher), "run", str(config), "--preset", "alpha", cwd=root)
        require(preset_result.returncode == 0, "explicit --preset should select deterministic binary")
        require(f"preset:{config}" in preset_result.stdout, "preset executable should receive config")
        require("[CHUI][LAUNCH]" in preset_result.stdout,
                "non-quiet launcher should report launch provenance")
        require(str(preset_exe) in preset_result.stdout and "selection=preset:alpha" in preset_result.stdout,
                "launch provenance should identify the selected executable and preset")

        auto_result = run(str(launcher), "run", str(config), cwd=root)
        require(auto_result.returncode == 0, "one auto-discovered executable should launch")
        require("selection=auto_discovery" in auto_result.stdout,
                "automatic selection should be visible in launch provenance")

        fake_mpi = root / "fake_mpiexec"
        make_executable(fake_mpi, "#!/bin/sh\nprintf '%s\\n' \"$@\"\n")
        cache_path = root / "build" / "alpha" / "CMakeCache.txt"
        cache_path.write_text(
            f"COSMOSIM_ENABLE_MPI:BOOL=OFF\nMPIEXEC_EXECUTABLE:FILEPATH={fake_mpi}\nMPIEXEC_NUMPROC_FLAG:STRING=-n\n",
            encoding="utf-8",
        )
        rejected_serial_mpi = run(
            str(launcher), "run", str(config), "--preset", "alpha", "--mpi", "2", cwd=root
        )
        require(rejected_serial_mpi.returncode == 2,
                "--mpi must reject a build that explicitly has MPI disabled")
        require("COSMOSIM_ENABLE_MPI=OFF" in rejected_serial_mpi.stderr,
                "serial-build MPI rejection should explain the cache mismatch")

        cache_path.write_text(
            f"COSMOSIM_ENABLE_MPI:BOOL=ON\nMPIEXEC_EXECUTABLE:FILEPATH={fake_mpi}\nMPIEXEC_NUMPROC_FLAG:STRING=-n\n",
            encoding="utf-8",
        )
        mpi_result = run(
            str(launcher), "run", str(config), "--preset", "alpha", "--mpi", "2",
            "--quiet", "--", "--bind-to", "core", cwd=root,
        )
        require(mpi_result.returncode == 0, "MPI command composition should execute configured launcher")
        mpi_args = mpi_result.stdout.splitlines()
        require(mpi_args[:4] == ["--bind-to", "core", "-n", "2"],
                "MPI passthrough and rank count should precede executable")
        require(str(preset_exe) in mpi_args, "MPI command should contain selected harness")
        require(str(config) in mpi_args, "MPI command should preserve config path")

        preset_exe.unlink()
        missing = run(str(launcher), "run", str(config), "--preset", "alpha", cwd=root)
        require(missing.returncode == 2, "missing preset executable should fail")
        require("cmake --preset alpha" in missing.stderr and "cmake --build --preset build-alpha" in missing.stderr,
                "missing build diagnostic should be actionable")
        no_auto = run(str(launcher), "run", str(config), cwd=root)
        require(no_auto.returncode == 2, "zero auto-discovered executables should fail")
        require("No built cosmosim_harness" in no_auto.stderr,
                "zero-candidate auto discovery should provide a build diagnostic")

        write_presets(root, ["alpha", "beta"])
        make_executable(root / "build" / "alpha" / "cosmosim_harness", "#!/bin/sh\nexit 0\n")
        make_executable(root / "build" / "beta" / "cosmosim_harness", "#!/bin/sh\nexit 0\n")
        ambiguous = run(str(launcher), "run", str(config), cwd=root)
        require(ambiguous.returncode == 2, "ambiguous automatic discovery should fail")
        require("refusing to guess" in ambiguous.stderr,
                "ambiguous automatic discovery should explain deterministic selection rule")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
