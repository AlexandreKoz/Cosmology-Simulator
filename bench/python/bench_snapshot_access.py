#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import tempfile
import time

import cosmosim


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark Python snapshot access overhead.")
    parser.add_argument("--particles", type=int, default=200000)
    parser.add_argument("--iterations", type=int, default=8)
    args = parser.parse_args()

    config = cosmosim.SimulationConfig()
    state = cosmosim.make_uniform_dark_matter_state(
        particle_count=args.particles,
        mass_code=1.0,
        box_size_mpc_comov=config.box_size_mpc_comov,
    )

    with tempfile.TemporaryDirectory(prefix="cosmosim_py_bench_") as tmp_dir:
        snapshot_path = pathlib.Path(tmp_dir) / "snapshot_bench_000.hdf5"
        cosmosim.write_snapshot(
            snapshot_path,
            state,
            config,
            normalized_config_text="schema_version = 1\n",
            git_sha="bench",
        )
        loaded = cosmosim.read_snapshot(snapshot_path, config)

        setup_start = time.perf_counter()
        mass_view = loaded.state.mass_code()
        _ = loaded.state.position_x_comoving()
        setup_seconds = time.perf_counter() - setup_start

        total = 0.0
        steady_start = time.perf_counter()
        for _ in range(args.iterations):
            total += float(mass_view.sum())
        steady_seconds = time.perf_counter() - steady_start

    element_count = args.particles * args.iterations
    throughput = element_count / steady_seconds if steady_seconds > 0.0 else 0.0
    print(
        "python_snapshot_access "
        f"particles={args.particles} iterations={args.iterations} "
        f"setup_s={setup_seconds:.6f} steady_s={steady_seconds:.6f} "
        f"elements_per_s={throughput:.3e} checksum={total:.6e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
