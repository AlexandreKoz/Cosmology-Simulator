#!/usr/bin/env python3
"""Synthetic clustered DMO particle box for TreePM performance diagnostics.

Writes the `id,x,y,z,mass,epsilon` CSV consumed by bench_tree_pm_sweep. The
distribution is a stand-in for an evolved low-redshift box: NFW-like halos
drawn from a power-law multiplicity function plus a uniform background. It is
a workload fixture for interaction-count and timing comparisons, not a
cosmological initial condition.

The `--order` option controls only the memory order of the rows: `random`
shuffles them (worst case for cache locality), `morton` sorts them along a
Morton curve (what a space-filling-curve reorder produces).
"""
import argparse

import numpy as np


def morton_order(positions, box, bits=10):
    cells = np.floor(positions / box * (1 << bits)).astype(np.uint64)
    cells = np.minimum(cells, np.uint64((1 << bits) - 1))

    def spread(values):
        out = np.zeros_like(values)
        for bit in range(bits):
            out |= ((values >> np.uint64(bit)) & np.uint64(1)) << np.uint64(3 * bit)
        return out

    keys = spread(cells[:, 0]) | (spread(cells[:, 1]) << np.uint64(1)) | (spread(cells[:, 2]) << np.uint64(2))
    return np.argsort(keys, kind="stable")


def generate(count, box, halo_fraction, halo_count, seed):
    rng = np.random.default_rng(seed)
    halo_particles = int(count * halo_fraction)
    weights = rng.pareto(0.9, halo_count) + 1.0
    weights /= weights.sum()
    members = rng.multinomial(halo_particles, weights)
    concentration = 8.0
    grid = np.linspace(1.0e-4, concentration, 4000)
    cdf = np.log1p(grid) - grid / (1.0 + grid)
    cdf /= cdf[-1]
    chunks = []
    for n in members:
        if n == 0:
            continue
        centre = rng.uniform(0.0, box, 3)
        virial_radius = 0.035 * box * (n / 2000.0) ** (1.0 / 3.0)
        radius = (virial_radius / concentration) * np.interp(rng.uniform(0.0, 1.0, n), cdf, grid)
        direction = rng.normal(size=(n, 3))
        direction /= np.linalg.norm(direction, axis=1)[:, None]
        chunks.append(centre + radius[:, None] * direction)
    chunks.append(rng.uniform(0.0, box, (count - halo_particles, 3)))
    return np.mod(np.vstack(chunks), box), rng


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output", help="CSV path to write")
    parser.add_argument("--count", type=int, default=64 ** 3)
    parser.add_argument("--box", type=float, default=10.0, help="periodic box length (code length)")
    parser.add_argument("--epsilon", type=float, default=0.005, help="uniform Plummer softening")
    parser.add_argument("--halo-fraction", type=float, default=0.55)
    parser.add_argument("--halo-count", type=int, default=600)
    parser.add_argument("--seed", type=int, default=1698)
    parser.add_argument("--order", choices=("random", "morton"), default="random")
    args = parser.parse_args()
    if args.count <= 0 or args.box <= 0.0 or args.epsilon < 0.0 or not 0.0 <= args.halo_fraction <= 1.0:
        parser.error("count and box must be positive, epsilon nonnegative, halo fraction in [0,1]")

    positions, rng = generate(args.count, args.box, args.halo_fraction, args.halo_count, args.seed)
    if args.order == "random":
        positions = positions[rng.permutation(len(positions))]
    else:
        positions = positions[morton_order(positions, args.box)]
    mass = 1.0 / len(positions)
    with open(args.output, "w") as out:
        out.write("id,x,y,z,mass,epsilon\n")
        for index, (x, y, z) in enumerate(positions):
            out.write(f"{index + 1},{x:.10f},{y:.10f},{z:.10f},{mass:.10e},{args.epsilon}\n")


if __name__ == "__main__":
    main()
