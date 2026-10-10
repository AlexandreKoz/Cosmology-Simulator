#!/usr/bin/env python3
"""Frozen-source TreePM diagnostics; never parses production .param.txt files."""
import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import statistics
import subprocess


def digest(path):
    if path is None:
        return None
    h = hashlib.sha256()
    with Path(path).open('rb') as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def export_snapshot(args):
    # Optional offline dependencies; no dependency installation is performed.
    import h5py
    import numpy as np
    with h5py.File(args.snapshot, 'r') as source:
        group = source[args.group]
        ids = np.asarray(group['ParticleIDs'], dtype=np.uint64)
        positions = np.asarray(group['Coordinates'], dtype=np.float64) * args.position_scale
        if 'Masses' in group:
            mass = np.asarray(group['Masses'], dtype=np.float64) * args.mass_scale
        else:
            species = int(args.group.removeprefix('PartType'))
            mass = np.full(len(ids), source['Header'].attrs['MassTable'][species] * args.mass_scale)
        order = np.argsort(ids, kind='stable')
        if len(ids) == 0 or len(np.unique(ids)) != len(ids):
            raise ValueError('snapshot needs nonempty unique particle IDs')
        if positions.shape != (len(ids), 3) or mass.shape != (len(ids),):
            raise ValueError('snapshot source dimensions disagree')
        if not np.isfinite(positions).all() or not np.isfinite(mass).all() or (mass <= 0).any():
            raise ValueError('snapshot has invalid position or mass')
        with args.output.open('w', newline='') as out:
            writer = csv.writer(out, lineterminator='\n')
            writer.writerow(['id', 'x', 'y', 'z', 'mass', 'epsilon'])
            # epsilon is an explicit physical convention supplied by the user,
            # in these comoving solver length units; it never scales with mesh.
            for row in order:
                writer.writerow([int(ids[row]), *[format(x, '.17g') for x in positions[row]],
                                 format(mass[row], '.17g'), format(args.epsilon, '.17g')])
        n_targets = min(args.targets, len(ids))
        target_rows = (np.arange(n_targets, dtype=np.uint64) * len(ids)) // n_targets
        target_ids = ids[order[target_rows]]
        args.target_ids.write_text(''.join(f'{int(i)}\n' for i in target_ids))
        # Preserve header information without inferring code units or G from it.
        header = {k: np.asarray(v).tolist() for k, v in source['Header'].attrs.items()
                  if k in ('Time', 'Redshift', 'BoxSize', 'Omega0', 'OmegaLambda', 'HubbleParam')}
    metadata = dict(snapshot=str(args.snapshot), snapshot_sha256=digest(args.snapshot),
                    particles_sha256=digest(args.output), target_ids_sha256=digest(args.target_ids),
                    group=args.group, position_scale=args.position_scale, mass_scale=args.mass_scale,
                    epsilon_comoving=args.epsilon, header=header)
    args.output.with_suffix('.metadata.json').write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n')


def comma_values(text, converter):
    values = [converter(value) for value in text.split(',')]
    if not values:
        raise ValueError('empty sweep axis')
    return values


def force_metrics(candidate, reference):
    def read(path):
        with Path(path).open(newline='') as source:
            rows = list(csv.DictReader(source))
        result = {int(row['id']): tuple(float(row[f'total_{axis}']) for axis in 'xyz') for row in rows}
        if len(result) != len(rows) or not result:
            raise ValueError('force reference has duplicate IDs or is empty')
        if not all(math.isfinite(x) for value in result.values() for x in value):
            raise ValueError('nonfinite total-force reference')
        return result
    a, b = read(candidate), read(reference)
    if a.keys() != b.keys():
        raise ValueError('total-force reference target IDs differ')
    errors = [math.dist(a[i], b[i]) for i in a]
    norm2 = sum(sum(x*x for x in b[i]) for i in b)
    floor = 1e-3 * math.sqrt(norm2 / len(b))
    normalized = sorted(e / max(math.hypot(*b[i]), floor, 1e-300) for i, e in zip(a, errors))
    absolute = sorted(errors)
    small_force_errors = [e for i, e in zip(a, errors) if math.hypot(*b[i]) <= floor]
    return dict(relative_l2=math.sqrt(sum(e*e for e in errors)/max(norm2, 1e-300)),
                absolute_max=max(errors), absolute_p95=absolute[int(.95*(len(a)-1))],
                absolute_p99=absolute[int(.99*(len(a)-1))],
                small_force_absolute_max=max(small_force_errors, default=0.0),
                p95=normalized[int(.95*(len(a)-1))],
                p99=normalized[int(.99*(len(a)-1))], maximum=normalized[-1], normalization_floor=floor)


def sweep(args):
    args.output.mkdir(parents=True, exist_ok=False)
    executable = args.exe.resolve()
    frozen_hash = digest(args.particles)
    targets_hash = digest(args.target_ids)
    binary_hash = digest(executable)
    axes = itertools.product(comma_values(args.meshes, int), comma_values(args.leaves, int),
                             comma_values(args.blocks, int), comma_values(args.threads, int),
                             comma_values(args.policies, str), comma_values(args.kernels, str))
    rows = []
    for mesh, leaf, block, threads, policy, kernel in axes:
        if min(mesh, leaf, block, threads) <= 0 or block > 4096:
            raise ValueError('positive mesh/leaf/block/thread axes required; block <=4096')
        name = f'm{mesh}_l{leaf}_b{block}_t{threads}_{policy}_{kernel}'
        forces = args.output / f'{name}.forces.csv'
        command = [str(executable), '--mesh', str(mesh), '--leaf', str(leaf), '--block', str(block),
                   '--policy', policy, '--kernel', kernel, '--history', args.history,
                   '--accounting', args.accounting, '--repeats', str(args.repeats),
                   '--warmups', str(args.warmups), '--count', str(args.count), '--targets', str(args.targets),
                   '--box-x', str(args.box_x), '--box-y', str(args.box_y), '--box-z', str(args.box_z),
                   '--asmth', str(args.asmth), '--rcut', str(args.rcut), '--epsilon', str(args.epsilon),
                   '--g', str(args.g), '--scale-factor', str(args.scale_factor),
                   '--ewald-level', str(args.ewald_level)]
        if args.particles:
            command += ['--particles', str(args.particles.resolve())]
        if args.target_ids:
            command += ['--target-ids', str(args.target_ids.resolve())]
        if args.oracle or args.total_reference:
            command += ['--forces', str(forces.resolve())]
        env = os.environ.copy()
        env.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE')
        result = subprocess.run(command, env=env, capture_output=True, text=True, check=False)
        (args.output / f'{name}.jsonl').write_text(result.stdout)
        (args.output / f'{name}.stderr.txt').write_text(result.stderr)
        manifest = dict(command=command, OMP_NUM_THREADS=threads, OMP_DYNAMIC='FALSE',
                        executable_sha256=binary_hash, particles_sha256=frozen_hash,
                        target_ids_sha256=targets_hash, reference_total_sha256=digest(args.total_reference),
                        returncode=result.returncode)
        (args.output / f'{name}.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        row = dict(candidate=name, mesh=mesh, leaf=leaf, block=block, threads=threads,
                   policy=policy, kernel=kernel, status='failed' if result.returncode else 'measured')
        if result.returncode == 0:
            records = [json.loads(line) for line in result.stdout.splitlines() if line.startswith('{')]
            timings = [r for r in records if r['record'] == 'timing']
            if len(timings) != args.repeats:
                raise ValueError(f'{name}: incomplete timing records')
            for field in ('wall_ms', 'pm_ms', 'tree_build_ms', 'traversal_ms', 'pairs', 'nodes', 'opened', 'multipoles'):
                row[field] = statistics.median(r[field] for r in timings)
            row['peak_rss_bytes'] = max((r['peak_rss_bytes'] or 0) for r in timings) or None
            row['solver_retained_bytes'] = max(r['solver_retained_bytes'] for r in timings)
            row['observed_workers'] = max(r['workers'] for r in timings)
            config = next(r for r in records if r['record'] == 'configuration')
            row.update(split=config['split'], cutoff=config['cutoff'])
            for record in records:
                if record['record'] == 'accuracy':
                    prefix = 'short' if record['kind'].startswith('matching_split') else 'ewald'
                    for key in ('relative_l2', 'absolute_max', 'absolute_p95', 'absolute_p99',
                                'small_force_absolute_max', 'normalization_floor', 'p95', 'p99', 'maximum'):
                        row[f'{prefix}_{key}'] = record[key]
            if args.total_reference:
                metrics = force_metrics(forces, args.total_reference)
                for key, value in metrics.items():
                    row[f'external_total_{key}'] = value
        rows.append(row)
        print(f'{name}: {row["status"]}', flush=True)
    # Compare wall times only at the same thread count, acceptance and history.
    for row in rows:
        baseline = next((b for b in rows if b['status']=='measured' and
                         b['mesh']==args.baseline_mesh and b['leaf']==args.baseline_leaf and
                         b['block']==args.baseline_block and b['kernel']=='analytic' and
                         b['threads']==row['threads'] and b['policy']==row['policy']), None)
        if baseline and row['status']=='measured':
            row['baseline'] = baseline['candidate']
            row['wall_speedup'] = baseline['wall_ms']/row['wall_ms']
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with (args.output / 'report.csv').open('w', newline='') as report:
        writer = csv.DictWriter(report, fieldnames=keys);writer.writeheader();writer.writerows(rows)
    lines = ['# TreePM frozen-source measurements', '',
             'Timing alone does not qualify a split. Each short reference uses the candidate split.',
             'Ewald and external total references require independent convergence/identity qualification.',
             'Source hashes and exact commands are in each manifest. PM is refreshed every measured iteration.', '',
             '| Candidate | wall ms | speedup | short L2 | total Ewald L2 | peak RSS bytes |',
             '|---|---:|---:|---:|---:|---:|']
    for r in rows:
        lines.append('| '+' | '.join(str(r.get(k, 'unmeasured')) for k in
                                     ('candidate','wall_ms','wall_speedup','short_relative_l2','ewald_relative_l2','peak_rss_bytes'))+' |')
    (args.output / 'report.md').write_text('\n'.join(lines)+'\n')
    if any(r['status']=='failed' for r in rows):
        raise SystemExit('One or more candidates failed; see saved stderr and manifests.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    export = commands.add_parser('export', help='export one frozen HDF5 particle group in explicit solver units')
    export.add_argument('snapshot', type=Path)
    export.add_argument('--group', default='PartType1')
    export.add_argument('--output', type=Path, required=True)
    export.add_argument('--target-ids', type=Path, required=True)
    export.add_argument('--epsilon', type=float, required=True)
    export.add_argument('--position-scale', type=float, default=1.0)
    export.add_argument('--mass-scale', type=float, default=1.0)
    export.add_argument('--targets', type=int, default=4096)
    export.set_defaults(action=export_snapshot)
    run = commands.add_parser('run', help='run the compiled serial-owner benchmark repeatedly')
    run.add_argument('--exe', type=Path, required=True)
    run.add_argument('--output', type=Path, required=True)
    run.add_argument('--particles', type=Path)
    run.add_argument('--target-ids', type=Path)
    run.add_argument('--meshes', default='64,96,128,192')
    run.add_argument('--leaves', default='4,8,16,32')
    run.add_argument('--blocks', default='64')
    run.add_argument('--threads', default='4')
    run.add_argument('--policies', default='strict,adaptive')
    run.add_argument('--kernels', default='analytic,lookup')
    run.add_argument('--history', choices=['fallback', 'bootstrap'], default='fallback')
    run.add_argument('--accounting', choices=['full', 'fast'], default='full')
    for key, default in [('count',4096),('targets',4096),('repeats',5),('warmups',1),
                         ('baseline-mesh',64),('baseline-leaf',16),('baseline-block',64),('ewald-level',0)]:
        run.add_argument('--'+key, type=int, default=default)
    for key, default in [('box-x',1.0),('box-y',1.0),('box-z',1.0),('asmth',1.25),('rcut',6.25),
                         ('epsilon',0.00035),('g',1.0),('scale-factor',1.0)]:
        run.add_argument('--'+key, type=float, default=default)
    run.add_argument('--oracle', action='store_true', help='generate a matching direct short reference for each split')
    run.add_argument('--total-reference', type=Path, help='independently qualified total CSV, id,total_x,total_y,total_z')
    run.set_defaults(action=sweep)
    args = parser.parse_args()
    if args.targets <= 0 or not math.isfinite(args.epsilon) or args.epsilon < 0:
        parser.error('positive target count and finite nonnegative softening required')
    if args.command == 'export':
        if not all(math.isfinite(x) and x > 0 for x in (args.position_scale, args.mass_scale)):
            parser.error('finite positive position/mass conversion scales required')
    elif args.count <= 0 or args.repeats <= 0 or args.warmups < 0 or not 0 <= args.ewald_level <= 32:
        parser.error('positive count/repeats, nonnegative warmups and Ewald level 0..32 required')
    args.action(args)


if __name__ == '__main__':
    main()
