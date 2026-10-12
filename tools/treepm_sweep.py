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

from chui_snapshot import add_snapshot_arguments, legacy_units, read_snapshot, species_index


def digest(path):
    if path is None:
        return None
    h = hashlib.sha256()
    with Path(path).open('rb') as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def select_target_ids(ids, count, fixed_path=None):
    """Preserve an existing oracle's target order; reject ambiguous identities."""
    import numpy as np
    if fixed_path is None:
        n_targets = min(count, len(ids))
        rows = (np.arange(n_targets, dtype=np.uint64)*len(ids))//n_targets
        return [int(ids[i]) for i in rows]
    selected = []
    seen = set()
    with Path(fixed_path).open() as source:
        for line in source:
            value = line.strip()
            if not value.isascii() or not value.isdecimal():
                raise ValueError('targets: each row must be a positive uint64 persistent ID')
            particle_id = int(value)
            if not 0 < particle_id <= 2**64-1 or particle_id in seen:
                raise ValueError('targets: invalid or duplicate persistent ID')
            index = int(np.searchsorted(ids, np.uint64(particle_id)))
            if index == len(ids) or int(ids[index]) != particle_id:
                raise ValueError('targets: persistent ID is absent from selected species')
            selected.append(particle_id)
            seen.add(particle_id)
    if not selected:
        raise ValueError('targets: fixed target file is empty')
    return selected


def export_snapshot(args):
    import numpy as np
    species_index(args.group)  # Validate even when an explicit Masses lane exists.
    state = read_snapshot(args.snapshot, args.group, args.legacy_dialect,
                          legacy_units(args), need_velocity=False)
    identity = state['provenance']['identity']
    if args.position_scale is None:
        args.position_scale = identity['stored_to_internal_position']
    if args.mass_scale is None:
        args.mass_scale = identity['stored_to_internal_mass']
    ids, positions, mass = state['ids'], state['x'], state['mass']
    targets = select_target_ids(ids, args.targets, getattr(args, 'fixed_target_ids', None))
    # Units remain explicit; conversion is applied in bounded row chunks.
    with args.output.open('w', newline='') as out:
        writer = csv.writer(out, lineterminator='\n')
        writer.writerow(['id', 'x', 'y', 'z', 'mass', 'epsilon'])
        for begin in range(0, len(ids), 65536):
            end = min(begin+65536, len(ids))
            xyz = positions[begin:end]*args.position_scale
            masses = mass[begin:end]*args.mass_scale
            if not np.isfinite(xyz).all() or not np.isfinite(masses).all() or (masses <= 0).any():
                raise ValueError('conversion: nonfinite coordinates or invalid mass')
            for row, i in enumerate(range(begin, end)):
                writer.writerow([int(ids[i]), *[format(x, '.17g') for x in xyz[row]],
                                 format(masses[row], '.17g'), format(args.epsilon, '.17g')])
    args.target_ids.write_text(''.join(str(i)+'\n' for i in targets))
    metadata = dict(logical_snapshot=state['provenance'], particles_sha256=digest(args.output),
                    target_ids_sha256=digest(args.target_ids), group=args.group,
                    position_scale=args.position_scale, mass_scale=args.mass_scale,
                    epsilon_comoving=args.epsilon, scale_factor=state['a'],
                    box_solver_units=(state['box']*args.position_scale).tolist())
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


def validate_total_reference(manifest, reference, particles, target_ids, epoch, box, epsilon, g_code):
    provenance = json.loads(Path(manifest).read_text())
    required = dict(reference_total_sha256=digest(reference), particles_sha256=digest(particles),
                    target_ids_sha256=digest(target_ids), scale_factor=epoch,
                    box_solver_units=box, epsilon_comoving=epsilon, g_code=g_code)
    if not isinstance(provenance, dict) or any(provenance.get(key) != value for key, value in required.items()):
        raise ValueError('reference: incompatible frozen source/epoch/units/softening or reference checksum')
    if (provenance.get('reference_kind') not in ('independent_periodic_ewald', 'converged_total') or
            not isinstance(provenance.get('qualification_evidence'), str) or
            not provenance['qualification_evidence'].strip()):
        raise ValueError('reference: independent/converged total qualification evidence must be declared')
    return provenance


def sweep(args):
    args.output.mkdir(parents=True, exist_ok=False)
    executable = args.exe.resolve()
    frozen_hash = digest(args.particles)
    targets_hash = digest(args.target_ids)
    source_metadata = None
    if args.particles:
        sidecar = args.particles.with_suffix('.metadata.json')
        if sidecar.is_file():
            source_metadata = json.loads(sidecar.read_text())
            if source_metadata.get('particles_sha256') != frozen_hash or source_metadata.get('target_ids_sha256') != targets_hash:
                raise ValueError('provenance: frozen CSV/target IDs differ from snapshot export')
            if source_metadata.get('scale_factor') != args.scale_factor:
                raise ValueError('provenance: benchmark scale factor differs from frozen snapshot')
            expected_box = source_metadata.get('box_solver_units')
            if expected_box != [args.box_x, args.box_y, args.box_z]:
                raise ValueError('provenance: benchmark box differs from frozen snapshot conversion')
            if source_metadata.get('epsilon_comoving') != args.epsilon:
                raise ValueError('provenance: benchmark softening differs from frozen snapshot export')
    total_reference_provenance = None
    if args.total_reference:
        if not args.particles or not args.target_ids or not args.total_reference_manifest:
            raise ValueError('reference: external total comparison requires frozen particles, target IDs and --total-reference-manifest')
        total_reference_provenance = validate_total_reference(args.total_reference_manifest,
            args.total_reference, args.particles, args.target_ids, args.scale_factor,
            [args.box_x, args.box_y, args.box_z], args.epsilon, args.g)
    elif args.total_reference_manifest:
        raise ValueError('reference: --total-reference-manifest requires --total-reference')
    binary_hash = digest(executable)
    axis_values = [comma_values(args.meshes, int), comma_values(args.leaves, int),
                   comma_values(args.blocks, int), comma_values(args.threads, int),
                   comma_values(args.policies, str), comma_values(args.kernels, str),
                   comma_values(args.asmths, float) if args.asmths else [args.asmth],
                   comma_values(args.rcuts, float) if args.rcuts else [args.rcut]]
    candidate_count = math.prod(len(values) for values in axis_values)
    if not 1 <= candidate_count <= args.maximum_candidates:
        raise ValueError(f'{candidate_count} candidates exceeds staged budget {args.maximum_candidates}')
    axes = itertools.product(*axis_values)
    if not math.isfinite(args.maximum_angle) or not 0 < args.maximum_angle <= 0.5:
        raise ValueError('maximum angle must be in (0,0.5]')
    rows = []
    for mesh, leaf, block, threads, policy, kernel, asmth, rcut in axes:
        if min(mesh, leaf, block, threads) <= 0 or block > 4096:
            raise ValueError('positive mesh/leaf/block/thread axes required; block <=4096')
        if not all(math.isfinite(x) and x > 0 for x in (asmth, rcut)):
            raise ValueError('split parameters must be finite positive')
        name = f'm{mesh}_l{leaf}_b{block}_t{threads}_{policy}_{kernel}_a{asmth}_c{rcut}'
        forces = args.output / f'{name}.forces.csv'
        command = [str(executable), '--mesh', str(mesh), '--leaf', str(leaf), '--block', str(block),
                   '--policy', policy, '--kernel', kernel, '--history', args.history,
                   '--maximum-angle', str(args.maximum_angle), '--reuse', args.reuse, '--refit', args.refit,
                   '--accounting', args.accounting, '--repeats', str(args.repeats),
                   '--warmups', str(args.warmups), '--count', str(args.count), '--targets', str(args.targets),
                   '--box-x', str(args.box_x), '--box-y', str(args.box_y), '--box-z', str(args.box_z),
                   '--asmth', str(asmth), '--rcut', str(rcut), '--epsilon', str(args.epsilon),
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
                        source_snapshot_export=source_metadata,
                        total_reference_provenance=total_reference_provenance,
                        total_reference_manifest_sha256=digest(args.total_reference_manifest),
                        source_provenance_status='snapshot_export_verified' if source_metadata else 'explicit_csv_or_generated_only',
                        returncode=result.returncode)
        (args.output / f'{name}.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        row = dict(candidate=name, mesh=mesh, leaf=leaf, block=block, threads=threads,
                   policy=policy, kernel=kernel, asmth_cells=asmth, rcut_cells=rcut, status='failed' if result.returncode else 'measured')
        if result.returncode == 0:
            records = [json.loads(line) for line in result.stdout.splitlines() if line.startswith('{')]
            timings = [r for r in records if r['record'] == 'timing']
            if len(timings) != args.repeats:
                raise ValueError(f'{name}: incomplete timing records')
            for field in ('wall_ms', 'pm_ms', 'tree_build_ms', 'traversal_ms', 'pairs', 'nodes', 'opened', 'multipoles',
                          'source_preprocess_ms', 'refit_ms', 'full_rebuild_ms', 'geometry_preparation_ms',
                          'residual_lookup_targets', 'simd_leaf_targets', 'skipped_geometry'):
                row[field] = statistics.median(r[field] for r in timings)
            row['peak_rss_bytes'] = max((r['peak_rss_bytes'] or 0) for r in timings) or None
            row['solver_retained_bytes'] = max(r['solver_retained_bytes'] for r in timings)
            row['observed_workers'] = max(r['workers'] for r in timings)
            config = next(r for r in records if r['record'] == 'configuration')
            row.update(split=config['split'], cutoff=config['cutoff'], maximum_angle=config['maximum_angle'],
                       softening=config['softening'], build=config['build'], cpu=config['cpu'],
                       reuse=config['reuse'], refit=config['refit'])
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
                         b['threads']==row['threads'] and b['policy']==row['policy'] and
                         b['asmth_cells']==row['asmth_cells'] and b['rcut_cells']==row['rcut_cells']), None)
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
    export = commands.add_parser('export', help='export a complete logical snapshot particle group in explicit solver units')
    export.add_argument('snapshot', type=Path)
    export.add_argument('--group', default='PartType1')
    export.add_argument('--output', type=Path, required=True)
    export.add_argument('--target-ids', type=Path, required=True)
    export.add_argument('--fixed-target-ids', type=Path, help='existing persistent-ID list; preserve its order instead of selecting new targets')
    export.add_argument('--epsilon', type=float, required=True)
    export.add_argument('--position-scale', type=float, help='stored-to-solver length factor; default: declared dialect to internal units')
    export.add_argument('--mass-scale', type=float, help='stored-to-solver mass factor; default: declared dialect to internal units')
    export.add_argument('--targets', type=int, default=4096)
    add_snapshot_arguments(export)
    export.set_defaults(action=export_snapshot)
    run = commands.add_parser('run', help='run the compiled serial-owner benchmark repeatedly')
    run.add_argument('--exe', type=Path, required=True)
    run.add_argument('--output', type=Path, required=True)
    run.add_argument('--particles', type=Path)
    run.add_argument('--target-ids', type=Path)
    run.add_argument('--meshes', default='64')
    run.add_argument('--leaves', default='16')
    run.add_argument('--blocks', default='64')
    run.add_argument('--threads', default='4')
    run.add_argument('--policies', default='strict')
    run.add_argument('--kernels', default='analytic')
    run.add_argument('--asmths', help='staged comma-separated smoothing-cell candidates')
    run.add_argument('--rcuts', help='staged comma-separated cutoff-cell candidates')
    run.add_argument('--maximum-angle', type=float, default=0.25)
    run.add_argument('--reuse', choices=['off','on'], default='off')
    run.add_argument('--refit', choices=['off','on'], default='off')
    run.add_argument('--maximum-candidates', type=int, default=40, help='fail before launch if a staged experiment exceeds this bound')
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
    run.add_argument('--total-reference-manifest', type=Path, help='qualified total-reference identity/evidence JSON sidecar')
    run.set_defaults(action=sweep)
    args = parser.parse_args()
    if args.targets <= 0 or not math.isfinite(args.epsilon) or args.epsilon < 0:
        parser.error('positive target count and finite nonnegative softening required')
    if args.command == 'export':
        if not all(x is None or (math.isfinite(x) and x > 0) for x in (args.position_scale, args.mass_scale)):
            parser.error('finite positive position/mass conversion scales required')
    elif args.count <= 0 or args.repeats <= 0 or args.warmups < 0 or not 0 <= args.ewald_level <= 32:
        parser.error('positive count/repeats, nonnegative warmups and Ewald level 0..32 required')
    args.action(args)


if __name__ == '__main__':
    main()
