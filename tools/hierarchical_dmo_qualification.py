#!/usr/bin/env python3
"""Offline assertions for matched global/hierarchical DMO snapshot convergence.

The requested k interval is an explicit qualification target, never a claim of
an already validated scale range. Dependencies: numpy and h5py, supplied by user.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    with h5py.File(path, 'r') as source:
        group = source['PartType1']
        ids = np.asarray(group['ParticleIDs'], dtype=np.uint64)
        positions = np.asarray(group['Coordinates'], dtype=np.float64)
        velocities = np.asarray(group['Velocities'], dtype=np.float64)
        header = source['Header'].attrs
        box = np.asarray(header['BoxSize'], dtype=np.float64)
        if box.ndim == 0:
            axis_keys = [f'CHUIBoxSize{axis}_MpcComoving' for axis in 'XYZ']
            if all(key in header for key in axis_keys):
                lengths = np.array([header[key] for key in axis_keys], dtype=np.float64)
                box = float(box) * lengths / lengths[0]
            else:
                box = np.repeat(box, 3)
        if box.shape != (3,) or not np.isfinite(box).all() or (box <= 0).any():
            raise ValueError('requires positive rectangular box lengths in Header/BoxSize')
        mass = (np.asarray(group['Masses'], dtype=np.float64) if 'Masses' in group else
                np.full(len(ids), float(header['MassTable'][1])))
        a = float(header['Time'])
        if not np.isfinite(a) or a <= 0:
            raise ValueError('requires a finite positive cosmological scale factor')
    if len(ids) == 0 or len(np.unique(ids)) != len(ids):
        raise ValueError('requires unique nonempty DMO particle IDs')
    if positions.shape != (len(ids), 3) or velocities.shape != positions.shape or mass.shape != (len(ids),):
        raise ValueError('inconsistent snapshot particle lanes')
    if not all(np.isfinite(x).all() for x in (positions, velocities, mass)) or (mass <= 0).any():
        raise ValueError('nonfinite state or nonpositive mass')
    order = np.argsort(ids, kind='stable')
    return dict(ids=ids[order], x=positions[order], u=velocities[order], mass=mass[order], box=box, a=a)


def power(state, mesh, k_min, k_max, bins):
    # Identical offline CIC estimator for every run. It is deliberately separate
    # from the PM force mesh, with no deconvolution/no shot-noise subtraction.
    x = (state['x'] % state['box']) / state['box'] * mesh
    base = np.floor(x).astype(np.int64)
    fraction = x-base
    density = np.zeros(mesh**3, dtype=np.float64)
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                offset = np.array([i, j, k])
                cell = (base+offset) % mesh
                weight = np.prod(np.where(offset, fraction, 1-fraction), axis=1) * state['mass']
                index = (cell[:, 0]*mesh+cell[:, 1])*mesh+cell[:, 2]
                density += np.bincount(index, weights=weight, minlength=mesh**3)
    density = density.reshape((mesh, mesh, mesh))
    delta = density/density.mean()-1.0
    transform = np.fft.fftn(delta)/mesh**3
    axes = [2*np.pi*np.fft.fftfreq(mesh, d=length/mesh) for length in state['box']]
    kval = np.sqrt(axes[0][:, None, None]**2+axes[1][None, :, None]**2+axes[2][None, None, :]**2)
    edges = np.linspace(k_min, k_max, bins+1)
    shell = np.searchsorted(edges, kval.ravel(), side='right')-1
    valid = (shell >= 0) & (shell < bins) & (kval.ravel() > 0)
    counts = np.bincount(shell[valid], minlength=bins)
    p = np.bincount(shell[valid], weights=np.abs(transform.ravel()[valid])**2, minlength=bins)
    if (counts < 8).any():
        raise ValueError('each requested k shell needs >=8 modes; reduce bins or widen the target range')
    return p/counts, counts, (edges[:-1]+edges[1:])/2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('initial', 'global-reference', 'hierarchical', 'refined'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--mesh', type=int, default=32)
    parser.add_argument('--bins', type=int, default=6)
    parser.add_argument('--k-min', type=float, required=True, help='in inverse snapshot coordinate units')
    parser.add_argument('--k-max', type=float, required=True)
    parser.add_argument('--max-power-relative-error', type=float, default=0.02)
    parser.add_argument('--max-growth-relative-error', type=float, default=0.01)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.mesh < 8 or args.bins < 1 or not (0 < args.k_min < args.k_max):
        parser.error('invalid estimator mesh/bins or k interval')
    if not (0 < args.max_power_relative_error < 1 and 0 < args.max_growth_relative_error < 1):
        parser.error('accuracy budgets must lie in (0,1)')
    states = {name: read(getattr(args, name.replace('-', '_'))) for name in
              ('initial', 'global-reference', 'hierarchical', 'refined')}
    initial, reference = states['initial'], states['global-reference']
    if args.k_max > 0.25*np.min(np.pi*args.mesh/reference['box']):
        raise ValueError('requested range exceeds the conservative quarter-Nyquist offline estimator limit')
    for name, state in states.items():
        if not (np.array_equal(state['ids'], initial['ids']) and np.array_equal(state['mass'], initial['mass'])
                and np.array_equal(state['box'], initial['box'])):
            raise ValueError(f'{name}: frozen identities, masses or geometry differ')
        if name != 'initial' and abs(state['a']-reference['a']) > 1e-12*max(1, reference['a']):
            raise ValueError('final cosmological epochs differ')
    if initial['a'] >= reference['a']:
        raise ValueError('initial snapshot must precede the final cosmological epoch')
    spectra = {name: power(state, args.mesh, args.k_min, args.k_max, args.bins)[0]
               for name, state in states.items()}
    p0, pref = spectra['initial'], spectra['global-reference']
    if (p0 <= 0).any() or (pref <= 0).any():
        raise ValueError('growth comparison requires resolved positive initial/final mode power')
    result = dict(status='qualification_failed', k_interval_requested=[args.k_min, args.k_max],
                  estimator='CIC/no deconvolution/no shot noise subtraction', mesh=args.mesh, bins=args.bins,
                  final_scale_factor=reference['a'], budgets=dict(power=args.max_power_relative_error,
                  growth=args.max_growth_relative_error), files={})
    for name in states:
        path = getattr(args, name.replace('-', '_'))
        result['files'][name] = dict(path=str(path), sha256=sha256(path))
    for name in ('hierarchical', 'refined'):
        p = spectra[name]
        power_error = np.abs(p/pref-1)
        # Mode growth is sqrt(P_final/P_initial), compared with global KDK.
        growth_error = np.abs(np.sqrt(p/p0)/np.sqrt(pref/p0)-1)
        dx = (states[name]['x']-reference['x'])/reference['box']
        dx -= np.rint(dx)
        result[name] = dict(power_relative_by_shell=power_error.tolist(),
                            growth_relative_by_shell=growth_error.tolist(),
                            position_rms_box_fraction=float(np.sqrt(np.mean(dx*dx))),
                            power_max=float(power_error.max()), growth_max=float(growth_error.max()))
    refined = result['refined']; coarse = result['hierarchical']
    passed = (refined['power_max'] <= args.max_power_relative_error and
              refined['growth_max'] <= args.max_growth_relative_error and
              refined['power_max'] <= coarse['power_max']+1e-12 and
              refined['position_rms_box_fraction'] <= coarse['position_rms_box_fraction']+1e-12)
    result['status'] = 'passed_requested_range' if passed else 'qualification_failed'
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    if not passed:
        raise SystemExit('hierarchical convergence budgets failed; inspect the saved report')


if __name__ == '__main__':
    main()
