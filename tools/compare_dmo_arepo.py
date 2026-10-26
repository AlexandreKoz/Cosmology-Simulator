#!/usr/bin/env python3
"""Matched-IC CHUI vs AREPO dark-matter-only comparison (2026-10-26).

AUTHORED, NOT RUN. Owner-executed validation tool.

Reads one CHUI science snapshot and one AREPO (or GADGET-4) snapshot of the
same initial conditions at the same scale factor through the validated
``tools/chui_snapshot.read_snapshot`` reader (IDs, units and completion
markers are checked there), then reports:

  * mass conservation (total mass of each code and their relative difference);
  * particle displacement statistics between codes (matched by ParticleIDs,
    minimum-image, in units of the mean interparticle spacing);
  * matter power spectra P(k) (CIC deposit, CIC window deconvolution,
    Poisson shot-noise subtraction) and their ratio;
  * large-scale phase agreement via the cross-correlation coefficient
    r(k) = P_xy / sqrt(P_xx P_yy);
  * FoF halo mass functions (linking length b times the mean spacing) and
    matched-halo position offsets / mass ratios (>= 50 % shared members).

Agreement with AREPO is NOT a correctness proof: both codes carry force,
time-integration and discretization error. Every threshold below is
PROVISIONAL, scientifically motivated but unqualified; ``--enforce`` turns the
provisional gates into a non-zero exit status, otherwise the verdicts are
informational only.

Example:
  python3 tools/compare_dmo_arepo.py \
      --chui output/run/snapshots/snap_010.complete \
      --arepo arepo_run/snap_010.hdf5 --arepo-dialect arepo_format3 \
      --legacy-length-unit kpc --legacy-mass-unit msun --legacy-velocity-unit km_s \
      --grid 256 --json comparison_snap_010.json
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chui_snapshot import SnapshotError, read_snapshot  # noqa: E402

# PROVISIONAL tolerances (owner qualification pending). Motivation:
#  * power ratio: TreePM codes at matched force resolution typically agree to
#    ~1 % on linear/quasi-linear scales; the band is k < k_Nyquist / 4.
#  * phase: r(k) > 0.99 on the same band tests large-scale phase coherence.
#  * mass function: Poisson-limited; compared only for bins with >= 20 halos.
#  * displacement: median offset below 10 % of the mean spacing at the
#    comparison epoch for an unconverged early-time run.
PROVISIONAL = {
    'mass_relative_difference': 1.0e-12,
    'power_ratio_band_max_abs_dev': 0.02,
    'cross_correlation_band_min': 0.99,
    'mass_function_bin_max_abs_dev': 0.10,
    'mass_function_min_halos_per_bin': 20,
    'median_displacement_over_spacing': 0.10,
}


def periodic_delta(a, b, box):
    d = a - b
    return d - box * np.round(d / box)


def cic_density(x, mass, box, n):
    """CIC deposit of positions (N,3) in [0,box) onto an n^3 periodic grid."""
    grid = np.zeros((n, n, n), dtype=np.float64)
    cell = box / n
    s = x / cell - 0.5
    i0 = np.floor(s).astype(np.int64)
    f = s - i0
    for dx in (0, 1):
        wx = f[:, 0] if dx else 1.0 - f[:, 0]
        ix = (i0[:, 0] + dx) % n
        for dy in (0, 1):
            wy = f[:, 1] if dy else 1.0 - f[:, 1]
            iy = (i0[:, 1] + dy) % n
            for dz in (0, 1):
                wz = f[:, 2] if dz else 1.0 - f[:, 2]
                iz = (i0[:, 2] + dz) % n
                np.add.at(grid, (ix, iy, iz), mass * wx * wy * wz)
    return grid


def overdensity_fft(x, mass, box, n):
    rho = cic_density(np.mod(x, box), mass, box, n)
    delta = rho / rho.mean() - 1.0
    dk = np.fft.rfftn(delta)
    k1 = 2.0 * np.pi * np.fft.fftfreq(n, d=box / n)
    kz = 2.0 * np.pi * np.fft.rfftfreq(n, d=box / n)
    kx, ky, kzz = np.meshgrid(k1, k1, kz, indexing='ij')
    # CIC window deconvolution W(k) = prod sinc^2(k_i dx / 2).
    half = box / n / 2.0

    def sinc(v):
        return np.where(v == 0.0, 1.0, np.sin(v) / np.where(v == 0.0, 1.0, v))
    window = (sinc(kx * half) * sinc(ky * half) * sinc(kzz * half)) ** 2
    dk = dk / window
    return dk, np.sqrt(kx * kx + ky * ky + kzz * kzz)


def binned_spectra(x_a, m_a, x_b, m_b, box, n, n_particles):
    da, kmag = overdensity_fft(x_a, m_a, box, n)
    db, _ = overdensity_fft(x_b, m_b, box, n)
    volume = box ** 3
    norm = volume / float(n) ** 6
    paa = (np.abs(da) ** 2) * norm
    pbb = (np.abs(db) ** 2) * norm
    pab = (da * np.conj(db)).real * norm
    # rfft halves the z axis: weight interior kz planes twice.
    weight = np.full(kmag.shape, 2.0)
    weight[..., 0] = 1.0
    if n % 2 == 0:
        weight[..., -1] = 1.0
    k_fund = 2.0 * np.pi / box
    k_nyq = np.pi * n / box
    edges = np.arange(k_fund * 0.5, k_nyq, k_fund)
    which = np.digitize(kmag.ravel(), edges)
    # Poisson shot noise for (possibly unequal) particle masses.
    shot_a = volume * float((m_a ** 2).sum()) / float(m_a.sum()) ** 2
    shot_b = volume * float((m_b ** 2).sum()) / float(m_b.sum()) ** 2
    rows = []
    for b in range(1, len(edges)):
        sel = which == b
        if not np.any(sel):
            continue
        w = weight.ravel()[sel]
        wsum = w.sum()
        k_mean = float((kmag.ravel()[sel] * w).sum() / wsum)
        p1 = float((paa.ravel()[sel] * w).sum() / wsum) - shot_a
        p2 = float((pbb.ravel()[sel] * w).sum() / wsum) - shot_b
        p12 = float((pab.ravel()[sel] * w).sum() / wsum)
        r = p12 / math.sqrt(max((p1 + shot_a) * (p2 + shot_b), 1e-300))
        rows.append({'k': k_mean, 'modes': float(wsum), 'p_chui': p1, 'p_arepo': p2,
                     'ratio': p1 / p2 if p2 != 0.0 else float('nan'), 'cross_correlation': r})
    return rows, k_nyq


def fof_groups(x, box, linking_length):
    try:
        from scipy.spatial import cKDTree
    except ImportError as error:  # pragma: no cover - optional dependency
        raise SystemExit('FoF halo comparison requires scipy (scipy.spatial.cKDTree)') from error
    tree = cKDTree(np.mod(x, box), boxsize=box)
    pairs = tree.query_pairs(linking_length, output_type='ndarray')
    parent = np.arange(len(x))

    def find(i):
        root = i
        while parent[root] != root:
            root = parent[root]
        while parent[i] != root:
            parent[i], i = root, parent[i]
        return root
    for i, j in pairs:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)
    return np.array([find(i) for i in range(len(x))])


def halo_catalogue(labels, x, mass, box, min_members):
    halos = {}
    order = np.argsort(labels, kind='stable')
    sorted_labels = labels[order]
    starts = np.flatnonzero(np.r_[True, sorted_labels[1:] != sorted_labels[:-1]])
    ends = np.r_[starts[1:], len(order)]
    for s, e in zip(starts, ends):
        members = order[s:e]
        if len(members) < min_members:
            continue
        ref = x[members[0]]
        rel = periodic_delta(x[members], ref, box)
        com = np.mod(ref + (rel * mass[members, None]).sum(axis=0) / mass[members].sum(), box)
        halos[int(sorted_labels[s])] = {'members': members, 'mass': float(mass[members].sum()), 'com': com}
    return halos


def mass_function(halos, edges):
    masses = np.array([h['mass'] for h in halos.values()])
    counts, _ = np.histogram(masses, bins=edges)
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--chui', required=True, help='CHUI snapshot (.complete marker, member or directory)')
    parser.add_argument('--arepo', required=True, help='AREPO/GADGET-4 snapshot')
    parser.add_argument('--arepo-dialect', choices=('arepo_format3', 'gadget4_hdf5'), default='arepo_format3')
    parser.add_argument('--legacy-length-unit', choices=('m', 'cm', 'kpc', 'mpc'))
    parser.add_argument('--legacy-mass-unit', choices=('kg', 'g', 'msun'))
    parser.add_argument('--legacy-velocity-unit', choices=('m_s', 'cm_s', 'km_s'))
    parser.add_argument('--grid', type=int, default=128, help='P(k) mesh per dimension')
    parser.add_argument('--fof-b', type=float, default=0.2, help='FoF linking length / mean spacing')
    parser.add_argument('--min-members', type=int, default=20)
    parser.add_argument('--skip-halos', action='store_true')
    parser.add_argument('--json', help='write the full comparison record here')
    parser.add_argument('--enforce', action='store_true', help='exit non-zero if a PROVISIONAL gate fails')
    args = parser.parse_args()

    units = (args.legacy_length_unit, args.legacy_mass_unit, args.legacy_velocity_unit)
    if any(units) and not all(units):
        raise SystemExit('legacy unit overrides require the complete length/mass/velocity triple')
    try:
        chui = read_snapshot(args.chui, need_velocity=False)
        arepo = read_snapshot(args.arepo, legacy_dialect=args.arepo_dialect,
                              legacy_units=units if all(units) else None, need_velocity=False)
    except SnapshotError as error:
        raise SystemExit('snapshot validation failed: ' + str(error))

    record = {'provisional_tolerances': PROVISIONAL, 'verdicts': {}}
    box = np.asarray(chui['box'], dtype=np.float64)
    if not np.allclose(box, np.asarray(arepo['box'], dtype=np.float64), rtol=1e-12):
        raise SystemExit('box sizes differ between snapshots')
    if not (box[0] == box[1] == box[2]):
        raise SystemExit('this tool assumes a cubic box')
    L = float(box[0])
    record['scale_factor'] = {'chui': float(chui['a']), 'arepo': float(arepo['a'])}
    if abs(chui['a'] - arepo['a']) > 1e-6 * max(1.0, abs(arepo['a'])):
        record['warning_scale_factor_mismatch'] = True

    ids_c, ids_a = np.asarray(chui['ids']), np.asarray(arepo['ids'])
    if len(ids_c) != len(ids_a) or not np.array_equal(np.sort(ids_c), np.sort(ids_a)):
        raise SystemExit('ParticleIDs differ: matched-IC comparison requires identical identity sets')
    oc, oa = np.argsort(ids_c, kind='stable'), np.argsort(ids_a, kind='stable')
    xc, xa = np.asarray(chui['x'])[oc], np.asarray(arepo['x'])[oa]
    mc, ma = np.asarray(chui['mass'])[oc], np.asarray(arepo['mass'])[oa]
    n = len(ids_c)
    spacing = L / n ** (1.0 / 3.0)

    total_c, total_a = float(mc.sum()), float(ma.sum())
    mass_dev = abs(total_c - total_a) / total_a
    record['mass'] = {'chui_total': total_c, 'arepo_total': total_a, 'relative_difference': mass_dev}
    record['verdicts']['mass_conservation'] = mass_dev <= PROVISIONAL['mass_relative_difference']

    d = periodic_delta(xc, xa, L)
    dist = np.sqrt((d * d).sum(axis=1)) / spacing
    record['displacement_over_mean_spacing'] = {
        'rms': float(np.sqrt((dist ** 2).mean())), 'median': float(np.median(dist)),
        'p90': float(np.percentile(dist, 90)), 'p99': float(np.percentile(dist, 99)), 'max': float(dist.max())}
    record['verdicts']['displacement'] = (
        record['displacement_over_mean_spacing']['median'] <= PROVISIONAL['median_displacement_over_spacing'])

    rows, k_nyq = binned_spectra(xc, mc, xa, ma, L, args.grid, n)
    record['power_spectrum'] = rows
    band = [r for r in rows if r['k'] < k_nyq / 4.0 and r['p_arepo'] > 0.0]
    if band:
        dev = max(abs(r['ratio'] - 1.0) for r in band)
        rmin = min(r['cross_correlation'] for r in band)
        record['power_band'] = {'k_max': k_nyq / 4.0, 'max_abs_ratio_deviation': dev, 'min_cross_correlation': rmin}
        record['verdicts']['power_ratio'] = dev <= PROVISIONAL['power_ratio_band_max_abs_dev']
        record['verdicts']['phase_cross_correlation'] = rmin >= PROVISIONAL['cross_correlation_band_min']

    if not args.skip_halos:
        ll = args.fof_b * spacing
        hc = halo_catalogue(fof_groups(xc, L, ll), xc, mc, L, args.min_members)
        ha = halo_catalogue(fof_groups(xa, L, ll), xa, ma, L, args.min_members)
        all_masses = [h['mass'] for h in list(hc.values()) + list(ha.values())]
        if all_masses:
            edges = np.logspace(np.log10(min(all_masses)), np.log10(max(all_masses) * 1.0001), 12)
            cc, ca = mass_function(hc, edges), mass_function(ha, edges)
            mf = []
            ok = True
            for i in range(len(edges) - 1):
                row = {'m_low': float(edges[i]), 'm_high': float(edges[i + 1]), 'n_chui': int(cc[i]), 'n_arepo': int(ca[i])}
                if ca[i] >= PROVISIONAL['mass_function_min_halos_per_bin']:
                    row['relative_deviation'] = (cc[i] - ca[i]) / ca[i]
                    ok = ok and abs(row['relative_deviation']) <= PROVISIONAL['mass_function_bin_max_abs_dev']
                mf.append(row)
            record['mass_function'] = mf
            record['verdicts']['mass_function'] = ok
            # Halo matching by shared members (>= 50 % of the AREPO halo).
            owner_c = np.full(n, -1, dtype=np.int64)
            for label, h in hc.items():
                owner_c[h['members']] = label
            matches = []
            for label, h in ha.items():
                owners = owner_c[h['members']]
                owners = owners[owners >= 0]
                if owners.size == 0:
                    continue
                values, counts = np.unique(owners, return_counts=True)
                best = int(values[np.argmax(counts)])
                if counts.max() >= 0.5 * len(h['members']):
                    offset = periodic_delta(hc[best]['com'][None, :], h['com'][None, :], L)[0]
                    matches.append({'mass_ratio': hc[best]['mass'] / h['mass'],
                                    'offset_over_spacing': float(np.sqrt((offset ** 2).sum()) / spacing)})
            record['halo_matching'] = {
                'arepo_halos': len(ha), 'chui_halos': len(hc), 'matched': len(matches),
                'median_mass_ratio': float(np.median([m['mass_ratio'] for m in matches])) if matches else None,
                'median_offset_over_spacing': float(np.median([m['offset_over_spacing'] for m in matches])) if matches else None}

    record['all_provisional_gates_pass'] = all(record['verdicts'].values())
    text = json.dumps(record, indent=2, default=float)
    if args.json:
        Path(args.json).write_text(text + '\n')
    print(json.dumps({'verdicts': record['verdicts'], 'power_band': record.get('power_band'),
                      'displacement_over_mean_spacing': record['displacement_over_mean_spacing'],
                      'halo_matching': record.get('halo_matching')}, indent=2, default=float))
    if args.enforce and not record['all_provisional_gates_pass']:
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
