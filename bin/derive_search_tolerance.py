#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import sys

import numpy as np

FLOOR_FROM = 0.75
START = 0.05
SUFFIX = '.fragment_mz.tolerance.tsv'


def read_counts(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    return np.array([float(r['upto_T']) for r in rows]), np.array([int(r['errors']) for r in rows])


def tolerance(upto, counts, share):
    width = upto[-1]
    floor = counts[upto > FLOOR_FROM * width + 1e-9].sum() / ((1 - FLOOR_FROM) * width)
    true = counts.sum() - floor * width
    if true <= 0:
        return None
    for u, c in zip(upto, np.cumsum(counts)):
        if u >= START - 1e-9 and (c - floor * u) / true >= share:
            return float(u)
    return None


def read_neighbours(path):
    out = {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            out[r['sample_id']] = [n for n in (r.get('neighbours') or '').split(',') if n]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--counts', required=True, help="folder of the runs' <sample>" + SUFFIX)
    ap.add_argument('--neighbours', help="CHOOSE_NEIGHBOURS's table: each run is pooled with its neighbours; "
                    "without it every run is pooled with all the others")
    ap.add_argument('--fragment-tolerance', type=float, required=True, help="T, the search arguments' -ftol")
    ap.add_argument('--fragment-unit', choices=('th', 'ppm'), required=True)
    ap.add_argument('--setting', required=True, help="search_tolerance: auto or off")
    ap.add_argument('--share', type=float, default=0.99,
                    help="the tolerance holds this share of the calibrants' true fragment matches")
    ap.add_argument('--where', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--runs-out', required=True)
    a = ap.parse_args()

    T, unit = a.fragment_tolerance, 'Th' if a.fragment_unit == 'th' else 'ppm'
    counts, upto = {}, None
    for f in sorted(glob.glob(os.path.join(a.counts, '*' + SUFFIX))):
        u, c = read_counts(f)
        if upto is None:
            upto = u
        if len(u) != len(upto) or not np.allclose(u, upto):
            sys.exit(f'derive_search_tolerance: {f} is not binned like the others')
        counts[os.path.basename(f)[:-len(SUFFIX)]] = c
    auto = a.setting == 'auto'
    neighbours = read_neighbours(a.neighbours) if a.neighbours else None

    def value(runs):
        if not auto or not runs:
            return None
        v = tolerance(upto, sum(counts[r] for r in runs), a.share)
        return None if v is None else round(v * T, 4)

    batch = value(sorted(counts))
    rows = []
    for sid in sorted(counts):
        own = tolerance(upto, counts[sid], a.share)
        pool = ([sid] + [n for n in neighbours.get(sid, []) if n in counts]) if neighbours is not None else sorted(counts)
        v = value(pool)
        rows.append([sid, int(counts[sid].sum()), '' if own is None else f'{round(own * T, 4):g}',
                     ','.join(r for r in pool if r != sid) if neighbours is not None else 'all runs',
                     f'{v if v is not None else T:g}', 'measured' if v is not None else '-ftol'])
    if not auto:
        source = 'the search arguments\' -ftol: search_tolerance off'
    elif batch is None:
        source = ('the search arguments\' -ftol: no run measured its fragment errors' if not counts else
                  'the search arguments\' -ftol: the runs\' fragment errors give none')
    else:
        source = (f'the width holding {100 * a.share:g} % of the calibrants\' true fragment matches after the m/z '
                  f'correction, the {int(sum(counts.values()).sum())} fragment errors of {len(counts)} runs together'
                  + ('; each run is searched at the value of itself and its neighbours (.runs.tsv), this one is '
                     'for the steps that read every run at once' if neighbours is not None else ''))
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['route', 'setting', 'search_tolerance', 'unit', 'runs_at_it', 'source'])
        w.writerow([a.where, a.setting, f'{batch if batch is not None else T:g}', unit,
                    len(counts) if batch is not None else 0, source])
    with open(a.runs_out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['sample_id', 'fragment_errors', 'own_search_tolerance', 'pooled_with', 'searched_at', 'source'])
        w.writerows(rows)
    print(f'{a.where}: {len(rows)} runs measured their fragment errors; '
          f'{sum(1 for r in rows if r[5] == "measured")} searched at the value of '
          f'{"themselves and their neighbours" if neighbours is not None else "all of them together"}; '
          f'the batch value {batch if batch is not None else T:g} {unit} ({source})', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
