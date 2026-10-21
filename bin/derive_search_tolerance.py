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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--counts', required=True, help="folder of the runs' <sample>" + SUFFIX)
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
    runs, total, upto = [], None, None
    for f in sorted(glob.glob(os.path.join(a.counts, '*' + SUFFIX))):
        u, c = read_counts(f)
        if upto is None:
            upto, total = u, np.zeros(len(c), int)
        if len(u) != len(upto) or not np.allclose(u, upto):
            sys.exit(f'derive_search_tolerance: {f} is not binned like the others')
        own = tolerance(u, c, a.share)
        runs.append((os.path.basename(f)[:-len(SUFFIX)], int(c.sum()), None if own is None else round(own * T, 4)))
        total += c
    pooled = tolerance(upto, total, a.share) if runs else None
    if a.setting != 'auto':
        value, at, source = T, 0, 'the search arguments\' -ftol: search_tolerance off'
    elif pooled is None:
        value, at, source = T, 0, ('the search arguments\' -ftol: no run measured its fragment errors' if not runs else
                                   'the search arguments\' -ftol: the runs\' fragment errors give none')
    else:
        value, at = round(pooled * T, 4), len(runs)
        source = (f'the width holding {100 * a.share:g} % of the calibrants\' true fragment matches after the m/z '
                  f'correction, the {int(total.sum())} fragment errors of {len(runs)} runs together')
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['route', 'setting', 'search_tolerance', 'unit', 'runs_at_it', 'source'])
        w.writerow([a.where, a.setting, f'{value:g}', unit, at, source])
    with open(a.runs_out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['sample_id', 'fragment_errors', 'own_search_tolerance', 'searched_at'])
        for sid, n, own in runs:
            w.writerow([sid, n, '' if own is None else f'{own:g}', f'{value:g}' if at else f'{T:g}'])
    print(f'{a.where}: search at a fragment tolerance of {value:g} {unit} ({source})', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
