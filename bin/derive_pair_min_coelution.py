#!/usr/bin/env python3
import argparse
import csv
import glob
import sys

import numpy as np

MIN_NULLS = 200
FALLBACK = 0.7
COLS = ['sample_id', 'setting', 'min_coelution', 'chance_pass', 'moved_traces', 'pooled_with', 'source']


def read_neighbours(path):
    out = {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            out[r['sample_id']] = [n for n in (r.get('neighbours') or '').split(',') if n]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tables', required=True, help='folder of the runs\' <sample>.standard_peak.tsv')
    ap.add_argument('--neighbours', help="CHOOSE_NEIGHBOURS's table: each run is pooled with its neighbours")
    ap.add_argument('--setting', required=True, help="pair_min_coelution: 'auto' or a number")
    ap.add_argument('--rate', type=float, default=0.05,
                    help='auto: the co-elution this share of the light traces where the light cannot be reach')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    setting = a.setting.strip()
    nulls = {}
    for f in sorted(glob.glob(f'{a.tables}/*.standard_peak.tsv')):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                v = [float(x) for x in (r.get('coelution_null') or '').split(',') if x]
                nulls.setdefault(r['sample_id'], []).extend(v)
    neighbours = read_neighbours(a.neighbours) if a.neighbours else {}
    fixed = None if setting.lower() == 'auto' else float(setting)

    rows = []
    for sid in sorted(nulls):
        pool = [sid] + [n for n in neighbours.get(sid, []) if n in nulls and n != sid]
        v = np.array([x for r in pool for x in nulls[r]])
        if fixed is not None:
            value, source = fixed, 'set by pair_min_coelution'
        elif len(v) >= MIN_NULLS:
            value = round(float(np.quantile(v, 1 - a.rate)), 3)
            source = (f'reached by {100 * a.rate:g} % of the light traces moved 2 and 3 median peak widths off '
                      f'their standard, where the light cannot be, in this run and its neighbours')
        else:
            value = FALLBACK
            source = f'{FALLBACK:g}: fewer than {MIN_NULLS} moved light traces in this run and its neighbours'
        rows.append(dict(sample_id=sid, setting=setting, min_coelution=f'{value:g}',
                         chance_pass=f'{float((v >= value).mean()):.3f}' if len(v) and value > 0 else '',
                         moved_traces=len(v), pooled_with=','.join(pool[1:]), source=source))
    other = fixed if fixed is not None else FALLBACK
    rows.append(dict(sample_id='*', setting=setting, min_coelution=f'{other:g}', chance_pass='', moved_traces=0,
                     pooled_with='', source='any other run' + (': set by pair_min_coelution' if fixed is not None else
                                                              f': {FALLBACK:g}')))
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    vals = [float(r['min_coelution']) for r in rows[:-1]]
    print(f'pair_min_coelution {setting}: ' + (f'{min(vals):g} to {max(vals):g} over {len(vals)} runs'
                                               if vals else 'no run measured'), file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
