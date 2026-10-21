#!/usr/bin/env python3
import argparse
import csv
import statistics
import sys

MIN_PEAKS = 20
FALLBACK = 9.0
COLS = ['sample_id', 'setting', 'tolerance_seconds', 'pooled_with', 'source']


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def read_neighbours(path):
    out = {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            out[r['sample_id']] = [n for n in (r.get('neighbours') or '').split(',') if n]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base', required=True, help='base_quant.long.tsv')
    ap.add_argument('--primary', required=True, help='primary_abundance')
    ap.add_argument('--fdr', type=float, required=True)
    ap.add_argument('--setting', required=True, help="pair_rt_tolerance_seconds: 'auto' or seconds")
    ap.add_argument('--neighbours', help="CHOOSE_NEIGHBOURS's table: each run is pooled with its neighbours")
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    setting = a.setting.strip()
    fixed = None if setting.lower() == 'auto' else float(setting)
    tools = ['encyclopedia', 'standard_peak'] if a.primary == 'standard_peak' else [a.primary]
    widths, runs = {t: {} for t in tools}, set()
    with open(a.base) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            runs.add(r['sample_id'])
            q, rank = num(r.get('id_qvalue')), num(r.get('peak_group_rank'))
            if r.get('channel') != 'heavy' or q is None or q > a.fdr or (rank is not None and rank > 1):
                continue
            for t in tools:
                w, ab = num(r.get(f'width_{t}_seconds')), num(r.get(f'abundance_{t}'))
                if w and w > 0 and ab and ab > 0:
                    widths[t].setdefault(r['sample_id'], []).append(w)
    neighbours = read_neighbours(a.neighbours) if a.neighbours else {}

    rows = []
    for sid in sorted(runs):
        pool = [sid] + [n for n in neighbours.get(sid, []) if n in runs and n != sid]
        if fixed is not None:
            tol, source = fixed, 'set by pair_rt_tolerance_seconds'
        else:
            tol, source = FALLBACK, f'{FALLBACK:g} s: fewer than {MIN_PEAKS} peaks of detected heavy standards with a width'
            for t in tools:
                w = [x for r in pool for x in widths[t].get(r, [])]
                if len(w) >= MIN_PEAKS:
                    median = statistics.median(w)
                    tol = round(median / 2, 2)
                    source = (f'half the median width, {median:.2f} s, of the {t} peaks of {len(w)} detected heavy '
                              f'standards in this run and its neighbours')
                    break
        rows.append(dict(sample_id=sid, setting=setting, tolerance_seconds=f'{tol:g}', pooled_with=','.join(pool[1:]),
                         source=source))
    other = fixed if fixed is not None else FALLBACK
    rows.append(dict(sample_id='*', setting=setting, tolerance_seconds=f'{other:g}', pooled_with='',
                     source='any other run'))
    with open(a.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    vals = [float(r['tolerance_seconds']) for r in rows[:-1]]
    print(f'pair_rt_tolerance_seconds {setting}: ' + (f'{min(vals):g} to {max(vals):g} s over {len(vals)} runs'
                                                      if vals else 'no run'), file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
