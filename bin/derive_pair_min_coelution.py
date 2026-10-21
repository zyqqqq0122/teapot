#!/usr/bin/env python3
import argparse
import csv
import glob
import sys

import numpy as np

MIN_NULLS = 200
FALLBACK = 0.7


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tables', required=True, help='folder of the runs\' <sample>.standard_peak.tsv')
    ap.add_argument('--setting', required=True, help="pair_min_coelution: 'auto' or a number")
    ap.add_argument('--rate', type=float, default=0.05,
                    help='auto: the co-elution this share of the light traces where the light cannot be reach')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    setting = a.setting.strip()
    nulls, pairs, runs = [], 0, set()
    for f in sorted(glob.glob(f'{a.tables}/*.standard_peak.tsv')):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                v = [float(x) for x in (r.get('coelution_null') or '').split(',') if x]
                if v:
                    nulls += v
                    pairs += 1
                    runs.add(r['sample_id'])
    if setting.lower() != 'auto':
        value, source = float(setting), 'set by pair_min_coelution'
    elif len(nulls) >= MIN_NULLS:
        value = round(float(np.quantile(nulls, 1 - a.rate)), 3)
        source = (f'reached by {100 * a.rate:g} % of {len(nulls)} light traces 2 and 3 median peak widths from '
                  f'their standard, where the light cannot be ({pairs} pairs, {len(runs)} runs)')
    else:
        value = FALLBACK
        source = f'{FALLBACK:g}: fewer than {MIN_NULLS} light traces away from their standard to learn from ({len(nulls)})'
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['setting', 'min_coelution', 'source'])
        w.writerow([setting, f'{value:g}', source])
    print(f'pair_min_coelution {setting}: {value:g} ({source})', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
