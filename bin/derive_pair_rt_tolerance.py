#!/usr/bin/env python3
import argparse
import csv
import statistics
import sys

MIN_PEAKS = 20
FALLBACK = 9.0


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base', required=True, help='base_quant.long.tsv')
    ap.add_argument('--primary', required=True, help='primary_abundance')
    ap.add_argument('--fdr', type=float, required=True)
    ap.add_argument('--setting', required=True, help="pair_rt_tolerance_seconds: 'auto' or seconds")
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    setting = a.setting.strip()
    if setting.lower() != 'auto':
        tol, source = float(setting), 'set by pair_rt_tolerance_seconds'
    else:
        tool = 'encyclopedia' if a.primary == 'standard_peak' else a.primary
        wcol, acol = f'width_{tool}_seconds', f'abundance_{tool}'
        widths = []
        with open(a.base) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                q, w, ab = num(r.get('id_qvalue')), num(r.get(wcol)), num(r.get(acol))
                rank = num(r.get('peak_group_rank'))
                if (r.get('channel') == 'heavy' and q is not None and q <= a.fdr
                        and (rank is None or rank <= 1)
                        and w and w > 0 and ab and ab > 0):
                    widths.append(w)
        if len(widths) >= MIN_PEAKS:
            median = statistics.median(widths)
            tol = round(median / 2, 2)
            source = (f'half the median width, {median:.2f} s, of the {tool} peaks of '
                      f'{len(widths)} detected heavy standards')
        else:
            tol = FALLBACK
            source = (f'{FALLBACK:g} s: fewer than {MIN_PEAKS} {tool} peaks of detected '
                      f'heavy standards with a width ({len(widths)})')

    with open(a.out, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['setting', 'tolerance_seconds', 'source'])
        w.writerow([setting, f'{tol:g}', source])
    print(f'pair_rt_tolerance_seconds {setting}: {tol:g} s ({source})', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
