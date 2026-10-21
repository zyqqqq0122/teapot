#!/usr/bin/env python3
import argparse
import csv
import glob
import json
import os
import sys
import warnings

import numpy as np

from align_retention_times import library_times, lowess, read_located, warp_onto

COLS = ['sample_id', 'anchors', 'neighbours', 'agreement_seconds', 'predicted_seconds', 'status']


def robust_sd(d, need):
    d = d[np.isfinite(d)]
    return 1.4826 * float(np.median(np.abs(d - np.median(d)))) if len(d) >= need else np.inf


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--samples', required=True, help='the route\'s runs, one sample_id per line')
    ap.add_argument('--located', required=True, help='folder of <sample>.fragment_mz.targets.tsv')
    ap.add_argument('--library', required=True, help='the search library, whose retention times every run is mapped onto')
    ap.add_argument('--labels', default='[]', help='JSON list of heavy labels with their delta; [] for label-free')
    ap.add_argument('--neighbours', type=int, default=4)
    ap.add_argument('--min-anchors', type=int, default=20)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    deltas = [float(r['delta']) for r in json.loads(a.labels)]
    samples = sorted({l.strip() for l in open(a.samples) if l.strip()})
    located = {}
    for p in glob.glob(os.path.join(a.located, '*.fragment_mz.targets.tsv')):
        sid = os.path.basename(p).split('.fragment_mz.targets.tsv')[0]
        if sid in samples:
            located[sid] = read_located(p, deltas)
    usable = [s for s in samples if len(located.get(s, {})) >= a.min_anchors]

    chosen, predicted, library_better = {}, {}, set()
    if len(usable) > 1:
        cons = library_times(a.library, deltas, set().union(*(located[r] for r in usable)))
        own_scale = len(set(cons.values())) < a.min_anchors
        if own_scale:
            cons = warp_onto(usable, located, None)
        keys = sorted(cons)
        col = {k: j for j, k in enumerate(keys)}
        lib = np.array([cons[k] for k in keys])
        M = np.full((len(usable), len(keys)), np.nan)
        for i, r in enumerate(usable):
            ks = [k for k in located[r] if k in cons]
            if len(ks) < 5:
                continue
            f = lowess([located[r][k] for k in ks], [cons[k] for k in ks])
            M[i, [col[k] for k in ks]] = f([located[r][k] for k in ks])
        for i, r in enumerate(usable):
            d = M - M[i]
            shared = np.isfinite(d).sum(axis=1)
            with np.errstate(all='ignore'), warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                sd = 1.4826 * np.nanmedian(np.abs(d - np.nanmedian(d, axis=1)[:, None]), axis=1)
                ranked = sorted((float(sd[j]), j) for j in range(len(usable))
                                if j != i and shared[j] >= a.min_anchors)[:a.neighbours]
                best, keep = (np.inf if own_scale else robust_sd(M[i] - lib, a.min_anchors)), 0
                for k in range(1, len(ranked) + 1):
                    v = robust_sd(M[i] - np.nanmedian(M[[j for _v, j in ranked[:k]]], axis=0), a.min_anchors)
                    if v < best:
                        best, keep = v, k
            chosen[r] = [(v, usable[j]) for v, j in ranked[:keep]]
            if np.isfinite(best):
                predicted[r] = best
            if ranked and not keep:
                library_better.add(r)

    rows = []
    for sid in samples:
        n = len(located.get(sid, {}))
        near = chosen.get(sid, [])
        status = ('no standards located' if sid not in located else
                  f'{n} anchors located, {a.min_anchors} needed' if sid not in usable else
                  "the library's retention times predict its own better than any of the runs most like it"
                  if sid in library_better else
                  'no other run shares enough anchors' if not near else 'chosen')
        rows.append(dict(sample_id=sid, anchors=n, neighbours=','.join(s for _v, s in near),
                         agreement_seconds=','.join('%.2f' % v for v, _s in near),
                         predicted_seconds='%.2f' % predicted[sid] if sid in predicted else '', status=status))
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print(f'{len(samples)} runs, {len(usable)} with at least {a.min_anchors} located anchors; each given, of the '
          f'{a.neighbours} whose anchors\' times, mapped onto the library\'s, agree best with its own, the closest '
          f'ones whose median predicts them best, or none where the library\'s times predict them better',
          file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
