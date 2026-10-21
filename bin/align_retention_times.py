#!/usr/bin/env python3
import argparse
import csv
import json
import os
import re
import sqlite3
import sys
from datetime import datetime

import numpy as np

MASS = re.compile(r'\[([+-]?\d+(?:\.\d+)?)\]')
FOLDS = 10
ROUNDS = 3
COLS = ['sample_id', 'status', 'anchors', 'neighbours', 'residual_sd_seconds', 'half_window_seconds',
        'precursors_with_window', 'windows_widened', 'acquired']


def is_label(text, deltas):
    tol = max(0.01, 0.5 * 10 ** -len(text.split('.')[1]) + 1e-9) if '.' in text else 0.5
    return any(abs(float(text) - d) <= tol for d in deltas)


def peptide_key(modseq, deltas):
    s = re.sub(r'^[A-Z-]?\.|\.[A-Z-]?$', '', modseq.strip())
    return MASS.sub(lambda m: '' if is_label(m.group(1), deltas) else '[%+.1f]' % float(m.group(1)), s)


def read_acquired(paths):
    out = {}
    for p in paths:
        with open(p, newline='') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                t = (r.get('acquired') or '').strip()
                if t:
                    out[r['sample_id']] = t
    return out


def seconds(stamp):
    return datetime.fromisoformat(stamp.replace('Z', '+00:00')).timestamp()


def lowess(x, y, frac=0.3, rounds=3):
    o = np.argsort(x)
    x, y = np.asarray(x, float)[o], np.asarray(y, float)[o]
    n = len(x)
    k = min(n, max(5, int(frac * n)))
    grid = np.linspace(x[0], x[-1], 60)
    w = np.ones(n)
    for _ in range(rounds):
        fit = np.empty(len(grid))
        for g, x0 in enumerate(grid):
            d = np.abs(x - x0)
            idx = np.argsort(d)[:k]
            h = d[idx].max() + 1e-9
            tw = (1 - (d[idx] / h) ** 3) ** 3 * w[idx]
            X = np.column_stack([np.ones(k), x[idx] - x0])
            beta = np.linalg.lstsq((X * tw[:, None]).T @ X, (X * tw[:, None]).T @ y[idx], rcond=None)[0]
            fit[g] = beta[0]
        r = y - np.interp(x, grid, fit)
        s = 6 * np.median(np.abs(r)) + 1e-9
        w = np.where(np.abs(r) < s, (1 - (r / s) ** 2) ** 2, 0.0)
    slope = np.polyfit(x, y, 1)[0] if n >= 2 else 1.0

    def f(v):
        v = np.asarray(v, float)
        out = np.interp(v, grid, fit)
        out = np.where(v < grid[0], fit[0] + slope * (v - grid[0]), out)
        return np.where(v > grid[-1], fit[-1] + slope * (v - grid[-1]), out)
    return f


def read_located(path, deltas, every=False):
    out = {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r.get('used', '').strip().lower() != 'true':
                continue
            if deltas and not any(is_label(m, deltas) for m in MASS.findall(r['modseq'])):
                continue
            try:
                t = float(r['apex_rt'])
            except (KeyError, ValueError):
                continue
            out.setdefault(peptide_key(r['modseq'], deltas), []).append(t)
    return out if every else {k: float(np.mean(v)) for k, v in out.items()}


def read_reference(path, deltas):
    lines = open(path, newline='').read().splitlines()
    rows = list(csv.reader(lines, delimiter='\t' if lines[0].count('\t') >= lines[0].count(',') else ','))
    header = [h.strip().strip('"').lower() for h in rows[0]]
    dcol = next((i for i, h in enumerate(header) if h in ('isdecoy', 'decoy', 'is_decoy')), None)
    keys = set()
    for r in rows[1:]:
        if len(r) < 5 or not r[0].strip():
            continue
        if dcol is not None and len(r) > dcol and r[dcol].strip().lower() in ('true', '1', 'yes', 'decoy'):
            continue
        keys.add(peptide_key(r[0], deltas))
    return keys


def read_library(path, deltas, keys):
    con = sqlite3.connect(f'file:{os.path.abspath(path)}?mode=ro', uri=True)
    out = {}
    for modseq, z in con.execute('SELECT DISTINCT PeptideModSeq, PrecursorCharge FROM entries'):
        k = peptide_key(modseq, deltas)
        if k in keys:
            out.setdefault(k, []).append((modseq, int(z)))
    con.close()
    return out


def warp_onto(runs, located, cons, rounds=ROUNDS):
    keys = sorted(set().union(*(located[r] for r in runs)))
    cur = {k: float(np.median([located[r][k] for r in runs if k in located[r]])) for k in keys} if cons is None else cons
    for _ in range(rounds):
        warped = {k: [] for k in keys}
        for r in runs:
            ks = [k for k in located[r] if k in cur]
            if len(ks) < 5:
                continue
            f = lowess([located[r][k] for k in ks], [cur[k] for k in ks])
            for k, v in zip(located[r], f([located[r][k] for k in located[r]])):
                warped[k].append(float(v))
        cur = {k: float(np.median(v)) for k, v in warped.items() if v}
    return cur


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--samples', required=True, help='the runs to align, one sample_id per line')
    ap.add_argument('--located', nargs='*', default=[], help='<sample>.fragment_mz.targets.tsv of the runs that have one')
    ap.add_argument('--acquired', nargs='*', default=[], help="the m/z step's summary rows (sample_id, acquired)")
    ap.add_argument('--library', required=True, help='the search library (.dlib or .elib)')
    ap.add_argument('--reference', required=True, help='reference list')
    ap.add_argument('--labels', default='[]', help='JSON list of heavy labels with their delta; [] for label-free')
    ap.add_argument('--neighbours', type=int, default=4)
    ap.add_argument('--half-window', default='auto', help='seconds, or auto')
    ap.add_argument('--min-anchors', type=int, default=20)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--summary', required=True, help='one QC row per run')
    ap.add_argument('--warnings', required=True)
    a = ap.parse_args()

    deltas = [float(r['delta']) for r in json.loads(a.labels)]
    os.makedirs(a.out_dir, exist_ok=True)
    samples = sorted({l.strip() for l in open(a.samples) if l.strip()})
    located, found = {}, {}
    for p in a.located:
        sid = os.path.basename(p).split('.fragment_mz.targets.tsv')[0]
        if sid in samples:
            located[sid] = read_located(p, deltas)
            found[sid] = read_located(p, deltas, every=True)
    acquired = read_acquired(a.acquired)
    stamps = {sid: seconds(acquired[sid]) if sid in acquired else None for sid in samples}
    ref = read_reference(a.reference, deltas)
    lib = read_library(a.library, deltas, ref)
    warnings, rows = [], []
    usable = [s for s in samples if s in located and len(located[s]) >= a.min_anchors and stamps[s] is not None]
    order = {s: i for i, s in enumerate(samples)}

    for sid in samples:
        row = dict(sample_id=sid, status='', anchors=len(located.get(sid, {})), neighbours='',
                   residual_sd_seconds='', half_window_seconds='', precursors_with_window=0,
                   windows_widened=0, acquired=acquired.get(sid, ''))
        out = os.path.join(a.out_dir, f'{sid}.apex_windows.tsv')
        with open(out, 'w') as fh:
            fh.write('PeptideModSeq\tPrecursorCharge\tMinApexRTInSeconds\tMaxApexRTInSeconds\n')
        others = [s for s in usable if s != sid]
        if sid not in usable:
            row['status'] = ('not aligned: no standards located (not a targeted mzML run)' if sid not in located
                             else 'not aligned: the run states no acquisition time' if stamps.get(sid) is None
                             else f'not aligned: {row["anchors"]} anchors located, {a.min_anchors} needed')
            warnings.append(f'{sid}: {row["status"]}; searched without retention-time windows')
            rows.append(row)
            continue
        if not others:
            row['status'] = 'not aligned: no other run to align to'
            warnings.append(f'{sid}: {row["status"]}; searched without retention-time windows')
            rows.append(row)
            continue
        near = sorted(others, key=lambda s: (abs(stamps[s] - stamps[sid]), order[s]))[:a.neighbours]
        cons = warp_onto(near, located, None)
        mine = located[sid]
        ks = [k for k in mine if k in cons]
        if len(ks) < a.min_anchors:
            row['status'] = f'not aligned: {len(ks)} anchors shared with the neighbours, {a.min_anchors} needed'
            warnings.append(f'{sid}: {row["status"]}; searched without retention-time windows')
            rows.append(row)
            continue
        x = np.array([cons[k] for k in ks])
        y = np.array([mine[k] for k in ks])
        oof = np.empty(len(ks))
        fold = np.arange(len(ks)) % FOLDS
        for f in range(FOLDS):
            te = fold == f
            oof[te] = lowess(x[~te], y[~te])(x[te])
        resid = y - oof
        rsd = 1.4826 * float(np.median(np.abs(resid - np.median(resid))))
        half = 4 * rsd if a.half_window == 'auto' else float(a.half_window)
        full = lowess(x, y)
        expected = {k: float(v) for k, v in zip(ks, oof)}
        for k in cons:
            if k not in expected:
                expected[k] = float(full(cons[k]))
        n = widened = 0
        with open(out, 'a') as fh:
            for k in sorted(lib):
                if k not in expected:
                    continue
                t = expected[k]
                lo, hi = t - half, t + half
                own = [v for v in found[sid].get(k, []) if not lo <= v <= hi]
                if own:
                    lo, hi = min([lo] + [v - half for v in own]), max([hi] + [v + half for v in own])
                    widened += 1
                for modseq, z in sorted(lib[k]):
                    fh.write(f'{modseq}\t{z}\t{lo:.2f}\t{hi:.2f}\n')
                    n += 1
        row.update(status='aligned', neighbours=','.join(near), residual_sd_seconds=round(rsd, 2),
                   half_window_seconds=round(half, 2), precursors_with_window=n, windows_widened=widened)
        if rsd > 5:
            warnings.append(f'{sid}: retention times agree with the neighbouring runs to {rsd:.1f} s (robust SD); '
                            f'windows of +/- {half:.0f} s')
        rows.append(row)
        print(f'{sid}: {len(ks)} anchors, neighbours {",".join(near)}; residual SD {rsd:.1f} s; '
              f'+/- {half:.1f} s windows for {n} precursors; {widened} peptides\' windows widened to their '
              f'standard\'s own apex in the run')

    with open(a.summary, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    with open(a.warnings, 'w') as fh:
        fh.write(''.join(f'{x}\n' for x in warnings))
    return 0


if __name__ == '__main__':
    sys.exit(main())
