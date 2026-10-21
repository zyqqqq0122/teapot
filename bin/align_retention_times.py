#!/usr/bin/env python3
import argparse
import csv
import glob
import json
import os
import re
import sqlite3
import sys

import numpy as np

MASS = re.compile(r'\[([+-]?\d+(?:\.\d+)?)\]')
FOLDS = 10
ROUNDS = 3
COLS = ['sample_id', 'status', 'anchors', 'reference', 'neighbours', 'residual_sd_seconds', 'half_window_seconds',
        'precursors_with_window', 'windows_widened']


def is_label(text, deltas):
    tol = max(0.01, 0.5 * 10 ** -len(text.split('.')[1]) + 1e-9) if '.' in text else 0.5
    return any(abs(float(text) - d) <= tol for d in deltas)


def peptide_key(modseq, deltas):
    s = re.sub(r'^[A-Z-]?\.|\.[A-Z-]?$', '', modseq.strip())
    return MASS.sub(lambda m: '' if is_label(m.group(1), deltas) else '[%+.1f]' % float(m.group(1)), s)


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


def read_neighbours(path):
    out = {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            out[r['sample_id']] = [n for n in (r.get('neighbours') or '').split(',') if n]
    return out


def library_times(path, deltas, keys):
    con = sqlite3.connect(f'file:{os.path.abspath(path)}?mode=ro', uri=True)
    out = {}
    for modseq, rt in con.execute('SELECT PeptideModSeq, RTInSeconds FROM entries'):
        k = peptide_key(modseq, deltas)
        if k in keys and rt is not None:
            out.setdefault(k, []).append(float(rt))
    con.close()
    return {k: float(np.median(v)) for k, v in out.items()}


def out_of_fold(x, y):
    oof = np.empty(len(x))
    fold = np.arange(len(x)) % FOLDS
    for f in range(FOLDS):
        te = fold == f
        oof[te] = lowess(x[~te], y[~te])(x[te])
    resid = y - oof
    return oof, 1.4826 * float(np.median(np.abs(resid - np.median(resid))))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--samples', required=True, help='the runs to align, one sample_id per line')
    ap.add_argument('--located', required=True, help='folder of <sample>.fragment_mz.targets.tsv')
    ap.add_argument('--neighbours', required=True, help="CHOOSE_NEIGHBOURS's table: each run's neighbours")
    ap.add_argument('--library', required=True, help='the search library (.dlib or .elib)')
    ap.add_argument('--reference', required=True, help='reference list')
    ap.add_argument('--labels', default='[]', help='JSON list of heavy labels with their delta; [] for label-free')
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
    for p in glob.glob(os.path.join(a.located, '*.fragment_mz.targets.tsv')):
        sid = os.path.basename(p).split('.fragment_mz.targets.tsv')[0]
        if sid in samples:
            located[sid] = read_located(p, deltas)
            found[sid] = read_located(p, deltas, every=True)
    neighbours = read_neighbours(a.neighbours)
    ref = read_reference(a.reference, deltas)
    lib = read_library(a.library, deltas, ref)
    lib_rt = library_times(a.library, deltas, set(lib))
    warnings, rows = [], []

    for sid in samples:
        mine = located.get(sid, {})
        row = dict(sample_id=sid, status='', anchors=len(mine), reference='', neighbours='',
                   residual_sd_seconds='', half_window_seconds='', precursors_with_window=0, windows_widened=0)
        out = os.path.join(a.out_dir, f'{sid}.apex_windows.tsv')
        with open(out, 'w') as fh:
            fh.write('PeptideModSeq\tPrecursorCharge\tMinApexRTInSeconds\tMaxApexRTInSeconds\n')
        if len(mine) < a.min_anchors:
            row['status'] = ('not aligned: no standards located (not a targeted mzML run)' if sid not in located
                             else f'not aligned: {len(mine)} anchors located, {a.min_anchors} needed')
            warnings.append(f'{sid}: {row["status"]}; searched without retention-time windows')
            rows.append(row)
            continue
        near = [s for s in neighbours.get(sid, []) if s in located and s != sid]
        candidates = [('library', lib_rt, '')]
        if near:
            candidates.insert(0, ('neighbours', warp_onto(near, located, None), ','.join(near)))
        best = None
        for name, cons, label in candidates:
            ks = [k for k in mine if k in cons]
            if len(ks) < a.min_anchors:
                continue
            x = np.array([cons[k] for k in ks])
            y = np.array([mine[k] for k in ks])
            oof, rsd = out_of_fold(x, y)
            if best is None or rsd < best[0]:
                best = (rsd, name, cons, label, ks, x, y, oof)
        if best is None:
            row['status'] = f'not aligned: fewer than {a.min_anchors} anchors shared with its neighbours or the library'
            warnings.append(f'{sid}: {row["status"]}; searched without retention-time windows')
            rows.append(row)
            continue
        rsd, name, cons, label, ks, x, y, oof = best
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
        row.update(status='aligned', reference=name, neighbours=label, residual_sd_seconds=round(rsd, 2),
                   half_window_seconds=round(half, 2), precursors_with_window=n, windows_widened=widened)
        if rsd > 5:
            warnings.append(f'{sid}: retention times agree with its {"neighbouring runs" if name == "neighbours" else "library"} '
                            f'to {rsd:.1f} s (robust SD); windows of +/- {half:.0f} s')
        rows.append(row)
        print(f'{sid}: {len(ks)} anchors, aligned to {"its neighbours " + label if name == "neighbours" else "the library"}; '
              f'residual SD {rsd:.1f} s; +/- {half:.1f} s windows for {n} precursors; {widened} peptides\' windows widened '
              f'to their standard\'s own apex in the run')

    with open(a.summary, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    with open(a.warnings, 'w') as fh:
        fh.write(''.join(f'{x}\n' for x in warnings))
    return 0


if __name__ == '__main__':
    sys.exit(main())
