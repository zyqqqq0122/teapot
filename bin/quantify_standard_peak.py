#!/usr/bin/env python3
import argparse
import base64
import csv
import json
import os
import re
import sqlite3
import sys
import zlib

import numpy as np

OWN_WINDOW_TOLERANCE = 0.01
FRAGMENTS = 12
MIN_FRAGMENTS = 3
PEAK_FLOOR = 0.1
MAX_HALF_WIDTHS = 2.5
MIN_SHARE = 0.05
MARGIN_WIDTHS = 1.0
NULL_SHIFTS = (-3, -2, 2, 3)
MASS = re.compile(r'\[([+-]?\d+(?:\.\d+)?)\]')
ARRAY = re.compile(rb'<binaryDataArray\b.*?</binaryDataArray>', re.S)
BINARY = re.compile(rb'<binary>(.*?)</binary>', re.S)
LVL = re.compile(rb'accession="MS:1000511"[^>]*value="([^"]+)"')
TGT = re.compile(rb'accession="MS:1000827"[^>]*value="([^"]+)"')
LO = re.compile(rb'accession="MS:1000828"[^>]*value="([^"]+)"')
UP = re.compile(rb'accession="MS:1000829"[^>]*value="([^"]+)"')
RT = re.compile(rb'accession="MS:1000016"[^>]*value="([^"]+)"[^>]*unitAccession="([^"]+)"')
SPECTRUM = re.compile(rb'<spectrum\s')
PSM = re.compile(r'.*:([\d.]+):(.*)\+(\d+)$')
COLS = ['sample_id', 'peptide', 'peptidoform', 'charge', 'heavy_mz', 'light_mz', 'status',
        'heavy_q_value', 'apex', 'start', 'stop', 'fragments', 'fragments_dropped', 'light_area',
        'heavy_area', 'ratio_L_H', 'coelution', 'rdotp', 'light_q_value', 'light_rt', 'coelution_null']


def is_label(text, deltas):
    tol = max(0.01, 0.5 * 10 ** -len(text.split('.')[1]) + 1e-9) if '.' in text else 0.5
    return any(abs(float(text) - d) <= tol for d in deltas)


def is_heavy(compound, deltas):
    return any(is_label(t, deltas) for t in MASS.findall(compound))


def light_form(compound, deltas):
    return MASS.sub(lambda m: '' if is_label(m.group(1), deltas) else '[%+.1f]' % float(m.group(1)), compound)


def stripped(compound):
    return re.sub(r'\[[^\]]*\]', '', compound)


def fragment_options(search_args):
    words = search_args.split()
    opt = {words[i]: words[i + 1] for i in range(len(words) - 1) if words[i].startswith('-')}
    tol = float(opt.get('-ftol', 10))
    unit = opt.get('-ftolunits', 'ppm').lower()
    offset = float(opt.get('-foffset', 0))
    if unit == 'resolution':
        return 0.5e6 / tol, True, offset
    if unit not in ('ppm', 'amu'):
        sys.exit(f'quantify_standard_peak: unknown -ftolunits {unit}')
    return tol, unit == 'ppm', offset


def read_reference(path, deltas):
    lines = open(path, newline='').read().splitlines()
    rows = list(csv.reader(lines, delimiter='\t' if lines[0].count('\t') >= lines[0].count(',') else ','))
    header = [h.strip().strip('"').lower() for h in rows[0]]
    dcol = next((i for i, h in enumerate(header) if h in ('isdecoy', 'decoy', 'is_decoy')), None)
    out = []
    for r in rows[1:]:
        if len(r) < 5 or not r[0].strip():
            continue
        if dcol is not None and len(r) > dcol and r[dcol].strip().lower() in ('true', '1', 'yes', 'decoy'):
            continue
        c = r[0].strip()
        out.append({'compound': c, 'mz': float(r[3]), 'z': int(float(r[4])), 'heavy': is_heavy(c, deltas),
                    'light_form': light_form(c, deltas)})
    return out


def read_library(path, deltas, peptides):
    con = sqlite3.connect(f'file:{os.path.abspath(path)}?mode=ro', uri=True)
    lib = {}
    peptides = sorted(peptides)
    for i in range(0, len(peptides), 500):
        chunk = peptides[i:i + 500]
        for modseq, z, mblob, iblob in con.execute(
                'SELECT PeptideModSeq, PrecursorCharge, MassArray, IntensityArray FROM entries '
                'WHERE PeptideSeq IN (%s)' % ','.join('?' * len(chunk)), chunk):
            mz = np.frombuffer(zlib.decompress(mblob), '>f8')
            it = np.frombuffer(zlib.decompress(iblob), '>f4').astype(float)
            key = (light_form(modseq, deltas), int(z), is_heavy(modseq, deltas))
            if key not in lib or it.sum() > lib[key][1].sum():
                lib[key] = (mz, it)
    con.close()
    return lib


def match_fragments(heavy, light, pep, deltas_by_site):
    hmz, hit = heavy
    lmz, _lit = light
    shifts = [0.0]
    for site, delta in deltas_by_site.items():
        n = pep.count(site) if len(site) == 1 else 1
        shifts = sorted({s + k * delta for s in shifts for k in range(n + 1)})
    out, used = [], set()
    for i in np.argsort(hit)[::-1][:FRAGMENTS]:
        best = None
        for j, l in enumerate(lmz):
            if j in used:
                continue
            for c in (1, 2, 3):
                for d in shifts:
                    e = abs((hmz[i] - l) - d / c)
                    if e < 0.02 and (best is None or e < best[0]):
                        best = (e, j)
        if best:
            used.add(best[1])
            out.append((float(hmz[i]), float(lmz[best[1]]), float(hit[i])))
    return out


def arrays(block):
    mz = it = np.zeros(0)
    for a in ARRAY.findall(block):
        raw = base64.b64decode(BINARY.search(a).group(1))
        if b'MS:1000574' in a and raw:
            raw = zlib.decompress(raw)
        elif b'MS:1000576' not in a and b'MS:1000574' not in a and raw:
            sys.exit('quantify_standard_peak: unsupported binary compression (numpress?); convert to zlib or none')
        values = np.frombuffer(raw, '<f8' if b'MS:1000523' in a else '<f4')
        if b'MS:1000514' in a:
            mz = values
        elif b'MS:1000515' in a:
            it = values
    return mz, it


def spectra(path):
    buf = b''
    with open(path, 'rb') as fh:
        while True:
            data = fh.read(64 << 20)
            if not data:
                break
            buf += data
            parts = buf.split(b'</spectrum>')
            buf = parts.pop()
            yield from parts


def ms2_header(block):
    head = block[:block.find(b'<binaryDataArrayList')] if b'<binaryDataArrayList' in block else block
    m = LVL.search(head)
    if not m or m.group(1) != b'2':
        return None
    t, lo, up = TGT.search(head), LO.search(head), UP.search(head)
    if not (t and lo and up):
        return None
    c = float(t.group(1))
    rt = RT.search(head)
    sec = float(rt.group(1)) * (60 if rt.group(2) == b'UO:0000031' else 1) if rt else float('nan')
    return c, c - float(lo.group(1)), c + float(up.group(1)), sec


def smooth(x):
    if len(x) < 3 or x.shape[1] == 0:
        return x
    k = np.ones(3) / 3
    return np.column_stack([np.convolve(x[:, j], k, mode='same') for j in range(x.shape[1])])


def cosine(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))


def area(rt, y, lo, hi):
    m = (rt >= lo) & (rt <= hi)
    if m.sum() < 2:
        return np.zeros(y.shape[1])
    r, v = rt[m], y[m]
    return ((v[1:] + v[:-1]) / 2 * np.diff(r)[:, None]).sum(0)


def peak(Hrt, Hx, apex_index, max_half):
    hs = smooth(Hx).sum(1)
    a = apex_index
    for _ in range(3):
        a = max([i for i in (a - 1, a, a + 1) if 0 <= i < len(hs)], key=lambda i: hs[i])
    lo = hi = a
    while lo > 0 and hs[lo - 1] >= PEAK_FLOOR * hs[a] and Hrt[a] - Hrt[lo - 1] <= max_half:
        lo -= 1
    while hi < len(hs) - 1 and hs[hi + 1] >= PEAK_FLOOR * hs[a] and Hrt[hi + 1] - Hrt[a] <= max_half:
        hi += 1
    return a, lo, hi


def measure(Hrt, Hx, Lrt, Lx, apex_index, max_half, width):
    a, lo, hi = peak(Hrt, Hx, apex_index, max_half)
    t0, t1 = Hrt[lo], Hrt[hi]
    Li = np.column_stack([np.interp(Hrt, Lrt, Lx[:, j]) for j in range(Lx.shape[1])])
    ha, la = area(Hrt, Hx, t0, t1), area(Hrt, Li, t0, t1)
    out = dict(apex=Hrt[a], start=t0, stop=t1)
    if ha.max() <= 0:
        return dict(out, status='heavy without signal')
    keep = ha >= MIN_SHARE * ha.max()
    if keep.sum() < MIN_FRAGMENTS:
        return dict(out, status='too few heavy fragments')
    r = np.log2((la[keep] + 1) / (ha[keep] + 1))
    good = np.abs(r - np.median(r)) <= 1
    idx = np.flatnonzero(keep)[good] if good.sum() >= MIN_FRAGMENTS else np.flatnonzero(keep)
    margin = MARGIN_WIDTHS * width
    w = (Hrt >= t0 - margin) & (Hrt <= t1 + margin)
    hs_w = smooth(Hx[:, idx]).sum(1)[w]
    ls = smooth(Lx[:, idx]).sum(1)
    lw = np.interp(Hrt[w], Lrt, ls)
    co = float(np.corrcoef(lw, hs_w)[0, 1]) if lw.std() > 0 and hs_w.std() > 0 else 0.0
    null = []
    for k in NULL_SHIFTS:
        t = Hrt[w] + k * width
        if len(t) and t[0] >= Lrt[0] and t[-1] <= Lrt[-1]:
            v = np.interp(t, Lrt, ls)
            if v.std() > 0 and hs_w.std() > 0:
                null.append(float(np.corrcoef(v, hs_w)[0, 1]))
    out.update(status='measured', fragments=len(idx), fragments_dropped=int(keep.sum() - len(idx)),
               light_area=float(la[idx].sum()), heavy_area=float(ha[idx].sum()),
               ratio_L_H=float(la[idx].sum() / ha[idx].sum()), coelution=co,
               rdotp=cosine(np.sqrt(la[keep]), np.sqrt(ha[keep])),
               coelution_null=','.join('%.4f' % v for v in null))
    return out


def search_results(path, deltas):
    out = {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            m = PSM.match(r.get('PSMId', ''))
            if not m or 'decoy' in m.group(2):
                continue
            pep = re.sub(r'^-\.|\.-$', '', m.group(2))
            key = (light_form(pep, deltas), int(m.group(3)), 'heavy' if is_heavy(pep, deltas) else 'light')
            q = float(r['q-value'])
            if key not in out or q < out[key][0]:
                out[key] = (q, float(m.group(1)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', required=True)
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--library', required=True, help='.dlib or .elib holding both forms of each pair')
    ap.add_argument('--reference', required=True, help='reference list with both forms')
    ap.add_argument('--labels', required=True, help='JSON list of heavy labels with site and delta')
    ap.add_argument('--search-args', default='', help='its -ftol, -ftolunits and -foffset are used')
    ap.add_argument('--search', required=True, help="the run's identifications (the search's peptide table)")
    ap.add_argument('--fdr', type=float, default=0.01)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    labels = json.loads(a.labels)
    deltas = [float(r['delta']) for r in labels]
    by_site = {r['site']: float(r['delta']) for r in labels}
    tol, ppm, offset = fragment_options(a.search_args)
    rows = []
    if not os.path.basename(a.run).lower().endswith('.mzml'):
        print(f'{a.run}: not mzML; not measured')
        return write(a.out, rows)

    ref = read_reference(a.reference, deltas)
    pairs = {}
    for r in ref:
        pairs.setdefault((r['light_form'], r['z']), {})['heavy' if r['heavy'] else 'light'] = r
    pairs = {k: v for k, v in pairs.items() if 'heavy' in v and 'light' in v}
    lib = read_library(a.library, deltas, {stripped(f) for f, _ in pairs})
    frags = {}
    for (form, z), p in pairs.items():
        h, l = lib.get((form, z, True)), lib.get((form, z, False))
        frags[(form, z)] = match_fragments(h, l, stripped(form), by_site) if h is not None and l is not None else []

    want = {}
    for key, p in pairs.items():
        for form in ('heavy', 'light'):
            want.setdefault(round(p[form]['mz'], 2), []).append((key, form))
    windows, traces = {}, {}
    for block in spectra(a.run):
        m = SPECTRUM.search(block)
        if m is None:
            continue
        block = block[m.start():]
        s = ms2_header(block)
        if not s:
            continue
        c, lo, hi, sec = s
        w = (round(lo, 4), round(hi, 4))
        windows[w] = (lo, hi, c)
        hits = [(k, f) for r in (round(c - 0.01, 2), round(c, 2), round(c + 0.01, 2)) for k, f in want.get(r, [])
                if abs(pairs[k][f]['mz'] - c) <= OWN_WINDOW_TOLERANCE]
        if not hits:
            continue
        mz, it = arrays(block)
        if offset and len(mz):
            mz = mz * (1 - offset * 1e-6)
        order = np.argsort(mz)
        mz, it = mz[order], it[order]
        for k, f in hits:
            fm = np.array([x[0] if f == 'heavy' else x[1] for x in frags[k]])
            v = np.zeros(len(fm))
            if len(fm) and len(mz):
                d = fm * tol * 1e-6 if ppm else np.full(len(fm), tol)
                lo_i, hi_i = np.searchsorted(mz, fm - d), np.searchsorted(mz, fm + d)
                for j in range(len(fm)):
                    if hi_i[j] > lo_i[j]:
                        v[j] = it[lo_i[j]:hi_i[j]].max()
            traces.setdefault((k, f, w), []).append((sec, v))

    def own(mz):
        best = None
        for k, (lo, hi, target) in windows.items():
            if abs(target - mz) <= OWN_WINDOW_TOLERANCE and lo <= mz <= hi:
                rank = (abs(target - mz), hi - lo)
                if best is None or rank < best[0]:
                    best = (rank, k)
        return best[1] if best else None

    owners = {}
    for r in ref:
        r['window'] = own(r['mz'])
        owners.setdefault(r['window'], set()).add((r['light_form'], r['z'], r['heavy']))
    found = search_results(a.search, deltas)

    todo = []
    for key in sorted(pairs):
        form, z = key
        p = pairs[key]
        pep = stripped(form)
        row = dict(sample_id=a.sample_id, peptide=pep, peptidoform=p['light']['compound'], charge=z,
                   heavy_mz=p['heavy']['mz'], light_mz=p['light']['mz'])
        hq, hrt = found.get((form, z, 'heavy'), (float('nan'), float('nan')))
        lq, lrt = found.get((form, z, 'light'), (float('nan'), float('nan')))
        row.update(heavy_q_value=hq, light_q_value=lq, light_rt=lrt)
        if not (hq == hq and hq <= a.fdr and hrt == hrt):
            rows.append(dict(row, status='heavy not identified'))
            continue
        hw, lw = p['heavy']['window'], p['light']['window']
        if hw is None or lw is None:
            rows.append(dict(row, status='no window of its own'))
            continue
        if len(owners[hw]) > 1 or len(owners[lw]) > 1:
            rows.append(dict(row, status='window shared with another target'))
            continue
        if len(frags[key]) < MIN_FRAGMENTS:
            rows.append(dict(row, status='fewer than 3 fragments in both forms'))
            continue
        H, L = sorted(traces.get((key, 'heavy', hw), [])), sorted(traces.get((key, 'light', lw), []))
        if len(H) < 3 or len(L) < 3:
            rows.append(dict(row, status='window without scans'))
            continue
        Hrt, Hx = np.array([t for t, _ in H]), np.array([v for _, v in H])
        Lrt, Lx = np.array([t for t, _ in L]), np.array([v for _, v in L])
        todo.append((row, Hrt, Hx, Lrt, Lx, int(np.argmin(np.abs(Hrt - hrt)))))
        rows.append(row)
    widths = []
    for _row, Hrt, Hx, _Lrt, _Lx, i in todo:
        _a, lo, hi = peak(Hrt, Hx, i, float('inf'))
        widths.append(Hrt[hi] - Hrt[lo])
    width = float(np.median(widths)) if widths else float('nan')
    max_half = MAX_HALF_WIDTHS * width if widths else float('inf')
    for row, Hrt, Hx, Lrt, Lx, i in todo:
        row.update(measure(Hrt, Hx, Lrt, Lx, i, max_half, width))
    write(a.out, rows)
    n = sum(1 for r in rows if r.get('status') == 'measured')
    print(f'{a.sample_id}: {len(pairs)} pairs; measured at the heavy\'s peak {n}; the heavies\' median peak '
          f'width {width:.1f} s, so a peak ends at most {max_half:.1f} s from its apex and co-elution is read '
          f'{MARGIN_WIDTHS * width:.1f} s either side of it')
    return 0


def write(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(COLS)
        for r in rows:
            w.writerow([('%.6g' % r[c]) if isinstance(r.get(c), float) else r.get(c, '') for c in COLS])
    return 0


if __name__ == '__main__':
    sys.exit(main())
