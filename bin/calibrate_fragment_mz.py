#!/usr/bin/env python3
import argparse, base64, hashlib, json, os, re, sqlite3, sys, time, zlib
import numpy as np

CHUNK = 64 << 20
PROTON = 1.007276467
H2O = 18.010564684
AA = {'G': 57.02146372, 'A': 71.03711379, 'S': 87.03202841, 'P': 97.05276385, 'V': 99.06841391,
      'T': 101.04767847, 'C': 103.00918478, 'L': 113.08406398, 'I': 113.08406398, 'N': 114.04292744,
      'D': 115.02694303, 'Q': 128.05857751, 'K': 128.09496302, 'E': 129.04259309, 'M': 131.04048491,
      'H': 137.05891186, 'F': 147.06841391, 'R': 156.10111103, 'Y': 163.06332853, 'W': 186.07931295,
      'U': 150.95363559, 'O': 237.14772686}
KNOWN_DELTAS = [57.021464, 15.994915, 8.014199, 10.008269, 42.010565, 79.966331, 0.984016,
                -17.026549, -18.010565, 6.020129, 4.025107, 6.013809, 43.005814, 14.01565,
                28.0313, 114.042927, 229.162932, 304.207146, 144.102063, 13.979265, 31.989829]

R_LVL = re.compile(rb'accession="MS:1000511"[^>]*?value="([^"]*)"')
R_TGT = re.compile(rb'accession="MS:1000827"[^>]*?value="([^"]*)"')
R_LO = re.compile(rb'accession="MS:1000828"[^>]*?value="([^"]*)"')
R_UP = re.compile(rb'accession="MS:1000829"[^>]*?value="([^"]*)"')
R_RT = re.compile(rb'accession="MS:1000016"[^>]*?value="([^"]*)"[^>]*?unitAccession="([^"]*)"')
R_TIC = re.compile(rb'accession="MS:1000285"[^>]*?value="([^"]*)"')
R_IIT = re.compile(rb'accession="MS:1000927"[^>]*?value="([^"]*)"')
R_SWL = re.compile(rb'accession="MS:1000501"[^>]*?value="([^"]*)"')
R_SWU = re.compile(rb'accession="MS:1000500"[^>]*?value="([^"]*)"')
R_FILTER = re.compile(rb'accession="MS:1000512"[^>]*?value="([^"]*)"')
R_PREC_IN_FILTER = re.compile(rb'\s[0-9.]+@([a-z]+)[0-9.]*')
R_RANGE_IN_FILTER = re.compile(rb'\s*\[[^\]]*\]')
R_INDEX = re.compile(rb'(<spectrum\s[^>]*?\bindex=")(\d+)(")')
R_COUNT = re.compile(rb'(<spectrumList\s[^>]*?\bcount=")(\d+)(")')
R_BPMZ = re.compile(rb'(accession="MS:1000504"[^>]*?value=")([^"]*)(")')
R_LOWMZ = re.compile(rb'(accession="MS:1000528"[^>]*?value=")([^"]*)(")')
R_HIGHMZ = re.compile(rb'(accession="MS:1000527"[^>]*?value=")([^"]*)(")')
R_BDA = re.compile(rb'(<binaryDataArray\b[^>]*>)(.*?)(</binaryDataArray>)', re.S)
R_BIN = re.compile(rb'(<binary>)(.*?)(</binary>)', re.S)
R_ENCLEN = re.compile(rb'(\bencodedLength=")(\d+)(")')
SPECTRUM = re.compile(rb'<spectrum\s')
TAG_ID = re.compile(rb'\bid="([^"]*)"')
CHROM = re.compile(rb'<chromatogram\s[^>]*')
NUMPRESS = (b'MS:1002312', b'MS:1002313', b'MS:1002314', b'MS:1002746', b'MS:1002747', b'MS:1002748')


def stream_spectra(path):
    buf = b''
    pos = 0
    with open(path, 'rb') as fh:
        while True:
            data = fh.read(CHUNK)
            if not data:
                break
            buf += data
            parts = buf.split(b'</spectrum>')
            buf = parts.pop()
            for p in parts:
                m = SPECTRUM.search(p)
                if m is None:
                    raise ValueError('%s: </spectrum> without <spectrum' % path)
                yield p[:m.start()], p[m.start():] + b'</spectrum>', pos + m.start()
                pos += len(p) + len(b'</spectrum>')
    yield buf, None, pos


def scan_meta(block):
    k = block.find(b'<binaryDataArrayList')
    if k > 0:
        block = block[:k]

    def f(rx, default=float('nan')):
        m = rx.search(block)
        return float(m.group(1)) if m and m.group(1) else default
    m = R_LVL.search(block)
    lvl = int(m.group(1)) if m else 0
    rt = R_RT.search(block)
    sec = float(rt.group(1)) * (60.0 if rt.group(2) == b'UO:0000031' else 1.0) if rt else float('nan')
    return dict(level=lvl, rt=sec, target=f(R_TGT), lo=f(R_LO), up=f(R_UP), tic=f(R_TIC), group=scan_group(block),
                iit=f(R_IIT), swl=f(R_SWL), swu=f(R_SWU))


def scan_group(block):
    m = R_FILTER.search(block)
    return R_RANGE_IN_FILTER.sub(b'', R_PREC_IN_FILTER.sub(rb' @\1', m.group(1))).decode() if m else ''


def window_bounds(block):
    t, lo, up = R_TGT.search(block), R_LO.search(block), R_UP.search(block)
    if not (t and lo and up):
        return None
    c = float(t.group(1))
    return round(c - float(lo.group(1)), 2), round(c + float(up.group(1)), 2)


def _array_kind(params):
    if any(n in params for n in NUMPRESS):
        raise ValueError('MS-Numpress compressed arrays are not supported')
    if b'MS:1000514' in params:
        kind = 'mz'
    elif b'MS:1000515' in params:
        kind = 'int'
    else:
        kind = None
    dt = '<f8' if b'MS:1000523' in params else '<f4' if b'MS:1000521' in params else \
         '<i8' if b'MS:1000522' in params else '<i4' if b'MS:1000519' in params else None
    return kind, dt, b'MS:1000574' in params


def decode_arrays(block):
    out = {}
    for m in R_BDA.finditer(block):
        body = m.group(2)
        b = R_BIN.search(body)
        params = body[:b.start()] if b else body
        kind, dt, z = _array_kind(params)
        if kind is None or dt is None:
            continue
        raw = base64.b64decode(b.group(2)) if b else b''
        if z and raw:
            raw = zlib.decompress(raw)
        out[kind] = np.frombuffer(raw, dtype=dt).astype(np.float64)
    return out.get('mz', np.zeros(0)), out.get('int', np.zeros(0))


def read_block(fh, offset):
    fh.seek(offset)
    buf = b''
    while True:
        data = fh.read(1 << 16)
        if not data:
            raise ValueError('no </spectrum> after offset %d' % offset)
        buf += data
        k = buf.find(b'</spectrum>')
        if k >= 0:
            return buf[:k + len(b'</spectrum>')]


def snap(delta):
    best = min(KNOWN_DELTAS, key=lambda k: abs(k - delta))
    tol = 0.0006 if len(('%r' % delta).split('.')[-1]) >= 4 else 0.051
    return best if abs(best - delta) <= tol else delta


TOKEN = re.compile(r'([A-Z])((?:\[[+-]?[0-9.]+\]|\([+-]?[0-9.]+\))*)')
DELTA = re.compile(r'[\[(]([+-]?[0-9.]+)[\])]')


def residue_masses(modseq):
    s = modseq.strip()
    nterm = 0.0
    m = re.match(r'^n?((?:\[[+-]?[0-9.]+\])+)-?', s)
    if m:
        nterm = sum(snap(float(d)) for d in DELTA.findall(m.group(1)))
        s = s[m.end():]
    res = []
    pos = 0
    for t in TOKEN.finditer(s):
        if t.start() != pos:
            raise ValueError('cannot parse %r' % modseq)
        pos = t.end()
        mass = AA[t.group(1)] + sum(snap(float(d)) for d in DELTA.findall(t.group(2)))
        res.append(mass)
    if pos != len(s) or not res:
        raise ValueError('cannot parse %r' % modseq)
    res[0] += nterm
    return res


def fragment_table(res, max_charge=2):
    n = len(res)
    cum = np.cumsum(res)
    total = cum[-1]
    out = []
    for i in range(1, n):
        b = cum[i - 1]
        y = total - cum[i - 1] + H2O
        for z in range(1, max_charge + 1):
            out.append(('b%d' % i, z, (b + z * PROTON) / z))
            out.append(('y%d' % (n - i), z, (y + z * PROTON) / z))
    return out


def precursor_mz(res, z):
    return (sum(res) + H2O + z * PROTON) / z


def _blob(b, n, dt):
    if b is None:
        return np.zeros(0)
    if n and len(b) == n * np.dtype(dt).itemsize:
        return np.frombuffer(b, dtype=dt).astype(np.float64)
    return np.frombuffer(zlib.decompress(b), dtype=dt).astype(np.float64)


def read_library(path):
    con = sqlite3.connect('file:%s?mode=ro' % os.path.abspath(path), uri=True)
    tables = {r[0] for r in con.execute("select name from sqlite_master where type='table'")}
    out = {}
    if 'RefSpectra' in tables:
        q = ('select r.peptideModSeq, r.precursorCharge, r.precursorMZ, r.numPeaks, p.peakMZ, p.peakIntensity '
             'from RefSpectra r join RefSpectraPeaks p on p.RefSpectraID = r.id')
        for modseq, z, pmz, n, mzb, itb in con.execute(q):
            key = (modseq, int(z))
            if key not in out:
                out[key] = (modseq, int(z), float(pmz), _blob(mzb, n, '<f8'), _blob(itb, n, '<f4'))
    elif 'entries' in tables:
        cols = {r[1] for r in con.execute('pragma table_info(entries)')}
        dec = ', IsDecoy' if 'IsDecoy' in cols else ''
        q = 'select PeptideModSeq, PrecursorCharge, PrecursorMz, MassArray, IntensityArray%s from entries' % dec
        for row in con.execute(q):
            if dec and row[5]:
                continue
            modseq, z, pmz, mzb, itb = row[:5]
            key = (modseq, int(z))
            if key not in out:
                out[key] = (modseq, int(z), float(pmz),
                            np.frombuffer(zlib.decompress(mzb), '>f8').astype(np.float64),
                            np.frombuffer(zlib.decompress(itb), '>f4').astype(np.float64))
    else:
        raise ValueError('%s: neither a .blib nor an EncyclopeDIA library' % path)
    con.close()
    return list(out.values())


def annotate(modseq, z, lib_mz, lib_int, annot_tol=0.02):
    try:
        res = residue_masses(modseq)
    except (KeyError, ValueError):
        return None
    ions = fragment_table(res, max_charge=min(2, max(1, z)))
    tmz = np.array([x[2] for x in ions])
    order = np.argsort(tmz)
    tmz = tmz[order]
    labels = [ions[i] for i in order]
    out = []
    for mz, it in zip(lib_mz, lib_int):
        if it <= 0:
            continue
        j = np.searchsorted(tmz, mz)
        cand = [k for k in (j - 1, j) if 0 <= k < len(tmz) and abs(tmz[k] - mz) <= annot_tol]
        if len(cand) != 1:
            continue
        k = cand[0]
        out.append((tmz[k], it, '%s%s' % (labels[k][0], '+' * labels[k][1])))
    if not out:
        return None
    best = {}
    for t, it, lab in out:
        if lab not in best or it > best[lab][1]:
            best[lab] = (t, it, lab)
    out = sorted(best.values(), key=lambda r: -r[1])
    return precursor_mz(res, z), out


def lowess_grid(x, y, grid, frac=0.3, iters=3):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    k = int(min(n, max(10, round(frac * n))))
    rw = np.ones(n)
    fit = np.zeros(len(grid))
    for it in range(iters + 1):
        for j, g in enumerate(grid):
            d = np.abs(x - g)
            h = max(np.partition(d, k - 1)[k - 1], 1e-9) * 1.0001
            u = d / h
            w = np.where(u < 1, (1 - u ** 3) ** 3, 0.0) * rw
            sw = w.sum()
            if sw <= 0:
                fit[j] = np.nan
                continue
            xm = (w * x).sum() / sw
            ym = (w * y).sum() / sw
            sxx = (w * (x - xm) ** 2).sum()
            b = (w * (x - xm) * (y - ym)).sum() / sxx if sxx > 1e-12 else 0.0
            fit[j] = ym + b * (g - xm)
        if np.isnan(fit).any():
            ok = ~np.isnan(fit)
            fit = np.interp(grid, grid[ok], fit[ok])
        if it == iters:
            break
        r = y - np.interp(x, grid, fit)
        s = np.median(np.abs(r))
        if s <= 0:
            break
        u = r / (6 * s)
        rw = np.where(np.abs(u) < 1, (1 - u ** 2) ** 2, 0.0)
    return fit, rw


def make_grid(x, step):
    lo, hi = np.percentile(x, [1, 99])
    lo, hi = np.floor(lo / step) * step, np.ceil(hi / step) * step
    return np.arange(lo, hi + step / 2, step)


class Curve:
    def __init__(self, mz_grid, mz_off, rt_grid=None, rt_off=None):
        self.mz_grid = np.asarray(mz_grid, float)
        self.mz_off = np.asarray(mz_off, float)
        self.rt_grid = None if rt_grid is None else np.asarray(rt_grid, float)
        self.rt_off = None if rt_off is None else np.asarray(rt_off, float)

    def __call__(self, mz, rt=None):
        off = np.interp(mz, self.mz_grid, self.mz_off)
        if self.rt_grid is not None and rt is not None:
            off = off + np.interp(rt, self.rt_grid, self.rt_off)
        return off

    def to_dict(self):
        d = dict(mz_grid=[float(v) for v in self.mz_grid], mz_offset_th=[float(v) for v in self.mz_off])
        if self.rt_grid is not None:
            d.update(rt_grid_s=[float(v) for v in self.rt_grid], rt_offset_th=[float(v) for v in self.rt_off])
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(d['mz_grid'], d['mz_offset_th'], d.get('rt_grid_s'), d.get('rt_offset_th'))


def strip_seq(s):
    return re.sub(r'[^A-Z]', '', re.sub(r'\[[^\]]*\]|\([^)]*\)', '', s))


def read_reference_list(path):
    import csv
    lines = open(path, newline='').read().splitlines()
    rows = list(csv.reader(lines, delimiter='\t' if lines[0].count('\t') >= lines[0].count(',') else ','))
    header = [h.strip().strip('"').lower() for h in rows[0]]
    dcol = next((i for i, h in enumerate(header) if h in ('isdecoy', 'decoy', 'is_decoy')), None)
    out = {}
    for r in rows[1:]:
        if len(r) < 5 or not r[0].strip():
            continue
        if dcol is not None and len(r) > dcol and r[dcol].strip().lower() in ('true', '1', 'yes', 'decoy'):
            continue
        out.setdefault(int(float(r[4])), []).append((float(r[3]), strip_seq(r[0])))
    return out


def in_reference(refs, modseq, z, pmz, tol=0.01):
    for mz, seq in refs.get(int(z), []):
        if abs(mz - pmz) <= tol and (not seq or seq == strip_seq(modseq)):
            return True
    return False


def th(a, f, mz):
    return round(f * a.fragment_tolerance, 12) * (np.asarray(mz, float) * 1e-6 if a.fragment_unit == 'ppm' else 1.0)


def tolerance_counts(u, a):
    edges = np.linspace(0, a.tol, int(round(a.tol / a.tolerance_bin)) + 1)
    d = np.abs(u[np.isfinite(u)])
    return edges[1:], np.histogram(d[d <= a.tol], edges)[0]


def tol_text(a, f):
    return '%g %s' % (round(f * a.fragment_tolerance, 12), 'Th' if a.fragment_unit == 'th' else 'ppm')


def window_key(target):
    return round(target, 4)


def find_calibrants(run, library, a, log):
    lib = read_library(library)
    log('library: %d entries' % len(lib))
    t0 = time.time()
    scans = []
    for _prefix, block, off in stream_spectra(run):
        if block is None:
            break
        m = scan_meta(block)
        if m['level'] != 2 or not np.isfinite(m['target']):
            continue
        width = (m['lo'] if np.isfinite(m['lo']) else 0) + (m['up'] if np.isfinite(m['up']) else 0)
        if width > a.max_window_width:
            continue
        m['offset'] = off
        scans.append(m)
    log('pass 1: %d narrow-window MS2 scans (%.0f s)' % (len(scans), time.time() - t0))
    if not scans:
        return [], dict(narrow_scans=0), []
    centres = np.array(sorted({window_key(s['target']) for s in scans}))
    refs = read_reference_list(a.reference_list) if a.reference_list else None
    if refs is not None:
        log('reference list: %d target precursors' % sum(len(v) for v in refs.values()))
    targets = []
    for modseq, z, pmz, lmz, lint in lib:
        if refs is not None and not in_reference(refs, modseq, z, pmz):
            continue
        j = np.searchsorted(centres, pmz)
        near = [centres[k] for k in (j - 1, j) if 0 <= k < len(centres) and abs(centres[k] - pmz) <= a.center_tol]
        if not near:
            continue
        ann = annotate(modseq, z, lmz, lint, a.annot_tol)
        if ann is None:
            continue
        theo_prec, frags = ann
        frags = [f for f in frags if abs(f[0] - pmz) > 1.5]
        if len(frags) < a.min_frags:
            continue
        targets.append(dict(modseq=modseq, z=z, prec=pmz, key=min(near, key=lambda c: abs(c - pmz)), frags=frags))
    claims = {}
    for i, t in enumerate(targets):
        claims.setdefault(t['key'], []).append(i)
    own = []
    for key, ids in claims.items():
        best = min(abs(targets[i]['prec'] - key) for i in ids)
        ids = [i for i in ids if abs(targets[i]['prec'] - key) <= best + 1e-4]
        if len({(strip_seq(targets[i]['modseq']), targets[i]['z']) for i in ids}) == 1:
            own.append(ids[0])
    left_out = len(targets) - len(own)
    targets = [targets[i] for i in sorted(own)]
    log('library targets at a window centre (+/-%.2f Th) with >= %d annotated fragments: %d%s'
        % (a.center_tol, a.min_frags, len(targets), '; %d left out: their window is centred on another '
           'target, or on an isobaric one' % left_out if left_out else ''))
    by_key = {}
    for i, t in enumerate(targets):
        by_key.setdefault(t['key'], []).append(i)
    t0 = time.time()
    per_target = {i: [] for i in range(len(targets))}
    by_offset = {s['offset']: si for si, s in enumerate(scans) if window_key(s['target']) in by_key}
    for _prefix, block, off in stream_spectra(run):
        if block is None:
            break
        si = by_offset.get(off)
        if si is None:
            continue
        s = scans[si]
        ids = by_key[window_key(s['target'])]
        mz, it = decode_arrays(block)
        if len(mz) == 0:
            continue
        o = np.argsort(mz, kind='stable')
        mz, it = mz[o], it[o]
        base = it.max()
        for i in ids:
            t = targets[i]
            fm = np.array([f[0] for f in t['frags'][:a.top_k]])
            d = th(a, a.tol, fm)
            lo = np.searchsorted(mz, fm - d)
            hi = np.searchsorted(mz, fm + d)
            inten = np.zeros(len(fm))
            for k in range(len(fm)):
                if hi[k] > lo[k]:
                    inten[k] = it[lo[k]:hi[k]].max()
            per_target[i].append((s['rt'], si, inten, base))
    log('scored %d targets in their windows (%.0f s)' % (len(targets), time.time() - t0))
    chosen, tstats = [], []
    for i, rows in per_target.items():
        t = targets[i]
        st = dict(modseq=t['modseq'], z=t['z'], prec=t['prec'], scans=len(rows), used=False)
        tstats.append(st)
        if len(rows) < 5:
            continue
        rows.sort(key=lambda r: r[0])
        k = min(a.top_k, len(t['frags']))
        L = np.sqrt(np.array([f[1] for f in t['frags'][:k]]))
        I = np.array([r[2] for r in rows])
        X = np.sqrt(I)
        Xc = X - X.mean(axis=1, keepdims=True)
        Lc = L - L.mean()
        den = np.linalg.norm(Xc, axis=1) * np.linalg.norm(Lc)
        r_lib = np.where(den > 0, Xc @ Lc / np.where(den > 0, den, 1), 0.0)
        nfr = (I > 0).sum(axis=1)
        tot = I.sum(axis=1)
        tot3 = np.convolve(tot, np.ones(3) / 3, mode='same')
        cand = (r_lib >= a.min_cos_neighbour) & (nfr >= min(a.min_frags, k))
        if not cand.any():
            continue
        j = int(np.argmax(np.where(cand, tot3, -1)))
        lo, hi = max(0, j - 3), min(len(rows), j + 4)
        W = I[lo:hi]
        T = W.sum(axis=1)
        cors = []
        for kk in range(k):
            w = W[:, kk]
            if w.std() > 0 and T.std() > 0:
                cors.append(np.corrcoef(w, T - w)[0, 1] if (T - w).std() > 0 else 0.0)
        xic = float(np.median(cors)) if cors else 0.0
        med = float(np.median(tot))
        peak_ratio = float(tot3[j] / med) if med > 0 else float('inf')
        st.update(apex_rt=round(rows[j][0], 2), r_lib=round(float(r_lib[j]), 3), xic_corr=round(xic, 3),
                  n_frags=int(nfr[j]), apex_intensity=float(tot[j]), peak_ratio=round(min(peak_ratio, 1e6), 2))
        if not (r_lib[j] >= a.min_cos and xic >= a.min_xic and peak_ratio >= a.min_peak_ratio and
                0 < j < len(rows) - 1):
            continue
        st['used'] = True
        chosen.append(dict(target=i, apex_rt=rows[j][0], cos=float(r_lib[j]),
                           scans=[rows[jj][1] for jj in (j - 1, j, j + 1)]))
    log('targets with a co-eluting apex (library correlation >= %.2f, fragment co-elution >= %.2f): %d of %d'
        % (a.min_cos, a.min_xic, len(chosen), len(targets)))
    pts = []
    with open(run, 'rb') as fh:
        for c in chosen:
            t = targets[c['target']]
            for rank, si in enumerate(c['scans']):
                s = scans[si]
                mz, it = decode_arrays(read_block(fh, s['offset']))
                if len(mz) == 0:
                    continue
                o = np.argsort(mz, kind='stable')
                mz, it = mz[o], it[o]
                strongest = 0.0
                rows = []
                for fmz, lint, lab in t['frags']:
                    if np.isfinite(s['swl']) and (fmz < s['swl'] + 0.5 or fmz > s['swu'] - 0.5):
                        continue
                    d = th(a, a.tol, fmz)
                    lo, hi = np.searchsorted(mz, fmz - d), np.searchsorted(mz, fmz + d)
                    if hi <= lo:
                        continue
                    jj = lo + int(np.argmax(it[lo:hi]))
                    rows.append((fmz, lint, lab, mz[lo:hi].copy(), it[lo:hi].copy(), mz[jj], it[jj]))
                    strongest = max(strongest, it[jj])
                for fmz, lint, lab, cmz, cit, omz, oit in rows:
                    if oit < a.min_rel * strongest:
                        continue
                    pts.append(dict(group=s['group'], modseq=t['modseq'], z=t['z'], prec=t['prec'], ion=lab, theo=fmz,
                                    lib_int=lint, rt=s['rt'], apex=rank == 1, cand_mz=cmz, cand_int=cit,
                                    obs=omz, inten=oit, rel=oit / strongest, tic=s['tic'], iit=s['iit']))
    info = dict(library_entries=len(lib), narrow_scans=len(scans), targets_at_window_centre=len(targets),
                targets_with_apex=len(chosen))
    return pts, info, tstats


def fit_curve(pts, a, log, use_rt):
    x = np.array([p['theo'] for p in pts])
    rt = np.array([p['rt'] for p in pts])
    e1 = np.array([p['obs'] - p['theo'] for p in pts])
    grid = make_grid(x, a.grid_step)
    f1, _ = lowess_grid(x, e1, grid, frac=a.frac)
    pred = np.interp(x, grid, f1)
    key = np.array(['%s/%d' % (p['modseq'], p['z']) for p in pts])
    r1 = np.abs(e1 - pred)
    per = {k: float(np.median(r1[key == k])) for k in np.unique(key)}
    bad = {k for k, v in per.items() if v > th(a, a.max_target_residual, float(np.median(x[key == k])))}
    stats = dict(targets=len(per), inconsistent=len(bad), share_inconsistent=round(len(bad) / max(1, len(per)), 3),
                 target_median_abs_residual=round(float(np.median(list(per.values()))), 4))
    log('consistency with the stage-1 curve: %d of %d targets have a median |residual| above %s (%.0f %%); '
        'median over targets %.3f Th' % (len(bad), len(per), tol_text(a, a.max_target_residual),
                                         100 * stats['share_inconsistent'], stats['target_median_abs_residual']))
    e2 = np.full(len(pts), np.nan)
    for i, p in enumerate(pts):
        if key[i] in bad:
            continue
        want = p['theo'] + pred[i]
        d = np.abs(p['cand_mz'] - want)
        okc = (d <= th(a, a.rematch_tol, p['theo'])) & (p['cand_int'] >= 0.2 * p['cand_int'].max())
        if okc.any():
            k = np.argmin(np.where(okc, d, np.inf))
            e2[i] = p['cand_mz'][k] - p['theo']
            p['obs2'] = p['cand_mz'][k]
    keep = np.isfinite(e2)
    log('stage 1 (most intense within +/-%s): %d errors; stage 2 (nearest to the stage-1 curve within '
        '+/-%s): %d' % (tol_text(a, a.tol), len(x), tol_text(a, a.rematch_tol), keep.sum()))
    xk, rk, ek = x[keep], rt[keep], e2[keep]
    f2, _ = lowess_grid(xk, ek, grid, frac=a.frac)
    if not use_rt:
        return Curve(grid, f2), e2, Curve(grid, f1), stats
    rgrid = make_grid(rk, a.rt_grid_step)
    g = np.zeros(len(rgrid))
    for _ in range(2):
        g, _w = lowess_grid(rk, ek - np.interp(xk, grid, f2), rgrid, frac=a.rt_frac)
        g -= np.median(np.interp(rk, rgrid, g))
        f2, _w = lowess_grid(xk, ek - np.interp(rk, rgrid, g), grid, frac=a.frac)
    return Curve(grid, f2, rgrid, g), e2, Curve(grid, f1), stats


def summarize(x, err, bins=(200, 400, 600, 800, 1000, 1200, 1500, 2000)):
    rows = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        s = err[(x > lo) & (x <= hi) & np.isfinite(err)]
        if len(s) == 0:
            continue
        q = np.percentile(s, [25, 50, 75])
        rows.append(dict(bin='%d-%d' % (lo, hi), n=int(len(s)), median=round(float(q[1]), 4),
                         q25=round(float(q[0]), 4), q75=round(float(q[2]), 4),
                         share_gt_0_2=round(float((np.abs(s) > 0.2).mean()), 3)))
    return rows


def _encode_like(params, values):
    kind, dt, z = _array_kind(params)
    raw = values.astype(dt).tobytes()
    if z:
        raw = zlib.compress(raw)
    return base64.b64encode(raw)


def recalibrate_block(block, curve, stats, groups):
    m = R_LVL.search(block)
    if not m or m.group(1) != b'2':
        return block
    if scan_group(block[:block.find(b'<binaryDataArrayList')]) not in groups:
        stats['other_ms2'] += 1
        return block
    rt = R_RT.search(block)
    sec = float(rt.group(1)) * (60.0 if rt.group(2) == b'UO:0000031' else 1.0) if rt else None
    out, pos, changed = [], 0, False
    for bda in R_BDA.finditer(block):
        body = bda.group(2)
        b = R_BIN.search(body)
        if b is None:
            continue
        params = body[:b.start()]
        kind, dt, z = _array_kind(params)
        if kind != 'mz' or dt is None:
            continue
        raw = base64.b64decode(b.group(2))
        if z and raw:
            raw = zlib.decompress(raw)
        mz = np.frombuffer(raw, dtype=dt).astype(np.float64)
        if len(mz) == 0:
            continue
        shift = curve(mz, sec)
        new = mz - shift
        if np.array_equal(mz.astype(np.float32).astype(np.float64), mz):
            new = new.astype(np.float32).astype(np.float64)
        stats['peaks'] += len(mz)
        stats['max_abs_shift'] = max(stats['max_abs_shift'], float(np.abs(shift).max()))
        enc = _encode_like(params, new)
        head = R_ENCLEN.sub(lambda x: x.group(1) + str(len(enc)).encode() + x.group(3), bda.group(1), count=1)
        newbody = body[:b.start()] + b'<binary>' + enc + b'</binary>' + body[b.end():]
        out.append(block[pos:bda.start()] + head + newbody + bda.group(3))
        pos = bda.end()
        changed = True
    if not changed:
        return block
    out.append(block[pos:])
    block = b''.join(out)

    def fix(rx, blk):
        mm = rx.search(blk)
        if not mm or not mm.group(2):
            return blk
        v = float(mm.group(2))
        nv = v - float(curve(np.array([v]), sec)[0])
        return blk[:mm.start()] + mm.group(1) + repr(nv).encode() + mm.group(3) + blk[mm.end():]
    for rx in (R_BPMZ, R_LOWMZ, R_HIGHMZ):
        block = fix(rx, block)
    stats['spectra'] += 1
    return block


def add_processing(header, note):
    sw = (b'<software id="teapot_fragment_mz_calibration" version="0.1">\n'
          b'        <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" '
          b'value="calibrate_fragment_mz.py"/>\n      </software>\n    ')
    m = re.search(rb'<softwareList\s+count="(\d+)">', header)
    if m:
        end = header.find(b'</softwareList>', m.end())
        header = (header[:m.start()] + b'<softwareList count="%d">' % (int(m.group(1)) + 1) +
                  header[m.end():end] + b'  ' + sw.rstrip() + b'\n    ' + header[end:])
    dp = re.search(rb'defaultDataProcessingRef="([^"]+)"', header)
    if dp:
        m2 = re.search(rb'<dataProcessing\s+id="%s">' % re.escape(dp.group(1)), header)
    else:
        m2 = re.search(rb'<dataProcessing\s+id="[^"]*">', header)
    if m2:
        end = header.find(b'</dataProcessing>', m2.end())
        orders = [int(o) for o in re.findall(rb'<processingMethod\s+order="(\d+)"', header[m2.end():end])]
        pm = (b'  <processingMethod order="%d" softwareRef="teapot_fragment_mz_calibration">\n'
              b'          <cvParam cvRef="MS" accession="MS:1001485" name="m/z calibration" value=""/>\n'
              b'          <userParam name="fragment m/z recalibration" value="%s"/>\n'
              b'        </processingMethod>\n      ' % (max(orders + [-1]) + 1, note.encode()))
        header = header[:end] + pm + header[end:]
    return header


def write_calibrated(src, dst, curve, note, groups, drop=frozenset(), n_drop=0):
    sha = hashlib.sha1()
    os.makedirs(os.path.dirname(dst) or '.', exist_ok=True)
    tmp = dst + '.part'
    out = open(tmp, 'wb')
    pos = 0

    def write(b):
        nonlocal pos
        out.write(b)
        sha.update(b)
        pos += len(b)

    stats = dict(spectra=0, peaks=0, max_abs_shift=0.0, other_ms2=0, dropped=0)
    offsets, first, tail, index = [], True, b'', 0
    for prefix, block, _off in stream_spectra(src):
        if block is None:
            tail = prefix
            break
        if first:
            if drop:
                prefix = R_COUNT.sub(lambda x: x.group(1) + str(int(x.group(2)) - n_drop).encode() + x.group(3),
                                     prefix, count=1)
            if curve is not None:
                prefix = add_processing(prefix, note)
            first = False
        write(prefix)
        if drop:
            head = block[:block.find(b'<binaryDataArrayList')]
            m = R_LVL.search(head)
            if m and m.group(1) == b'2' and window_bounds(head) in drop:
                stats['dropped'] += 1
                continue
            block = R_INDEX.sub(lambda x: x.group(1) + str(index).encode() + x.group(3), block, count=1)
            index += 1
        open_tag = block[:block.index(b'>') + 1]
        offsets.append((TAG_ID.search(open_tag).group(1), pos))
        write(recalibrate_block(block, curve, stats, groups) if curve is not None else block)
    cut = tail.find(b'<indexList')
    if cut < 0:
        write(tail)
        out.close()
        os.replace(tmp, dst)
        return stats
    body = tail[:cut]
    chrom = [(TAG_ID.search(t.group(0)).group(1), pos + t.start()) for t in CHROM.finditer(body)]
    write(body)
    index_at = pos
    lines = [b'<indexList count="%d">\n' % (2 if chrom else 1), b'    <index name="spectrum">\n']
    lines += [b'      <offset idRef="%s">%d</offset>\n' % (i, o) for i, o in offsets]
    lines.append(b'    </index>\n')
    if chrom:
        lines.append(b'    <index name="chromatogram">\n')
        lines += [b'      <offset idRef="%s">%d</offset>\n' % (i, o) for i, o in chrom]
        lines.append(b'    </index>\n')
    lines.append(b'  </indexList>\n  <indexListOffset>%d</indexListOffset>\n  <fileChecksum>' % index_at)
    write(b''.join(lines))
    out.write(sha.hexdigest().encode() + b'</fileChecksum>\n</indexedmzML>\n')
    out.close()
    os.replace(tmp, dst)
    return stats


def survey_windows(path):
    import csv
    if not path:
        return frozenset(), 0
    with open(path, newline='') as fh:
        rows = [r for r in csv.DictReader(fh, delimiter='\t') if r.get('survey') == 'yes']
    return (frozenset((round(float(r['isolation_lower']), 2), round(float(r['isolation_upper']), 2)) for r in rows),
            sum(int(r['scans']) for r in rows))


def as_acquired(a, survey, log=print):
    drop, n_drop = survey
    if os.path.lexists(a.out):
        os.remove(a.out)
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    if drop and a.run.lower().endswith('.mzml'):
        st = write_calibrated(a.run, a.out, None, '', set(), drop, n_drop)
        log('%s: the run without its %d survey scans' % (a.out, st['dropped']))
        return st
    os.symlink(os.path.realpath(a.run), a.out)
    log('%s: a link to the run as it is' % a.out)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', required=True)
    ap.add_argument('--library', required=True, help='.blib, .dlib or .elib holding the targets')
    ap.add_argument('--report', required=True, help='output prefix for .json, .txt, .calibrants.tsv')
    ap.add_argument('--out', help='calibrated mzML (omit to measure only)')
    ap.add_argument('--fragment-tolerance', type=float, default=0.4, help="the search's fragment tolerance")
    ap.add_argument('--fragment-unit', choices=('th', 'ppm'), default='th')
    ap.add_argument('--tol', type=float, default=1.5,
                    help='fragment tolerance for the search of calibrants, in fragment tolerances')
    ap.add_argument('--rematch-tol', type=float, default=0.5,
                    help='tolerance around the stage-1 curve, in fragment tolerances')
    ap.add_argument('--center-tol', type=float, default=0.05, help='precursor to window centre (Th)')
    ap.add_argument('--max-window-width', type=float, default=4.0, help='wider isolation windows are ignored (Th)')
    ap.add_argument('--annot-tol', type=float, default=0.02, help='library peak to b/y ion (Th)')
    ap.add_argument('--top-k', type=int, default=8, help='strongest library fragments used for the apex')
    ap.add_argument('--min-frags', type=int, default=4)
    ap.add_argument('--min-cos', type=float, default=0.7, help='library correlation at the apex')
    ap.add_argument('--min-cos-neighbour', type=float, default=0.5, help='library correlation of apex candidates')
    ap.add_argument('--min-xic', type=float, default=0.7, help='fragment co-elution around the apex')
    ap.add_argument('--max-target-residual', type=float, default=0.375,
                    help='a target whose fragments are off the stage-1 curve by more than this (median, in fragment '
                         'tolerances) is dropped')
    ap.add_argument('--max-share-inconsistent', type=float, default=0.3,
                    help='more dropped targets than this share: the run is not calibrated')
    ap.add_argument('--min-peak-ratio', type=float, default=3.0,
                    help='apex (3-scan mean of the summed fragments) over the median of the window')
    ap.add_argument('--reference-list', help='only these precursors are calibrants (EncyclopeDIA/Skyline '
                    'reference list: Compound, m/z, z[, isDecoy]); default: every library entry at a window centre')
    ap.add_argument('--min-rel', type=float, default=0.05, help='fragment at least this share of the strongest one')
    ap.add_argument('--frac', type=float, default=0.2, help='LOWESS span (share of points)')
    ap.add_argument('--grid-step', type=float, default=20.0)
    ap.add_argument('--no-rt-term', action='store_true', help='m/z term only')
    ap.add_argument('--rt-min-targets', type=int, default=50, help='fewer targets: m/z term only')
    ap.add_argument('--rt-frac', type=float, default=0.3)
    ap.add_argument('--rt-grid-step', type=float, default=30.0)
    ap.add_argument('--min-offset', default='auto',
                    help='correct only when the largest offset of the m/z curve is above this (Th); 0: always; '
                         'auto: a quarter of the fragment tolerance')
    ap.add_argument('--measure-only', action='store_true',
                    help='measure and report, locate the targets, but never correct (the run passes as a link)')
    ap.add_argument('--tolerance-counts', help="the calibrants' fragment errors the search will see, as counts by "
                    "absolute error in fragment tolerances, for the search tolerance")
    ap.add_argument('--tolerance-bin', type=float, default=0.0125, help='bin of those counts, in fragment tolerances')
    ap.add_argument('--min-targets', type=int, default=20)
    ap.add_argument('--min-points', type=int, default=150)
    ap.add_argument('--windows', help="REMOVE_SURVEY_SCANS's isolation-window table: a run without narrow "
                    "targeted windows passes through without being read, and the scans of the windows it "
                    "marks as survey windows are left out of the file written")
    ap.add_argument('--sample-id', help='for the summary row and the warnings')
    ap.add_argument('--summary', help='one-row TSV for the per-run table')
    ap.add_argument('--warnings', help='lines for warnings.txt')
    a = ap.parse_args()
    return run(a)


SUMMARY = ['sample_id', 'status', 'corrected', 'calibrants', 'share_disagreeing', 'offset_300', 'offset_500',
           'offset_700', 'offset_900', 'offset_1100', 'offset_1300', 'largest_offset', 'rt_term_min',
           'rt_term_max', 'median_abs_residual_after', 'scan_types', 'acquired']
STAMP = re.compile(rb'startTimeStamp="([^"]+)"')


def acquired(path):
    if not path.lower().endswith('.mzml'):
        return ''
    with open(path, 'rb') as fh:
        m = STAMP.search(fh.read(1 << 20))
    return m.group(1).decode() if m else ''


def finish(a, report, status, corrected):
    sid = a.sample_id or os.path.basename(a.run)
    if a.summary:
        off = report.get('offset_at', {})
        rt = report.get('curve', {}).get('rt_offset_th') if corrected else None
        row = dict(sample_id=sid, status=status, corrected='yes' if corrected else 'no',
                   calibrants=report.get('targets_with_apex', ''),
                   share_disagreeing=report.get('consistency', {}).get('share_inconsistent', ''),
                   largest_offset=report.get('largest_offset_th', ''),
                   rt_term_min=round(min(rt), 4) if rt else '', rt_term_max=round(max(rt), 4) if rt else '',
                   median_abs_residual_after=report.get('median_abs_residual_th', ''),
                   scan_types='; '.join(report.get('scan_types_calibrated', [])),
                   acquired=acquired(a.run))
        for v in (300, 500, 700, 900, 1100, 1300):
            row['offset_%d' % v] = off.get(str(v), '')
        with open(a.summary, 'w') as fh:
            fh.write('\t'.join(SUMMARY) + '\n' + '\t'.join(str(row[c]) for c in SUMMARY) + '\n')
    if a.warnings:
        with open(a.warnings, 'w') as fh:
            if corrected:
                fh.write('%s: fragment m/z read up to %+.2f Th off; corrected (qc/%s/). '
                         'The instrument\'s calibration for this scan type may need attention\n'
                         % (sid, report['largest_offset_signed_th'], sid))
            elif status.startswith('not calibrated: the calibrants disagree'):
                fh.write('%s: fragment m/z not calibrated, %s\n' % (sid, status[len('not calibrated: '):]))


def passthrough(a, status, report):
    print(status, flush=True)
    if a.out:
        as_acquired(a, survey_windows(a.windows))
    report['status'] = status
    with open(a.report + '.json', 'w') as fh:
        json.dump(report, fh, indent=1, default=float)
    with open(a.report + '.txt', 'w') as fh:
        fh.write(status + '\n')
    finish(a, report, status, False)
    return 0


def narrow_windows(path, max_width):
    import csv
    with open(path, newline='') as fh:
        return any(float(r['width']) <= max_width and r.get('survey', 'no') != 'yes'
                   for r in csv.DictReader(fh, delimiter='\t'))


def run(a):

    lines = []

    def log(msg):
        lines.append(msg)
        print(msg, flush=True)

    name = os.path.basename(a.run)
    if not name.lower().endswith('.mzml'):
        return passthrough(a, 'not applicable: not mzML', dict(run=name))
    if a.windows and not narrow_windows(a.windows, a.max_window_width):
        return passthrough(a, 'not applicable: no isolation window narrower than %g Th (DIA)'
                           % a.max_window_width, dict(run=name))
    log('%s: fragment m/z calibration' % name)
    t0 = time.time()
    pts, info, tstats = find_calibrants(a.run, a.library, a, log)
    if not info['narrow_scans']:
        return passthrough(a, 'not applicable: no isolation window narrower than %g Th (DIA)'
                           % a.max_window_width, dict(run=name))
    if tstats:
        cols = ['modseq', 'z', 'prec', 'scans', 'used', 'apex_rt', 'r_lib', 'xic_corr', 'n_frags', 'apex_intensity',
                'peak_ratio']
        with open(a.report + '.targets.tsv', 'w') as fh:
            fh.write('\t'.join(cols) + '\n')
            for st in tstats:
                fh.write('\t'.join(str(st.get(c, 'NA')) for c in cols) + '\n')
    report = dict(run=name, library=os.path.basename(a.library), settings={k: v for k, v in vars(a).items()},
                  **info, calibrant_errors=len(pts))
    status = 'calibrated'
    curve = fitted = None
    if len(pts) == 0 or info.get('targets_with_apex', 0) < a.min_targets or len(pts) < a.min_points:
        status = 'not calibrated: too few calibrants (%d targets, %d fragment errors; need %d and %d)' % (
            info.get('targets_with_apex', 0), len(pts), a.min_targets, a.min_points)
        log(status)
    else:
        use_rt = not a.no_rt_term and info['targets_with_apex'] >= a.rt_min_targets
        curve, e2, c1, cstats = fit_curve(pts, a, log, use_rt)
        fitted = curve
        report['consistency'] = cstats
        x = np.array([p['theo'] for p in pts])
        rt = np.array([p['rt'] for p in pts])
        resid = e2 - curve(x, rt)
        groups = sorted({p['group'] for p in pts})
        report['scan_types_calibrated'] = groups
        report['rt_term'] = use_rt
        log('scan type(s) of the calibrants: %s' % ('; '.join(groups) or '(no filter strings)'))
        report['curve'] = curve.to_dict()
        report['curve_stage1'] = c1.to_dict()
        report['before'] = summarize(x, e2)
        report['after_fit'] = summarize(x, resid)
        grid = curve.mz_grid
        report['offset_at'] = {str(v): round(float(curve(np.array([float(v)]))[0]), 4)
                               for v in (300, 500, 700, 900, 1100, 1300)}
        report['median_abs_residual_th'] = round(float(np.nanmedian(np.abs(resid))), 4)
        log('offset curve (Th) at m/z 300/500/700/900/1100/1300: %s' % ' / '.join(
            '%+.3f' % v for v in report['offset_at'].values()))
        log('calibrant m/z range %.0f-%.0f (constant outside); median |residual| %.3f Th'
            % (grid[0], grid[-1], report['median_abs_residual_th']))
        if use_rt:
            log('retention-time term: %+.3f to %+.3f Th (%.1f-%.1f min, constant outside)' % (
                curve.rt_off.min(), curve.rt_off.max(), curve.rt_grid[0] / 60, curve.rt_grid[-1] / 60))
        log('before (median error by fragment m/z): ' + '; '.join(
            '%s %+.3f' % (r['bin'], r['median']) for r in report['before']))
        log('after the fit (residual): ' + '; '.join('%s %+.3f' % (r['bin'], r['median']) for r in report['after_fit']))
        with open(a.report + '.calibrants.tsv', 'w') as fh:
            fh.write('modseq\tz\tprecursor_mz\tion\ttheo_mz\tlib_intensity\trt_s\tapex\tobs_mz_stage1\tobs_mz\t'
                     'intensity\trel\ttic\tinjection_ms\terror_th\tresidual_th\n')
            for p, e, r in zip(pts, e2, resid):
                fh.write('%s\t%d\t%.4f\t%s\t%.4f\t%.4g\t%.2f\t%d\t%.4f\t%s\t%.4g\t%.3f\t%.4g\t%.4g\t%s\t%s\n' % (
                    p['modseq'], p['z'], p['prec'], p['ion'], p['theo'], p['lib_int'], p['rt'], p['apex'],
                    p['obs'], '%.4f' % p['obs2'] if 'obs2' in p else 'NA', p['inten'], p['rel'], p['tic'],
                    p['iit'], '%.4f' % e if np.isfinite(e) else 'NA', '%.4f' % r if np.isfinite(r) else 'NA'))
    if curve is not None and report['consistency']['targets'] - report['consistency']['inconsistent'] < a.min_targets:
        status = 'not calibrated: too few consistent calibrants (%d; need %d)' % (
            report['consistency']['targets'] - report['consistency']['inconsistent'], a.min_targets)
        log(status)
        curve = None
    if curve is not None and report['consistency']['share_inconsistent'] > a.max_share_inconsistent:
        status = ('not calibrated: the calibrants disagree (%.0f %% of the targets off the common curve by more '
                  'than %s; at most %.0f %% allowed)' % (100 * report['consistency']['share_inconsistent'],
                                                        tol_text(a, a.max_target_residual),
                                                        100 * a.max_share_inconsistent))
        log(status)
        curve = None
    if curve is not None:
        largest = float(np.abs(curve.mz_off).max())
        report['largest_offset_th'] = round(largest, 4)
        report['largest_offset_signed_th'] = round(float(curve.mz_off[np.argmax(np.abs(curve.mz_off))]), 4)
        if a.min_offset == 'auto':
            small = bool((np.abs(curve.mz_off) < th(a, 0.25, curve.mz_grid)).all())
            limit = '%s, a quarter of the fragment tolerance' % tol_text(a, 0.25)
        else:
            small = largest < float(a.min_offset)
            limit = '%g Th' % float(a.min_offset)
        if small:
            status = ('not corrected: the largest offset, %.3f Th, is below %s (the run is well calibrated)'
                      % (largest, limit))
            log(status)
            curve = None
        elif a.measure_only:
            status = 'not corrected: measured only (largest offset %.3f Th)' % largest
            log(status)
            curve = None
    if fitted is not None and not status.startswith('not calibrated'):
        seen = np.array([p['obs'] - p['theo'] for p in pts])
        if curve is not None:
            seen = seen - fitted(x, rt)
        if a.tolerance_counts:
            upto, n = tolerance_counts(seen / th(a, 1.0, x), a)
            with open(a.tolerance_counts, 'w') as fh:
                fh.write('upto_T\terrors\n' + ''.join('%.4f\t%d\n' % (u, c) for u, c in zip(upto, n)))
            log('fragment errors the search will see, for its tolerance: %d within +/-%s' % (n.sum(), tol_text(a, a.tol)))
    report['status'] = status
    if a.out:
        os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
        survey = survey_windows(a.windows)
        if curve is None:
            as_acquired(a, survey, log)
        else:
            note = 'MS2 m/z minus a smooth function of m/z%s fitted to %d fragment errors of %d library targets' % (
                ' and retention time' if use_rt else '', len(pts), info['targets_with_apex'])
            st = write_calibrated(a.run, a.out, curve, note, set(groups), *survey)
            report['written'] = st
            log('wrote %s: %d MS2 spectra corrected (%d peaks, largest shift %.3f Th); %d MS2 spectra of other '
                'scan types left as they were; %d survey scans left out' % (
                    a.out, st['spectra'], st['peaks'], st['max_abs_shift'], st['other_ms2'], st['dropped']))
    report['seconds'] = round(time.time() - t0, 1)
    with open(a.report + '.json', 'w') as fh:
        json.dump(report, fh, indent=1, default=float)
    with open(a.report + '.txt', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    finish(a, report, status, curve is not None)
    return 0


if __name__ == '__main__':
    sys.exit(main())
