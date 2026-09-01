#!/usr/bin/env python3
import argparse
import csv
import json
import re
import shutil
import sqlite3
import struct
import sys
import zlib
from bisect import bisect_left
from collections import defaultdict


AA_MASS = {
    'A': 71.037114, 'R': 156.101111, 'N': 114.042927, 'D': 115.026943,
    'C': 103.009185, 'E': 129.042593, 'Q': 128.058578, 'G': 57.021464,
    'H': 137.058912, 'I': 113.084064, 'L': 113.084064, 'K': 128.094963,
    'M': 131.040485, 'F': 147.068414, 'P': 97.052764, 'S': 87.032028,
    'T': 101.047679, 'W': 186.079313, 'Y': 163.063329, 'V': 99.068414,
    'U': 150.953636, 'O': 237.147727,
}
PROTON = 1.007276466812
WATER = 18.010565


MOD_RE = re.compile(r'\[([+\-][\d.]+)\]')


def parse_modseq(modseq):
    s = modseq
    nterm = 0.0
    m = re.match(r'^\[([+\-][\d.]+)\]', s)
    if m:
        nterm = float(m.group(1))
        s = s[m.end():]

    stripped_chars = []
    mods = {}
    i = 0
    while i < len(s):
        ch = s[i]
        if not ch.isalpha():
            raise ValueError(f"unexpected char {ch!r} at pos {i} in {modseq!r}")
        stripped_chars.append(ch)
        pos = len(stripped_chars)
        i += 1
        if i < len(s) and s[i] == '[':
            m2 = re.match(r'\[([+\-][\d.]+)\]', s[i:])
            if not m2:
                raise ValueError(f"malformed mod at pos {i} in {modseq!r}")
            mods[pos] = mods.get(pos, 0.0) + float(m2.group(1))
            i += m2.end()
    return ''.join(stripped_chars), mods, nterm, 0.0


def build_modseq(stripped, mods, nterm=0.0, cterm=0.0, decimals=4):
    def fmt(d):
        return f"[{d:+.{decimals}f}]"

    parts = []
    if abs(nterm) > 0:
        parts.append(fmt(nterm))
    for i, ch in enumerate(stripped, start=1):
        parts.append(ch)
        d = mods.get(i, 0.0)
        if abs(d) > 0:
            parts.append(fmt(d))
    if abs(cterm) > 0:
        parts.append(fmt(cterm))
    return ''.join(parts)


def entry_key(stripped, mods, nterm, cterm, charge, tol=0.01):
    def q(d):
        return round(d / tol) * tol

    positions = []
    if abs(nterm) > tol:
        positions.append(('n_term', q(nterm)))
    for pos in sorted(mods.keys()):
        d = mods[pos]
        if abs(d) > tol:
            positions.append((pos, q(d)))
    if abs(cterm) > tol:
        positions.append(('c_term', q(cterm)))
    return (stripped, charge, tuple(positions))


def resolve_label_positions(stripped, label_spec):
    out = []
    for rule in label_spec:
        site = rule['site']
        delta = float(rule['delta'])
        if site == 'n_term':
            out.append(('n_term', delta, site))
        elif site == 'c_term':
            out.append(('c_term', delta, site))
        elif len(site) == 1 and site.isalpha():
            for i, ch in enumerate(stripped, start=1):
                if ch == site:
                    out.append((i, delta, site))
        else:
            raise ValueError(f"unrecognised label site {site!r}")
    return out


def classify_entry(stripped, mods, nterm, cterm, label_spec, tol):
    positions = resolve_label_positions(stripped, label_spec)
    if not positions:
        return 'unlabelable'
    hits = 0
    for (pos, delta, _s) in positions:
        if pos == 'n_term':
            observed = nterm
        elif pos == 'c_term':
            observed = cterm
        else:
            observed = mods.get(pos, 0.0)
        if abs(observed - delta) <= tol:
            hits += 1
    if hits == len(positions):
        return 'heavy'
    if hits == 0:
        return 'light'
    return 'mixed'


def apply_channel_shift(stripped, mods, nterm, cterm, label_spec, direction, tol):
    new_mods = dict(mods)
    new_nterm = nterm
    new_cterm = cterm
    shifts = []
    for (pos, delta, _s) in resolve_label_positions(stripped, label_spec):
        if direction == 'light_to_heavy':
            if pos == 'n_term':
                if abs(new_nterm - delta) > tol:
                    new_nterm += delta
                    shifts.append(('n_term', delta))
            elif pos == 'c_term':
                if abs(new_cterm - delta) > tol:
                    new_cterm += delta
                    shifts.append(('c_term', delta))
            else:
                if abs(new_mods.get(pos, 0.0) - delta) > tol:
                    new_mods[pos] = new_mods.get(pos, 0.0) + delta
                    shifts.append((pos, delta))
        elif direction == 'heavy_to_light':
            if pos == 'n_term':
                if abs(new_nterm - delta) <= tol:
                    new_nterm -= delta
                    if abs(new_nterm) <= tol:
                        new_nterm = 0.0
                    shifts.append(('n_term', -delta))
            elif pos == 'c_term':
                if abs(new_cterm - delta) <= tol:
                    new_cterm -= delta
                    if abs(new_cterm) <= tol:
                        new_cterm = 0.0
                    shifts.append(('c_term', -delta))
            else:
                have = new_mods.get(pos, 0.0)
                if abs(have - delta) <= tol:
                    new_mods[pos] = have - delta
                    if abs(new_mods[pos]) <= tol:
                        del new_mods[pos]
                    shifts.append((pos, -delta))
        else:
            raise ValueError(direction)
    return new_mods, new_nterm, new_cterm, shifts


def theoretical_ions(stripped, mods, nterm=0.0, cterm=0.0,
                     ion_types=('b', 'y'), max_charge=2):
    n = len(stripped)

    if 'b' in ion_types:
        cum = nterm
        for i in range(1, n):
            cum += AA_MASS[stripped[i - 1]] + mods.get(i, 0.0)
            for z in range(1, max_charge + 1):
                yield ('b', i, z, (cum + z * PROTON) / z)

    if 'y' in ion_types:
        cum = WATER + cterm
        for i in range(1, n):
            rpos = n - i + 1
            cum += AA_MASS[stripped[rpos - 1]] + mods.get(rpos, 0.0)
            for z in range(1, max_charge + 1):
                yield ('y', i, z, (cum + z * PROTON) / z)


def ion_covers_position(series, index, pos, n):
    if series == 'b':
        return 1 <= pos <= index
    if series == 'y':
        return (n - index + 1) <= pos <= n
    return False


def ion_covers_terminus(series, index, terminus, n):
    if terminus == 'n_term':
        return series == 'b'
    if terminus == 'c_term':
        return series == 'y'
    return False


def match_peaks(obs_mzs, obs_intens, theo, tol_ppm):
    theo_sorted = sorted(theo, key=lambda t: t[3])
    theo_mzs = [t[3] for t in theo_sorted]
    ann = []
    matched_tic = 0.0
    total_tic = 0.0
    for m, it in zip(obs_mzs, obs_intens):
        total_tic += it
        best = None
        best_err = None
        idx = bisect_left(theo_mzs, m)
        for j in (idx - 1, idx, idx + 1):
            if 0 <= j < len(theo_sorted):
                cand_mz = theo_mzs[j]
                err = abs(m - cand_mz) / m * 1e6
                if err <= tol_ppm and (best_err is None or err < best_err):
                    best_err = err
                    best = theo_sorted[j]
        if best is None:
            ann.append((m, it, None, None, None))
        else:
            ann.append((m, it, best[0], best[1], best[2]))
            matched_tic += it
    coverage = 0.0 if total_tic == 0 else matched_tic / total_tic
    return ann, coverage


def shifted_mz(obs_mz, series, index, charge, shifts, stripped_len):
    total_shift = 0.0
    for (pos, signed_delta) in shifts:
        if pos in ('n_term', 'c_term'):
            covers = ion_covers_terminus(series, index, pos, stripped_len)
        else:
            covers = ion_covers_position(series, index, pos, stripped_len)
        if covers:
            total_shift += signed_delta
    return obs_mz + total_shift / charge


def detect_format(con):
    tables = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    if 'RefSpectra' in tables and 'RefSpectraPeaks' in tables:
        return 'blib'
    if 'entries' in tables:
        return 'dlib'
    raise SystemExit("complete_channels: input SQLite is neither .blib "
                     "(RefSpectra) nor .dlib (entries).")


def blib_has_fileid(con):
    return 'fileid' in {r[1].lower()
                        for r in con.execute("PRAGMA table_info(RefSpectra)")}


def read_blib_entries(con):
    has_fileid = blib_has_fileid(con)
    fileid_col = "rs.fileID" if has_fileid else "NULL"
    cur = con.execute(f"""
        SELECT rs.id, rs.peptideSeq, rs.peptideModSeq, rs.precursorCharge,
               rs.precursorMZ, rs.numPeaks, rs.retentionTime,
               rp.peakMZ, rp.peakIntensity, {fileid_col}
        FROM RefSpectra rs
        JOIN RefSpectraPeaks rp ON rp.RefSpectraId = rs.id
    """)
    for (rid, seq, modseq, z, pmz, npeaks, rt, mz_blob, it_blob, file_id) in cur:
        mzs = list(struct.unpack(f'<{npeaks}d', mz_blob))
        its = list(struct.unpack(f'<{npeaks}f', it_blob))
        yield {
            'id': rid, 'peptide_seq': seq, 'peptide_mod_seq': modseq,
            'charge': z, 'precursor_mz': pmz, 'num_peaks': npeaks,
            'rt': rt, 'peak_mzs': mzs, 'peak_intens': its,
            'file_id': file_id,
        }


def append_blib_entry(con, entry, proteins_map, has_fileid=True):
    cur = con.cursor()
    npeaks = len(entry['peak_mzs'])
    mz_blob = struct.pack(f'<{npeaks}d', *entry['peak_mzs'])
    it_blob = struct.pack(f'<{npeaks}f', *entry['peak_intens'])
    fileid_col = ", fileID" if has_fileid else ""
    fileid_val = ", ?" if has_fileid else ""
    params = [entry['peptide_seq'], entry['peptide_mod_seq'], entry['charge'],
              entry['precursor_mz'], npeaks, entry.get('rt', 0.0)]
    if has_fileid:
        params.append(entry.get('file_id'))
    cur.execute(f"""
        INSERT INTO RefSpectra
          (peptideSeq, peptideModSeq, precursorCharge, precursorAdduct,
           precursorMZ, moleculeName, chemicalFormula, inchiKey, otherKeys,
           prevAA, nextAA, copies, numPeaks, retentionTime,
           SpecIDinFile, score, scoreType{fileid_col})
        VALUES (?, ?, ?, '', ?, '', '', '', '', '-', '-', 1, ?, ?, '', 0, 0{fileid_val})
    """, params)
    new_id = cur.lastrowid
    cur.execute(
        "INSERT INTO RefSpectraPeaks (RefSpectraId, peakMZ, peakIntensity) "
        "VALUES (?, ?, ?)", (new_id, mz_blob, it_blob))
    for (pos, delta) in entry.get('mods', {}).items():
        if abs(delta) > 0:
            cur.execute(
                "INSERT INTO Modifications (RefSpectraId, position, mass) "
                "VALUES (?, ?, ?)", (new_id, pos, delta))
    for pid in proteins_map.get(entry['source_id'], ()):
        cur.execute(
            "INSERT INTO RefSpectraProteins (RefSpectraId, ProteinId) "
            "VALUES (?, ?)", (new_id, pid))

    rt_cols = [r[1] for r in con.execute("PRAGMA table_info(RetentionTimes)")]
    if rt_cols:
        copyable = [c for c in rt_cols
                    if c.lower() not in ('id', 'refspectraid')]
        src_rows = cur.execute(
            f"SELECT {', '.join(copyable)} FROM RetentionTimes "
            "WHERE RefSpectraId = ?", (entry['source_id'],)).fetchall()
        if src_rows:
            placeholders = ', '.join('?' * (len(copyable) + 1))
            for row in src_rows:
                cur.execute(
                    f"INSERT INTO RetentionTimes "
                    f"(RefSpectraId, {', '.join(copyable)}) "
                    f"VALUES ({placeholders})", (new_id, *row))
        else:
            defaults = {'RefSpectraId': new_id,
                        'retentionTime': entry.get('rt', 0.0)}
            if 'RedundantRefSpectraId' in rt_cols:
                defaults['RedundantRefSpectraId'] = -1
            if 'SpectrumSourceId' in rt_cols:
                defaults['SpectrumSourceId'] = entry.get('file_id') or 1
            if 'bestSpectrum' in rt_cols:
                defaults['bestSpectrum'] = 1
            cols = [c for c in rt_cols if c in defaults]
            cur.execute(
                f"INSERT INTO RetentionTimes ({', '.join(cols)}) "
                f"VALUES ({', '.join('?' * len(cols))})",
                tuple(defaults[c] for c in cols))
    return new_id


def load_blib_proteins_map(con):
    m = defaultdict(list)
    for rid, pid in con.execute(
            "SELECT RefSpectraId, ProteinId FROM RefSpectraProteins"):
        m[rid].append(pid)
    return m


def load_blib_accessions(con):
    """RefSpectraId -> set of protein accession strings."""
    try:
        rows = con.execute(
            "SELECT rp.RefSpectraId, p.accession "
            "FROM RefSpectraProteins rp JOIN Proteins p ON p.id = rp.ProteinId")
    except sqlite3.Error:
        return {}
    m = defaultdict(set)
    for rid, acc in rows:
        if acc:
            m[rid].add(str(acc))
    return m


def read_dlib_entries(con):
    cur = con.execute("""
        SELECT PrecursorMz, PrecursorCharge, PeptideModSeq, PeptideSeq,
               RTInSeconds, MassEncodedLength, MassArray,
               IntensityEncodedLength, IntensityArray
        FROM entries
    """)
    for (pmz, z, mseq, seq, rt, mlen, mblob, ilen, iblob) in cur:
        n = mlen // 8
        mzs = list(struct.unpack(f'>{n}d', zlib.decompress(mblob)))
        its = list(struct.unpack(f'>{ilen // 4}f', zlib.decompress(iblob)))
        yield {
            'id': None, 'peptide_seq': seq, 'peptide_mod_seq': mseq,
            'charge': z, 'precursor_mz': pmz, 'num_peaks': n,
            'rt': rt, 'peak_mzs': mzs, 'peak_intens': its,
        }


def append_dlib_entry(con, entry):
    cur = con.cursor()
    n = len(entry['peak_mzs'])
    mblob = zlib.compress(struct.pack(f'>{n}d', *entry['peak_mzs']))
    iblob = zlib.compress(struct.pack(f'>{n}f', *entry['peak_intens']))
    cur.execute("""
        INSERT INTO entries
          (PrecursorMz, PrecursorCharge, PeptideModSeq, PeptideSeq,
           Copies, RTInSeconds, Score,
           MassEncodedLength, MassArray,
           IntensityEncodedLength, IntensityArray, SourceFile)
        VALUES (?, ?, ?, ?, 1, ?, 0.0, ?, ?, ?, ?, 'complete_channels')
    """, (entry['precursor_mz'], entry['charge'],
          entry['peptide_mod_seq'], entry['peptide_seq'],
          entry.get('rt', 0.0),
          n * 8, mblob, n * 4, iblob))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--labels', required=True,
                    help="JSON list: [{\"site\":\"K\",\"delta\":8.014199}, ...]")
    ap.add_argument('--mod-tol', type=float, default=0.01,
                    help="Tolerance for identifying label-matching deltas.")
    ap.add_argument('--ppm', type=float, default=20.0,
                    help="ppm tolerance for peak <-> theoretical matching.")
    ap.add_argument('--ion-types', default='b,y')
    ap.add_argument('--max-frag-charge', type=int, default=2)
    ap.add_argument('--coverage-hard-floor', type=float, default=0.5,
                    help="Below this the process fails, naming the peptide.")
    ap.add_argument('--coverage-warn-below', type=float, default=0.9,
                    help="Below this we drop unannotated peaks (producing "
                         "an annotated-ions-only twin spectrum, which is "
                         "exactly what a Prosit-predicted .dlib is).")
    ap.add_argument('--qc', default='complete_channels_qc.tsv')
    ap.add_argument('--skip-proteins', default='',
                    help="Comma-separated protein accessions (substring match) "
                         "whose peptides get NO synthesised twin. For a "
                         "spike-in kit whose constructs share a purification / "
                         "RT fusion tag, the tag's peptides exist ONLY as "
                         "heavy standard -- there is no endogenous light form "
                         "and the acquisition does not target one. Synthesising "
                         "light twins for them adds precursors the instrument "
                         "never isolated; they then pick up signal from "
                         "whatever unrelated peptide shares their isolation "
                         "window and can be reported as identifications. "
                         "Example: --skip-proteins iRT_Tag")
    args = ap.parse_args()

    label_spec = json.loads(args.labels)
    ion_types = tuple(x.strip() for x in args.ion_types.split(','))

    shutil.copy(args.inp, args.out)
    con = sqlite3.connect(args.out)

    fmt = detect_format(con)
    print(f"format: {fmt}", file=sys.stderr)

    skip_proteins = [s.strip() for s in args.skip_proteins.split(',') if s.strip()]

    if fmt == 'blib':
        proteins_map = load_blib_proteins_map(con)
        accessions = load_blib_accessions(con)
        has_fileid = blib_has_fileid(con)
        entry_iter = read_blib_entries(con)

        def append(entry):
            return append_blib_entry(con, entry, proteins_map, has_fileid)

        existing_key_query = ("SELECT peptideSeq, peptideModSeq, "
                              "precursorCharge FROM RefSpectra")
    else:
        proteins_map = None
        accessions = {}
        entry_iter = read_dlib_entries(con)

        def append(entry):
            return append_dlib_entry(con, entry)

        existing_key_query = ("SELECT PeptideSeq, PeptideModSeq, "
                              "PrecursorCharge FROM entries")

    if skip_proteins and not accessions:
        print(f"WARNING --skip-proteins {skip_proteins} given but this library "
              "carries no protein accessions; nothing will be skipped.",
              file=sys.stderr)

    existing_keys = set()
    for row in con.execute(existing_key_query):
        try:
            st, mds, nt, ct = parse_modseq(row[1])
            existing_keys.add(entry_key(st, mds, nt, ct, row[2], args.mod_tol))
        except Exception:
            pass

    entries = list(entry_iter)
    qc_rows = []
    counters = defaultdict(int)

    for e in entries:
        try:
            stripped, mods, nterm, cterm = parse_modseq(e['peptide_mod_seq'])
        except Exception as exc:
            counters['parse_fail'] += 1
            qc_rows.append((e['peptide_mod_seq'], e['charge'],
                            'parse_fail', 0.0, 0, str(exc)))
            continue

        entry_accs = accessions.get(e['id'], set())
        hit = next((pat for pat in skip_proteins
                    if any(pat in acc for acc in entry_accs)), None)
        if hit:
            counters['skipped_protein'] += 1
            qc_rows.append((e['peptide_mod_seq'], e['charge'], 'n/a', 0.0,
                            e['num_peaks'], f'skipped_protein:{hit}'))
            continue

        cls = classify_entry(stripped, mods, nterm, cterm, label_spec, args.mod_tol)
        counters[f'input_{cls}'] += 1

        if cls == 'unlabelable':
            qc_rows.append((e['peptide_mod_seq'], e['charge'],
                            'unlabelable', 0.0, e['num_peaks'], 'skipped'))
            continue
        if cls == 'mixed':
            qc_rows.append((e['peptide_mod_seq'], e['charge'],
                            'mixed', 0.0, 0, 'skipped'))
            continue

        direction = 'light_to_heavy' if cls == 'light' else 'heavy_to_light'
        new_mods, new_nterm, new_cterm, shifts = apply_channel_shift(
            stripped, mods, nterm, cterm, label_spec, direction, args.mod_tol)

        twin_key = entry_key(stripped, new_mods, new_nterm, new_cterm,
                             e['charge'], args.mod_tol)
        if twin_key in existing_keys:
            counters['already_present'] += 1
            qc_rows.append((e['peptide_mod_seq'], e['charge'], cls,
                            1.0, e['num_peaks'], 'already_present'))
            continue

        if not e['peak_mzs']:
            counters['no_peaks'] = counters.get('no_peaks', 0) + 1
            qc_rows.append((e['peptide_mod_seq'], e['charge'], cls,
                            0.0, 0, 'no_peaks'))
            continue

        theo = list(theoretical_ions(stripped, mods, nterm, cterm,
                                     ion_types=ion_types,
                                     max_charge=args.max_frag_charge))
        ann, coverage = match_peaks(e['peak_mzs'], e['peak_intens'],
                                    theo, args.ppm)

        if coverage < args.coverage_hard_floor:
            sys.exit(
                f"HARD FAIL: annotation coverage {coverage:.3f} < "
                f"hard_floor {args.coverage_hard_floor} for peptide "
                f"{e['peptide_mod_seq']} charge {e['charge']}. Check "
                f"ion_types / ppm tolerance / label spec.")

        new_mzs, new_intens = [], []
        n_matched = 0
        for (obs_mz, obs_it, series, idx, z) in ann:
            if series is None:
                continue
            new_mzs.append(shifted_mz(obs_mz, series, idx, z, shifts, len(stripped)))
            new_intens.append(obs_it)
            n_matched += 1
        low_coverage = coverage < args.coverage_warn_below

        pshift = sum(sd for (_p, sd) in shifts) / e['charge']
        twin = {
            'peptide_seq': stripped,
            'peptide_mod_seq': build_modseq(stripped, new_mods, new_nterm, new_cterm),
            'charge': e['charge'],
            'precursor_mz': e['precursor_mz'] + pshift,
            'rt': e['rt'],
            'peak_mzs': new_mzs,
            'peak_intens': new_intens,
            'mods': new_mods,
            'source_id': e['id'],
            'file_id': e.get('file_id'),
        }
        append(twin)
        existing_keys.add(twin_key)
        counters['twins_added'] += 1
        if low_coverage:
            counters['low_coverage_warn'] += 1
        status = 'added_low_coverage' if low_coverage else 'added'
        qc_rows.append((e['peptide_mod_seq'], e['charge'], cls,
                        coverage, n_matched, status))

    if fmt == 'blib':
        try:
            n = con.execute("SELECT COUNT(*) FROM RefSpectra").fetchone()[0]
            con.execute("UPDATE LibInfo SET numSpecs = ?", (n,))
        except sqlite3.Error as exc:
            print(f"note: could not update LibInfo.numSpecs: {exc}",
                  file=sys.stderr)

    con.commit()
    con.close()

    with open(args.qc, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['peptide_mod_seq', 'charge', 'input_channel',
                    'annotation_coverage', 'n_kept_peaks', 'status'])
        for r in qc_rows:
            w.writerow(r)

    print(f"input entries:      {len(entries)}", file=sys.stderr)
    for k in ('input_light', 'input_heavy', 'input_mixed', 'input_unlabelable',
              'parse_fail', 'skipped_protein', 'twins_added',
              'low_coverage_warn', 'already_present', 'no_peaks'):
        print(f"  {k:26} {counters[k]}", file=sys.stderr)


if __name__ == '__main__':
    main()
