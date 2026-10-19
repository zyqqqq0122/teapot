#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sqlite3
import sys

DECOY = ('true', '1', 'yes', 'decoy')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference', required=True, help='reference list, Skyline layout')
    ap.add_argument('--library', required=True, help='.dlib or .elib')
    ap.add_argument('--out', required=True)
    ap.add_argument('--warnings', required=True)
    a = ap.parse_args()

    con = sqlite3.connect(f'file:{a.library}?mode=ro', uri=True)
    charges, sequences = {}, set()
    for modseq, seq, z in con.execute(
            'SELECT DISTINCT PeptideModSeq, PeptideSeq, PrecursorCharge FROM entries'):
        charges.setdefault(modseq, set()).add(int(z))
        sequences.add(seq)
    con.close()

    lines = open(a.reference, newline='').read().splitlines()
    delim = '\t' if lines[0].count('\t') >= lines[0].count(',') else ','
    rows = list(csv.reader(lines, delimiter=delim))
    header = [h.strip().strip('"').lower() for h in rows[0]]
    dcol = next((i for i, h in enumerate(header) if h in ('isdecoy', 'decoy', 'is_decoy')), None)

    out = []
    for r in rows[1:]:
        if not r or not r[0].strip():
            continue
        if dcol is not None and len(r) > dcol and r[dcol].strip().lower() in DECOY:
            continue
        peptide, z = r[0].strip(), int(float(r[4]))
        if z in charges.get(peptide, ()):
            status = 'in_library'
        elif peptide in charges:
            status = 'other_charge_only'
        elif re.sub(r'\[[^\]]*\]', '', peptide) in sequences:
            status = 'other_modifications_only'
        else:
            status = 'not_in_library'
        out.append((peptide, z, status))

    with open(a.out, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['peptide', 'charge', 'library'])
        w.writerows(out)

    n = {s: sum(1 for _p, _z, t in out if t == s) for s in
         ('in_library', 'not_in_library', 'other_charge_only', 'other_modifications_only')}
    missing = len(out) - n['in_library']
    print(f'reference precursors: {len(out)}; ' + ', '.join(f'{k} {v}' for k, v in n.items()),
          file=sys.stderr)
    with open(a.warnings, 'w') as f:
        if missing:
            f.write(f'{missing} of {len(out)} reference precursors cannot be detected, because '
                    f'the library does not hold them: {n["not_in_library"]} not at all, '
                    f'{n["other_charge_only"]} only at another charge, '
                    f'{n["other_modifications_only"]} only with other modifications. They are '
                    f'reported as not detected; they are listed in qc/{os.path.basename(a.out)}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
