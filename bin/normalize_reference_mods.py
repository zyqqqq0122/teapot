#!/usr/bin/env python3
import argparse
import csv
import re
import sqlite3
import sys


MOD = re.compile(r"\[([+-]?\d+(?:\.\d+)?)\]")


def key(modseq):
    seq, mods, pos, i = [], [], 0, 0
    while i < len(modseq):
        c = modseq[i]
        if c == "[":
            m = MOD.match(modseq, i)
            if not m:
                raise ValueError(f"bad modification in {modseq!r}")
            mods.append((pos, round(float(m.group(1)), 1)))
            i = m.end()
        elif c.isalpha():
            seq.append(c)
            pos += 1
            i += 1
        else:
            i += 1
    return "".join(seq), tuple(sorted(mods))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--library", required=True, help=".dlib or .elib")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    con = sqlite3.connect(f"file:{a.library}?mode=ro", uri=True)
    exact, lib = set(), {}
    for (modseq,) in con.execute("SELECT DISTINCT PeptideModSeq FROM entries"):
        exact.add(modseq)
        lib.setdefault(key(modseq), set()).add(modseq)
    con.close()

    lines = open(a.reference, newline="").read().splitlines()
    if not lines:
        sys.exit(f"normalize_reference_mods: {a.reference} is empty")
    delim = "\t" if lines[0].count("\t") >= lines[0].count(",") else ","
    rows = list(csv.reader(lines, delimiter=delim))
    changed = unmatched = ambiguous = 0
    for row in rows[1:]:
        if not row or not row[0].strip() or row[0].strip() in exact:
            continue
        hits = lib.get(key(row[0].strip()), set())
        if len(hits) == 1:
            row[0] = next(iter(hits))
            changed += 1
        elif hits:
            ambiguous += 1
        else:
            unmatched += 1
    with open(a.out, "w", newline="") as f:
        csv.writer(f, delimiter=delim, lineterminator="\n").writerows(rows)
    print(f"normalize_reference_mods: {len(rows) - 1} rows, {changed} rewritten to the "
          f"library notation, {unmatched} without a library entry, {ambiguous} with "
          f"several library spellings (both kept as is)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
