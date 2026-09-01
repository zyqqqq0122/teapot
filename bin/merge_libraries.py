#!/usr/bin/env python3
import argparse
import shutil
import sqlite3
import sys


ENTRIES_COLUMNS = (
    "PrecursorMz", "PrecursorCharge", "PeptideModSeq", "PeptideSeq",
    "Copies", "RTInSeconds", "Score",
    "MassEncodedLength", "MassArray",
    "IntensityEncodedLength", "IntensityArray",
    "CorrelationEncodedLength", "CorrelationArray",
    "QuantifiedIonsArray",
    "RTInSecondsStart", "RTInSecondsStop", "IonMobility",
    "MedianChromatogramEncodedLength", "MedianChromatogramArray",
    "SourceFile",
)

PROVENANCE_COLUMN = "provenance"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primary", required=True,
                    help="Primary .dlib (kept in full; wins on collision).")
    ap.add_argument("--secondary", required=True,
                    help="Secondary .dlib (colliding entries dropped).")
    ap.add_argument("--out", required=True, help="Output .dlib path.")
    args = ap.parse_args()

    shutil.copy(args.primary, args.out)

    con = sqlite3.connect(args.out)
    cur = con.cursor()

    def _table_columns(cur, table):
        return {r[1] for r in cur.execute(f"PRAGMA table_info({table})")}

    if PROVENANCE_COLUMN not in _table_columns(cur, "entries"):
        cur.execute(
            f"ALTER TABLE entries ADD COLUMN {PROVENANCE_COLUMN} TEXT")
        cur.execute(
            f"UPDATE entries SET {PROVENANCE_COLUMN} = 'target' "
            f"WHERE {PROVENANCE_COLUMN} IS NULL")

    primary_keys = set()
    for row in cur.execute("SELECT PeptideSeq, PrecursorCharge FROM entries"):
        primary_keys.add((row[0], row[1]))
    n_primary = len(primary_keys)

    src = sqlite3.connect(args.secondary)
    src_cur = src.cursor()

    src_has_provenance = PROVENANCE_COLUMN in _table_columns(src_cur, "entries")
    entries_cols = list(ENTRIES_COLUMNS) + [PROVENANCE_COLUMN]

    col_sql = ", ".join(entries_cols)
    placeholders = ", ".join(["?"] * len(entries_cols))
    insert_sql = f"INSERT INTO entries ({col_sql}) VALUES ({placeholders})"

    if src_has_provenance:
        src_select = f"SELECT {col_sql} FROM entries"
    else:
        src_select = (
            f"SELECT {', '.join(ENTRIES_COLUMNS)}, 'background' "
            f"FROM entries"
        )

    n_secondary_seen = 0
    n_secondary_added = 0
    n_secondary_dropped_collision = 0

    for row in src_cur.execute(src_select):
        n_secondary_seen += 1
        key = (row[3], row[1])  # PeptideSeq, PrecursorCharge
        if key in primary_keys:
            n_secondary_dropped_collision += 1
            continue
        cur.execute(insert_sql, row)
        n_secondary_added += 1

    existing_ptp = set()
    for row in cur.execute(
            "SELECT PeptideSeq, isDecoy, ProteinAccession FROM peptidetoprotein"):
        existing_ptp.add(row)
    n_ptp_added = 0
    for row in src_cur.execute(
            "SELECT PeptideSeq, isDecoy, ProteinAccession FROM peptidetoprotein"):
        if row in existing_ptp:
            continue
        cur.execute(
            "INSERT INTO peptidetoprotein (PeptideSeq, isDecoy, ProteinAccession) "
            "VALUES (?, ?, ?)", row)
        existing_ptp.add(row)
        n_ptp_added += 1

    con.commit()
    con.close()
    src.close()

    con = sqlite3.connect(args.out)
    n_total = con.execute("SELECT COUNT(*) FROM entries").fetchone()[0]
    n_prot = con.execute("SELECT COUNT(*) FROM peptidetoprotein").fetchone()[0]
    con.close()

    print(f"primary entries kept:        {n_primary} unique (PeptideSeq, "
          f"PrecursorCharge) keys", file=sys.stderr)
    print(f"secondary entries seen:      {n_secondary_seen}", file=sys.stderr)
    print(f"secondary added:             {n_secondary_added}", file=sys.stderr)
    print(f"secondary dropped collision: {n_secondary_dropped_collision}",
          file=sys.stderr)
    print(f"peptidetoprotein rows added: {n_ptp_added}", file=sys.stderr)
    print(f"union entries total:         {n_total}", file=sys.stderr)
    print(f"union peptidetoprotein rows: {n_prot}", file=sys.stderr)


if __name__ == "__main__":
    main()
