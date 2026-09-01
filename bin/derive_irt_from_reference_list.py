#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sqlite3
import struct
import sys
import zlib


PEP_COLS   = {"compound", "peptidesequence", "sequence", "peptide",
              "modifiedpeptidesequence", "peptidemodifiedsequence"}
RT_COLS    = {"rt", "rt time (min)", "retention time", "retentiontime",
              "time", "explicitretentiontime", "rt time",
              "normalizedretentiontime"}
CHARGE_COLS = {"z", "charge", "precursor charge", "precursorcharge"}
MZ_COLS     = {"m/z", "mz", "precursor m/z", "precursor mz",
               "precursormz"}


def strip_mods(seq):
    if seq is None:
        return ""
    s = re.sub(r"\[[^\]]*\]", "", str(seq))
    s = re.sub(r"\([^)]*\)", "", s)
    return s.strip()


def autodetect_delim(path):
    with open(path) as f:
        head = f.readline()
    return "\t" if head.count("\t") > head.count(",") else ","


def find_column(header, candidates):
    for i, c in enumerate(header):
        if c.strip().lower() in candidates:
            return i
    return -1


def load_reference_list(path):
    delim = autodetect_delim(path)
    accum = {}
    with open(path) as f:
        reader = csv.reader(f, delimiter=delim)
        header = next(reader)
        pep_idx    = find_column(header, PEP_COLS)
        rt_idx     = find_column(header, RT_COLS)
        charge_idx = find_column(header, CHARGE_COLS)
        mz_idx     = find_column(header, MZ_COLS)
        if pep_idx < 0:
            sys.exit(
                f"derive_irt_from_reference_list: no peptide-sequence "
                f"column found in reference_list header {header!r}. "
                f"Looked for any of {sorted(PEP_COLS)}."
            )
        if rt_idx < 0:
            sys.exit(
                f"derive_irt_from_reference_list: no RT column found in "
                f"reference_list header {header!r}. Looked for any of "
                f"{sorted(RT_COLS)}. Without an RT column this route "
                f"cannot derive anchors. Supply --irt_library and "
                f"--irt_anchors explicitly, or accept the full-window "
                f"fallback."
            )
        for row in reader:
            if len(row) <= max(pep_idx, rt_idx):
                continue
            stripped = strip_mods(row[pep_idx])
            if not stripped:
                continue
            try:
                rt = float(row[rt_idx])
            except ValueError:
                continue
            charge = 0
            if charge_idx >= 0 and charge_idx < len(row):
                try:
                    charge = int(float(row[charge_idx]))
                except ValueError:
                    pass
            mz = 0.0
            if mz_idx >= 0 and mz_idx < len(row):
                try:
                    mz = float(row[mz_idx])
                except ValueError:
                    pass
            accum.setdefault(stripped, []).append((rt, charge, mz))

    out = {}
    for stripped, rows in accum.items():
        rts = sorted(r[0] for r in rows)
        median_rt = rts[len(rts) // 2]
        spread = rts[-1] - rts[0]
        if spread > 1.0:
            print(
                f"derive_irt_from_reference_list WARNING: peptide "
                f"{stripped} appears {len(rows)} times with RT spread "
                f"{spread:.2f} min ({rts[0]:.2f}..{rts[-1]:.2f}). "
                f"Using median {median_rt:.2f}.",
                file=sys.stderr,
            )
        charge = 0
        mz = 0.0
        for r in rows:
            if r[1] > 0 or r[2] > 0:
                charge = r[1]
                mz = r[2]
                break
        out[stripped] = {"rt": median_rt, "charge": charge, "mz": mz}
    return out


def decode_peak_arrays(mz_blob, int_blob,
                       compressed):
    mz_bytes = mz_blob
    int_bytes = int_blob
    if compressed:
        try:
            mz_bytes = zlib.decompress(mz_bytes)
            int_bytes = zlib.decompress(int_bytes)
        except zlib.error:
            pass
    n = len(mz_bytes) // 8
    mzs = list(struct.unpack(f"<{n}d", mz_bytes[: n * 8]))
    n_int = len(int_bytes) // 4
    ints = list(struct.unpack(f"<{n_int}f", int_bytes[: n_int * 4]))
    if len(ints) != len(mzs):
        n_int8 = len(int_blob) // 8
        try:
            ints = list(struct.unpack(f"<{n_int8}d", int_blob[: n_int8 * 8]))
        except struct.error:
            pass
    return mzs, ints


def top_n(mzs, ints, n
          ):
    if not mzs:
        return []
    pairs = sorted(zip(mzs, ints), key=lambda p: -p[1])[:n]
    return sorted(pairs, key=lambda p: p[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-list", required=True,
                    help="Skyline-style assay CSV/TSV with RT column.")
    ap.add_argument("--koina-dlib", required=True,
                    help="Koina .dlib to source transitions from (for -tr_irt).")
    ap.add_argument("--anchors-tsv", required=True)
    ap.add_argument("--osw-tsv", required=True)
    ap.add_argument("--top-transitions", type=int, default=6)
    args = ap.parse_args()

    if not os.path.exists(args.reference_list):
        sys.exit(f"derive_irt_from_reference_list: reference_list "
                 f"{args.reference_list} not found.")
    if not os.path.exists(args.koina_dlib):
        sys.exit(f"derive_irt_from_reference_list: koina dlib "
                 f"{args.koina_dlib} not found.")

    ref = load_reference_list(args.reference_list)
    print(
        f"derive_irt_from_reference_list: loaded {len(ref)} anchor "
        f"peptides from reference_list.",
        file=sys.stderr,
    )
    if not ref:
        sys.exit(
            "derive_irt_from_reference_list: reference_list contained no "
            "usable peptide+RT rows. Supply --irt_library / --irt_anchors "
            "or set rt_harmonize_fallback='full_window'."
        )

    with open(args.anchors_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["PeptideModSeq", "Irt", "Standard"])
        for stripped, meta in ref.items():
            w.writerow([stripped, f"{meta['rt']:.4f}", 1])
    print(
        f"derive_irt_from_reference_list: wrote {len(ref)} anchors to "
        f"{args.anchors_tsv} (Irt = run-minutes; Standard=1).",
        file=sys.stderr,
    )

    con = sqlite3.connect(args.koina_dlib)
    have_peaks_table = con.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='entries'"
    ).fetchone() is not None
    if not have_peaks_table:
        sys.exit("derive_irt_from_reference_list: koina .dlib has no "
                 "`entries` table.")

    koina_entries = {}
    for (pep_seq, charge, precmz, mz_blob, int_blob) in con.execute(
            "SELECT PeptideSeq, PrecursorCharge, PrecursorMz, "
            "       MassArray, IntensityArray "
            "FROM entries"):
        stripped = strip_mods(pep_seq)
        if not stripped:
            continue
        koina_entries.setdefault(stripped, []).append(
            (int(charge or 0), float(precmz or 0.0), mz_blob, int_blob)
        )
    con.close()

    header = [
        "PrecursorMz", "ProductMz", "LibraryIntensity",
        "NormalizedRetentionTime",
        "PeptideSequence", "ModifiedPeptideSequence",
        "ProteinId", "PrecursorCharge",
        "TransitionGroupId", "TransitionId",
        "Decoy", "PeptideGroupLabel", "DetectingTransition",
    ]

    n_matched = 0
    n_transitions = 0
    with open(args.osw_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for stripped, meta in ref.items():
            options = koina_entries.get(stripped, [])
            if not options:
                continue
            preferred = None
            if meta["charge"] > 0:
                for (charge, precmz, mz_blob, int_blob) in options:
                    if charge == meta["charge"]:
                        preferred = (charge, precmz, mz_blob, int_blob)
                        break
            if preferred is None:
                preferred = options[0]
            charge, precmz, mz_blob, int_blob = preferred
            mzs, ints = decode_peak_arrays(mz_blob, int_blob, 1)
            transitions = top_n(mzs, ints, args.top_transitions)
            if not transitions:
                continue
            n_matched += 1
            tg_id = f"{stripped}_{charge}"
            rt = meta["rt"]
            eff_precmz = meta["mz"] if meta["mz"] > 0 else precmz
            for k, (mz, intensity) in enumerate(transitions, 1):
                w.writerow([
                    f"{eff_precmz:.5f}",
                    f"{mz:.5f}",
                    f"{intensity:.2f}",
                    f"{rt:.4f}",
                    stripped,
                    stripped,
                    "reference_list_anchor",
                    charge,
                    tg_id,
                    f"{tg_id}_t{k}",
                    0,
                    stripped,
                    1,
                ])
                n_transitions += 1

    print(
        f"derive_irt_from_reference_list: wrote {n_transitions} "
        f"transitions across {n_matched} anchor spectra to "
        f"{args.osw_tsv} (matched {n_matched} of {len(ref)} "
        f"reference_list peptides against the Koina library).",
        file=sys.stderr,
    )
    if n_matched == 0:
        sys.exit(
            "derive_irt_from_reference_list HARD ERROR: 0 reference_list "
            "peptides matched Koina library entries. Check that --fasta "
            "actually contains the reference_list peptides (a proteome-wide "
            "Koina library should contain them; a small vendor FASTA may "
            "not). Alternative: supply --irt_library / --irt_anchors "
            "explicitly."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
