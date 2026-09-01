#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sqlite3
import struct
import sys
import zlib


def strip_mods(seq):
    s = re.sub(r"\[[^\]]*\]", "", str(seq or ""))
    s = re.sub(r"\([^)]*\)", "", s)
    return s.strip()


def blib_has_table(con, name):
    row = con.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
        (name,),
    ).fetchone()
    return row is not None


def _inflate_if_needed(blob, *expected_sizes):
    if any(len(blob) == s for s in expected_sizes if s):
        return blob
    try:
        return zlib.decompress(blob)
    except zlib.error:
        return blob


def decode_peak_arrays(mz_blob, int_blob,
                       num_peaks = 0):
    mz_bytes = _inflate_if_needed(mz_blob, num_peaks * 8)
    int_bytes = _inflate_if_needed(int_blob, num_peaks * 4, num_peaks * 8)
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
    ap.add_argument("--in", dest="blib", required=True)
    ap.add_argument("--anchors-tsv", required=True,
                    help="Two-column TSV: PeptideModSeq (stripped) \\t Irt")
    ap.add_argument("--osw-tsv", required=True,
                    help="OpenSwath -tr_irt input TSV with transitions.")
    ap.add_argument("--top-transitions", type=int, default=6)
    ap.add_argument(
        "--require-label-deltas", default="",
        help="Comma-separated monoisotopic deltas (e.g. '8.0142,10.0083' for "
             "K:UniMod:259 + R:UniMod:267). When set, only RefSpectra rows "
             "whose peptideModSeq contains at least one of these delta tokens "
             "are emitted to --osw-tsv. Anchors that fail this filter are "
             "skipped from -tr_irt but still contribute to the anchors TSV "
             "(HARMONIZE_RT's regression is label-agnostic since it fits on "
             "stripped sequence). The purpose is to prevent -tr_irt from "
             "asking OpenSwathWorkflow to find light iRT peptides that were "
             "never spiked. Those exist in the library but not in the "
             "sample, and OSW's RT normalisation would then fail with a "
             "misleading 'not enough iRT peptides' error."
    )
    ap.add_argument("--label-tolerance", type=float, default=0.02,
                    help="ppm/Da tolerance for delta-token matching.")
    args = ap.parse_args()

    required_deltas: list[float] = []
    for tok in args.require_label_deltas.split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            required_deltas.append(float(tok))
        except ValueError:
            sys.exit(f"emit_irt_library: cannot parse --require-label-deltas "
                     f"token '{tok}' as float.")

    if not os.path.exists(args.blib):
        sys.exit(f"emit_irt_library: input blib {args.blib} not found.")

    con = sqlite3.connect(args.blib)
    if not blib_has_table(con, "IrtLibrary"):
        con.close()
        sys.exit(
            "emit_irt_library: blib has no IrtLibrary table. "
            "Supply --irt_library explicitly to the pipeline."
        )

    ir_cols = {r[1] for r in con.execute("PRAGMA table_info(IrtLibrary)")}
    has_standard = "Standard" in ir_cols
    if not has_standard:
        print(
            "emit_irt_library WARNING: IrtLibrary has no `Standard` column. "
            "The regression anchors (--anchors-tsv) will draw from ALL rows: "
            "fine, more anchors mean better RT coverage, but the "
            "run-anchor set (--osw-tsv → OSW -tr_irt) cannot be split into "
            "reliable spike-ins vs endogenous variable-abundance peptides. "
            "Emitting ALL rows into -tr_irt as a fallback, which may cause "
            "OSW to hunt for endogenous peptides that are not always "
            "detectable in every run and produce a 'not enough iRT peptides' "
            "warning at extraction. If that becomes a problem, override with "
            "--irt_library / --irt_anchors supplying a curated run-anchor "
            "set (e.g. Biognosys iRT kit transitions).",
            file=sys.stderr,
        )

    select_cols = "PeptideModSeq, Irt" + (", Standard" if has_standard else "")
    accum: dict[str, list[tuple[float, int]]] = {}
    for row in con.execute(
            f"SELECT {select_cols} FROM IrtLibrary WHERE Irt IS NOT NULL"):
        stripped = strip_mods(row[0])
        if not stripped:
            continue
        try:
            irt = float(row[1])
        except (TypeError, ValueError):
            continue
        standard = int(row[2]) if has_standard else 1
        accum.setdefault(stripped, []).append((irt, standard))

    def _collapse(items, key
                  ):
        vs = sorted(v for v, _ in items)
        spread = vs[-1] - vs[0]
        std1 = [v for v, s in items if s == 1]
        chosen = (sorted(std1)[len(std1) // 2] if std1
                  else vs[len(vs) // 2])
        std_out = 1 if std1 else 0
        if spread > 1.0 and len(items) > 1:
            print(
                f"emit_irt_library WARNING: anchor {key} has {len(items)} "
                f"rows disagreeing by {spread:.2f} iRT "
                f"({vs[0]:.2f}..{vs[-1]:.2f}); using "
                f"{'Standard=1 median' if std1 else 'median'} {chosen:.2f}.",
                file=sys.stderr,
            )
        return chosen, std_out

    anchors_written = 0
    anchor_seqs_all: dict[str, float] = {}
    anchor_seqs_standard1: dict[str, float] = {}
    with open(args.anchors_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["PeptideModSeq", "Irt", "Standard"])
        for stripped, items in accum.items():
            irt, standard = _collapse(items, stripped)
            w.writerow([stripped, f"{irt:.4f}", standard])
            anchor_seqs_all[stripped] = irt
            if standard == 1:
                anchor_seqs_standard1[stripped] = irt
            anchors_written += 1
    print(
        f"emit_irt_library: wrote {anchors_written} anchors to "
        f"{args.anchors_tsv} "
        f"(Standard=1: {len(anchor_seqs_standard1)}, "
        f"Standard=0: {len(anchor_seqs_all) - len(anchor_seqs_standard1)}).",
        file=sys.stderr,
    )

    header = [
        "PrecursorMz", "ProductMz", "LibraryIntensity",
        "NormalizedRetentionTime",
        "PeptideSequence", "ModifiedPeptideSequence",
        "ProteinId", "PrecursorCharge",
        "TransitionGroupId", "TransitionId",
        "Decoy", "PeptideGroupLabel", "DetectingTransition",
    ]

    ref_rows = con.execute(
        "SELECT id, peptideSeq, peptideModSeq, precursorCharge, "
        "       precursorMZ, numPeaks "
        "FROM RefSpectra"
    ).fetchall()

    have_peaks_table = blib_has_table(con, "RefSpectraPeaks")
    if not have_peaks_table:
        con.close()
        sys.exit("emit_irt_library: RefSpectraPeaks table missing.")

    def modseq_carries_required_delta(modseq):
        if not required_deltas:
            return True
        deltas_in_seq = [float(t) for t in re.findall(r"[-+]?\d+\.\d+",
                                                       str(modseq or ""))]
        for req in required_deltas:
            for got in deltas_in_seq:
                if abs(got - req) <= args.label_tolerance:
                    return True
        return False

    n_transitions_written = 0
    n_anchor_spectra = 0
    n_skipped_label = 0
    n_skipped_not_standard1 = 0
    with open(args.osw_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for (spec_id, pep_seq, pep_modseq, charge, precmz, npeaks) in ref_rows:
            stripped = strip_mods(pep_seq or pep_modseq)
            if stripped not in anchor_seqs_all:
                continue
            if stripped not in anchor_seqs_standard1:
                n_skipped_not_standard1 += 1
                continue
            if required_deltas and not modseq_carries_required_delta(pep_modseq):
                n_skipped_label += 1
                continue
            row = con.execute(
                "SELECT peakMZ, peakIntensity "
                "FROM RefSpectraPeaks WHERE RefSpectraID = ?",
                (spec_id,),
            ).fetchone()
            if row is None:
                continue
            mz_blob, int_blob = row
            mzs, ints = decode_peak_arrays(mz_blob, int_blob, int(npeaks or 0))
            transitions = top_n(mzs, ints, args.top_transitions)
            if not transitions:
                continue
            n_anchor_spectra += 1
            irt = anchor_seqs_standard1[stripped]
            transition_group_id = f"{stripped}_{int(charge)}"
            for k, (mz, intensity) in enumerate(transitions, 1):
                w.writerow([
                    f"{precmz:.5f}",
                    f"{mz:.5f}",
                    f"{intensity:.2f}",
                    f"{irt:.4f}",
                    stripped,
                    pep_modseq or stripped,
                    "iRT_anchor",
                    int(charge),
                    transition_group_id,
                    f"{transition_group_id}_t{k}",
                    0,
                    stripped,
                    1,
                ])
                n_transitions_written += 1
    con.close()

    print(
        f"emit_irt_library: wrote {n_transitions_written} transitions "
        f"across {n_anchor_spectra} anchor spectra to {args.osw_tsv} "
        f"(-tr_irt: Standard=1 heavy only).",
        file=sys.stderr,
    )
    if n_skipped_not_standard1:
        print(
            f"emit_irt_library: skipped {n_skipped_not_standard1} RefSpectra "
            f"whose IrtLibrary rows were Standard=0 (endogenous anchors "
            f"kept in the HARMONIZE_RT anchors TSV, excluded from -tr_irt "
            f"because their abundance in any given run is not guaranteed).",
            file=sys.stderr,
        )
    if required_deltas:
        print(
            f"emit_irt_library: skipped {n_skipped_label} anchor spectra "
            f"whose peptideModSeq lacked required label deltas "
            f"{required_deltas} (only heavy anchors are useful for "
            f"-tr_irt when the sample is spiked heavy-only).",
            file=sys.stderr,
        )
    if n_anchor_spectra == 0:
        sys.exit(
            "emit_irt_library HARD ERROR: 0 anchor spectra passed all "
            "filters. Either IrtLibrary has anchor peptide names but no "
            "spectra for them exist in RefSpectra, or --require-label-deltas "
            "removed every Standard=1 candidate, or the blib has no "
            "Standard=1 rows at all. OSW cannot use this for -tr_irt. "
            "Supply --irt_library explicitly, run EMIT_IRT_LIBRARY on the "
            "pre-COMPLETE_CHANNELS vendor blib, or fall back to full-window "
            "extraction (loud warning; expensive)."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
