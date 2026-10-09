#!/usr/bin/env python3
import argparse
import csv
import glob
import math
import os
import re
import sqlite3
import statistics
import sys


UNIMOD = re.compile(r"\(UniMod:(\d+)\)")

# Basis: the window covers how far apart the runs put the same precursor, taken over the
# precursors identified in two or more runs. The 99.9th percentile of that, plus 30 s for
# the peak, and at least 60 s.
# Single runs and sparse data: with one run it’s 60 s. If fewer than 20 precursors are
# found in more than one run, it stays at 180 s.
# Likely wrong identifications are left out first: an identification counts as out of line
# when it sits more than 10 robust standard deviations from its run’s typical shift
RT_QUANTILE = 0.999
RT_OUTLIER = 10
RT_PEAK = 30
RT_MIN = 60
RT_MIN_SHARED = 20
RT_UNMEASURED = 180


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def from_encyclopedia(paths, fdr):
    hits = {}
    for p in paths:
        n = 0
        for r in read_tsv(p):
            q = num(r.get("q-value"))
            if q is None or q > fdr:
                continue
            head, z = r["PSMId"].rsplit("+", 1)
            _file, rt, pf = head.rsplit(":", 2)
            hits.setdefault((pf, int(z)), {})[p] = float(rt)
            n += 1
        print(f"{os.path.basename(p)}: {n} precursors at q <= {fdr}", file=sys.stderr)
    return hits


def unimod_to_delta(pf, deltas):
    s = pf
    nterm = ""
    if s.startswith(".("):
        m = UNIMOD.match(s, 1)
        if not m or m.group(1) not in deltas:
            return None
        nterm = "[%+.6f]" % deltas[m.group(1)]
        s = s[m.end():]
    if ".(" in s:
        return None

    def sub(m):
        d = deltas.get(m.group(1))
        if d is None:
            raise KeyError(m.group(1))
        return "[%+.6f]" % d

    try:
        return nterm + UNIMOD.sub(sub, s)
    except KeyError:
        return None


def from_pyprophet(paths, fdr, deltas):
    hits, mz = {}, {}
    skipped = set()
    for p in paths:
        n = 0
        for r in read_tsv(p):
            q = num(r.get("m_score"))
            if q is None or q > fdr or r.get("decoy", "0") not in ("0", "0.0"):
                continue
            if r.get("peak_group_rank", "1") not in ("1", "1.0"):
                continue
            pf = unimod_to_delta(r["FullPeptideName"], deltas)
            if pf is None:
                skipped.add(r["FullPeptideName"])
                continue
            key = (pf, int(float(r["Charge"])))
            hits.setdefault(key, {})[p] = float(r["RT"])
            mz.setdefault(key, float(r["mz"]))
            n += 1
        print(f"{os.path.basename(p)}: {n} precursors at q <= {fdr}", file=sys.stderr)
    if skipped:
        print(f"skipped {len(skipped)} peptidoform(s) whose modification has no "
              f"[+delta] form diathem can read, e.g. {sorted(skipped)[:3]}", file=sys.stderr)
    return hits, mz


def rt_window(hits, n_runs):
    if n_runs == 1:
        return RT_MIN, "one run: each target is read around the RT it was identified at"
    multi = {k: rts for k, rts in hits.items() if len(rts) > 1}
    shifts = {}
    for k, rts in multi.items():
        c = statistics.median(rts.values())
        for run, rt in rts.items():
            shifts.setdefault(run, []).append((k, rt - c))
    out = set()
    for run, ks in shifts.items():
        if len(ks) < RT_MIN_SHARED:
            continue
        d = [x for _, x in ks]
        m = statistics.median(d)
        sd = max(1.0, 1.4826 * statistics.median(abs(x - m) for x in d))   # MAD as an SD
        out.update((k, run) for k, x in ks if abs(x - m) > RT_OUTLIER * sd)
    spread = []
    for k, rts in multi.items():
        kept = [rt for run, rt in rts.items() if (k, run) not in out]
        if len(kept) > 1:
            spread.append(max(kept) - min(kept))
    if len(spread) < RT_MIN_SHARED:
        return RT_UNMEASURED, (f"only {len(spread)} precursor(s) identified in two or more "
                               f"runs, too few to measure how far RTs move between runs")
    spread.sort()
    p = spread[int(RT_QUANTILE * (len(spread) - 1))]
    return max(RT_MIN, math.ceil(p + RT_PEAK)), (
        f"99.9 % of the {len(spread)} precursors identified in two or more runs lie within "
        f"{p:.1f} s across runs, leaving out {len(out)} identification(s) more than "
        f"{RT_OUTLIER} robust SDs from their run's median shift; plus {RT_PEAK} s for the "
        f"peak, at least {RT_MIN} s")


def bare(pf):
    return re.sub(r"\[[^\]]*\]", "", pf)


def is_encyclopedia_library(path):
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        return "entries" in {r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
    finally:
        con.close()


def library_mz(path):
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        return {(pf, int(z)): float(mz) for pf, z, mz in con.execute(
            "SELECT PeptideModSeq, PrecursorCharge, PrecursorMz FROM entries")}
    finally:
        con.close()


def subset_library(src, out, keys):
    if os.path.exists(out):
        os.remove(out)
    con = sqlite3.connect(out)
    try:
        con.execute("ATTACH DATABASE ? AS src", (src,))
        con.execute("CREATE TEMP TABLE keys (seq TEXT, z INTEGER)")
        con.executemany("INSERT INTO keys VALUES (?, ?)", sorted(keys))
        con.execute("CREATE TABLE entries AS SELECT e.* FROM src.entries e "
                    "JOIN keys k ON e.PeptideSeq = k.seq AND e.PrecursorCharge = k.z "
                    "ORDER BY e.rowid")
        n = con.execute("SELECT COUNT(*) FROM entries").fetchone()[0]
        con.commit()
    finally:
        con.close()
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", required=True, help="directory of per-run identification tables")
    ap.add_argument("--format", required=True, choices=("encyclopedia", "pyprophet"))
    ap.add_argument("--fdr", type=float, required=True)
    ap.add_argument("--library", help=".dlib/.elib the search used, for precursor m/z "
                                      "(encyclopedia format)")
    ap.add_argument("--prior", help="library diathem takes its fragment-ratio prior from")
    ap.add_argument("--unimod", required=True, help="assets/unimod_delta_masses.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--prior-out-dir", required=True)
    ap.add_argument("--rt-window-out", required=True, help="file for diathem's --rt-window")
    a = ap.parse_args()

    paths = sorted(p for p in glob.glob(os.path.join(a.ids, "*")) if os.path.isfile(p))
    if not paths:
        sys.exit(f"select_diathem_targets: no identification tables in {a.ids}")

    if a.format == "encyclopedia":
        hits = from_encyclopedia(paths, a.fdr)
        mz = library_mz(a.library)
        missing = [k for k in hits if k not in mz]
        if missing:
            print(f"{len(missing)} identified precursor(s) not in {a.library}, "
                  f"dropped, e.g. {missing[:3]}", file=sys.stderr)
        hits = {k: v for k, v in hits.items() if k in mz}
    else:
        with open(a.unimod) as f:
            deltas = {r["accession"].split(":")[-1]: float(r["monoisotopic_delta"])
                      for r in csv.DictReader(f, delimiter="\t")}
        hits, mz = from_pyprophet(paths, a.fdr, deltas)

    if not hits:
        sys.exit(f"select_diathem_targets: no precursor passed q <= {a.fdr} in any of "
                 f"{len(paths)} run(s), so diathem has nothing to quantify. Check the "
                 f"search, or set run_diathem = false.")

    with open(a.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Compound", "m/z", "z", "RT Time (min)"])
        for (pf, z) in sorted(hits):
            w.writerow([pf, "%.6f" % mz[(pf, z)], z,
                        "%.4f" % (statistics.median(hits[(pf, z)].values()) / 60.0)])
    print(f"wrote {len(hits)} target precursor(s) identified in at least one of "
          f"{len(paths)} run(s) to {a.out}", file=sys.stderr)

    rtw, why = rt_window(hits, len(paths))
    with open(a.rt_window_out, "w") as f:
        f.write(f"{rtw}\n")
    print(f"RT window: +/-{rtw} s ({why})", file=sys.stderr)

    if a.prior and os.path.basename(a.prior) != "NO_FILE":
        os.makedirs(a.prior_out_dir, exist_ok=True)
        if is_encyclopedia_library(a.prior):
            out = os.path.join(a.prior_out_dir, "diathem_prior.dlib")
            n = subset_library(a.prior, out, {(bare(pf), z) for pf, z in hits})
            print(f"prior: {n} entries of {a.prior} for the targets -> {out}", file=sys.stderr)
        else:
            out = os.path.join(a.prior_out_dir, os.path.basename(a.prior))
            os.symlink(os.path.realpath(a.prior), out)
            print(f"prior: {a.prior} used whole (not an EncyclopeDIA library)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
