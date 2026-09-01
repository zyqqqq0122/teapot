#!/usr/bin/env python3
import argparse
import csv
import re
import shutil
import sqlite3
import statistics
import sys


def strip_mods(seq):
    if seq is None:
        return ""
    s = re.sub(r"\[[^\]]*\]", "", str(seq))
    s = re.sub(r"\([^)]*\)", "", s)
    return s.strip()


def load_anchors(path, conflict_tolerance = 1.0,
                 ):
    accum = {}
    has_standard_col = False
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if header and len(header) >= 3 and header[2].strip().lower() == "standard":
            has_standard_col = True
        for row in reader:
            if len(row) < 2:
                continue
            stripped = strip_mods(row[0])
            try:
                irt = float(row[1])
            except ValueError:
                continue
            standard = 0
            if has_standard_col and len(row) >= 3:
                try:
                    standard = int(row[2])
                except ValueError:
                    pass
            accum.setdefault(stripped, []).append((irt, standard))

    out = {}
    n_conflicts = 0
    for stripped, items in accum.items():
        vs = sorted(v for v, _ in items)
        median = vs[len(vs) // 2]
        spread = vs[-1] - vs[0]
        std1_vals = [v for v, s in items if s == 1]
        chosen = statistics.median(std1_vals) if std1_vals else median
        if spread > conflict_tolerance:
            n_conflicts += 1
            print(
                f"harmonize_rt WARNING: anchor {stripped} has {len(items)} "
                f"rows disagreeing by {spread:.2f} "
                f"({vs[0]:.2f}..{vs[-1]:.2f}); using "
                f"{'Standard=1 median' if std1_vals else 'median'} "
                f"{chosen:.2f}.",
                file=sys.stderr,
            )
        out[stripped] = chosen
    if n_conflicts:
        print(
            f"harmonize_rt: {n_conflicts} anchor(s) had conflicting values "
            f"across mod-strip collapses; see warnings above.",
            file=sys.stderr,
        )
    return out


def load_dlib_rts(dlib_path):
    con = sqlite3.connect(dlib_path)
    rows = con.execute(
        "SELECT rowid, PeptideSeq, RTInSeconds FROM entries"
    ).fetchall()
    con.close()
    return [(r[0], r[1], float(r[2]) / 60.0 if r[2] is not None else None)
            for r in rows]


def ols_fit(xs, ys):
    n = len(xs)
    if n < 2:
        raise ValueError("Need >= 2 anchors for OLS fit.")
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0:
        raise ValueError("Anchor x-values have zero variance.")
    slope = sxy / sxx
    intercept = my - slope * mx
    r2 = (sxy ** 2) / (sxx * syy) if syy > 0 else 0.0
    return slope, intercept, r2


def robust_fit(xs, ys, target_r2 = 0.98,
               n_mad = 3.0, min_keep_frac = 0.3,
               min_keep = 8, max_iter = 40
               ):
    keep_x, keep_y = list(xs), list(ys)
    slope, intercept, r2_raw = ols_fit(keep_x, keep_y)
    r2 = r2_raw
    floor = max(min_keep, int(round(min_keep_frac * len(xs))))

    for _ in range(max_iter):
        resid = [abs(y - (slope * x + intercept)) for x, y in zip(keep_x, keep_y)]
        srt = sorted(resid)
        med = srt[len(srt) // 2]
        devs = sorted(abs(r - med) for r in resid)
        mad = devs[len(devs) // 2] * 1.4826
        if mad <= 0:
            break
        cutoff = med + n_mad * mad
        nx = [x for x, r in zip(keep_x, resid) if r <= cutoff]
        ny = [y for y, r in zip(keep_y, resid) if r <= cutoff]
        if len(nx) < floor or len(nx) == len(keep_x):
            break
        keep_x, keep_y = nx, ny
        slope, intercept, r2 = ols_fit(keep_x, keep_y)

    while r2 < target_r2 and len(keep_x) > floor:
        resid = [abs(y - (slope * x + intercept)) for x, y in zip(keep_x, keep_y)]
        cut = sorted(resid)[max(1, int(0.90 * len(resid))) - 1]
        nx = [x for x, r in zip(keep_x, resid) if r <= cut]
        ny = [y for y, r in zip(keep_y, resid) if r <= cut]
        if len(nx) < floor or len(nx) == len(keep_x):
            break
        keep_x, keep_y = nx, ny
        slope, intercept, r2 = ols_fit(keep_x, keep_y)

    return slope, intercept, r2, len(keep_x), r2_raw


def median_iqr(values):
    if not values:
        return float("nan"), float("nan"), float("nan")
    vs = sorted(values)
    med = statistics.median(vs)
    q1 = statistics.median(vs[: len(vs) // 2]) if len(vs) >= 4 else vs[0]
    q3 = statistics.median(vs[(len(vs) + 1) // 2 :]) if len(vs) >= 4 else vs[-1]
    return med, q1, q3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="input_dlib", required=True)
    ap.add_argument("--out", dest="output_dlib", required=True)
    ap.add_argument("--anchors", required=False,
                    help="TSV of stripped_seq<TAB>iRT (from emit_irt_library.py). "
                         "Omit for --assume-koina-passthrough.")
    ap.add_argument("--provenance", required=True,
                    help="'target' | 'background', written into a new "
                         "column so ASSERT_RT_COHERENT can compare strata "
                         "post-union.")
    ap.add_argument("--fit-report", required=True,
                    help="Output path for <basename>.rt_fit.tsv.")
    ap.add_argument("--min-anchors", type=int, default=10,
                    help="Minimum overlapping anchor peptides for a fit.")
    ap.add_argument("--no-robust-fit", dest="robust_fit", action="store_false",
                    default=True,
                    help="Disable iteratively-trimmed fitting and use plain OLS "
                         "over every matched anchor. Robust fitting is on by "
                         "default: a clean anchor set converges immediately and "
                         "drops nothing, while a contaminated one would "
                         "otherwise hard-fail on --min-r2.")
    ap.add_argument("--robust-min-keep", type=float, default=0.3,
                    help="Robust fitting will never trim below this fraction of "
                         "the matched anchors (default 0.3). The floor is what "
                         "stops it manufacturing a perfect fit from a handful "
                         "of points.")
    ap.add_argument("--robust-mad", type=float, default=3.0,
                    help="Outlier cutoff for robust fitting, in MADs of the "
                         "residual (default 3.0).")
    ap.add_argument("--min-r2", type=float, default=0.98,
                    help="Reject fits with R^2 below this, which usually means "
                         "IrtLibrary anchors are not appropriate for this "
                         "library.")
    ap.add_argument("--min-anchor-coverage", type=float, default=0.8,
                    help="Fraction of the library's RT range that must be "
                         "covered by the anchor peptides' RT range. A "
                         "near-perfect R^2 over the middle of the gradient "
                         "tells you nothing about the ends. LC non-linearity "
                         "is worst exactly where there are no anchors. Below "
                         "this coverage a loud WARNING is emitted and a "
                         "per-entry `extrapolated` column is added to the fit "
                         "report; the run is not failed, since anchor-coverage "
                         "gaps are common with 10-anchor kits.")
    ap.add_argument("--fallback", choices=("fail", "full_window"),
                    default="fail",
                    help="What to do when fit is not possible.")
    ap.add_argument("--assume-koina-passthrough", action="store_true",
                    help="Skip fit; assume input RT is already iRT. Still "
                         "sanity-checks the range.")
    ap.add_argument("--koina-min-rt", type=float, default=-150.0)
    ap.add_argument("--koina-max-rt", type=float, default=400.0)
    ap.add_argument("--emit-full-window-sentinel", action="store_true",
                    help="Also write the `rt_extraction_full_window` sentinel "
                         "unconditionally. Used by the Koina-only OSW route "
                         "when no --irt_library is supplied. The library RT "
                         "is iRT but we have no -tr_irt to fit "
                         "library_iRT -> run_RT, so OPENSWATH_WORKFLOW must "
                         "fall back to full-gradient extraction.")
    ap.add_argument("--full-window-warn-size", type=int, default=20000,
                    help="Emit a loud warning when the full-window sentinel "
                         "is engaged AND the library has more than this many "
                         "precursor entries. -rt_extraction_window -1 extracts "
                         "across the entire gradient for every library entry. "
                         "On proteome-wide predicted libraries that is a very "
                         "large amount of chromatogram work per sample.")
    args = ap.parse_args()

    shutil.copy(args.input_dlib, args.output_dlib)

    con = sqlite3.connect(args.output_dlib)
    cur = con.cursor()

    cols = [r[1] for r in cur.execute("PRAGMA table_info(entries)").fetchall()]
    if "provenance" not in cols:
        cur.execute("ALTER TABLE entries ADD COLUMN provenance TEXT")
    cur.execute("UPDATE entries SET provenance = ?", (args.provenance,))
    con.commit()

    rows = load_dlib_rts(args.input_dlib)
    all_rt_min = [r[2] for r in rows if r[2] is not None]
    if not all_rt_min:
        con.close()
        sys.exit("HARMONIZE_RT: dlib has no RTInSeconds values to harmonize.")

    med, q1, q3 = median_iqr(all_rt_min)
    print(
        f"HARMONIZE_RT: input dlib RT stats (min): "
        f"n={len(all_rt_min)}  median={med:.2f}  IQR=[{q1:.2f}, {q3:.2f}]  "
        f"range=[{min(all_rt_min):.2f}, {max(all_rt_min):.2f}]",
        file=sys.stderr,
    )

    if args.assume_koina_passthrough:
        undo_rows = [(r[0], (r[2] * 60.0) if r[2] is not None else None) for r in rows]
        pseudo_range = [r[1] for r in undo_rows if r[1] is not None]
        if pseudo_range and not (args.koina_min_rt <= min(pseudo_range) and
                                  max(pseudo_range) <= args.koina_max_rt):
            print(
                f"HARMONIZE_RT WARNING: passthrough assumed iRT but values "
                f"[{min(pseudo_range):.2f}, {max(pseudo_range):.2f}] fall "
                f"outside expected Koina iRT range "
                f"[{args.koina_min_rt}, {args.koina_max_rt}]. Continuing "
                f"but check upstream.",
                file=sys.stderr,
            )
        with open(args.fit_report, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["mode", "passthrough"])
            w.writerow(["n_entries", len(rows)])
            w.writerow(["median_rt_stored", med])
        if args.emit_full_window_sentinel:
            with open("rt_extraction_full_window", "w") as sf:
                sf.write("koina_passthrough_no_tr_irt\n")
            print(
                "HARMONIZE_RT: emitted full-window sentinel, Koina-only OSW "
                "route without --irt_library; OpenSwathWorkflow will use "
                "-rt_extraction_window -1.",
                file=sys.stderr,
            )
            if len(rows) > args.full_window_warn_size:
                print(
                    f"HARMONIZE_RT WARNING: full-window fallback engaged on a "
                    f"large library ({len(rows)} precursors > "
                    f"{args.full_window_warn_size}). "
                    f"-rt_extraction_window -1 will extract across the entire "
                    f"gradient for every one of them, per sample. This can "
                    f"be prohibitive on proteome-wide predicted libraries. "
                    f"Consider supplying --irt_library + --irt_anchors so "
                    f"HARMONIZE_RT can regress the library onto a proper iRT "
                    f"scale instead. Extraction then uses a narrow window "
                    f"around each entry's fitted retention time.",
                    file=sys.stderr,
                )
        con.close()
        return 0

    if not args.anchors:
        con.close()
        sys.exit("HARMONIZE_RT: --anchors is required unless "
                 "--assume-koina-passthrough is set.")

    anchors = load_anchors(args.anchors)
    n_anchors_loaded = len(anchors)
    print(f"HARMONIZE_RT: loaded {n_anchors_loaded} anchor peptides.",
          file=sys.stderr)

    matches = []
    for _rowid, seq, rt_min in rows:
        if rt_min is None:
            continue
        stripped = strip_mods(seq)
        if stripped in anchors:
            matches.append((stripped, rt_min, anchors[stripped]))
    by_pep = {}
    for stripped, x, y in matches:
        by_pep.setdefault(stripped, []).append((x, y))
    xs, ys = [], []
    anchor_report_rows = []
    for stripped, pts in by_pep.items():
        xm = statistics.median(p[0] for p in pts)
        ym = statistics.median(p[1] for p in pts)
        xs.append(xm)
        ys.append(ym)
        anchor_report_rows.append((stripped, xm, ym))

    n_anchors = len(xs)
    matched_frac = (n_anchors / n_anchors_loaded) if n_anchors_loaded > 0 else 0.0
    print(
        f"HARMONIZE_RT: {n_anchors} of {n_anchors_loaded} anchor peptides "
        f"overlap the library ({matched_frac:.1%}).",
        file=sys.stderr,
    )
    if n_anchors_loaded > 0 and matched_frac < 0.05:
        print(
            f"HARMONIZE_RT WARNING: only {matched_frac:.1%} of anchors matched "
            f"the library. That is unexpectedly low. Check that the anchors "
            f"TSV was generated from a compatible library (matching organism, "
            f"digestion, mod conventions). Stripped-sequence match is "
            f"mod-agnostic, so the drop is not from mod-notation differences.",
            file=sys.stderr,
        )

    if n_anchors < args.min_anchors:
        with open(args.fit_report, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["status", "insufficient_anchors"])
            w.writerow(["n_anchors_loaded", n_anchors_loaded])
            w.writerow(["n_anchors_matched", n_anchors])
            w.writerow(["n_anchors_matched_fraction",
                        f"{n_anchors / n_anchors_loaded:.4f}"
                        if n_anchors_loaded > 0 else "nan"])
            w.writerow(["min_anchors_required", args.min_anchors])
        if args.fallback == "fail":
            con.close()
            sys.exit(
                f"HARMONIZE_RT HARD ERROR: {n_anchors} anchor peptides "
                f"overlap this library (need >= {args.min_anchors}). "
                f"Cannot fit vendor_RT -> iRT regression. "
                f"Either supply --irt_library with anchors that overlap, "
                f"or set rt_harmonize_fallback='full_window' to fall back "
                f"to `-rt_extraction_window -1` in OpenSwathWorkflow (loud "
                f"warning; correct but expensive)."
            )
        with open("rt_extraction_full_window", "w") as f:
            f.write(f"insufficient_anchors n_matched={n_anchors}\n")
        print(
            f"HARMONIZE_RT: fallback=full_window engaged; "
            f"library RT column left unchanged. "
            f"OpenSwathWorkflow must set `-rt_extraction_window -1` for "
            f"this library.",
            file=sys.stderr,
        )
        con.close()
        return 0

    if args.robust_fit:
        slope, intercept, r2, n_kept, r2_raw = robust_fit(
            xs, ys, target_r2=args.min_r2, n_mad=args.robust_mad,
            min_keep_frac=args.robust_min_keep, min_keep=args.min_anchors)
        print(
            f"HARMONIZE_RT: fit iRT = {slope:.4f} * RT_min + {intercept:.4f}  "
            f"R^2 = {r2:.4f}  (robust: kept {n_kept} of {len(xs)} anchors; "
            f"R^2 over all anchors was {r2_raw:.4f})",
            file=sys.stderr,
        )
        if n_kept < len(xs):
            dropped = len(xs) - n_kept
            print(
                f"HARMONIZE_RT: dropped {dropped} anchor(s) as outliers "
                f"({100.0 * dropped / len(xs):.0f}%). A large fraction usually "
                f"means the anchor set does not suit this library. Check that "
                f"the anchors and the library come from comparable "
                f"chromatography.",
                file=sys.stderr,
            )
    else:
        slope, intercept, r2 = ols_fit(xs, ys)
        n_kept, r2_raw = len(xs), r2
        print(
            f"HARMONIZE_RT: fit iRT = {slope:.4f} * RT_min + {intercept:.4f}  "
            f"R^2 = {r2:.4f}",
            file=sys.stderr,
        )

    if r2 < args.min_r2:
        with open(args.fit_report, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["status", "low_r2"])
            w.writerow(["slope", slope])
            w.writerow(["intercept", intercept])
            w.writerow(["r2", r2])
            w.writerow(["min_r2", args.min_r2])
        if args.fallback == "fail":
            con.close()
            sys.exit(
                f"HARMONIZE_RT HARD ERROR: fit R^2={r2:.4f} below "
                f"threshold {args.min_r2}. Anchor peptides do not fit a "
                f"linear model. Inspect {args.fit_report}."
            )
        with open("rt_extraction_full_window", "w") as f:
            f.write(f"low_r2 r2={r2:.4f}\n")
        con.close()
        return 0

    cols = [r[1] for r in cur.execute("PRAGMA table_info(entries)").fetchall()]
    if "RTInSeconds_orig" not in cols:
        cur.execute("ALTER TABLE entries ADD COLUMN RTInSeconds_orig REAL")
    cur.execute("UPDATE entries SET RTInSeconds_orig = RTInSeconds")

    cur.execute(
        "UPDATE entries SET RTInSeconds = "
        "  (? * (RTInSeconds / 60.0) + ?)",
        (slope, intercept),
    )
    con.commit()
    con.close()

    residuals = [(y - (slope * x + intercept)) for x, y in zip(xs, ys)]

    anchor_rt_min, anchor_rt_max = min(xs), max(xs)
    lib_rt_min, lib_rt_max       = min(all_rt_min), max(all_rt_min)
    lib_range     = lib_rt_max - lib_rt_min
    anchor_range  = anchor_rt_max - anchor_rt_min
    coverage      = (anchor_range / lib_range) if lib_range > 0 else 1.0

    n_extrapolated = sum(
        1 for r in all_rt_min if r < anchor_rt_min or r > anchor_rt_max
    )

    coverage_warning = coverage < args.min_anchor_coverage
    if coverage_warning:
        print(
            f"HARMONIZE_RT WARNING: anchor RT coverage {coverage:.2%} is "
            f"below --min-anchor-coverage {args.min_anchor_coverage:.2%}. "
            f"Anchors span RT=[{anchor_rt_min:.2f}, {anchor_rt_max:.2f}] "
            f"vs library RT=[{lib_rt_min:.2f}, {lib_rt_max:.2f}] "
            f"(min). {n_extrapolated} of {len(all_rt_min)} library entries "
            f"({100 * n_extrapolated / len(all_rt_min):.1f}%) fall outside "
            f"the anchor range and are being mapped by EXTRAPOLATION. "
            f"LC non-linearity is worst at gradient ends. Extrapolated "
            f"entries carry the largest RT uncertainty. See fit report "
            f"for the per-entry `extrapolated` flag.",
            file=sys.stderr,
        )

    with open(args.fit_report, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["status", "ok" if not coverage_warning else "ok_low_coverage"])
        w.writerow(["slope", slope])
        w.writerow(["intercept", intercept])
        w.writerow(["r2", r2])
        w.writerow(["n_anchors_loaded", n_anchors_loaded])
        w.writerow(["n_anchors_matched", n_anchors])
        w.writerow(["n_anchors_matched_fraction",
                    f"{n_anchors / n_anchors_loaded:.4f}"
                    if n_anchors_loaded > 0 else "nan"])
        w.writerow(["anchor_rt_min", anchor_rt_min])
        w.writerow(["anchor_rt_max", anchor_rt_max])
        w.writerow(["library_rt_min", lib_rt_min])
        w.writerow(["library_rt_max", lib_rt_max])
        w.writerow(["anchor_coverage", f"{coverage:.4f}"])
        w.writerow(["min_anchor_coverage", f"{args.min_anchor_coverage:.4f}"])
        w.writerow(["n_extrapolated_entries", n_extrapolated])
        w.writerow(["residual_median", statistics.median(residuals)])
        w.writerow(["residual_iqr",
                    (max(residuals) - min(residuals))])
        w.writerow(["--", "--"])
        w.writerow(["stripped_seq", "library_rt_min", "iRT_anchor", "extrapolated"])
        for row in anchor_report_rows:
            extrapolated = (
                row[1] < anchor_rt_min or row[1] > anchor_rt_max
            )
            w.writerow(list(row) + [str(extrapolated).lower()])

    return 0


if __name__ == "__main__":
    sys.exit(main())
