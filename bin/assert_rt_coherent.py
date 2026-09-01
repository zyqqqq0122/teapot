#!/usr/bin/env python3
import argparse
import csv
import sqlite3
import statistics
import sys


def median_iqr(values):
    if not values:
        return (float("nan"),) * 5
    vs = sorted(values)
    med = statistics.median(vs)
    half = len(vs) // 2
    lower = vs[:half]
    upper = vs[-half:] if half > 0 else vs
    q1 = statistics.median(lower) if lower else vs[0]
    q3 = statistics.median(upper) if upper else vs[-1]
    return med, q1, q3, min(vs), max(vs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--union", required=True,
                    help="Path to the union .dlib after MERGE_LIBRARIES.")
    ap.add_argument("--report", required=True,
                    help="Output TSV with per-stratum stats.")
    ap.add_argument("--max-range-ratio", type=float, default=3.0,
                    help="Fail if max(range(target), range(background)) / "
                         "min(...) exceeds this.")
    ap.add_argument("--require-iqr-overlap", type=int, default=1,
                    help="If 1, require the IQRs of the two strata to "
                         "overlap.")
    ap.add_argument("--allow-single-stratum", type=int, default=0,
                    help="If 1, pass through when only one provenance is "
                         "present (i.e. no union happened).")
    args = ap.parse_args()

    con = sqlite3.connect(args.union)
    cols = [r[1] for r in con.execute("PRAGMA table_info(entries)").fetchall()]
    if "provenance" not in cols:
        if args.allow_single_stratum:
            con.close()
            with open(args.report, "w", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(["status", "no_provenance_column"])
                w.writerow(["note", "HARMONIZE_RT not detected on inputs; "
                                    "assuming intentional passthrough."])
            return 0
        con.close()
        sys.exit(
            "ASSERT_RT_COHERENT HARD ERROR: union dlib has no `provenance` "
            "column: HARMONIZE_RT did not run on either sub-library. "
            "Either configure HARMONIZE_RT correctly or set "
            "--allow-single-stratum if you know only one library was in the "
            "union."
        )

    strata = {}
    for prov, rt in con.execute(
            "SELECT provenance, RTInSeconds FROM entries "
            "WHERE RTInSeconds IS NOT NULL"):
        if prov is None:
            continue
        strata.setdefault(prov, []).append(float(rt) / 60.0)
    con.close()

    rows = []
    for prov, vals in strata.items():
        med, q1, q3, lo, hi = median_iqr(vals)
        rows.append({
            "provenance": prov, "n": len(vals),
            "median": med, "q1": q1, "q3": q3,
            "min": lo, "max": hi, "range": hi - lo,
        })

    with open(args.report, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["provenance", "n", "median", "q1", "q3",
                    "min", "max", "range"])
        for r in rows:
            w.writerow([r["provenance"], r["n"],
                        f'{r["median"]:.4f}', f'{r["q1"]:.4f}',
                        f'{r["q3"]:.4f}', f'{r["min"]:.4f}',
                        f'{r["max"]:.4f}', f'{r["range"]:.4f}'])

    for r in rows:
        print(f"  {r['provenance']:>12}: n={r['n']}  "
              f"median={r['median']:.2f}  "
              f"IQR=[{r['q1']:.2f}, {r['q3']:.2f}]  "
              f"range=[{r['min']:.2f}, {r['max']:.2f}]",
              file=sys.stderr)

    if len(rows) < 2:
        if args.allow_single_stratum:
            return 0
        sys.exit(
            "ASSERT_RT_COHERENT HARD ERROR: only one provenance stratum "
            "found in union. MERGE_LIBRARIES did not add both target and "
            "background entries, or one sub-library had no RTInSeconds."
        )

    tgt = next((r for r in rows if r["provenance"] == "target"), None)
    bg  = next((r for r in rows if r["provenance"] == "background"), None)
    if tgt is None or bg is None:
        sys.exit(
            "ASSERT_RT_COHERENT HARD ERROR: expected 'target' and "
            "'background' provenance strata, saw: "
            f"{[r['provenance'] for r in rows]}"
        )

    ranges = [tgt["range"], bg["range"]]
    if min(ranges) > 0 and max(ranges) / min(ranges) > args.max_range_ratio:
        sys.exit(
            f"ASSERT_RT_COHERENT HARD ERROR: range ratio exceeds "
            f"{args.max_range_ratio}. "
            f"target range = {tgt['range']:.2f}, "
            f"background range = {bg['range']:.2f}. "
            f"HARMONIZE_RT likely did not place both sub-libraries on the "
            f"same iRT scale. Inspect fit reports."
        )

    if args.require_iqr_overlap:
        max_q1 = max(tgt["q1"], bg["q1"])
        min_q3 = min(tgt["q3"], bg["q3"])
        if max_q1 > min_q3:
            sys.exit(
                f"ASSERT_RT_COHERENT HARD ERROR: IQR ranges disjoint. "
                f"target IQR=[{tgt['q1']:.2f}, {tgt['q3']:.2f}], "
                f"background IQR=[{bg['q1']:.2f}, {bg['q3']:.2f}]. "
                f"HARMONIZE_RT likely did not put both sub-libraries on "
                f"the same iRT scale: one appears to still be in minutes "
                f"while the other is in iRT."
            )

    print("ASSERT_RT_COHERENT: OK", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
