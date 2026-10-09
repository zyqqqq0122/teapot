#!/usr/bin/env python3
import argparse, base64, re, struct, sys, zlib


R = {
    'lvl': re.compile(rb'accession="MS:1000511"[^>]*value="([^"]+)"'),
    'lo':  re.compile(rb'accession="MS:1000828"[^>]*value="([^"]+)"'),
    'up':  re.compile(rb'accession="MS:1000829"[^>]*value="([^"]+)"'),
    'bin': re.compile(rb'<binaryDataArray.*?</binaryDataArray>', re.S),
    'b64': re.compile(rb'<binary>([^<]*)</binary>', re.S),
}
IS_MZ, IS_INT = re.compile(rb'accession="MS:1000514"'), re.compile(rb'accession="MS:1000515"')
IS64,  ISZ    = re.compile(rb'accession="MS:1000523"'), re.compile(rb'accession="MS:1000574"')


def decode(block):
    m = R['b64'].search(block)
    if not m:
        return []
    raw = base64.b64decode(m.group(1))
    if ISZ.search(block):
        raw = zlib.decompress(raw)
    fmt = 'd' if IS64.search(block) else 'f'
    n = len(raw) // struct.calcsize(fmt)
    return struct.unpack('<%d%s' % (n, fmt), raw[:n * struct.calcsize(fmt)])


def occupancy(mzs, tol, lo=300.0, hi=1500.0):
    inr = sorted(m for m in mzs if lo <= m <= hi)
    if not inr:
        return 0.0
    total = 0.0
    cs, ce = inr[0] - tol, inr[0] + tol
    for m in inr[1:]:
        a, b = m - tol, m + tol
        if a > ce:
            total += ce - cs
            cs, ce = a, b
        else:
            ce = max(ce, b)
    total += ce - cs
    return min(total / (hi - lo), 1.0)


def measure(path, tol, max_scans, width=None, chunk=32 << 20):
    widths, buf = {}, b''
    kept, pending = [], {}
    with open(path, 'rb') as fh:
        while True:
            data = fh.read(chunk)
            if not data:
                break
            buf += data
            parts = buf.split(b'</spectrum>')
            buf = parts.pop()
            for p in parts:
                lvl = R['lvl'].search(p)
                if not lvl or lvl.group(1) != b'2':
                    continue
                lo, up = R['lo'].search(p), R['up'].search(p)
                if not (lo and up):
                    continue
                w = round(float(lo.group(1)) + float(up.group(1)), 2)
                widths[w] = widths.get(w, 0) + 1
                if width is None:
                    b = pending.setdefault(w, [])
                    if len(b) < max_scans:
                        b.append(p)
                    for extra in sorted(pending)[3:]:
                        del pending[extra]
                elif abs(w - width) < 1e-6 and len(kept) < max_scans:
                    kept.append(p)

    if width is None:
        real = [w for w, n in widths.items() if n >= 5] or list(widths)
        width = min(real) if real else None
        kept = pending.get(width, [])[:max_scans]

    occs = []
    for p in kept:
        mzs = None
        for b in R['bin'].findall(p):
            if IS_MZ.search(b):
                mzs = decode(b)
                break
        if mzs:
            occs.append(occupancy(mzs, tol))
    return occs, widths, width


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mzml', required=True)
    ap.add_argument('--tolerance', type=float, required=True,
                    help='extraction tolerance in Th, the same one OpenSWATH will use')
    ap.add_argument('--max-occupancy', type=float, default=0.35,
                    help='fail above this median occupancy (default 0.35)')
    ap.add_argument('--scans', type=int, default=1200,
                    help='MS2 scans to sample (default 1200)')
    ap.add_argument('--warn-only', action='store_true',
                    help='report and pass regardless')
    ap.add_argument('--window-width', type=float, default=None,
                    help='isolation width to measure, in Th; default is the '
                         'narrowest present -- the scheduled/targeted windows')
    ap.add_argument('--report', default=None)
    a = ap.parse_args()

    occs, widths, width = measure(a.mzml, a.tolerance, a.scans, a.window_width)
    lines = ["OpenSWATH feasibility -- %s" % a.mzml.split('/')[-1], ""]

    if not occs:
        lines.append("No MS2 scans with a decodable m/z array were found.")
        verdict, ok = "INDETERMINATE", True   # not our call to make; let it run
    else:
        occs.sort()
        med = occs[len(occs) // 2]
        lines += [
            "  scans sampled        : %d at %.2f Th windows" % (len(occs), width),
            "  extraction tolerance : +/-%.3f Th" % a.tolerance,
            "  median occupancy     : %.1f %%" % (100 * med),
            "  threshold            : %.1f %%" % (100 * a.max_occupancy),
            "",
            "  a 10-transition decoy is fully supported by chance with",
            "  probability %.4g" % (med ** 10),
            "",
        ]
        ok = med < a.max_occupancy
        verdict = "OK" if ok else "TOO CROWDED"

    if widths:
        total = sum(widths.values())
        lines.append("  isolation windows:")
        for w in sorted(widths):
            lines.append("    %7.2f Th : %6d scans (%4.1f %%)"
                         % (w, widths[w], 100 * widths[w] / total))
    lines += ["", "VERDICT: %s" % verdict]

    report = "\n".join(lines)
    print(report)
    if a.report:
        with open(a.report, 'w') as fh:
            fh.write(report + "\n")

    if not ok and not a.warn_only:
        sys.stderr.write(
            "\nOpenSWATH would not produce a trustworthy FDR on this acquisition:\n"
            "its decoys are extracted from the same scans as their targets and, at\n"
            "%.1f %% occupancy, land on real signal often enough to look genuine.\n\n"
            "Either run this sample with engine=encyclopedia, which selects scans by\n"
            "isolation window and does not depend on decoy transitions missing, or\n"
            "override with --openswath_max_occupancy / --skip_openswath_feasibility\n"
            "if you know why this acquisition is an exception.\n"
            % (100 * (occs[len(occs) // 2] if occs else 0)))
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
