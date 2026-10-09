#!/usr/bin/env python3
import argparse, bisect, csv, glob, re, sys


CHUNK = 64 << 20
OWN = 0.05        # a target's own window is centred within this of its m/z
LVL = re.compile(rb'accession="MS:1000511"[^>]*value="([^"]+)"')
TGT = re.compile(rb'accession="MS:1000827"[^>]*value="([^"]+)"')
LO = re.compile(rb'accession="MS:1000828"[^>]*value="([^"]+)"')
UP = re.compile(rb'accession="MS:1000829"[^>]*value="([^"]+)"')


def spectra(path):
    buf = b''
    with open(path, 'rb') as fh:
        while True:
            data = fh.read(CHUNK)
            if not data:
                break
            buf += data
            parts = buf.split(b'</spectrum>')
            buf = parts.pop()
            yield from parts
    yield buf


def isolation_targets(path):
    for block in spectra(path):
        m = LVL.search(block)
        if not m or m.group(1) != b'2':
            continue
        t = TGT.search(block)
        if t and LO.search(block) and UP.search(block):
            yield float(t.group(1))


def target_mzs(path):
    try:
        with open(path, newline='') as fh:
            first = fh.readline()
            fh.seek(0)
            rows = list(csv.reader(fh, delimiter='\t' if first.count('\t') > first.count(',') else ','))
    except (OSError, UnicodeDecodeError):
        return []
    if not rows:
        return []
    h = [x.strip().lower() for x in rows[0]]
    if 'm/z' not in h:
        return []
    i = h.index('m/z')
    d = h.index('isdecoy') if 'isdecoy' in h else None
    out = []
    for r in rows[1:]:
        if len(r) <= i or (d is not None and len(r) > d and r[d].strip().lower() == 'true'):
            continue
        try:
            out.append(float(r[i]))
        except ValueError:
            pass
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--targets', required=True, help='target list with an m/z column')
    ap.add_argument('--fallback', type=float, required=True)
    ap.add_argument('runs', nargs='+', help='mzML files')
    a = ap.parse_args()
    runs = [r for p in a.runs for r in (sorted(glob.glob(p)) or [p])]

    mzs = target_mzs(a.targets)
    centres = set()
    for run in runs:
        centres.update(round(c, 5) for c in isolation_targets(run))
    centres = sorted(centres)
    say = lambda s: sys.stderr.write(s + '\n')
    if not mzs or not centres:
        say('precursor tolerance: no target m/z or no MS2 scans; using %g' % a.fallback)
        print(a.fallback)
        return 0
    offs = []
    for mz in mzs:
        i = bisect.bisect_left(centres, mz)
        offs.append(min(abs(centres[j] - mz) for j in (i - 1, i) if 0 <= j < len(centres)))
    own = [o for o in offs if o <= OWN]
    if len(own) < 0.5 * len(mzs):
        say('precursor tolerance: %d of %d targets have an isolation window centred on them; '
            'not a targeted acquisition, using %g' % (len(own), len(mzs), a.fallback))
        print(a.fallback)
        return 0
    # robust to targets the method did not acquire that sit near another window
    own.sort()
    tol = round(min(OWN, 2 * own[int(0.95 * (len(own) - 1))] + 0.0005), 4)
    shared = 0
    for mz in mzs:
        i = bisect.bisect_left(centres, mz - tol)
        if sum(1 for c in centres[i:i + 3] if abs(c - mz) <= tol) > 1:
            shared += 1
    beyond = sum(1 for o in offs if o > tol)
    say('precursor tolerance: %d of %d targets have an isolation window centred on them '
        '(95 %% within %.4f Th of the target m/z); using %g, which selects each target\'s own '
        'window%s%s' % (len(own), len(mzs), own[int(0.95 * (len(own) - 1))], tol,
                        '' if not shared else '; %d targets have a second window that close' % shared,
                        '' if not beyond else '; %d targets have no window that close and get no scans' % beyond))
    print(tol)
    return 0


if __name__ == '__main__':
    sys.exit(main())
