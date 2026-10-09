#!/usr/bin/env python3
import argparse, hashlib, os, re, shutil, sqlite3, sys


CHUNK = 64 << 20
EPS = 0.01
R = dict(lvl=re.compile(rb'accession="MS:1000511"[^>]*value="([^"]+)"'),
         tgt=re.compile(rb'accession="MS:1000827"[^>]*value="([^"]+)"'),
         lo=re.compile(rb'accession="MS:1000828"[^>]*value="([^"]+)"'),
         up=re.compile(rb'accession="MS:1000829"[^>]*value="([^"]+)"'),
         rt=re.compile(rb'accession="MS:1000016"[^>]*value="([^"]+)"[^>]*unitAccession="([^"]+)"'))
SPECTRUM = re.compile(rb'<spectrum\s')
INDEX = re.compile(rb'(<spectrum\s[^>]*?\bindex=")(\d+)(")')
COUNT = re.compile(rb'(<spectrumList\s[^>]*?\bcount=")(\d+)(")')
TAG_ID = re.compile(rb'\bid="([^"]*)"')
CHROM = re.compile(rb'<chromatogram\s[^>]*')
HEADER = 'isolation_lower\tisolation_upper\twidth\tscans\trt_start_min\trt_end_min\tsurvey\n'


def key(lo, hi):
    return round(lo, 2), round(hi, 2)


def ms2_scans(block):
    m = R['lvl'].search(block)
    if not m or m.group(1) != b'2':
        return None
    t, lo, up = R['tgt'].search(block), R['lo'].search(block), R['up'].search(block)
    if not (t and lo and up):
        return None
    c = float(t.group(1))
    rt = R['rt'].search(block)
    sec = float(rt.group(1)) * (60 if rt.group(2) == b'UO:0000031' else 1) if rt else float('nan')
    return c, c - float(lo.group(1)), c + float(up.group(1)), sec


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


def mzml_scans(path):
    for p in spectra(path):
        s = ms2_scans(p)
        if s:
            yield s


def dia_scans(path):
    con = sqlite3.connect(path)
    for lo, c, hi, rt in con.execute('SELECT IsolationWindowLower, IsolationWindowCenter, '
                                     'IsolationWindowUpper, ScanStartTime FROM spectra'):
        yield c, lo, hi, rt
    con.close()


def windows(scans):
    win = {}
    for c, lo, hi, rt in scans:
        w = win.setdefault(key(lo, hi), [0, rt, rt])
        w[0] += 1
        w[1], w[2] = min(w[1], rt), max(w[2], rt)
    return win


def survey_windows(win):
    keys = sorted(win, key=lambda k: k[1] - k[0])
    out = set()
    for i, (lo, hi) in enumerate(keys):
        w = hi - lo
        for lo2, hi2 in keys[:i]:
            if (hi2 - lo2) * 4 <= w and lo - EPS <= lo2 and hi2 <= hi + EPS:
                out.add((lo, hi))
                break
    return out


def write_table(path, win, survey):
    with open(path, 'w') as fh:
        fh.write(HEADER)
        for (lo, hi), (n, a, b) in sorted(win.items()):
            fh.write('%.2f\t%.2f\t%.2f\t%d\t%.2f\t%.2f\t%s\n'
                     % (lo, hi, hi - lo, n, a / 60, b / 60, 'yes' if (lo, hi) in survey else 'no'))


def filter_mzml(src, dst, drop, n_keep):
    sha = hashlib.sha1()
    out = open(dst, 'wb')
    pos = 0

    def write(b):
        nonlocal pos
        out.write(b)
        sha.update(b)
        pos += len(b)

    offsets, index, first = [], 0, True
    parts = spectra(src)
    tail = b''
    for p in parts:
        m = SPECTRUM.search(p)
        if m is None:
            tail = p
            break
        head, block = p[:m.start()], p[m.start():] + b'</spectrum>'
        if first:
            head = COUNT.sub(lambda x: x.group(1) + str(n_keep).encode() + x.group(3), head, count=1)
            first = False
        write(head)
        s = ms2_scans(block)
        if s and key(s[1], s[2]) in drop:
            continue
        open_tag = block[:block.index(b'>') + 1]
        offsets.append((TAG_ID.search(open_tag).group(1), pos))
        block = INDEX.sub(lambda x: x.group(1) + str(index).encode() + x.group(3), block, count=1)
        index += 1
        write(block)

    cut = tail.find(b'<indexList')
    if cut < 0:
        write(tail)
        out.close()
        return
    body = tail[:cut]
    chrom = [(TAG_ID.search(t.group(0)).group(1), pos + t.start()) for t in CHROM.finditer(body)]
    write(body)
    index_at = pos
    lines = [b'<indexList count="%d">\n' % (2 if chrom else 1), b'    <index name="spectrum">\n']
    lines += [b'      <offset idRef="%s">%d</offset>\n' % (i, o) for i, o in offsets]
    lines.append(b'    </index>\n')
    if chrom:
        lines.append(b'    <index name="chromatogram">\n')
        lines += [b'      <offset idRef="%s">%d</offset>\n' % (i, o) for i, o in chrom]
        lines.append(b'    </index>\n')
    lines.append(b'  </indexList>\n  <indexListOffset>%d</indexListOffset>\n  <fileChecksum>' % index_at)
    write(b''.join(lines))
    out.write(sha.hexdigest().encode() + b'</fileChecksum>\n</indexedmzML>\n')
    out.close()


def filter_dia(src, dst, drop):
    shutil.copyfile(src, dst)
    con = sqlite3.connect(dst)
    for lo, hi in drop:
        for table, a, b in (('spectra', 'IsolationWindowLower', 'IsolationWindowUpper'),
                            ('ranges', 'Start', 'Stop')):
            con.execute('DELETE FROM %s WHERE abs(%s - ?) < ? AND abs(%s - ?) < ?' % (table, a, b),
                        (lo, EPS, hi, EPS))
    con.commit()
    con.execute('VACUUM')
    con.close()


def count_spectra(path):
    with open(path, 'rb') as fh:
        head = fh.read(1 << 20)
    m = COUNT.search(head)
    if m:
        return int(m.group(2))
    return sum(1 for p in spectra(path) if SPECTRUM.search(p))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', required=True, help='mzML or EncyclopeDIA .dia')
    ap.add_argument('--out', required=True, help='the run without its survey scans, or a link to it')
    ap.add_argument('--table', required=True, help='isolation windows, TSV')
    a = ap.parse_args()

    name = os.path.basename(a.run)
    low = name.lower()
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    if low.endswith('.mzml'):
        win = windows(mzml_scans(a.run))
    elif low.endswith('.dia'):
        win = windows(dia_scans(a.run))
    else:
        print('%s: not mzML or .dia; isolation windows not inspected' % name)
        open(a.table, 'w').write(HEADER)
        os.symlink(os.path.realpath(a.run), a.out)
        return 0
    survey = survey_windows(win)
    write_table(a.table, win, survey)
    total = sum(n for n, _, _ in win.values())
    print('%s: %d MS2 scans in %d isolation windows' % (name, total, len(win)))
    if not survey:
        print('no survey windows; the run is used as it is')
        os.symlink(os.path.realpath(a.run), a.out)
        return 0
    widths = sorted(set(round(hi - lo, 2) for lo, hi in survey))
    removed = sum(win[k][0] for k in survey)
    print('%d survey windows (%s Th wide) enclose narrower ones; removing the %d survey scans (%.1f %%)'
          % (len(survey), ', '.join('%g' % w for w in widths), removed, 100 * removed / total))
    if low.endswith('.dia'):
        filter_dia(a.run, a.out, survey)
    else:
        filter_mzml(a.run, a.out, survey, count_spectra(a.run) - removed)
    return 0


if __name__ == '__main__':
    sys.exit(main())
