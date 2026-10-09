#!/usr/bin/env python3
import argparse, csv, io, os, re, shutil, sqlite3, sys, zipfile


csv.field_size_limit(sys.maxsize)
ATLAS = ('peptide.tsv', 'peptide_instance.tsv')
HEADERS = {'peptide', 'peptides', 'sequence', 'peptide_sequence', 'peptidesequence', 'stripped_sequence',
           'modified_sequence', 'peptidoform', 'compound'}


def is_atlas(path):
    if os.path.isdir(path):
        return all(os.path.exists(os.path.join(path, n)) for n in ATLAS)
    if zipfile.is_zipfile(path):
        names = set(zipfile.ZipFile(path).namelist())
        if not all(n in names for n in ATLAS):
            sys.exit('%s is a zip without %s: not the tables of a PeptideAtlas build'
                     % (path, ' and '.join(ATLAS)))
        return True
    return False


def table(atlas, name):
    if os.path.isdir(atlas):
        return open(os.path.join(atlas, name), newline='')
    zf = zipfile.ZipFile(atlas)
    return io.TextIOWrapper(zf.open(name), newline='')


def most_observed(atlas, top):
    seq = {r['peptide_id']: r['peptide_sequence']
           for r in csv.DictReader(table(atlas, 'peptide.tsv'), delimiter='\t')}
    obs = {}
    for r in csv.DictReader(table(atlas, 'peptide_instance.tsv'), delimiter='\t'):
        s = seq.get(r['peptide_id'])
        if s is None:
            continue
        try:
            n = float(r['n_adjusted_observations'])
        except ValueError:
            n = 0.0
        obs[s] = max(n, obs.get(s, 0.0))
    ranked = sorted(obs, key=lambda s: -obs[s])
    return set(ranked[:top]), len(obs), obs[ranked[min(top, len(ranked)) - 1]]


def bare(s):
    """the residues of a peptide written with flanks or modifications."""
    s = re.sub(r'^[A-Za-z\-]\.|\.[A-Za-z\-]$', '', s.strip().strip('"'))
    s = re.sub(r'\[[^\]]*\]|\([^)]*\)', '', s)
    return re.sub(r'[^A-Z]', '', s.upper())


def listed(path):
    with open(path, newline='') as fh:
        lines = fh.read().splitlines()
    if any(l.startswith('>') for l in lines):
        peps, cur = [], None
        for l in lines:
            if l.startswith('>'):
                cur = []
                peps.append(cur)
            elif cur is not None:
                cur.append(l.strip())
        peps = [''.join(p) for p in peps]
    else:
        peps = [re.split(r'[\t,; ]', l.strip(), 1)[0] for l in lines if l.strip()]
        if peps and peps[0].strip('"').lower() in HEADERS:
            peps = peps[1:]
    return {b for b in map(bare, peps) if b}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--peptides', required=True, help='a PeptideAtlas build\'s tables, or a list of peptides')
    ap.add_argument('--dlib', required=True, help='the predicted background')
    ap.add_argument('--top', type=int, default=5000, help='of a PeptideAtlas build, how many of the most observed (5000)')
    ap.add_argument('--out', required=True, help='the restricted background .dlib')
    a = ap.parse_args()

    if is_atlas(a.peptides):
        keep, n_atlas, floor = most_observed(a.peptides, a.top)
        what = ('%d peptides in the build; the %d most observed (at least %d observations)'
                % (n_atlas, len(keep), floor))
    else:
        keep = listed(a.peptides)
        what = '%d peptides listed in %s' % (len(keep), os.path.basename(a.peptides))
    shutil.copy(a.dlib, a.out)
    con = sqlite3.connect(a.out)
    con.execute('CREATE TEMP TABLE keep (s TEXT PRIMARY KEY)')
    con.executemany('INSERT INTO keep VALUES (?)', [(s,) for s in keep])
    con.execute('DELETE FROM entries WHERE PeptideSeq NOT IN (SELECT s FROM keep)')
    con.execute('DELETE FROM peptidetoprotein WHERE PeptideSeq NOT IN (SELECT s FROM keep)')
    con.commit()
    con.execute('VACUUM')
    n_entries, n_peptides = con.execute('SELECT COUNT(*), COUNT(DISTINCT PeptideSeq) FROM entries').fetchone()
    con.close()
    print('%s: %d of them in %s, %d precursors -> %s'
          % (what, n_peptides, os.path.basename(a.dlib), n_entries, a.out))
    if n_peptides == 0:
        sys.exit('none of these peptides is in %s; are they of the same organism as background_fasta?'
                 % os.path.basename(a.dlib))


if __name__ == '__main__':
    main()
