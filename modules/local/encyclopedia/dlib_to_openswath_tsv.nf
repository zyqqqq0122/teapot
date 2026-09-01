process DLIB_TO_OPENSWATH_TSV {
    label 'python'
    tag   "${dlib.baseName}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/prepare_library_openswath" },
               mode: params.publish_mode

    input:
    path dlib
    path unimod_table

    output:
    path "${dlib.baseName}.osw.tsv",             emit: tsv
    path "${dlib.baseName}_unmapped_mods.tsv",   emit: unmapped
    path "${dlib.baseName}.dlib_to_tsv.log",     emit: log
    path "versions.yml",                         emit: versions

    script:
    """
    python3 - <<'PY' > ${dlib.baseName}.dlib_to_tsv.log 2>&1
    import sqlite3, struct, zlib, re, sys, csv, os

    DLIB          = "${dlib}"
    UNIMOD_TABLE  = "${unimod_table}"
    OUT           = "${dlib.baseName}.osw.tsv"
    UNMAPPED_OUT  = "${dlib.baseName}_unmapped_mods.tsv"
    MASS_TOL      = float("${params.mod_mass_tolerance}")
    MAX_UNMAPPED  = float("${params.unmapped_mod_max_fraction}")

    unimod = []
    with open(UNIMOD_TABLE) as f:
        r = csv.DictReader(f, delimiter="\\t")
        for row in r:
            unimod.append({
                "accession": row["accession"],
                "delta":     float(row["monoisotopic_delta"]),
                "residues":  set(row["target_residues"].split(",")),
            })

    def match_unimod(delta, residue):
        cands = [u for u in unimod if abs(u["delta"] - delta) <= MASS_TOL]
        if not cands:
            return None
        residue_hits = [u for u in cands if residue in u["residues"]]
        if residue_hits:
            residue_hits.sort(key=lambda u: abs(u["delta"] - delta))
            return residue_hits[0]["accession"]
        cands.sort(key=lambda u: abs(u["delta"] - delta))
        return cands[0]["accession"]

    MOD_RE = re.compile(r'\\[([+\\-][0-9.]+)\\]')

    def to_unimod(mod_seq):
        s = mod_seq
        out = []
        i = 0
        pending_nterm = None
        m = re.match(r'^\\[([+\\-][0-9.]+)\\]', s)
        if m:
            pending_nterm = float(m.group(1))
            s = s[m.end():]
        while i < len(s):
            ch = s[i]
            if s[i:i+1] == '[':
                return None, [("MID", 0.0)]
            i += 1
            if i < len(s) and s[i] == '[':
                m = re.match(r'\\[([+\\-][0-9.]+)\\]', s[i:])
                if not m:
                    return None, [("MID", 0.0)]
                delta = float(m.group(1))
                acc = match_unimod(delta, ch)
                if acc is None:
                    return None, [(ch, delta)]
                out.append(f"{ch}({acc})")
                i += m.end()
            else:
                out.append(ch)
        prefix = ""
        if pending_nterm is not None:
            acc = match_unimod(pending_nterm, "N-term")
            if acc is None:
                return None, [("N-term", pending_nterm)]
            prefix = f".({acc})"
        return prefix + "".join(out), None

    def decode_doubles(blob, n):
        return struct.unpack('>%dd' % n, zlib.decompress(blob))
    def decode_floats(blob, n):
        return struct.unpack('>%df' % n, zlib.decompress(blob))

    con = sqlite3.connect(DLIB)
    cur = con.cursor()
    cur.execute("SELECT PrecursorMz, PrecursorCharge, PeptideModSeq, "
                "PeptideSeq, RTInSeconds, MassEncodedLength, MassArray, "
                "IntensityEncodedLength, IntensityArray FROM entries")

    pep2prot = {}
    try:
        cur2 = con.cursor()
        cur2.execute("SELECT PeptideSeq, ProteinAccession FROM peptidetoprotein")
        for pep, prot in cur2:
            pep2prot.setdefault(pep, set()).add(prot)
    except sqlite3.OperationalError:
        pass

    n_seen = n_ok = 0
    unmapped_rows = []
    best = {}
    for row in cur:
        n_seen += 1
        pmz, pz, mod_seq, seq, rt, mlen, mblob, ilen, iblob = row
        unimod_seq, err = to_unimod(mod_seq or seq)
        if unimod_seq is None:
            res, delta = err[0]
            unmapped_rows.append((seq, mod_seq, res, delta))
            continue
        n_ok += 1
        mzs    = decode_doubles(mblob, mlen // 8)
        intens = decode_floats (iblob, ilen // 4)
        total  = sum(i for i in intens if i > 0)
        key    = (unimod_seq, pz)
        if key not in best or total > best[key][0]:
            best[key] = (total, (pmz, pz, mod_seq, seq, rt, mzs, intens, unimod_seq))

    n_dup = n_ok - len(best)
    with open(OUT, 'w', newline='') as fout:
        w = csv.writer(fout, delimiter='\\t')
        w.writerow([
            'PrecursorMz', 'ProductMz', 'LibraryIntensity',
            'NormalizedRetentionTime', 'ProteinId',
            'PeptideSequence', 'ModifiedPeptideSequence',
            'PrecursorCharge', 'ProductCharge',
            'FragmentType', 'FragmentSeriesNumber',
            'TransitionGroupId', 'TransitionId', 'Decoy',
            'PeptideGroupLabel',
        ])
        for _total, (pmz, pz, mod_seq, seq, rt, mzs, intens, unimod_seq) in best.values():
            prot  = ','.join(sorted(pep2prot.get(seq, [seq]))) or seq
            tg_id = f"{unimod_seq}_{pz}"
            for k, (mz, inten) in enumerate(zip(mzs, intens)):
                if inten <= 0:
                    continue
                w.writerow([
                    f'{pmz:.6f}', f'{mz:.6f}', f'{inten:.4f}',
                    f'{(rt or 0):.4f}',
                    prot, seq, unimod_seq, pz, 1,
                    'y', 1, tg_id, f'{tg_id}_{k}', 0, tg_id,
                ])
    con.close()

    with open(UNMAPPED_OUT, 'w', newline='') as fum:
        w = csv.writer(fum, delimiter='\\t')
        w.writerow(['PeptideSequence', 'PeptideModSeq', 'unmapped_residue', 'delta_mass'])
        for row in unmapped_rows:
            w.writerow(row)

    frac_unmapped = 0.0 if n_seen == 0 else len(unmapped_rows) / n_seen
    print(f"entries seen:    {n_seen}", file=sys.stderr)
    print(f"entries mapped:  {n_ok}", file=sys.stderr)
    print(f"transition groups written: {len(best)} "
          f"({n_dup} duplicate per-source entries collapsed)", file=sys.stderr)
    print(f"entries dropped: {len(unmapped_rows)} ({100*frac_unmapped:.2f}%)", file=sys.stderr)
    if frac_unmapped > MAX_UNMAPPED:
        print(f"ERROR: unmapped fraction {frac_unmapped:.3f} exceeds "
              f"unmapped_mod_max_fraction={MAX_UNMAPPED}", file=sys.stderr)
        sys.exit(2)

    residual_re = re.compile(r'\\[[+\\-][0-9.]+\\]')
    with open(OUT) as f:
        for line_no, line in enumerate(f, 1):
            if residual_re.search(line):
                print(f"ERROR: residual delta-mass survived in {OUT}:{line_no}",
                      file=sys.stderr)
                sys.exit(3)
    print("residual delta-mass check: clean", file=sys.stderr)
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    touch ${dlib.baseName}.osw.tsv ${dlib.baseName}_unmapped_mods.tsv ${dlib.baseName}.dlib_to_tsv.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
