process DLIB_TO_REFERENCE_LIST {
    label 'python'
    tag   "${library.baseName}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/prepare_library" },
               mode: params.publish_mode

    input:
    path library
    val  source_kind

    output:
    path "${library.baseName}.reference_list.csv", emit: reference_list
    path "${library.baseName}.reference_list.log", emit: log
    path "versions.yml",                           emit: versions

    script:
    def keep_rt = (source_kind == 'elib') ? 1 : 0
    """
    python3 - <<'PY' > ${library.baseName}.reference_list.log 2>&1
    import sqlite3, csv, sys

    LIB     = "${library}"
    OUT     = "${library.baseName}.reference_list.csv"
    KEEP_RT = bool(${keep_rt})
    KIND    = "${source_kind}"

    con = sqlite3.connect(LIB)
    have = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    if 'entries' not in have:
        sys.exit("DLIB_TO_REFERENCE_LIST: %s has no 'entries' table; it is not a "
                 "readable EncyclopeDIA .dlib/.elib." % LIB)

    rows = list(con.execute(
        "SELECT PeptideModSeq, PrecursorCharge, PrecursorMz, RTInSeconds "
        "FROM entries ORDER BY PeptideModSeq, PrecursorCharge"))
    con.close()

    if not rows:
        sys.exit("DLIB_TO_REFERENCE_LIST: %s has an empty 'entries' table." % LIB)

    with open(OUT, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Compound', 'Formula', 'Adduct', 'm/z', 'z',
                    'RT Time (min)', 'Window (min)'])
        seen = set()
        for seq, z, mz, rt in rows:
            key = (seq, int(z))
            if key in seen:
                continue
            seen.add(key)
            rt_min = ''
            if KEEP_RT and rt is not None:
                rt_min = '%.4f' % (float(rt) / 60.0)
            w.writerow([seq, '', '(no adduct)', '%.6f' % float(mz), int(z),
                        rt_min, 5.0])

    print("extracted %d peptidoform-charge targets from %s (%s)"
          % (len(seen), LIB, KIND), file=sys.stderr)
    if KEEP_RT:
        print("  RT written as minutes from RTInSeconds", file=sys.stderr)
    else:
        print("  RT column left empty: a .dlib's RTInSeconds may be iRT rather "
              "than gradient time. diathem then uses every scan in the "
              "isolation window, which is correct but less selective.",
              file=sys.stderr)
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      source_kind: '${source_kind}'
    EOF
    """

    stub:
    """
    touch ${library.baseName}.reference_list.csv ${library.baseName}.reference_list.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
