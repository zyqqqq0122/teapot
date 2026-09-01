process BLIB_TO_REFERENCE_LIST {
    label 'python'
    tag   "${blib.baseName}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/prepare_library_openswath" },
               mode: params.publish_mode

    input:
    path blib

    output:
    path "${blib.baseName}.reference_list.csv",  emit: reference_list
    path "${blib.baseName}.reference_list.log",  emit: log
    path "versions.yml",                    emit: versions

    script:
    """
    python3 - <<'PY' > ${blib.baseName}.reference_list.log 2>&1
    import sqlite3, csv, sys

    BLIB = "${blib}"
    OUT  = "${blib.baseName}.reference_list.csv"

    con = sqlite3.connect(BLIB)
    cur = con.cursor()

    rows = list(con.execute(
        "SELECT peptideModSeq, precursorCharge, precursorMZ, retentionTime "
        "FROM RefSpectra"))
    con.close()

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
            w.writerow([seq, '', '(no adduct)', f'{float(mz):.6f}', int(z),
                        f'{float(rt):.4f}' if rt is not None else '', 5.0])

    print(f"extracted {len(seen)} peptidoform-charge entries from {BLIB}",
          file=sys.stderr)
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    touch ${blib.baseName}.reference_list.csv ${blib.baseName}.reference_list.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
