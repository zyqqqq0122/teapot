process EXTRACT_OPENSWATH_INTENSITY {
    label 'python'
    tag   "${meta.id}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/extract_openswath_intensity/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(reference_list), path(osw)

    output:
    tuple val(meta), path("${meta.id}.osw_intensity.tsv"), emit: intensity
    tuple val(meta), path("${meta.id}.extract.log"),       emit: log
    path  "versions.yml", emit: versions

    script:
    """
    python3 - <<'PY' > ${meta.id}.extract.log 2>&1
    import sqlite3, sys, csv

    OSW = "${osw}"
    OUT = "${meta.id}.osw_intensity.tsv"

    con = sqlite3.connect(OSW)
    cur = con.cursor()

    def has_column(table, col):
        for r in con.execute(f"PRAGMA table_info({table})"):
            if r[1].lower() == col.lower():
                return True
        return False

    order = ""
    if has_column('SCORE_MS2', 'SCORE'):
        joins = "LEFT JOIN SCORE_MS2 SC ON SC.FEATURE_ID = F.ID"
        order = "ORDER BY SC.SCORE DESC NULLS LAST, FMS2.AREA_INTENSITY DESC"
    else:
        joins = ""
        order = "ORDER BY FMS2.AREA_INTENSITY DESC"

    sql = f'''
      SELECT P.UNMODIFIED_SEQUENCE AS stripped_seq,
             P.MODIFIED_SEQUENCE   AS peptidoform,
             PRE.CHARGE            AS charge,
             P.DECOY               AS decoy,
             F.EXP_RT              AS rt_apex_seconds,
             FMS2.AREA_INTENSITY   AS abundance_openswath
      FROM FEATURE F
      JOIN PRECURSOR_PEPTIDE_MAPPING PPM ON PPM.PRECURSOR_ID = F.PRECURSOR_ID
      JOIN PEPTIDE P                     ON P.ID              = PPM.PEPTIDE_ID
      JOIN PRECURSOR PRE                 ON PRE.ID            = F.PRECURSOR_ID
      JOIN FEATURE_MS2 FMS2              ON FMS2.FEATURE_ID   = F.ID
      {joins}
      {order}
    '''
    seen = set()
    rows = []
    for r in con.execute(sql):
        key = (r[1], int(r[2]), int(r[3] or 0))
        if key in seen:
            continue
        seen.add(key)
        rows.append(r)
    con.close()

    with open(OUT, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\\t')
        w.writerow(['stripped_seq', 'peptidoform', 'charge', 'decoy',
                    'rt_apex_seconds', 'abundance_openswath'])
        w.writerows(rows)

    print(f"rows written: {len(rows)}", file=sys.stderr)
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    touch ${meta.id}.osw_intensity.tsv ${meta.id}.extract.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
