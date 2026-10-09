process EXTRACT_OPENSWATH_INTENSITY {
    label 'python'
    tag   "${meta.id}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/extract_openswath_intensity/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(reference_list), path(osw), path(peptide_table)

    output:
    tuple val(meta), path("${meta.id}.osw_intensity.tsv"), emit: intensity
    tuple val(meta), path("${meta.id}.extract.log"),       emit: log
    path  "versions.yml", emit: versions

    script:
    """
    python3 - <<'PY' > ${meta.id}.extract.log 2>&1
    import sqlite3, sys, csv, re

    OSW = "${osw}"
    IDS = "${peptide_table}"
    OUT = "${meta.id}.osw_intensity.tsv"

    chosen = set()
    with open(IDS) as f:
        for r in csv.DictReader(f, delimiter='\\t'):
            m = re.match(r'f(-?\\d+)_', r.get('SpecId') or '')
            if m: chosen.add(int(m.group(1)))

    con = sqlite3.connect(OSW)
    sql = '''
      SELECT F.ID                  AS feature_id,
             P.UNMODIFIED_SEQUENCE AS stripped_seq,
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
    '''
    best = {}
    for r in con.execute(sql):
        key = (r[2], int(r[3]), int(r[4] or 0))
        rank = (r[0] in chosen, r[6] or 0.0)
        if key not in best or rank > best[key][0]:
            best[key] = (rank, r[1:])
    con.close()

    with open(OUT, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\\t')
        w.writerow(['stripped_seq', 'peptidoform', 'charge', 'decoy',
                    'rt_apex_seconds', 'abundance_openswath'])
        w.writerows(row for _rank, row in best.values())

    n_id = sum(1 for rank, _row in best.values() if rank[0])
    print(f"precursors written: {len(best)}; {n_id} from the feature the FDR step "
          f"identified, the rest from their largest peak", file=sys.stderr)
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
