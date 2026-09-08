process STRIP_REFERENCE_DECOYS {
    label 'python'
    tag   "${reference_list.baseName}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/prepare_library" },
               mode: params.publish_mode

    input:
    path reference_list

    output:
    path "${reference_list.baseName}.targets_only.${reference_list.extension}", emit: reference_list
    path "strip_reference_decoys.log",                                          emit: log
    path "versions.yml",                                                        emit: versions

    script:
    def out_name = "${reference_list.baseName}.targets_only.${reference_list.extension}"
    """
    set -o pipefail

    python3 - <<'PY' 2>&1 | tee strip_reference_decoys.log
    import csv, sys

    SRC = "${reference_list}"
    OUT = "${out_name}"
    NAMES = ('isdecoy', 'decoy', 'is_decoy')
    TRUE  = ('true', '1', 'yes', 'decoy')

    rows = [l for l in open(SRC).read().splitlines() if l.strip()]
    if not rows:
        sys.exit("strip_reference_decoys: %s is empty." % SRC)

    delim = '\\t' if rows[0].count('\\t') > rows[0].count(',') else ','
    header = next(csv.reader([rows[0]], delimiter=delim))
    idx = next((i for i, c in enumerate(header)
                if c.strip().lower().strip('"') in NAMES), -1)

    with open(OUT, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter=delim)
        if idx < 0:
            for r in rows:
                fh.write(r + '\\n')
            print("no decoy column in %s; copied %d data row(s) unchanged"
                  % (SRC, len(rows) - 1))
            sys.exit(0)
        w.writerow(header)
        kept = dropped = 0
        for line in rows[1:]:
            cells = next(csv.reader([line], delimiter=delim))
            v = cells[idx].strip().lower().strip('"') if len(cells) > idx else ''
            if v in TRUE:
                dropped += 1
            else:
                w.writerow(cells); kept += 1

    if kept == 0:
        sys.exit("strip_reference_decoys HARD ERROR: every row in %s is a decoy "
                 "(%d dropped, 0 kept). diathem would have nothing to quantify."
                 % (SRC, dropped))
    print("kept %d target row(s), dropped %d decoy row(s)" % (kept, dropped))
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    def stub_out = "${reference_list.baseName}.targets_only.${reference_list.extension}"
    """
    touch ${stub_out} strip_reference_decoys.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
