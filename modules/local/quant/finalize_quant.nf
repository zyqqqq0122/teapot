process FINALIZE_QUANT {
    label 'python'
    tag   "${route}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/quant/${route}" },
               mode: params.publish_mode

    input:
    path  base_long
    path  reference_list, stageAs: 'ref/*'
    val   heavy_label
    val   route

    output:
    path "quant.full.tsv",       emit: full
    path "quant.confident.tsv",  emit: confident
    path "quant.reference.tsv",  emit: reference
    path "quant_summary.tsv",    emit: summary
    path "finalize_quant.log",   emit: log
    path "versions.yml",         emit: versions

    script:
    def _rows_uni = file("${projectDir}/assets/unimod_delta_masses.tsv").readLines()
    def _delta_by_acc = [:]
    _rows_uni.drop(1).each { l ->
        def c = l.split('\t')
        if (c.size() >= 3) _delta_by_acc[c[0]] = c[2] as Double
    }
    def _aug = (heavy_label ?: []).collect { r ->
        def d = r.containsKey('delta') ? (r.delta as Double) : _delta_by_acc[r.unimod]
        [ site: (r.site ?: r.residue), unimod: r.unimod, delta: d ]
    }
    def label_json = groovy.json.JsonOutput.toJson(_aug)
    def mod_tol    = params.mod_mass_tolerance
    """
    python3 - <<'PY' > finalize_quant.log 2>&1
    import csv, os, sys, glob, itertools, json, re

    BASE  = "${base_long}"
    ROUTE = "${route}"
    FDR   = float("${params.fdr}")
    LABEL = json.loads(r'''${label_json}''')
    MTOL  = float("${mod_tol}")
    refs  = [p for p in glob.glob('ref/*') if os.path.basename(p) != 'NO_FILE']
    REF   = refs[0] if refs else None

    def sniff(path):
        with open(path) as f:
            head = f.readline()
        return '\\t' if head.count('\\t') >= head.count(',') else ','

    def read(path):
        with open(path) as f:
            return list(csv.DictReader(f, delimiter=sniff(path)))

    rows = read(BASE)
    if not rows:
        sys.exit("FINALIZE_QUANT [%s]: %s is empty." % (ROUTE, BASE))

    def strip_mods(s):
        s = re.sub(r'\\.\\([^)]+\\)|\\([^)]+\\)|\\[[^\\]]*\\]', '', str(s))
        return re.sub(r'^[A-Z\\-]\\.|\\.[A-Z\\-]\$', '', s).strip('.')

    _accs   = [r["unimod"] for r in LABEL if r.get("unimod")]
    _deltas = [float(r["delta"]) for r in LABEL if r.get("delta") is not None]
    _acc_re   = re.compile("|".join(re.escape(a) for a in _accs)) if _accs else None
    _delta_re = re.compile(r'[\\[\\(]([+-]?\\d+(?:\\.\\d+)?)[\\]\\)]')
    HAS_HEAVY = bool(LABEL)

    def is_heavy(pf):
        if not HAS_HEAVY: return False
        s = str(pf)
        if _acc_re and _acc_re.search(s): return True
        for m in _delta_re.finditer(s):
            try: v = float(m.group(1))
            except ValueError: continue
            if any(abs(v - d) <= MTOL for d in _deltas): return True
        return False

    channel_of = lambda pf: 'heavy' if is_heavy(pf) else 'light'

    targets = {}
    if REF:
        for r in read(REF):
            cmp_ = (r.get('Compound') or r.get('compound') or '').strip()
            if not cmp_: continue
            is_dec = str(r.get('isDecoy', '')).strip().lower() in ('true','1','yes')
            targets.setdefault((strip_mods(cmp_), channel_of(cmp_), is_dec), cmp_)
        n_dec = sum(1 for k in targets if k[2])
        print("reference list: %s -> %d (sequence, channel) targets + %d decoys"
              % (REF, len(targets) - n_dec, n_dec), file=sys.stderr)
    else:
        print("no reference list supplied; reporting only what the search produced",
              file=sys.stderr)

    samples = sorted({r['sample_id'] for r in rows if r.get('sample_id')})

    def num(v):
        v = (v or '').strip()
        if v in ('', 'NA', 'nan', 'None'): return None
        try: return float(v)
        except ValueError: return None

    ABUND = [c for c in rows[0] if c.startswith('abundance_')]

    def annotate(r):
        q      = num(r.get('id_qvalue'))
        primary= num(r.get('abundance_primary'))
        any_ab = any(num(r.get(c)) not in (None, 0.0) for c in ABUND)
        ident  = q is not None and q <= FDR
        quant  = any_ab
        r['passes_fdr']  = 'true' if ident else 'false'
        r['status']      = ('identified_quantified' if ident and quant else
                            'identified_only'       if ident else
                            'quantified_only'       if quant else
                            'not_detected')
        pf  = r.get('peptidoform') or r.get('peptide') or ''
        key = strip_mods(pf)
        r['stripped_seq'] = key
        ch = (r.get('channel') or '').strip() or channel_of(pf)
        r['channel'] = ch
        known_t = (key, ch, False) in targets
        known_d = (key, ch, True)  in targets
        r['in_reference_list'] = 'true' if (not targets or known_t or known_d) else 'false'
        if not (r.get('is_decoy') or '').strip():
            r['is_decoy'] = 'true' if (known_d and not known_t) else 'false'
        return r

    out = [annotate(dict(r)) for r in rows]

    if targets:
        seen = {(r['sample_id'], r['stripped_seq'], r['channel']) for r in out}
        blank = {c: '' for c in out[0]}
        added = 0
        for sid, key in itertools.product(samples, sorted(targets)):
            seq, ch, is_dec = key
            if (sid, seq, ch) in seen: continue
            r = dict(blank)
            r.update(sample_id=sid, peptide=targets[key], peptidoform=targets[key],
                     stripped_seq=seq, channel=ch, protein='', passes_fdr='false',
                     status='no_feature', in_reference_list='true',
                     is_decoy='true' if is_dec else 'false')
            out.append(r); added += 1
        print("grid completion: added %d (target, sample) rows with no feature" % added,
              file=sys.stderr)

    cols = [c for c in rows[0]]
    for c in ('stripped_seq','in_reference_list','is_decoy','passes_fdr','status'):
        if c not in cols: cols.append(c)

    def write(path, data):
        with open(path,'w',newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter='\\t', extrasaction='ignore')
            w.writeheader()
            for r in sorted(data, key=lambda x: (x.get('sample_id',''),
                                                 x.get('stripped_seq',''),
                                                 x.get('channel',''))):
                w.writerow(r)

    write('quant.full.tsv', out)

    confident = [r for r in out if r['passes_fdr'] == 'true' and r['is_decoy'] != 'true']
    write('quant.confident.tsv', confident)
    reference = [r for r in out if r['in_reference_list'] == 'true' and r['is_decoy'] != 'true']
    write('quant.reference.tsv', reference)

    with open('quant_summary.tsv','w',newline='') as f:
        w = csv.writer(f, delimiter='\\t')
        w.writerow(['route','sample_id','rows','targeted','targeted_peptides',
                    'decoys','identified','identified_peptides',
                    'quantified','quantified_peptides',
                    'identified_and_quantified','no_signal'])
        def npep(rows):
            return len({r.get('stripped_seq') or r.get('peptide') for r in rows})
        for sid in samples:
            g    = [r for r in out if r['sample_id'] == sid]
            real = [r for r in g if r['is_decoy'] != 'true']
            tgt  = [r for r in real if r['in_reference_list'] == 'true']
            idn  = [r for r in real if r['passes_fdr'] == 'true']
            qnt  = [r for r in real if r['status'] in ('quantified_only','identified_quantified')]
            w.writerow([ROUTE, sid, len(g),
                        len(tgt), npep(tgt),
                        sum(1 for r in g if r['is_decoy'] == 'true'),
                        len(idn), npep(idn),
                        len(qnt), npep(qnt),
                        sum(1 for r in real if r['status'] == 'identified_quantified'),
                        sum(1 for r in real if r['status'] in ('not_detected','no_feature'))])

    print("wrote quant.full.tsv (%d) / quant.confident.tsv (%d) / quant.reference.tsv (%d) "
          "over %d sample(s)" % (len(out), len(confident), len(reference), len(samples)),
          file=sys.stderr)
    for line in open('quant_summary.tsv'):
        print("  " + line.rstrip(), file=sys.stderr)

    if not confident:
        print("WARNING: FINALIZE_QUANT [%s]: no row passes q<=%g. The full table is "
              "still written; check ASSERT_SEARCH_PRODUCTIVE and ASSERT_QUANT_PRODUCTIVE."
              % (ROUTE, FDR), file=sys.stderr)
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      route:  ${route}
      fdr:    ${params.fdr}
    EOF
    """

    stub:
    """
    for f in quant.full.tsv quant.confident.tsv quant.reference.tsv quant_summary.tsv; do touch \$f; done
    touch finalize_quant.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
