process FINALIZE_QUANT {
    label 'python'
    tag   "${route}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/quant/${route}" },
               mode: params.publish_mode,
               saveAs: { fn -> fn == 'versions.yml' ? 'finalize_quant.versions.yml' :
                               fn.endsWith('.warnings.txt') ? null : fn }

    input:
    path  base_long
    path  reference_list, stageAs: 'ref/*'
    path  fasta,          stageAs: 'fasta/*'
    val   heavy_label
    path  tolerance,      stageAs: 'tolerance/*'
    path  min_coelution,  stageAs: 'coelution/*'
    path  sample_map,     stageAs: 'runs/*'
    val   route

    output:
    path "quant.full.tsv",                emit: full
    path "quant.reference.confident.tsv", emit: confident
    path "quant.reference.tsv",           emit: reference
    path "quant_summary.tsv",             emit: summary
    path "finalize_quant.log",            emit: log
    path "finalize_quant.warnings.txt",   emit: warnings
    path "versions.yml",                  emit: versions

    script:
    def _rows_uni = file("${projectDir}/assets/unimod_delta_masses.tsv").readLines()
    def _delta_by_acc = [:]
    _rows_uni.drop(1).each { l ->
        def c = l.split('\t')
        if (c.size() >= 3) _delta_by_acc[c[0]] = c[2] as Double
    }
    def _aug = (heavy_label ?: []).collect { r ->
        def d = r.containsKey('delta') ? (r.delta as Double) : _delta_by_acc[r.unimod]
        [ site: r.site, unimod: r.unimod, delta: d ]
    }
    def label_json = groovy.json.JsonOutput.toJson(_aug)
    def mod_tol    = params.mod_mass_tolerance
    """
    set -o pipefail

    python3 - <<'PY' 2>&1 | tee finalize_quant.log
    import csv, os, sys, glob, json, re

    BASE    = "${base_long}"
    ROUTE   = "${route}"
    FDR     = float("${params.fdr}")
    PRIMARY = "${params.primary_abundance}"
    LABEL   = json.loads(r'''${label_json}''')
    MTOL    = float("${mod_tol}")
    def per_run(pattern, column, default):
        files = [p for p in glob.glob(pattern) if os.path.basename(p) != 'NO_FILE']
        table = {r['sample_id']: float(r[column]) for r in csv.DictReader(open(files[0]), delimiter='\\t')} if files else {}
        other = table.pop('*', default)
        return lambda sid: table.get(sid, other)
    TOL      = per_run('tolerance/*', 'tolerance_seconds', 0.0)
    MIN_COEL = per_run('coelution/*', 'min_coelution',
                       float("${params.pair_min_coelution}") if "${params.pair_min_coelution}" != 'auto' else 0.7)
    RUNS     = [l.split('\\t')[1].strip() for p in glob.glob('runs/*') if os.path.basename(p) != 'NO_FILE'
                for l in open(p) if l.strip()]
    MIN_IDS = int("${params.search_min_reference_ids}")
    warnings = []
    refs    = [p for p in glob.glob('ref/*') if os.path.basename(p) != 'NO_FILE']
    REF     = refs[0] if refs else None
    FASTA   = "${fasta}"
    TOOLS   = ['encyclopedia', 'openswath', 'diathem', 'tric', 'standard_peak']
    SHARED  = 'abundance_%s_shared_ions' % PRIMARY

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

    _flank = re.compile(r'^[A-Z\\-]\\.|\\.[A-Z\\-]\$')
    def unflank(s):
        return _flank.sub('', str(s or '').strip())

    def strip_mods(s):
        s = re.sub(r'\\.\\([^)]+\\)|\\([^)]+\\)|\\[[^\\]]*\\]', '', unflank(s))
        return s.strip('.')

    _accs   = [r["unimod"] for r in LABEL if r.get("unimod")]
    _deltas = [float(r["delta"]) for r in LABEL if r.get("delta") is not None]
    _acc_re   = re.compile("|".join(re.escape(a) for a in _accs)) if _accs else None
    _delta_re = re.compile(r'[\\[\\(]([+-]?\\d+(?:\\.\\d+)?)[\\]\\)]')

    def channel_of(pf):
        s = str(pf)
        if _acc_re and _acc_re.search(s): return 'heavy'
        for m in _delta_re.finditer(s):
            try: v = float(m.group(1))
            except ValueError: continue
            if any(abs(v - d) <= MTOL for d in _deltas): return 'heavy'
        return 'light'

    _mod_re = re.compile(r'\\[([^\\]]*)\\]|\\(([^)]*)\\)')
    def form_of(pf):
        def one(m):
            text = (m.group(1) if m.group(1) is not None else m.group(2)).strip()
            try:
                v = float(text)
            except ValueError:
                return '' if _acc_re and _acc_re.fullmatch(text) else m.group(0)
            return '' if any(abs(v - d) <= MTOL for d in _deltas) else '[%+.1f]' % v
        return _mod_re.sub(one, unflank(pf)).strip('.')

    def num(v):
        v = (v or '').strip()
        if v in ('', 'NA', 'nan', 'None'): return None
        try: return float(v)
        except ValueError: return None

    def in_fasta(path, peptides):
        peptides = {p for p in peptides if p}
        short = {p for p in peptides if len(p) < 5}
        by_start = {}
        for p in peptides - short:
            by_start.setdefault(p[:5], []).append(p)
        found = {}

        def scan(name, s):
            for i in range(len(s) - 4):
                for p in by_start.get(s[i:i + 5], ()):
                    if s.startswith(p, i):
                        found.setdefault(p, set()).add(name)
            for p in short:
                if p in s:
                    found.setdefault(p, set()).add(name)
        name, seq = None, []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if name is not None:
                        scan(name, ''.join(seq).upper())
                    name, seq = line[1:].split(' ', 1)[0], []
                else:
                    seq.append(line)
        if name is not None:
            scan(name, ''.join(seq).upper())
        return {p: ';'.join(sorted(n)) for p, n in found.items()}

    def col(row, names):
        low = {k.strip().lower(): k for k in row}
        return next((low[n] for n in names if n in low), None)

    def charge(v):
        z = num(v)
        return str(int(z)) if z is not None else ''

    targets, decoys = {}, set()
    if REF:
        ref_rows = read(REF)
        pcol = col(ref_rows[0], ('compound', 'sequence', 'peptide', 'peptidesequence',
                                 'modifiedpeptidesequence')) if ref_rows else None
        zcol = col(ref_rows[0], ('z', 'charge', 'precursorcharge')) if ref_rows else None
        dcol = col(ref_rows[0], ('isdecoy', 'decoy', 'is_decoy')) if ref_rows else None
        if pcol is None:
            sys.exit("FINALIZE_QUANT [%s]: no peptide column in %s" % (ROUTE, REF))
        for r in ref_rows:
            cmp_ = (r.get(pcol) or '').strip()
            if not cmp_: continue
            key = (form_of(cmp_), channel_of(cmp_))
            if dcol and (r.get(dcol) or '').strip().lower() in ('true', '1', 'yes', 'decoy'):
                decoys.add(key)
            else:
                targets.setdefault(key, {}).setdefault(charge(r.get(zcol)) if zcol else '', unflank(cmp_))
        print("reference list: %s -> %d reference peptides (peptide, channel), %d decoys"
              % (REF, len(targets), len(decoys)), file=sys.stderr)
    from_list = bool(targets)

    samples = sorted({r['sample_id'] for r in rows if r.get('sample_id')} | set(RUNS))
    no_results = sorted(set(RUNS) - {r['sample_id'] for r in rows if r.get('sample_id')})
    abund = {t: 'abundance_' + t for t in TOOLS}
    ran = [t for t in TOOLS if any(num(r.get(abund[t])) not in (None, 0.0) for r in rows)]
    print("tools with abundances on this route: %s; abundance = abundance_%s"
          % (', '.join(ran) or 'none', PRIMARY), file=sys.stderr)

    def best_peak_group(r):
        v = num(r.get('peak_group_rank'))
        return v is None or int(v) <= 1

    out = []
    for r in rows:
        pf   = r.get('peptidoform') or r.get('peptide') or ''
        pep  = strip_mods(pf)
        form = form_of(pf)
        ch   = (r.get('channel') or '').strip() or channel_of(pf)
        q    = num(r.get('id_qvalue'))
        key  = (form, ch)
        dec  = (r.get('is_decoy') or '').strip().lower() in ('true', '1') or \\
               (key in decoys and key not in targets)
        o = {
            'sample_id': r['sample_id'], 'protein': r.get('protein') or '',
            'peptide': pep, 'peptidoform': unflank(pf), 'charge': charge(r.get('charge')),
            'channel': ch,
            'detected': 'true' if (q is not None and q <= FDR and not dec) else 'false',
            'q_value': r.get('id_qvalue') or '', 'pep': r.get('id_pep') or '',
            'abundance': r.get(abund.get(PRIMARY, '')) or '',
        }
        for t in ran:
            o[abund[t]] = r.get(abund[t]) or ''
        for c in ('consistency', 'n_effective_transitions', 'rt_apex_seconds',
                  'peak_group_rank') + tuple('rt_%s_seconds' % t for t in TOOLS):
            o[c] = r.get(c) or ''
        o['is_decoy'] = 'true' if dec else 'false'
        o['_best'] = best_peak_group(r)
        o['_shared'] = r.get(SHARED) or ''
        o['_coel'] = r.get('coelution_standard_peak') or ''
        o['coelution_standard_peak'] = o['_coel']
        o['_form'] = form
        out.append(o)

    del rows
    RTQ = 'rt_%s_seconds' % PRIMARY
    if not any(o.get(RTQ) for o in out):
        RTQ = 'rt_apex_seconds'

    if not targets:
        for o in out:
            if o['detected'] == 'true':
                targets.setdefault((o['_form'], o['channel']), {}).setdefault(o['charge'], o['peptidoform'])
        print("no reference list: the reference set is the %d peptide(s) detected at "
              "q <= %g in at least one run" % (len(targets), FDR), file=sys.stderr)
    for o in out:
        o['in_reference_list'] = 'true' if (o['_form'], o['channel']) in targets else 'false'

    seen = {(o['sample_id'], o['_form'], o['charge'], o['channel']) for o in out}
    added = 0
    for sid in samples:
        for (form, ch), charges in sorted(targets.items()):
            for z, pf in sorted(charges.items()):
                if (sid, form, z, ch) in seen: continue
                out.append({'sample_id': sid, 'peptide': strip_mods(form), 'peptidoform': pf, 'charge': z,
                            'channel': ch, 'detected': 'false', 'is_decoy': 'false',
                            'in_reference_list': 'true', '_best': True, '_form': form})
                added += 1
    print("added %d reference peptide x sample row(s) that no tool reported" % added,
          file=sys.stderr)

    protein_of = {}
    for o in out:
        if o.get('protein'): protein_of.setdefault(o['peptide'], o['protein'])
    for o in out:
        if not o.get('protein'): o['protein'] = protein_of.get(o['peptide'], '')
    blank = {o['peptide'] for o in out
             if not o.get('protein') and o['peptide'] and o['is_decoy'] != 'true'}
    if blank:
        found = in_fasta(FASTA, blank)
        for o in out:
            if not o.get('protein'): o['protein'] = found.get(o['peptide'], '')
        print("protein: from the FASTA for %d of %d peptide(s) no tool gave one"
              % (len(found), len(blank)), file=sys.stderr)

    has = lambda c: any((o.get(c) or '') != '' for o in out)
    cols = ['sample_id', 'protein', 'peptide', 'peptidoform', 'charge', 'channel',
            'detected', 'q_value', 'pep', 'abundance'] + [abund[t] for t in ran]
    if 'diathem' in ran:
        cols += ['consistency', 'n_effective_transitions']
    cols += [c for c in ['rt_apex_seconds'] + ['rt_%s_seconds' % t for t in ran] if has(c)]
    if 'standard_peak' in ran and has('coelution_standard_peak'):
        cols += ['coelution_standard_peak']
    full_cols = cols + ['in_reference_list', 'is_decoy'] + \\
                (['peak_group_rank'] if has('peak_group_rank') else [])

    def write(path, data, fields):
        with open(path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter='\\t', extrasaction='ignore',
                               lineterminator='\\n')
            w.writeheader()
            for r in sorted(data, key=lambda x: (x['sample_id'], x['peptide'],
                                                 str(x.get('charge', '')), x['channel'])):
                w.writerow(r)

    paired = {p for p, c in targets if c == 'light'} & {p for p, c in targets if c == 'heavy'}
    heavy_only = ({p for p, c in targets if c == 'heavy'} - paired) if (from_list and paired) else set()
    if heavy_only:
        print("left out of quant.reference.tsv and quant.reference.confident.tsv: %d reference peptide(s) "
              "with no light form in the reference list" % len(heavy_only), file=sys.stderr)

    write('quant.full.tsv', out, full_cols)
    shown = [o for o in out if o['is_decoy'] != 'true' and o['_best']
             and o['_form'] not in heavy_only]
    reference = [o for o in shown if o['in_reference_list'] == 'true']
    confident = [o for o in reference if o['detected'] == 'true']
    write('quant.reference.tsv', reference, cols)
    write('quant.reference.confident.tsv', confident, cols)

    by_sample = {}
    for o in reference:
        if o['detected'] == 'true':
            by_sample.setdefault(o['sample_id'], []).append(o)

    def detected_in(sid, peptides):
        det = {}
        for o in by_sample.get(sid, ()):
            if o['_form'] in peptides:
                d = det.setdefault((o['_form'], o['charge'], o['channel']), [False, None, None, None])
                if not d[0] and num(o.get('abundance')) not in (None, 0.0):
                    d[0], d[1], d[2], d[3] = True, num(o.get(RTQ)), num(o.get('_shared')), num(o.get('_coel'))
        return det

    def co_elute(sid, det, p, z):
        if PRIMARY == 'standard_peak':
            c, m = det[(p, z, 'light')][3], MIN_COEL(sid)
            return m <= 0 or c is None or c >= m
        a, b, t = det[(p, z, 'light')][1], det[(p, z, 'heavy')][1], TOL(sid)
        return t <= 0 or a is None or b is None or abs(a - b) <= t

    n_ids = []
    with open('quant_summary.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\\t', lineterminator='\\n')
        if paired:
            w.writerow(['route', 'sample_id', 'reference_pairs',
                        'heavy_detected', 'heavy_quantified', 'light_detected', 'light_quantified',
                        'pairs_detected', 'pairs_quantified'])
            for sid in samples:
                det = detected_in(sid, paired)
                pairs = [(p, z) for p, z, c in det if c == 'light' and (p, z, 'heavy') in det]
                counts = []
                for ch in ('heavy', 'light'):
                    counts += [len({p for p, _z, c in det if c == ch}),
                               len({p for (p, _z, c), q in det.items() if c == ch and q[0]})]
                quant = [(p, z) for p, z in pairs
                         if det[(p, z, 'light')][0] and det[(p, z, 'heavy')][0]
                         and 0.0 not in (det[(p, z, 'light')][2], det[(p, z, 'heavy')][2])]
                n_ids.append((sid, counts[0] + counts[2]))
                w.writerow([ROUTE, sid, len(paired)] + counts +
                           [len({p for p, _z in pairs}),
                            len({p for p, z in quant if co_elute(sid, det, p, z)})])
        else:
            w.writerow(['route', 'sample_id', 'reference_peptides', 'detected', 'quantified'])
            for sid in samples:
                det = detected_in(sid, {p for p, _c in targets})
                hit = {(p, c) for p, _z, c in det}
                quant = {(p, c) for (p, _z, c), q in det.items() if q[0]}
                n_ids.append((sid, len(hit)))
                w.writerow([ROUTE, sid, len(targets), len(hit), len(quant)])

    if paired:
        print("detected: identified at q <= %g; a pair is quantified when its light also %s" % (
              FDR, "follows its standard over the standard's peak (co-elution at least the run's pair_min_coelution)"
              if PRIMARY == 'standard_peak' else "peaks within the run's pair_rt_tolerance_seconds of its standard"),
              file=sys.stderr)
    print("wrote quant.full.tsv (%d rows), quant.reference.tsv (%d), quant.reference.confident.tsv (%d) "
          "over %d sample(s)" % (len(out), len(reference), len(confident), len(samples)),
          file=sys.stderr)
    for line in open('quant_summary.tsv'):
        print("  " + line.rstrip(), file=sys.stderr)
    if not confident:
        warnings.append("no reference peptide was detected at q <= %g in any run; the tables are "
                        "still written. Check the search settings for your instrument." % FDR)
    else:
        warnings += ["run %s: %d reference peptides detected at q <= %g, fewer than "
                     "search_min_reference_ids (%d)" % (sid, n, FDR, MIN_IDS)
                     for sid, n in n_ids if n < MIN_IDS and sid not in no_results]
    warnings += ["run %s: no results; its search reported no peptide, or a step failed for it (the "
                 "pipeline's log and pipeline_info/trace.txt say which)" % sid for sid in no_results]
    for line in warnings:
        print("WARNING: " + line, file=sys.stderr)
    with open('finalize_quant.warnings.txt', 'w') as f:
        f.writelines(line + '\\n' for line in warnings)
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
    for f in quant.full.tsv quant.reference.confident.tsv quant.reference.tsv quant_summary.tsv; do touch \$f; done
    touch finalize_quant.log finalize_quant.warnings.txt
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
