#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPARE_LIBRARY_ENCYCLOPEDIA } from './subworkflows/local/prepare_library_encyclopedia'
include { PREPARE_LIBRARY_OPENSWATH    } from './subworkflows/local/prepare_library_openswath'
include { ENCYCLOPEDIA_DIA             } from './subworkflows/local/encyclopedia_dia'
include { ENCYCLOPEDIA_PRM             } from './subworkflows/local/encyclopedia_prm'
include { OPENSWATH_DIA                } from './subworkflows/local/openswath_dia'
include { OPENSWATH_PRM                } from './subworkflows/local/openswath_prm'
include { MSCONVERT                    } from './modules/local/msconvert/msconvert'
include { REMOVE_SURVEY_SCANS          } from './modules/local/runs/remove_survey_scans'
include { CALIBRATE_FRAGMENT_MZ        } from './modules/local/runs/calibrate_fragment_mz'
include { ADD_DECOYS_TO_REFERENCE_LIST } from './modules/local/encyclopedia/add_decoys_to_reference_list'
include { NORMALIZE_REFERENCE_MODS     } from './modules/local/library/normalize_reference_mods'
include { CHECK_REFERENCE_COVERAGE     } from './modules/local/library/check_reference_coverage'


def helpMessage() {
    log.info """
    Usage:
      nextflow run main.nf -profile apptainer -params-file my_run.yml
      nextflow run main.nf -profile apptainer --input samplesheet.csv --fasta proteins.fasta [settings]

    Start from templates/: a samplesheet (.csv) and settings (.yml) for your kind of
    data, and a launcher, run_local.sh or run_slurm.sh. README.md says which to use.
    A setting goes in the .yml as name: value, or on the command line as --name value.

    The samplesheet (--input) is a CSV with one row per run and the header
    sample_id,file,mode,reference_list,engine:
      sample_id       your name for the run, a different one on every row
      file            the run: .mzML, .mzXML (engine openswath only), .dia (engine
                      encyclopedia only), or a vendor file the pipeline converts:
                      .raw, .wiff, .wiff2, .d, .lcd. Each with a name of its own
      mode            PRM or DIA
      reference_list  your reference peptides, in the isolation-list layout Skyline
                      exports, the same file for every run of an engine and mode.
                      Left empty, it is derived from the library when that is a
                      .blib. PRM needs one or the other; DIA without either reports
                      the peptides it detects. openswath does not use it on DIA.
      engine          encyclopedia (the default) or openswath

    Some settings (all of them in PARAMETERS.md):
      --primary_abundance   whose amounts fill the abundance column, and the SIL
                            amounts: encyclopedia (the default), openswath, tric,
                            diathem, or standard_peak (SIL runs, EncyclopeDIA:
                            both forms measured over the heavy standard's peak)
      --run_diathem         also quantify with diathem: true (the default) or false
      --diathem_library     a library for diathem to take each peptide's expected
                            fragment intensities from, instead of the assay library
      --heavy_label         SIL runs: the heavy labels, as in [{site: K, unimod:
                            'UniMod:259'}, {site: R, unimod: 'UniMod:267'}]. Empty,
                            the default, for label-free runs
      --standard_amounts    SIL runs: a CSV of the known amount (pmol) of each heavy
                            standard, by peptide or by protein
      --pair_min_coelution  SIL runs measured with standard_peak: a light/heavy pair
                            counts only when the light's signal follows its
                            standard's over the standard's peak, their correlation
                            at least this. auto (the default): the correlation 1 in
                            20 light signals reach when moved away from their
                            standard; a number from 0 to 1; 0 = off
      --pair_rt_tolerance_seconds
                            SIL runs, the other tools: a light/heavy pair whose two
                            peaks lie further apart than this many seconds does not
                            count. auto (the default): half the median peak width of
                            your heavy standards; 0 = off
      --fragment_mz_calibration
                            targeted runs: measure each run's fragment m/z error on
                            its own targets and correct it when larger than
                            fragment_mz_calibration_min_offset (auto: a quarter of
                            the search's fragment tolerance): auto (the default) or off
      --search_tolerance    targeted runs on the EncyclopeDIA routes: search at the
                            fragment tolerance the runs' own targets show after the
                            correction, all runs together, instead of -ftol: auto
                            (the default) or off
      --rt_alignment        SIL runs on the EncyclopeDIA routes: align the runs on
                            their heavy standards and search each peptide only
                            where it elutes in that run: auto (the default) or off
      --min_consistency     SIL runs: leave out of the amounts every measurement
                            whose diathem consistency is below this (default 0: none)

    INPUTS.md describes every input file, OUTPUTS.md every result.
    """.stripIndent()
}

def reference_list_has_decoys(rl) {
    if (!rl || !rl.exists()) return false
    def lines = rl.readLines().findAll { it.trim() }
    if (lines.size() < 2) return false
    def delim = (lines[0].count('\t') > lines[0].count(',')) ? '\t' : ','
    def header = lines[0].split(delim, -1).collect { it.trim().toLowerCase().replaceAll('"', '') }
    def idx = header.findIndexOf { it in ['isdecoy', 'decoy', 'is_decoy'] }
    if (idx < 0) return false
    return lines.drop(1).any { line ->
        def cells = line.split(delim, -1)
        if (cells.size() <= idx) return false
        def v = cells[idx].trim().toLowerCase().replaceAll('"', '')
        v in ['true', '1', 'yes', 'decoy']
    }
}

workflow {
    if (!params.input || !params.fasta) { helpMessage(); exit 1 }

    def missing_assets = ['NO_FILE', 'unimod_delta_masses.tsv', 'masslist-decoys.jar', 'wine_cpus_online']
        .findAll { !file("${projectDir}/assets/${it}").exists() }
    if (missing_assets) {
        error "Pre-flight: missing shipped asset(s) in ${projectDir}/assets: " +
              missing_assets.join(', ')
    }

    def hasPlaceholder
    hasPlaceholder = { x ->
        if (x instanceof CharSequence) return (x =~ /<[\w .\/-]+>/) as Boolean
        if (x instanceof Map)          return x.values().any { hasPlaceholder(it) }
        if (x instanceof Collection)   return x.any { hasPlaceholder(it) }
        return false
    }
    def unfilled = params.findAll { k, v -> hasPlaceholder(v) }
    if (unfilled) {
        error "Pre-flight: unfilled template placeholder(s):\n" +
              unfilled.collect { k, v -> "  ${k} = '${v}'" }.join('\n')
    }

    def rows = file(params.input).splitCsv(header: true)
    def bad_mode = rows.findAll { !(it.mode?.toUpperCase() in ['PRM', 'DIA']) }
    if (bad_mode) {
        error "Pre-flight: samplesheet mode must be PRM or DIA, found " +
              bad_mode.collect { "'${it.mode}' (${it.sample_id})" }.join(', ')
    }
    def bad_id = rows.findAll { !(it.sample_id ==~ /[A-Za-z0-9._+-]+/) }
    if (bad_id) {
        error "Pre-flight: a sample_id may hold letters, digits and . _ + - only, found " +
              bad_id.collect { "'${it.sample_id}'" }.join(', ')
    }
    def repeated = { names -> names.countBy { it }.findAll { _k, n -> n > 1 }.keySet() }
    def dup_id = repeated(rows.collect { it.sample_id })
    if (dup_id) {
        error "Pre-flight: every samplesheet row needs a sample_id of its own, repeated: ${dup_id.join(', ')}"
    }
    def dup_run = repeated(rows.collect { file(it.file).getBaseName() })
    if (dup_run) {
        error "Pre-flight: run files must differ in name, whatever their folder or extension, " +
              "repeated: ${dup_run.join(', ')}. Rename or link them under distinct names."
    }
    def ref_of = { r -> r.reference_list?.trim() ?: '' }
    def sig  = { r -> [ (r.engine?.trim() ?: 'encyclopedia').toLowerCase(),
                        r.mode.toUpperCase() ] }
    def routes = rows.collect(sig) as Set
    def has_enc     = routes.any { it[0] == 'encyclopedia' }
    def has_osw     = routes.any { it[0] == 'openswath' }
    def has_osw_dia = routes.any { it[0] == 'openswath' && it[1] == 'DIA' }
    def has_osw_prm = routes.any { it[0] == 'openswath' && it[1] == 'PRM' }
    def has_enc_dia = routes.any { it[0] == 'encyclopedia' && it[1] == 'DIA' }
    def has_enc_prm = routes.any { it[0] == 'encyclopedia' && it[1] == 'PRM' }
    def osw_labels = has_osw ? (params.heavy_label ?: []) : []

    def bg_source = params.background_source
    if (!(bg_source in ['none', 'measured', 'predicted', 'library'])) {
        error "background_source must be one of none, measured, predicted, library, got: '${bg_source}'"
    }
    def bg_files = ['background_library', 'background_fasta', 'background_peptides',
                    'background_sheet', 'background_search_library']
    def bg_blib  = params.background_library &&
                   params.background_library.toString().toLowerCase().endsWith('.blib')
    def bg_needs = [none:      [],
                    library:   ['background_library'] + (bg_blib ? ['background_fasta'] : []),
                    predicted: ['background_fasta'],
                    measured:  ['background_sheet', 'background_fasta']]
    def bg_takes = [none:      [],
                    library:   bg_needs.library,
                    predicted: ['background_fasta', 'background_peptides'],
                    measured:  ['background_sheet', 'background_fasta', 'background_search_library']]
    if (bg_source == 'measured' && params.background_library) {
        error "Pre-flight: background_library is a background you already have, for background_source: " +
              "library. To search the background_sheet runs with a library, give it as " +
              "background_search_library."
    }
    if (bg_source == 'library' && params.background_search_library) {
        error "Pre-flight: background_search_library is what the runs of a measured background are " +
              "searched with. A background you already have is background_library."
    }
    def bg_missing = bg_needs[bg_source].findAll { !params[it] }
    if (bg_missing) {
        error "Pre-flight: background_source '${bg_source}' needs ${bg_missing.join(' and ')}. " +
              ('background_fasta' in bg_missing && bg_blib
                  ? " A .blib is placed in proteins by a FASTA: entries no FASTA entry holds are dropped." : '')
    }
    def bg_unused = bg_files.findAll { params[it] && !(it in bg_takes[bg_source]) }
    if (bg_unused) {
        error "Pre-flight: ${bg_unused.join(' and ')} ${bg_unused.size() > 1 ? 'are' : 'is'} set, but background_source '${bg_source}' " +
              "does not use ${bg_unused.size() > 1 ? 'them' : 'it'}."
    }
    bg_files.findAll { params[it] && !file(params[it]).exists() }.each {
        error "Pre-flight: ${it}: ${params[it]} does not exist."
    }
    ['background_library', 'background_search_library'].findAll { params[it] }.each {
        if (!(params[it].toString().toLowerCase() ==~ /.*\.(elib|dlib|blib)$/)) {
            error "Pre-flight: ${it} must be an .elib, .dlib or .blib, got: ${params[it]}"
        }
    }
    if (Args.mode(params.fragment_mz_calibration) == null) {
        error "fragment_mz_calibration must be auto or off, got '${params.fragment_mz_calibration}'."
    }
    if (Args.mode(params.rt_alignment) == null) {
        error "rt_alignment must be auto or off, got '${params.rt_alignment}'."
    }
    if (Args.mode(params.search_tolerance) == null) {
        error "search_tolerance must be auto or off, got '${params.search_tolerance}'."
    }
    def half = params.rt_alignment_half_window_seconds as String
    if (!(half == 'auto' || (half.isNumber() && (half as double) > 0))) {
        error "rt_alignment_half_window_seconds must be auto or a number of seconds above 0, got '${half}'."
    }
    def cal_min = params.fragment_mz_calibration_min_offset as String
    if (!(cal_min == 'auto' || (cal_min.isNumber() && (cal_min as double) >= 0))) {
        error "fragment_mz_calibration_min_offset must be auto or a number of Th (0 or more), got '${cal_min}'."
    }
    def enc_quant_mode = params.encyclopedia_quant in [true, 'true', 'on'] ? 'on' :
                         params.encyclopedia_quant in [false, 'false', 'off'] ? 'off' :
                         params.encyclopedia_quant == 'auto' ? 'auto' : null
    if (enc_quant_mode == null) {
        error "encyclopedia_quant must be auto, on or off, got '${params.encyclopedia_quant}'."
    }
    def enc_quant = enc_quant_mode == 'on' || (enc_quant_mode == 'auto' && params.primary_abundance == 'encyclopedia')
    def tol = params.pair_rt_tolerance_seconds as String
    if (!(tol == 'auto' || (tol.isNumber() && (tol as double) >= 0))) {
        error "pair_rt_tolerance_seconds must be auto or a number of seconds (0 = off), got '${tol}'."
    }
    def coel = params.pair_min_coelution as String
    if (!(coel == 'auto' || (coel.isNumber() && (coel as double) >= 0 && (coel as double) <= 1))) {
        error "pair_min_coelution must be auto or a number from 0 to 1 (0 = off), got '${coel}'."
    }
    if (!((params.background_peptides_top as String).isInteger() && (params.background_peptides_top as int) > 0)) {
        error "background_peptides_top must be a positive whole number, got '${params.background_peptides_top}'."
    }

    ['library_sheet', 'background_sheet'].findAll { params[it] }.each { k ->
        if (!file(params[k]).exists()) {
            error "Pre-flight: ${k}: ${params[k]} does not exist."
        }
        def sheet_rows = file(params[k]).splitCsv(header: true)
        if (!sheet_rows || !['condition', 'fraction_file'].every { sheet_rows[0].containsKey(it) }) {
            error "Pre-flight: ${k} ${params[k]} needs the columns condition and fraction_file, " +
                  "and a row per run file (templates/library_sheet.csv)."
        }
        sheet_rows.each { r ->
            def f = r.fraction_file?.trim()
            if (!r.condition?.trim() || !f) {
                error "Pre-flight: ${k}: every row needs a condition and a fraction_file."
            }
            if (!file(f).exists()) {
                error "Pre-flight: ${k}: ${f} does not exist."
            }
            if (!(f.toLowerCase() ==~ /.*\.(mzml|dia|raw|wiff2?|d|lcd)$/)) {
                error "Pre-flight: ${k}: ${f}: the runs must be mzML, EncyclopeDIA .dia, or a vendor " +
                      "format the pipeline converts (.raw, .wiff, .wiff2, .d, .lcd)."
            }
        }
    }

    def enc_prm_uses_derived = has_enc_prm && params.blib && rows.any {
        (it.engine?.trim() ?: 'encyclopedia').toLowerCase() == 'encyclopedia' &&
        it.mode.toUpperCase() == 'PRM' && !ref_of(it)
    }
    def openswath_prm_uses_derived = has_osw_prm && (params.openswath_blib || params.blib) && rows.any {
        (it.engine?.trim() ?: 'encyclopedia').toLowerCase() == 'openswath' &&
        it.mode.toUpperCase() == 'PRM' && !ref_of(it)
    }
    def openswath_prm_needs_bg = openswath_prm_uses_derived &&
                           (params.openswath_prm_fdr ?: 'context') != 'pyprophet'
    if ((enc_prm_uses_derived || openswath_prm_needs_bg) && bg_source == 'none') {
        error """
            Pre-flight: PRM sample(s) will use a reference_list derived from the
            vendor .blib (samplesheet reference_list column empty). Without a
            background library the reference covers the entire target library
            and Context's mProphet LDA has zero background targets to train on.

            For the OpenSWATH PRM path you can instead set
              openswath_prm_fdr = 'pyprophet'
            which scores target vs decoy within the search and needs no
            background at all.

            Configure a background:
              background_source = 'measured'  + background_sheet (a few runs of
                                               samples like yours) + background_fasta
                                               (the proteins of your samples)
              background_source = 'predicted' + background_fasta, best with
                                               background_peptides
              background_source = 'library'   + background_library (.elib, .dlib, .blib)
        """.stripIndent().trim()
    }
    def has_vendor_blib = params.blib || params.openswath_blib
    if (has_vendor_blib && params.library_sheet) {
        error """
            Pre-flight: a .blib and --library_sheet were both supplied. The
            library_sheet runs would be searched with the .blib, and the chromatogram
            library they give would replace it: every standard the runs do not show
            would be lost.

            To build a background from runs, list them in background_sheet instead:
              --blib / --openswath_blib <vendor.blib>   (the assay library)
              background_source: measured
              background_sheet:  <runs.csv>             (laid out as library_sheet)
              background_fasta:  <the proteins of your samples>
        """.stripIndent().trim()
    }
    if (has_vendor_blib && params.elib) {
        error """
            Pre-flight: a .blib and --elib were both supplied as the assay library.
            The .elib would be searched and the .blib dropped, with its heavy channel
            and the reference list derived from it.

            To use the .elib as the background:
              background_source:  library
              background_library: <chromatogram.elib>
        """.stripIndent().trim()
    }
    def other_lib = ['dlib', 'use_koina'].findAll { params[it] }
    if (params.blib && other_lib) {
        error "Pre-flight: --blib and --${other_lib.join(' and --')} were both supplied as the assay " +
              "library, and only one would be searched. Supply one. To add peptides to a .blib as a " +
              "background, use background_source."
    }
    if (params.background_fasta && has_vendor_blib && file(params.background_fasta) == file(params.fasta)) {
        log.warn "background_fasta is --fasta. With a vendor .blib, --fasta usually holds only the " +
                 "standards, and a background of their peptides gives Context nothing to learn from: " +
                 "give the proteins of your samples."
    }
    if (bg_source != 'none') {
        def osw_native = params.openswath_pqp || params.openswath_tsv || params.openswath_traml
        def why = "On EncyclopeDIA a background is added only to a --blib library, and on OpenSWATH " +
                  "to any library but --openswath_pqp, --openswath_tsv and --openswath_traml. Other " +
                  "libraries are searched as they are: their entries outside your reference list are " +
                  "the background."
        def ignored = (has_enc && !params.blib ? ['EncyclopeDIA'] : []) +
                      (has_osw && osw_native   ? ['OpenSWATH']    : [])
        if (!((has_enc && params.blib) || (has_osw && !osw_native))) {
            error "Pre-flight: background_source is '${bg_source}', but no library here takes a background. ${why}"
        }
        if (ignored) {
            log.warn "background_source '${bg_source}': the ${ignored.join(' and ')} routes do not use it. ${why}"
        }
        if (!(has_enc_prm || (has_osw_prm && (params.openswath_prm_fdr ?: 'context') != 'pyprophet'))) {
            log.warn "background_source '${bg_source}': no route here learns from a background (EncyclopeDIA " +
                     "PRM, or OpenSWATH PRM with openswath_prm_fdr: context). It only adds to the search " +
                     "library, where it competes with your peptides."
        }
    }
    if (bg_source == 'none' && (params.heavy_label ?: []) && has_enc_prm && params.blib && !enc_prm_uses_derived) {
        log.warn "EncyclopeDIA PRM: the library is --blib, with no background. Context learns from the " +
                 "library's peptides outside your reference list; if the list covers the whole .blib, as a " +
                 "vendor kit's does, there are none. Set background_source."
    }
    if (bg_source == 'none' && (params.heavy_label ?: []) && has_osw_prm && (params.openswath_blib || params.blib) &&
        (params.openswath_prm_fdr ?: 'context') != 'pyprophet' && !openswath_prm_uses_derived) {
        log.warn "OpenSWATH PRM: the library is the .blib, with no background. Context-MS learns from the " +
                 "library's peptides outside your reference list; if the list covers the whole .blib, as a " +
                 "vendor kit's does, there are none. Set background_source, or " +
                 "openswath_prm_fdr: pyprophet."
    }

    def diathem_weights = ['uniform', 'poisson', 'beta_kl', 'beta_is', 'sqrt_model', 'top3', 'top5',
                           'library_prior', 'nmf2_varexp', 'nmf3_varexp', 'nmf5_varexp', 'nmf5_share',
                           'gated_nmf5_varexp', 'gated_nmf5_varexp_beta_kl']
    if (!(params.diathem_weight in diathem_weights)) {
        error "diathem_weight must be one of ${diathem_weights.join(', ')}, got: '${params.diathem_weight}'"
    }
    if (!(params.diathem_max_fragments instanceof Integer) || params.diathem_max_fragments < 0) {
        error "diathem_max_fragments must be a whole number, 0 for every fragment, got: '${params.diathem_max_fragments}'"
    }
    if (!(params.diathem_baseline instanceof Boolean)) {
        error "diathem_baseline must be true or false, got: '${params.diathem_baseline}'"
    }
    if (!(params.diathem_targets_source in ['reference', 'library'])) {
        error "diathem_targets_source must be 'reference' or 'library', got: " +
              "'${params.diathem_targets_source}'"
    }
    if (params.run_diathem && has_osw_prm && params.diathem_targets_source == 'library' &&
        (params.openswath_pqp || params.openswath_tsv || params.openswath_traml)) {
        error """
            Pre-flight: diathem_targets_source = 'library' quantifies the spectral
            library the OpenSWATH library is built from. --openswath_pqp,
            --openswath_tsv and --openswath_traml are used as supplied, with no
            such library. Use diathem_targets_source = 'reference'.
        """.stripIndent().trim()
    }

    def not_mzml = { r -> !(file(r.file).getName().toLowerCase() ==~ /.*\.(mzml|raw|wiff2?|d|lcd)$/) }
    if (params.run_diathem && rows.any(not_mzml)) {
        log.warn "run_diathem: diathem reads mzML only, so it does not quantify " +
                 rows.findAll(not_mzml).collect { it.sample_id }.join(', ')
    }
    def diathemable = { engine, mode -> rows.any { r ->
        (r.engine?.trim() ?: 'encyclopedia').toLowerCase() == engine &&
        r.mode.toUpperCase() == mode && !not_mzml(r)
    } }
    def diathem_ok_enc_dia = params.run_diathem && diathemable('encyclopedia', 'DIA')
    def diathem_ok_enc_prm = params.run_diathem && diathemable('encyclopedia', 'PRM')
    def diathem_ok_osw_dia = params.run_diathem && diathemable('openswath', 'DIA')
    def diathem_ok_osw_prm = params.run_diathem && diathemable('openswath', 'PRM')

    rows.findAll { sig(it) != ['openswath', 'DIA'] }.groupBy(sig).each { route, rs ->
        def lists = rs.collect { ref_of(it) }.unique()
        if (lists.size() > 1) {
            error """
                Pre-flight: the ${route[0]} ${route[1]} rows name different reference lists:
                  ${lists.collect { it ?: '(blank)' }.join('\n                  ')}
                One list serves every run of an engine and mode. Name the same file on
                all their rows, or leave it blank on all of them, or run the pipeline
                once per list.
            """.stripIndent().trim()
        }
    }

    def diathem_on =[ ['encyclopedia', 'DIA']: diathem_ok_enc_dia, ['encyclopedia', 'PRM']: diathem_ok_enc_prm,
                       ['openswath', 'DIA']: diathem_ok_osw_dia,    ['openswath', 'PRM']: diathem_ok_osw_prm ]
    rows.findAll { ref_of(it) }.unique { [ref_of(it), sig(it)] }.each { r ->
        def first  = file(ref_of(r)).withReader { it.readLine() } ?: ''
        def header = first.split(first.count('\t') >= first.count(',') ? '\t' : ',', -1)
                          .collect { it.trim().replaceAll('"', '') }
        def low    = header.collect { it.toLowerCase() }
        def problem = null
        if (sig(r) == ['encyclopedia', 'PRM'] &&
            !(low.size() >= 7 && low[0] == 'compound' && low[3] == 'm/z' && low[4] == 'z')) {
            problem = "on the EncyclopeDIA PRM route it is read by column position"
        } else if (diathem_on[sig(r)] && !['Compound', 'm/z', 'z'].every { it in header }) {
            problem = "diathem reads it by the column names Compound, m/z and z"
        }
        if (problem) {
            error """
                Pre-flight: reference_list ${ref_of(r)} (sample ${r.sample_id}):
                ${problem}. Supply it in the isolation-list layout Skyline exports,
                with a numeric RT and window on every row:
                  Compound,Formula,Adduct,m/z,z,RT Time (min),Window (min)[,isDecoy]
                found: ${header.join(',')}
            """.stripIndent().trim()
        }
    }

    def gives_primary = { engine, mode ->
        switch (params.primary_abundance) {
            case 'encyclopedia': return engine == 'encyclopedia'
            case 'openswath':    return engine == 'openswath'
            case 'tric':         return engine == 'openswath' && mode == 'DIA' && params.run_tric
            case 'diathem':      return params.run_diathem && diathemable(engine, mode)
            case 'standard_peak': return engine == 'encyclopedia' && (params.heavy_label ?: []).size() > 0
            default:             return false
        }
    }
    if (has_enc && !enc_quant && params.primary_abundance == 'encyclopedia') {
        error "Pre-flight: encyclopedia_quant is off, but primary_abundance is encyclopedia, whose amounts it " +
              "gives. Set encyclopedia_quant to auto, or choose another primary_abundance."
    }
    def no_primary = routes.findAll { r -> !gives_primary(r[0], r[1]) }
    if (no_primary) {
        error """
            Pre-flight: primary_abundance = '${params.primary_abundance}' gives no
            abundance for ${no_primary.collect { "${it[0]} ${it[1]}" }.join(', ')}, so the
            abundance column would be empty. Choose one that runs there:
              encyclopedia  engine encyclopedia
              openswath     engine openswath
              diathem       run_diathem = true, runs given as or converted to mzML
              tric          engine openswath, mode DIA, run_tric = true
              standard_peak engine encyclopedia, SIL runs (heavy_label)
        """.stripIndent().trim()
    }

    def runs_by_route = rows
        .groupBy { r -> [ (r.engine?.trim() ?: 'encyclopedia').toLowerCase(),
                          r.mode.toUpperCase() ] }
        .collectEntries { k, v -> [ (k): v.size() ] }
    runs_by_route.each { route, n ->
        if (n < params.quant_min_runs) {
            def msg = "quantification cohort for ${route[0]}/${route[1]} is ${n} " +
                      "run(s), below quant_min_runs=${params.quant_min_runs}. " +
                      "Each SIL run's retention times are aligned to the other runs there " +
                      "are instead of the 4 acquired nearest in time, a single run is " +
                      "searched without retention-time windows, and cross-run quantification " +
                      "(EncyclopeDIA, diathem consistency, pyprophet, TRIC) is less reliable " +
                      "at this cohort size; identifications and per-run amounts are still reported."
            if (params.quant_require_min_runs) {
                error """
                    Pre-flight: ${msg}

                    Either supply more runs, lower quant_min_runs, or set
                    quant_require_min_runs = false to downgrade this to a warning.
                """.stripIndent().trim()
            }
            log.warn "Pre-flight: ${msg} Treat abundances as within-run values."
        }
    }

    samples = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def engine = (row.engine?.trim() ?: 'encyclopedia').toLowerCase()
            if (!(engine in ['encyclopedia', 'openswath'])) {
                error "samplesheet: unknown engine '${engine}' for ${row.sample_id}"
            }
            def meta = [id: row.sample_id, mode: row.mode.toUpperCase(), engine: engine]
            def f    = file(row.file)
            def rl_val = ref_of(row)
            def rl   = rl_val ? file(rl_val) : file("${projectDir}/assets/NO_FILE")
            if (engine == 'openswath' && f.getName().toLowerCase().endsWith('.dia')) {
                error ".dia is encyclopedia-specific; sample '${meta.id}' " +
                      "uses engine=openswath. Provide the .raw or a .mzML instead."
            }
            if (engine == 'encyclopedia' && f.getName().toLowerCase().endsWith('.mzxml')) {
                error "sample '${meta.id}': EncyclopeDIA does not read mzXML. Provide the " +
                      "vendor file or an .mzML instead."
            }
            tuple(meta, f, rl)
        }

    fasta            = file(params.fasta)
    no_file          = file("${projectDir}/assets/NO_FILE")
    standard_amts_ch = params.standard_amounts ? file(params.standard_amounts) : no_file
    diathem_lib_ch   = params.diathem_library ? file(params.diathem_library) : no_file

    def VENDOR_EXT = ['.raw', '.wiff', '.wiff2', '.d', '.lcd']
    def is_vendor  = { f -> VENDOR_EXT.any { f.getName().toLowerCase().endsWith(it) } }

    def sidecars_of = { f ->
        ['.scan', '2.scan', '.timeseries.data'].collect { file("${f}${it}") }
                                               .findAll { it.exists() }
    }

    ext_split = samples.branch { meta, f, _ml ->
        raw:   is_vendor(f)
        ready: true
    }

    MSCONVERT(ext_split.raw.map { meta, f, rl -> tuple(meta, f, sidecars_of(f), rl) })

    converted = MSCONVERT.out.mzml.flatMap { meta, mzml, rl ->
        def files = (mzml instanceof List) ? mzml : [mzml]
        if (files.size() == 1) return [ tuple(meta, files[0], rl) ]
        log.warn "MSCONVERT: ${meta.id} is a multi-sample acquisition -> " +
                 "${files.size()} runs: ${files.collect{ it.baseName }.join(', ')}"
        files.collect { f -> tuple(meta + [id: f.baseName], f, rl) }
    }

    normalized = ext_split.ready.mix(converted)

    REMOVE_SURVEY_SCANS(normalized)

    def warning_chs = []
    def runs_ready = REMOVE_SURVEY_SCANS.out.run
    def calibration_library = params.blib ?: params.openswath_blib ?: params.elib ?: params.dlib
    def mz_on = Args.mode(params.fragment_mz_calibration) == 'auto'
    def rt_on = Args.mode(params.rt_alignment) == 'auto'
    def located  = Channel.empty()
    def acquired = Channel.empty()
    def fragment_errors = Channel.empty()
    if ((mz_on || rt_on) && calibration_library) {
        def fragment_of = { meta -> meta.engine == 'openswath' ? [] :
            Args.fragmentTolerance(meta.mode == 'DIA' ? params.dia_search_args : params.context_search_args) }
        CALIBRATE_FRAGMENT_MZ(REMOVE_SURVEY_SCANS.out.run.join(REMOVE_SURVEY_SCANS.out.windows)
                                  .map { meta, f, rl, w -> tuple(meta, f, rl, w, fragment_of(meta)) },
                              file(calibration_library), mz_on)
        runs_ready = CALIBRATE_FRAGMENT_MZ.out.run
        located    = CALIBRATE_FRAGMENT_MZ.out.located
        acquired   = CALIBRATE_FRAGMENT_MZ.out.row
        fragment_errors = CALIBRATE_FRAGMENT_MZ.out.tolerance
        CALIBRATE_FRAGMENT_MZ.out.row
            .collectFile(name: 'fragment_mz_calibration.tsv', storeDir: "${params.outdir}/qc",
                         keepHeader: true, skip: 1, sort: true)
        warning_chs << CALIBRATE_FRAGMENT_MZ.out.warnings.flatMap { Warnings.inFile(it, 'calibration') }
    } else if (mz_on || rt_on) {
        log.info "Fragment m/z calibration and retention-time alignment: no spectral library given " +
                 "(blib, openswath_blib, elib or dlib), so the runs are used as acquired and searched " +
                 "without retention-time windows."
    }
    def rt_ready = rt_on && calibration_library as boolean

    by_route = runs_ready.branch { meta, _f, _ml ->
        enc_dia: meta.engine == 'encyclopedia' && meta.mode == 'DIA'
        enc_prm: meta.engine == 'encyclopedia' && meta.mode == 'PRM'
        osw_dia: meta.engine == 'openswath'    && meta.mode == 'DIA'
        osw_prm: meta.engine == 'openswath'    && meta.mode == 'PRM'
    }

    def enc_prm_row = has_enc_prm
        ? rows.find {
            (it.engine?.trim() ?: 'encyclopedia').toLowerCase() == 'encyclopedia' &&
            it.mode.toUpperCase() == 'PRM' }
        : null
    def _enc_prm_ref = enc_prm_row ? ref_of(enc_prm_row) : ''
    def enc_prm_sample_rl = _enc_prm_ref ? file(_enc_prm_ref) : null

    def enc_dia_row = has_enc_dia
        ? rows.find {
            (it.engine?.trim() ?: 'encyclopedia').toLowerCase() == 'encyclopedia' &&
            it.mode.toUpperCase() == 'DIA' && ref_of(it) }
        : null
    def _enc_dia_ref = enc_dia_row ? ref_of(enc_dia_row) : ''
    def enc_dia_sample_rl = _enc_dia_ref ? file(_enc_dia_ref) : null

    def enc_labels = has_enc ? (params.heavy_label ?: []) : []

    if (has_enc) {
        def enc_dia_has_reference_list =
                has_enc_dia && (enc_dia_sample_rl != null || params.blib)

        PREPARE_LIBRARY_ENCYCLOPEDIA(fasta, enc_labels, bg_source,
                        params.background_library, params.background_fasta,
                        params.background_min_targets, no_file)
        warning_chs << PREPARE_LIBRARY_ENCYCLOPEDIA.out.warnings

        if (has_enc_prm) {
            def enc_prm_samples
            def enc_prm_targets
            if (enc_prm_sample_rl != null) {
                NORMALIZE_REFERENCE_MODS(Channel.value(enc_prm_sample_rl),
                                         PREPARE_LIBRARY_ENCYCLOPEDIA.out.library)
                def enc_prm_rl = NORMALIZE_REFERENCE_MODS.out.reference_list
                CHECK_REFERENCE_COVERAGE(enc_prm_rl, PREPARE_LIBRARY_ENCYCLOPEDIA.out.library, 'enc_prm')
                warning_chs << CHECK_REFERENCE_COVERAGE.out.warnings.flatMap { Warnings.inFile(it, 'enc_prm') }
                if (params.add_decoys_to_reference_list &&
                    !reference_list_has_decoys(enc_prm_sample_rl)) {
                    log.warn "encyclopedia+PRM: the supplied reference_list " +
                             "'${enc_prm_sample_rl.name}' carries no decoy rows. " +
                             "Generating them, because mProphet cannot estimate a " +
                             "reference-side FDR without a null distribution. Set " +
                             "add_decoys_to_reference_list = false to keep the list as-is."
                    ADD_DECOYS_TO_REFERENCE_LIST(
                        enc_prm_rl,
                        Channel.value(file(params.reference_list_decoy_jar)))
                    enc_prm_targets = ADD_DECOYS_TO_REFERENCE_LIST.out.reference_list
                } else {
                    enc_prm_targets = enc_prm_rl
                }
                enc_prm_samples = by_route.enc_prm
                    .combine(enc_prm_targets)
                    .map { meta, f, _orig_rl, rl -> tuple(meta, f, rl) }
            } else if (params.blib) {
                enc_prm_samples = by_route.enc_prm
                    .combine(PREPARE_LIBRARY_ENCYCLOPEDIA.out.reference_list_derived)
                    .map { meta, f, _orig_rl, derived -> tuple(meta, f, derived) }
                enc_prm_targets = PREPARE_LIBRARY_ENCYCLOPEDIA.out.reference_list_derived
            } else {
                error "encyclopedia+PRM samples need a reference_list, either " +
                      "in the samplesheet's reference_list column, or via " +
                      "--blib for automatic " +
                      "derivation from a vendor .blib."
            }
            ENCYCLOPEDIA_PRM(enc_prm_samples, PREPARE_LIBRARY_ENCYCLOPEDIA.out.library, fasta,
                         diathem_lib_ch, standard_amts_ch, params.heavy_label,
                         enc_prm_targets, no_file, diathem_ok_enc_prm, located, acquired, fragment_errors,
                         rt_ready, enc_quant)
            warning_chs << ENCYCLOPEDIA_PRM.out.warnings
        }

        if (has_enc_dia) {
            def enc_dia_targets = no_file
            if (enc_dia_has_reference_list) {
                enc_dia_targets = (enc_dia_sample_rl != null)
                    ? Channel.value(enc_dia_sample_rl)
                    : PREPARE_LIBRARY_ENCYCLOPEDIA.out.reference_list_derived
            }
            ENCYCLOPEDIA_DIA(by_route.enc_dia, PREPARE_LIBRARY_ENCYCLOPEDIA.out.library, fasta,
                         diathem_lib_ch, standard_amts_ch, params.heavy_label,
                         enc_dia_targets, enc_dia_has_reference_list, no_file,
                         diathem_ok_enc_dia, located, acquired, fragment_errors, rt_ready, enc_quant)
            warning_chs << ENCYCLOPEDIA_DIA.out.warnings
        }
    }

    if (has_osw) {
        def openswath_prm_refs = has_osw_prm ? rows.findAll {
            (it.engine?.trim() ?: 'encyclopedia').toLowerCase() == 'openswath' &&
            it.mode.toUpperCase() == 'PRM' && ref_of(it)
        }.collect { ref_of(it) }.unique() : []
        def openswath_reflist_for_anchors = openswath_prm_refs ?
                                          file(openswath_prm_refs[0]) : no_file

        PREPARE_LIBRARY_OPENSWATH(fasta, osw_labels, bg_source,
                            params.background_library, params.background_fasta,
                            params.background_min_targets,
                            openswath_reflist_for_anchors, no_file)
        warning_chs << PREPARE_LIBRARY_OPENSWATH.out.warnings

        if (has_osw_prm) {
            def openswath_prm_row = rows.find {
                (it.engine?.trim() ?: 'encyclopedia').toLowerCase() == 'openswath' &&
                it.mode.toUpperCase() == 'PRM'
            }
            def _openswath_ref = ref_of(openswath_prm_row)
            def sample_rl = _openswath_ref ? file(_openswath_ref) : null

            def openswath_prm_samples
            def openswath_prm_targets
            if (sample_rl != null) {
                openswath_prm_samples = by_route.osw_prm
                openswath_prm_targets = Channel.value(sample_rl)
            } else if (params.openswath_blib || params.blib) {
                openswath_prm_samples = by_route.osw_prm
                    .combine(PREPARE_LIBRARY_OPENSWATH.out.reference_list_derived)
                    .map { meta, f, _orig_rl, derived -> tuple(meta, f, derived) }
                openswath_prm_targets = PREPARE_LIBRARY_OPENSWATH.out.reference_list_derived
            } else {
                error "openswath+PRM samples need a reference_list, either " +
                      "in the samplesheet's reference_list column, or via " +
                      "--openswath_blib or --blib for automatic " +
                      "derivation from a vendor .blib."
            }

            OPENSWATH_PRM(openswath_prm_samples,
                                   PREPARE_LIBRARY_OPENSWATH.out.library,
                                   PREPARE_LIBRARY_OPENSWATH.out.irt_library,
                                   PREPARE_LIBRARY_OPENSWATH.out.full_window_flag,
                                   PREPARE_LIBRARY_OPENSWATH.out.spectral_library,
                                   diathem_lib_ch, standard_amts_ch,
                                   params.heavy_label,
                                   openswath_prm_targets, no_file, diathem_ok_osw_prm, fasta)
            warning_chs << OPENSWATH_PRM.out.warnings
        }
        if (has_osw_dia) {
            OPENSWATH_DIA(by_route.osw_dia,
                                   PREPARE_LIBRARY_OPENSWATH.out.library,
                                   PREPARE_LIBRARY_OPENSWATH.out.irt_library,
                                   PREPARE_LIBRARY_OPENSWATH.out.full_window_flag,
                                   PREPARE_LIBRARY_OPENSWATH.out.spectral_library,
                                   diathem_lib_ch, standard_amts_ch,
                                   params.heavy_label,
                                   no_file, diathem_ok_osw_dia, fasta)
            warning_chs << OPENSWATH_DIA.out.warnings
        }
    }

    if (warning_chs) {
        def warnings = Channel.empty().mix(*warning_chs)
        warnings.subscribe { log.warn it }
        warnings.collectFile(name: 'warnings.txt', storeDir: params.outdir, newLine: true, sort: true)
    }
}
