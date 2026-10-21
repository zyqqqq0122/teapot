# TEAPOT parameters

Every setting, its default, and when to change it. Set them in a YAML params file,
best a copy of one of the `templates/`, or as `--flags` on the command line:

```yaml
input:  my_samplesheet.csv
fasta:  my_proteins.fasta
fdr:    0.01
```
This page lists every setting you might want to change. `nextflow.config` holds a few
more that the pipeline needs internally; you can ignore them.

## Settings the pipeline measures from your runs

Some settings depend on the state of your runs (how far off their fragment m/z are, how their
peptides elute, how wide their peaks are), so the pipeline measures them from your runs each time
instead of using a fixed value. What it found is written to a table (see OUTPUTS.md); give a value
for the setting to use your own instead.

| setting | what is measured from your runs | to set it yourself |
|---|---|---|
| `fragment_mz_calibration`, `fragment_mz_calibration_min_offset` | how far each run's fragment m/z are off, and whether to correct them ([section 10](#10-acquisition-checks)) | `off`, or a minimum offset in Th |
| `search_tolerance` | the fragment tolerance the searches need after the correction ([section 10](#10-acquisition-checks)) | `off`: the `-ftol` of your search arguments |
| `rt_alignment`, `rt_alignment_half_window_seconds` | where each peptide elutes in each run, and how far from there its peak may lie ([section 7](#7-retention-time)) | `off`, or a number of seconds |
| `pair_min_coelution` | how closely a light peptide's signal must follow its standard's to count ([section 8](#8-heavylight-channels-and-absolute-quantification)) | a number from 0 to 1 |
| `pair_rt_tolerance_seconds` | how far apart a light's and its standard's peaks may lie, for the tools other than `standard_peak` ([section 8](#8-heavylight-channels-and-absolute-quantification)) | a number of seconds |
| `encyclopedia_quant` | whether EncyclopeDIA's cross-run quantification is needed ([section 11](#11-quantification)) | `on` or `off` |

---

## 1. Required

| setting | default | |
|---|---|---|
| `input` | required | your samplesheet |
| `fasta` | required | protein FASTA |
| `outdir` | `results` | where results go |

## 2. Run files

Vendor formats are converted for you. These two settings are passed to
msconvert.

| setting | default | |
|---|---|---|
| `peakpicking` | `--filter "peakPicking vendor msLevel=1-"` | centroids the spectra, using the instrument vendor's own algorithm |
| `otherfilters` | `''` | further msconvert options, usually more `--filter` arguments |

## 3. The assay library

Supply exactly **one**.

| setting | default | |
|---|---|---|
| `elib` | `null` | an EncyclopeDIA chromatogram library. Best quality |
| `dlib` | `null` | an EncyclopeDIA spectrum library |
| `blib` | `null` | a Skyline / BiblioSpec library, such as a vendor kit's |
| `library_sheet` | `null` | GPF runs to build a chromatogram library from, searched with `dlib` or with a predicted library (`use_koina`) |
| `use_koina` | `false` | predict a library from `fasta` |
| `koina_models` | `Prosit_2020_HCD;IM2Deep_CCS;Prosit_2019_iRT` | prediction models |
| `koina_url` | `https://koina.wilhelmlab.org:443/` | prediction server |
| `openswath_pqp` | `null` | OpenSWATH only: a ready `.pqp`, used as it is |
| `openswath_tsv` | `null` | OpenSWATH only: a transition list |
| `openswath_traml` | `null` | OpenSWATH only: a TraML |
| `openswath_blib` | `null` | OpenSWATH only: equivalent to `blib` |

## 4. The background library

Only Context and Context-MS need one, only when your library holds nothing but
your reference peptides, as a vendor kit's library does.

| `background_source` | the background | you supply |
|---|---|---|
| `library` | a library you already have | `background_library`, and `background_fasta` if it is a `.blib` |
| `predicted` | predicted from the proteins of your samples | `background_fasta`, best with `background_peptides` |
| `measured` | the peptides found in a few of your runs, such as QC injections | `background_sheet`, `background_fasta`, and if you like `background_search_library` |

A setting the chosen source does not use stops the run.

| setting | default | |
|---|---|---|
| `background_source` | `none` | `none`, `library`, `predicted` or `measured` |
| `background_library` | `null` | an `.elib`, `.dlib` or `.blib` |
| `background_fasta` | `null` | the proteins that may be in your samples. Not the kit's FASTA |
| `background_peptides` | `null` | `predicted`: which peptides to keep: a PeptideAtlas build's tables as downloaded, or your own list ([INPUTS.md](INPUTS.md#peptides-to-keep-background_peptides)) |
| `background_peptides_top` | `5000` | how many of a PeptideAtlas build's most often seen peptides to keep |
| `background_sheet` | `null` | `measured`: the runs, listed like `library_sheet` with each run under a `condition` of its own ([template](templates/background_sheet.csv)). Your PRM runs can be used as they are |
| `background_search_library` | `null` | `measured`: a library to search the runs with, instead of a prediction from `background_fasta` |
| `background_min_targets` | `500` | stop if the background has fewer peptides than this |

With `measured`, also set `sheet_search_args` for your instrument (§6).

## 5. Modifications

How a mass shift in your library is matched to a named UniMod modification.

| setting | default | |
|---|---|---|
| `mod_mass_tolerance` | `0.01` | tolerance in Da for recognizing a modification in your library by its mass |
| `unmapped_mod_max_fraction` | `0.05` | stop if more than this fraction of library entries carries a modification that is not recognized |

## 6. Instrument settings

The defaults suit an Orbitrap. For other instruments, set the fragment tolerance:

| your instrument | EncyclopeDIA | OpenSWATH | diathem |
|---|---|---|---|
| ion trap (Thermo Stellar, LTQ) | `-ftol 0.4 -ftolunits AMU -lftol 0.4 -lftolunits AMU -frag HCD` | `-mz_extraction_window 0.8 -mz_extraction_window_unit Th` | `diathem_fragment_tol_th: 0.4` |
| TOF (SCIEX TripleTOF) | `-ftol 20 -ftolunits PPM -frag CID` | `-mz_extraction_window 50 -mz_extraction_window_unit ppm -mz_extraction_window_ms1 20 -mz_extraction_window_ms1_unit ppm` | `diathem_fragment_ppm: 20` |
| Orbitrap | the default | the default | the default |

OpenSWATH's `-mz_extraction_window` is the full width of the window.

Give every EncyclopeDIA search your run uses the same values. For targeted runs, the
pipeline then measures the tolerance your runs need and searches at it
([section 10](#10-acquisition-checks)); your value is where it starts. EncyclopeDIA's `-frag`
chooses the fragment ions it uses (`HCD`: y ions; `CID`: b and y ions), not your
instrument's fragmentation.

| setting | default | |
|---|---|---|
| `dia_search_args` | `''` | EncyclopeDIA, `mode: DIA` runs |
| `context_search_args` | `''` | EncyclopeDIA, `mode: PRM` runs |
| `sheet_search_args` | `''` | EncyclopeDIA, the runs in `library_sheet` or `background_sheet` |
| `libexport_args` | `''` | extra options for EncyclopeDIA's quantification. It already inherits the search's tolerances and `-frag` |
| `openswath_search_args` | Orbitrap settings (30 ppm extraction), see `nextflow.config` | [OpenSwathWorkflow](https://openms.de/documentation/html/TOPP_OpenSwathWorkflow.html) options. Setting it replaces the whole default, so start from a template's list |
| `openswath_assaygen_args` | `''` | extra [OpenSwathAssayGenerator](https://openms.de/documentation/html/TOPP_OpenSwathAssayGenerator.html) options |
| `openswath_decoygen_args` | `''` | extra [OpenSwathDecoyGenerator](https://openms.de/documentation/html/TOPP_OpenSwathDecoyGenerator.html) options |
| `java_mem` | `32g` | memory for EncyclopeDIA |

## 7. Retention time

The pipeline puts every library on one retention time (RT) scale, using anchor
peptides ([INPUTS.md](INPUTS.md#rt-anchors)).

On SIL runs of the EncyclopeDIA routes (in both modes, PRM and DIA), the runs are also
aligned with each other on their heavy standards, which are in every run. Before any
search, the pipeline finds where each standard elutes in each run and compares the run with
the runs acquired just before and after it. From that it knows where every peptide of your
reference list should elute in that run, and the search looks for the peptide, its light
form and its decoy only within a short window of time around that point. When a standard is
found in the run itself outside that window, the window is widened to include where it was
found. A light peptide is therefore always looked for where its standard elutes. Each run's
alignment is in `qc/rt_alignment.<route>.tsv` and its windows in
`qc/rt_alignment/<route>/`. A run that cannot be aligned (too few standards found, or no
acquisition time in the file) is searched without windows, with a warning.

| setting | default | |
|---|---|---|
| `rt_alignment` | `auto` | `auto`: as above. `off`: every peptide is looked for over the whole time its isolation window was acquired. Needs a spectral library (`blib`, `elib` or `dlib`) |
| `rt_alignment_half_window_seconds` | `auto` | how far, in seconds, a peptide's peak may lie from where it is expected. `auto`: set for each run from how closely its standards' times agree with its neighbours' |

| setting | default | |
|---|---|---|
| `irt_library` | `null` | your own iRT library (TraML, TSV or PQP) |
| `irt_anchors` | `null` | your own anchor peptides. **Required together with `irt_library`** |
| `irt_top_transitions` | `6` | transitions kept per anchor peptide |
| `rt_harmonize_fallback` | `fail` | `fail` or `full_window`: search the whole gradient instead when the RT scales cannot be matched. Correct, but much slower |
| `rt_harmonize_min_anchors` | `5` | fewest anchors accepted |
| `rt_harmonize_min_r2` | `0.98` | lowest fit quality accepted |
| `rt_harmonize_min_anchor_coverage` | `0.8` | warn if the anchors span less of the library than this |
| `rt_coherence_max_range_ratio` | `3.0` | check that a library and its background ended up on one scale |
| `rt_coherence_require_iqr_overlap` | `true` | the same check |

## 8. Heavy/light channels and absolute quantification

Leave `heavy_label` empty for label-free runs.

| setting | default | |
|---|---|---|
| `heavy_label` | `[]` | the heavy labels, e.g. `[{site: K, unimod: 'UniMod:259'}, {site: R, unimod: 'UniMod:267'}]`. `site` is a residue, `n_term` or `c_term` |
| `standard_amounts` | `null` | the known pmol of your SIL internal standards ([INPUTS.md](INPUTS.md#known-amounts---standard_amounts)) |
| `pair_min_coelution` | `auto` | a light peptide and its standard elute together. With `primary_abundance: standard_peak` (the SIL template), a pair gets an amount only when the light's signal rises and falls with its standard's over the standard's peak: this is the lowest correlation of the two accepted. `auto` learns it from your runs: the pipeline also reads each light's signal moved away from its standard, where the light cannot be, and takes the correlation only 1 in 20 of those reach, so that a light that does not elute with its standard rarely passes. Or give a number from 0 to 1 (1: they follow each other exactly); `0` turns the check off. The value used is in `quant/<route>/pair_min_coelution.tsv` |
| `pair_rt_tolerance_seconds` | `auto` | the same check for the other tools, which give only where each peak lies: a pair gets an amount only when its two peaks lie within this many seconds. `auto` takes half the typical width of your heavy standards' peaks in these runs; give a number of seconds to set it yourself, or `0` to turn the check off |
| `complete_channels_skip_proteins` | `[]` | proteins that exist only as heavy, such as `['iRT_Tag']`: their peptides give no amount |
| `min_consistency` | `0.0` | leave out measurements whose diathem `consistency` is lower |
| `complete_channels_ion_types` | `['b','y']` | fragment ions of the light or heavy forms the pipeline adds to your library. Only `b` and `y` are supported |
| `complete_channels_max_frag_charge` | `2` | highest fragment charge of those forms |
| `complete_channels_tolerance_ppm` | `20.0` | tolerance for matching the library's fragment peaks to their ions |
| `complete_channels_coverage_hard_floor` | `0.5` | stop the run if any one spectrum in your library has less than this share of its signal in `b` and `y` ions: its other form could not be made reliably |
| `complete_channels_coverage_warn_below` | `0.9` | a spectrum below this still gets its other form, made from its `b` and `y` ions only, and is marked in `complete_channels/complete_channels_qc.tsv` |

## 9. Confidence estimation

| setting | default | |
|---|---|---|
| `fdr` | `0.01` | the $q$ value threshold: a peptide is detected at or below it |
| `openswath_prm_fdr` | `context` | OpenSWATH PRM: `context` (needs a background) or `pyprophet` (does not) |
| `add_decoys_to_reference_list` | `true` | add decoys to a reference list that has none |
| `search_min_reference_ids` | `10` | warn when a run detects fewer reference peptides (on screen and in `warnings.txt`). The pipeline stops only when no run detects any |
| `pyprophet_level` | `ms2` | `ms2` or `ms1ms2` |
| `pyprophet_learn_args` | XGBoost classifier | add `--pi0_lambda 0 0 0` if the run stops with `pi0 <= 0` |
| `pyprophet_apply_args` | `--classifier=XGBoost` | |
| `pyprophet_export_args` | `''` | |
| `context_engine` | `mprophet` | Context-MS: `mprophet` or `percolator` |
| `context_input_profile` | `auto` | Context-MS `mprophet`: which features it uses: `auto`, `pin` or `encyclopedia` |
| `context_seed_coefficients` | `encyclopedia` | Context-MS `mprophet`: its starting weights: `encyclopedia`, `none`, or a JSON file |
| `context_args` | `''` | extra Context-MS options |

## 10. Acquisition checks

Wide survey scans are left out of the runs in your samplesheet, and each peptide is
measured only in its own isolation window; there is nothing to set.

Before any search, the pipeline also checks how far off each targeted run's fragment m/z
are, using the run's own targets, each read in its own isolation window, and your library, and
corrects runs that are clearly off. Two different peptides with the same precursor m/z share a
window and are left out of this check.
Runs with wide isolation windows (DIA) are used as they are. What it found for every run is
in `qc/fragment_mz_calibration.tsv`. The runs of `library_sheet` are checked and corrected the
same way, with the library they are searched with (`qc/fragment_mz_calibration.chromatogram.tsv`).
The runs of `background_sheet` are used as acquired: a background only needs to be roughly
right.
The targets it finds in each run are also what the
retention-time alignment uses ([section 7](#7-retention-time)), so the runs are still read
for that when the correction is off.

After the correction, the same targets show how far each run's fragments still lie from
where they should, and so how wide the search's fragment tolerance needs to be. With
`search_tolerance: auto`, the EncyclopeDIA searches of these runs use the value measured on
all your runs together instead of the `-ftol` you gave; a run it could not be measured in keeps
your `-ftol`. Your `-ftol` is where the check starts looking, so give one wide enough for your
instrument before any correction ([section 6](#6-instrument-settings)). The value used is in
`qc/fragment_tolerance.<route>.tsv` (and `qc/fragment_tolerance.chromatogram.tsv` for the runs
of `library_sheet`).

| setting | default | |
|---|---|---|
| `fragment_mz_calibration` | `auto` | `auto`: measure every targeted run, and correct the runs whose fragments are off by more than `fragment_mz_calibration_min_offset`. `off`: use the runs as acquired. Needs a spectral library (`blib`, `openswath_blib`, `elib` or `dlib`) |
| `fragment_mz_calibration_min_offset` | `auto` | a run whose fragments are off by less than this anywhere is left as it is: correcting a well-calibrated run gains nothing. `auto`: a quarter of your search's fragment tolerance (`-ftol` in `dia_search_args`, `context_search_args` or `sheet_search_args`); or a number, in Th (m/z) |
| `search_tolerance` | `auto` | `auto`: search the targeted runs at the fragment tolerance their own targets show after the correction (above). `off`: search at your `-ftol` |

## 11. Quantification

| setting | default | |
|---|---|---|
| `primary_abundance` | `encyclopedia` | the tool whose abundance fills `abundance` and gives the amounts. It has to run on your data: `encyclopedia` or `openswath` with that engine, `tric` on OpenSWATH DIA runs with `run_tric: true`, `diathem` with `run_diathem: true`, `standard_peak` on SIL runs of the EncyclopeDIA routes: each light peptide and its heavy standard measured by the pipeline over the standard's peak, with the fragments both share (the SIL template sets it) |
| `run_diathem` | `true` | run diathem interference-corrected quantification |
| `run_tric` | `false` | align features across runs. OpenSWATH DIA only |
| `tric_args` | `--method LocalMST --realign_method lowess_cython --max_rt_diff 60 --mst:useRTCorrection True --mst:Stdev_multiplier 3.0 --fdr_cutoff 0.01 --max_fdr_quality 0.05` | TRIC options. Setting it replaces the whole default |
| `quant_min_runs` | `5` | warn when a route has fewer runs: each SIL run is then aligned to the other runs there are instead of the 4 acquired nearest it, a single run is searched without retention-time windows, and quantification across runs is less reliable. The run goes on |
| `quant_require_min_runs` | `false` | `true`: stop instead of warning |
| `encyclopedia_quant` | `auto` | EncyclopeDIA's own quantification across all your runs (`abundance_encyclopedia`). `auto`: only when its amounts are the ones reported (`primary_abundance: encyclopedia`). `on`: also with another `primary_abundance`, for its amounts alongside. `off`: never. It is one job over every run, the slowest step of a large batch |
| `quant_max_zero_ion_fraction` | `0.5` | warn when more than this fraction of measurements matched no fragment ion (on screen and in `warnings.txt`). If none matched any, the run stops |
| `augment_context_with_diathem` | `false` | OpenSWATH PRM: give Context-MS diathem's statistics as extra features |

### diathem

Left at `null`, `diathem_precursor_tol` and `diathem_rt_window` are chosen for you.
Without a reference list, diathem quantifies the peptides detected in at least
one of your runs.

| setting | default | |
|---|---|---|
| `diathem_weight` | `nmf5_varexp` | how diathem weighs a peptide's fragments. `nmf5_varexp`, `uniform`, `poisson`, `beta_kl`, `beta_is`, `sqrt_model`, `top3`, `top5`, `library_prior`, `nmf2_varexp`, `nmf3_varexp`, `nmf5_share`, `gated_nmf5_varexp`, or `gated_nmf5_varexp_beta_kl` |
| `diathem_max_fragments` | `8` | fit only this many of each peptide's brightest fragments, from 200 m/z. `0` fits them all |
| `diathem_baseline` | `true` | fit a flat baseline under each peak, so other peptides' signal is not counted. Unrelated to background peptides |
| `diathem_precursor_tol` | `null` | the scans that count for a peptide: those isolating within this many Th of its m/z |
| `diathem_rt_window` | `null` | half-width in seconds of the retention time window read for each peptide. `0` reads the whole run |
| `diathem_fragment_tol_th` | `null` | fragment tolerance in Th, for ion traps. Overrides `diathem_fragment_ppm` |
| `diathem_fragment_ppm` | `null` | fragment tolerance in ppm. diathem's default, 20, suits an Orbitrap |
| `diathem_library_ppm` | `null` | tolerance in ppm for matching fragments to the library |
| `diathem_library` | `null` | the library of expected fragment ratios. Defaults to the assay library |
| `diathem_targets_source` | `reference` | OpenSWATH PRM: quantify the `reference` list or the whole `library` (not with `openswath_pqp`, `openswath_tsv` or `openswath_traml`) |
| `diathem_args` | `''` | extra diathem options |

## 12. Containers

Every tool runs in a container, pulled for you. The exact versions are in
`nextflow.config`.

| setting | default | |
|---|---|---|
| `encyclopedia_container` | `ghcr.io/shannon225/encyclopedia`, a fixed version | EncyclopeDIA |
| `openswath_container` | `ghcr.io/openswath/openswath`, a fixed version | OpenMS |
| `pyprophet_container` | `ghcr.io/pyprophet/pyprophet`, a fixed version | PyProphet |
| `context_container` | `ghcr.io/shannon225/context-ms:main` | Context-MS |
| `diathem_container` | `ghcr.io/statisticalbiotechnology/diathem`, a fixed version | diathem |
| `python_container` | `quay.io/biocontainers/pandas:2.2.1` | Python |
| `msconvert_container` | `proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses`, a fixed version | msconvert |
| `container_cache_dir` | `<pipeline>/apptainer_cache` | where images are kept, to be reused by later runs |
| `container_tmp_dir` | `<pipeline>/apptainer_tmp` | scratch space for unpacking images |

## 13. Execution

| setting | default | |
|---|---|---|
| `max_cpus` | detected cores | caps a run on one machine |
| `cluster_options` | `''` | added to every cluster job, e.g. `-A <account>` |
| `work_dir` | `<outdir>/pipeline_info/work` | where intermediate files go. They can be deleted once the run has finished |
| `publish_mode` | `copy` | `copy` or `symlink`: how results are placed in `outdir`. Symlinks save space, but break once `work_dir` is deleted |

## 14. CPU, memory and time are not parameters

They are Nextflow *process directives*, set in `conf/base.config` (default) and
`conf/slurm.config` (loaded automatically by the SLURM launcher). If a step runs
out of CPUs, memory or time, raise it there.
