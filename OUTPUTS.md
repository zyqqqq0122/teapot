# TEAPOT results

What a run writes, and what the columns mean.

Your results are in `results/quant/<route>/`. `<route>` is the engine and the
acquisition mode of your runs: `enc_prm`, `enc_dia`, `osw_prm` or `osw_dia`, one
folder for each in your samplesheet. The rest of `results/` is what the individual
steps wrote, kept so you can check them.

## Which file to read

| your runs | read | it holds |
|---|---|---|
| label-free | `quant.reference.confident.tsv` | the reference peptides each run detected, with their abundances |
| with SIL internal standards | `peptide.abs_quant.confident.tsv` | the light/heavy peptide pairs that give an amount, with the amount in pmol |

`quant_summary.tsv` counts, for each run, the reference peptides (or light/heavy
pairs) it detected and quantified.

A peptide is **detected** in a run when its $q$ value is ≤ `fdr`. On SIL runs a
light peptide also has to elute with its heavy standard, as a confident pair's does
([below](#absolute-amounts)): one the search identified elsewhere is not detected.
A peptide is **quantified** when the `primary_abundance` tool gave
it an abundance.

---

## Peptide tables

| file | rows |
|---|---|
| `quant.reference.tsv` | every reference peptide in every run, detected or not |
| `quant.reference.confident.tsv` | the same, only where the run detected them |
| `quant.full.tsv` | everything the searches reported, also background peptides not on your list. For troubleshooting |

Your **reference peptides** are the peptides on your reference list. Without a
list, as on a DIA run, they are the peptides detected in at least one of your runs.

Peptides from iRT tags are listed only in `quant.full.tsv`.

| column | meaning |
|---|---|
| `sample_id` | the run, as named in your samplesheet |
| `protein` | the protein, as your library names it |
| `peptide`, `peptidoform` | the sequence, without and with its modifications |
| `charge` | precursor charge |
| `channel` | `light` or `heavy`. Label-free runs are all `light` |
| `detected` | `true` when `q_value` is at most `fdr` and, for a light peptide on SIL runs, it elutes with its standard |
| `q_value`, `pep` | the identification's $q$ value and posterior error probability |
| `abundance` | the abundance to use, from the tool `primary_abundance` names |
| `abundance_encyclopedia`, `abundance_openswath`, `abundance_diathem`, `abundance_tric` | each tool's own abundance, for the tools that ran (EncyclopeDIA's: see `encyclopedia_quant` in [PARAMETERS.md](PARAMETERS.md)) |
| `abundance_standard_peak` | SIL runs, EncyclopeDIA routes: measured by the pipeline over the heavy standard's peak, with the fragments both forms share ([below](#the-light-measured-over-its-standards-peak-standard_peak)) |
| `rt_apex_seconds` | retention time of the peak the search identified |
| `rt_encyclopedia_seconds`, `rt_diathem_seconds` | retention time where that tool measured the abundance |
| `coelution_standard_peak` | SIL runs with `standard_peak`, light peptides: how closely the light's signal follows its standard's over the standard's peak (1 = exactly) |
| `consistency`, `n_effective_transitions` | diathem only: how well the measurement fits the expected fragments, from 0 to 1, and how many fragments support it. Leave out measurements below `quant_min_runs` |

`quant.full.tsv` also has `in_reference_list` (`true` for your reference peptides),
`is_decoy`, and on OpenSWATH DIA `peak_group_rank` (the other tables keep rank 1).

## `quant_summary.tsv`

One row per run.

| column | label-free runs |
|---|---|
| `reference_peptides` | your reference peptides |
| `detected` | how many of them this run detected, in at least one charge state |
| `quantified` | how many of those were also quantified |

| column | runs with SIL internal standards |
|---|---|
| `reference_pairs` | your reference peptides that have a light and a heavy form |
| `heavy_detected`, `light_detected` | how many had that form detected |
| `heavy_quantified`, `light_quantified` | how many of those were also quantified |
| `pairs_detected` | how many had both forms detected at the same charge |
| `pairs_quantified` | how many of those were quantified in both forms, the two eluting together (see the confident pairs below). The number to read |

---

## Absolute amounts

Written when your runs have SIL internal standards (`heavy_label` and
`standard_amounts` set).

| file | holds |
|---|---|
| `peptide.abs_quant.tsv` | every light/heavy pair in every run, from each tool that measured it |
| `peptide.abs_quant.confident.tsv` | the confident pairs, with their amounts |
| `peptide.abs_quant.matrix.tsv` | their amounts, one row per peptide, form and charge, one column per run |
| `protein.abs_quant.tsv` | each protein's amount in each run: the median of its peptides' amounts, with their number (`n_peptides`) |

A pair is **confident** when the `primary_abundance` tool quantified both forms,
both were detected, and the two elute together: with `standard_peak`, the light's
signal rises and falls with its standard's over the standard's peak (`coelution` at
least `pair_min_coelution`); with the other tools, their peaks lie within
`pair_rt_tolerance_seconds` of each other.

| column | meaning |
|---|---|
| `method` | the tool that measured the pair |
| `peptide`, `peptidoform` | the sequence, and its form: the modifications other than the label. Each form is a pair of its own |
| `status` | which forms have an abundance: `both_quantified`, `light_missing`, `heavy_missing` or `both_missing` |
| `confident` | the pair is in `peptide.abs_quant.confident.tsv` |
| `light_abundance`, `heavy_abundance` | the two forms' abundances. EncyclopeDIA's are summed over the same fragments in both forms, so they can be lower than the `abundance` in `quant.reference.tsv`; `standard_peak` measures both over the same fragments and the same time |
| `ratio_L_H`, `ratio_H_L` | light / heavy, and heavy / light |
| `heavy_pmol_known` | the standard's amount, from `standard_amounts` |
| `abs_pmol_light` | the endogenous peptide's amount: `ratio_L_H` × `heavy_pmol_known` |
| `light_q_value`, `heavy_q_value` | each form's $q$ value |
| `light_detected`, `heavy_detected` | that form's $q$ value is at most `fdr`; for the light, it also elutes with its standard (`co_eluting` not `false`) |
| `light_rt_seconds`, `heavy_rt_seconds` | where each form was measured; for `standard_peak`, where the search identified each form |
| `coelution` | `standard_peak`: how closely the light's signal follows its standard's over the standard's peak, from -1 to 1 (1 = exactly) |
| `co_eluting` | the two elute together: for `standard_peak`, `coelution` is at least `pair_min_coelution`; for the other tools, their peaks lie within `pair_rt_tolerance_seconds`. Blank when there was nothing to compare, or the check is off |
| `light_consistency`, `heavy_consistency` | when diathem ran: each form's `consistency` |

Every `both_quantified` pair has a ratio, but only those marked `confident` are
reliable. Peptides without a standard, and standards that exist only as heavy, are
not listed.

### The light measured over its standard's peak (`standard_peak`)

SIL runs, EncyclopeDIA routes. The heavy standard is in every run, so its peak shows where
its light peptide elutes. The pipeline measures both forms itself over the standard's peak,
with the fragments both share; a fragment whose light/heavy ratio is far off the others'
(an interference) is left out. These are the amounts of `method: standard_peak` above, and
the main amounts when `primary_abundance` is `standard_peak`, as in the SIL template.
EncyclopeDIA's own amounts stay in `method: encyclopedia`. Which peptides were identified,
and their $q$ values, are the search's: a pair is measured when the search identified its
standard, and it is confident when the search also identified the light and, over the
standard's peak, the light's signal follows the standard's (`coelution` at least
`pair_min_coelution`).

Each run's measurements are in `quant/<route>/standard_peak/<sample>.standard_peak.tsv`:

| column | meaning |
|---|---|
| `status` | `measured`, or why not: `heavy not identified`, `window shared with another target`, too few fragments, and so on |
| `heavy_q_value`, `light_q_value`, `light_rt` | the search's $q$ values, and where it identified the light |
| `apex`, `start`, `stop` | the standard's peak, in seconds; both forms are measured over `start` to `stop` |
| `light_area`, `heavy_area`, `ratio_L_H` | the two forms' signal over the shared fragments, and light / heavy |
| `coelution` | how closely the light's signal follows the standard's over the peak (1 = exactly). The pair counts from `pair_min_coelution` |
| `coelution_null` | the same, with the light's signal moved away from its standard, where the light cannot be: what `pair_min_coelution: auto` learns from |
| `rdotp` | how well the light's fragment intensities match the standard's (1 = perfectly) |
| `fragments`, `fragments_dropped` | fragments used, and fragments left out because their light/heavy ratio was more than 2-fold off the others' |

---

## Checks and warnings

The pipeline checks itself as it runs. A problem it cannot work around stops the run.
Anything else you should look at is printed as a warning while the pipeline runs and
collected in `warnings.txt`, which is written only when there is a warning: for
example a run that detected fewer than `search_min_reference_ids` reference peptides,
or reference peptides your library does not hold.

| file | route | what it holds |
|---|---|---|
| `warnings.txt` | all | every warning of the run |
| `qc/<sample>/<sample>.<engine>.search_qc.tsv` | PRM routes | each run's identifications. The run stops when a report is empty or every $q$ value came out 1 |
| `qc/encyclopedia_prm.quant_qc.tsv`, `qc/encyclopedia_dia.quant_qc.tsv` | EncyclopeDIA routes | how many measurements matched no fragment ion. The run stops when none matched any |
| `qc/enc_prm.reference_coverage.tsv` | EncyclopeDIA PRM, with your own reference list | each reference peptide and charge, and whether your library holds it (`in_library`), holds it only at another charge or only with other modifications, or does not hold it. A peptide your library does not hold cannot be detected |
| `quant/<route>/pair_rt_tolerance.tsv` | SIL runs | the `pair_rt_tolerance_seconds` used for the tools that give only where each peak lies, and how `auto` arrived at it |
| `quant/<route>/pair_min_coelution.tsv` | SIL runs measured with `standard_peak` | the `pair_min_coelution` used, and how `auto` arrived at it |
| `qc/fragment_tolerance.<route>.tsv` | EncyclopeDIA routes | the fragment tolerance the searches used (`search_tolerance`, in `unit`), how many runs it was measured on (`runs_at_it`) and where it came from (`source`); `qc/fragment_tolerance.chromatogram.tsv` for the runs of `library_sheet` |
| `qc/fragment_tolerance.<route>.runs.tsv` | the same | per run: how many fragment errors it gave (`fragment_errors`), the tolerance it would need on its own (`own_search_tolerance`) and the one it was searched at (`searched_at`) |
| `qc/fragment_mz_calibration.tsv` | targeted runs (`fragment_mz_calibration`, `rt_alignment`) | per run: whether its fragment m/z were corrected (`corrected`) and why not (`status`), how far off they were at m/z 300 to 1300 (`offset_300` ... `offset_1300`, in Th), how many of its targets served as reference points (`calibrants`) and the share that disagreed, and when the run was acquired (`acquired`). A run corrected is also a warning: the instrument's calibration may need attention |
| `qc/fragment_mz_calibration.chromatogram.tsv` | runs of `library_sheet` | the same for those runs |
| `qc/<sample>/<sample>.fragment_mz.*` | the same | each run's details: the correction curve, the reference points, the error before and after |
| `qc/rt_alignment.<route>.tsv` | SIL runs, EncyclopeDIA routes (`rt_alignment`) | per run: whether it was aligned (`status`, and why not), its standards found (`anchors`), the runs it was compared with (`neighbours`), how closely its standards' times agree with theirs (`residual_sd_seconds`), how far a peak may lie from where it is expected (`half_window_seconds`), the peptides given a window, how many windows were widened to include where the standard itself was found in the run (`windows_widened`), and when it was acquired. A run whose times agree poorly, or that could not be aligned, is also a warning |
| `qc/rt_alignment/<route>/<sample>.apex_windows.tsv` | the same | each reference peptide's window in that run, in seconds: the search took its peak only there |

## Other folders

| folder | holds |
|---|---|
| `encyclopedia_libexport/background/background.elib` | with `background_source: measured`: the background built from your runs. Later runs can reuse it with `background_source: library` |
| `encyclopedia_libexport/analysis/` | EncyclopeDIA's own quantification across your runs, when it ran (`encyclopedia_quant`) |
| `msconvert/` | the mzML converted from your vendor files |
| `pipeline_info/` | the Nextflow report, timeline and log |
| everything else | each step's own output and logs |

---

## Two things worth knowing

### Are the $q$ values per run or pooled?

| route | rescored by | scope |
|---|---|---|
| `enc_prm` | mProphet, inside Context | per run |
| `enc_dia` | Percolator, inside EncyclopeDIA | per run |
| `osw_prm` | Context-MS, or pyprophet | per run |
| `osw_dia` | pyprophet | **pooled across runs** |

### How many runs you need

EncyclopeDIA, diathem and TRIC rely on multiple runs for quantification,
and their accuracy improves as more runs are added to a cohort; each SIL run's
retention times are aligned to the 4 runs acquired nearest it. With fewer than
`quant_min_runs` (default 5) runs the pipeline warns and goes on: each run is aligned
to the other runs there are, and a single run is searched without retention-time
windows, so the results can be less accurate. `quant_require_min_runs: true` stops it
instead.
