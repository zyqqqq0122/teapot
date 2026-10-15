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

A peptide is **detected** in a run when its $q$ value is ≤ `fdr`.
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
| `detected` | `true` when `q_value` is at most `fdr` |
| `q_value`, `pep` | the identification's $q$ value and posterior error probability |
| `abundance` | the abundance to use, from the tool `primary_abundance` names |
| `abundance_encyclopedia`, `abundance_openswath`, `abundance_diathem`, `abundance_tric` | each tool's own abundance, for the tools that ran |
| `rt_apex_seconds` | retention time of the peak the search identified |
| `rt_encyclopedia_seconds`, `rt_diathem_seconds` | retention time where that tool measured the abundance |
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
| `pairs_quantified` | how many of those were quantified in both forms, with their peaks within `pair_rt_tolerance_seconds`. The number to read |

---

## Absolute amounts

Written when your runs have SIL internal standards (`heavy_label` and
`standard_amounts` set).

| file | holds |
|---|---|
| `peptide.abs_quant.tsv` | every light/heavy pair in every run, from each tool that measured it |
| `peptide.abs_quant.confident.tsv` | the confident pairs, with their amounts |
| `peptide.abs_quant.matrix.tsv` | their amounts, one row per peptide and charge, one column per run |
| `protein.abs_quant.tsv` | each protein's amount in each run: the median of its peptides' amounts, with their number (`n_peptides`) |

A pair is **confident** when the `primary_abundance` tool quantified both forms,
both were detected, and their peaks lie within `pair_rt_tolerance_seconds` of each
other.

| column | meaning |
|---|---|
| `method` | the tool that measured the pair |
| `status` | which forms have an abundance: `both_quantified`, `light_missing`, `heavy_missing` or `both_missing` |
| `confident` | the pair is in `peptide.abs_quant.confident.tsv` |
| `light_abundance`, `heavy_abundance` | the two forms' abundances. EncyclopeDIA's are summed over the same fragments in both forms, so they can be lower than the `abundance` in `quant.reference.tsv` |
| `ratio_L_H`, `ratio_H_L` | light / heavy, and heavy / light |
| `heavy_pmol_known` | the standard's amount, from `standard_amounts` |
| `abs_pmol_light` | the endogenous peptide's amount: `ratio_L_H` × `heavy_pmol_known` |
| `light_q_value`, `heavy_q_value` | each form's $q$ value |
| `light_detected`, `heavy_detected` | that form's $q$ value is at most `fdr` |
| `light_rt_seconds`, `heavy_rt_seconds` | where each form was measured |
| `co_eluting` | the two lie within `pair_rt_tolerance_seconds`. Blank when there was nothing to compare, or the check is off |
| `light_consistency`, `heavy_consistency` | when diathem ran: each form's `consistency` |

Every `both_quantified` pair has a ratio, but only those marked `confident` are
reliable. Peptides without a standard, and standards that exist only as heavy, are
not listed.

---

## Checks

The pipeline checks itself as it runs. If your run finished, the checks passed.

| file | route | what it catches |
|---|---|---|
| `qc/<sample>/<sample>.<engine>.search_qc.tsv` | PRM routes | the search identified nothing, or every $q$ value came out 1 |
| `qc/<engine>.quant_qc.tsv` | EncyclopeDIA routes | the quantification pass matched no fragment ions |

## Other folders

| folder | holds |
|---|---|
| `encyclopedia_libexport/background/background.elib` | with `background_source: measured`: the background built from your runs. Later runs can reuse it with `background_source: library` |
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
and their accuracy improves as more runs are added to a cohort. `quant_min_runs`
(default 2) warns below that; `quant_require_min_runs: true` makes it stop the
run.
