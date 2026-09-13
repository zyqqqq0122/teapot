# TEAPOT results

What a run writes, and what every column means.

Everything you probably want is in `results/quant/<route>/`. Everything
else in `results/` is intermediate output from the individual steps, kept so
you can check them if you need to.

`<route>` names the analysis engine and the acquisition mode: `enc_prm`, `enc_dia`,
`osw_prm` or `osw_dia`. One route directory is written per combination present
in your samplesheet.

```
results/
  quant/<route>/       the answer: what was found, and how much
  qc/                  checks the pipeline ran on itself
  pipeline_info/       report, timeline, logs, and the work directory
  ...                  output from every individual step
```

---

## The main tables: `results/quant/<route>/`

| file | rows | read it for |
|---|---|---|
| `quant_summary.tsv` | one per sample | how many peptides were targeted, identified and quantified in each run |
| `quant.reference.tsv` | every reference peptide in every sample, whatever happened to it, decoys removed | the answer for your reference list: what you monitored, and what you found |
| `quant.confident.tsv` | only peptides that passed your `fdr` threshold, decoys removed | your analysis set |
| `quant.full.tsv` | everything the run saw, including decoys and background | checking why a peptide is missing |

`quant.reference.tsv` and `quant.confident.tsv` are both filtered subsets of
`quant.full.tsv`.

`finalize_quant.log` alongside them records how the tables were built.

### `quant_summary.tsv`

One row per sample, counted from `quant.full.tsv`. Decoys are counted only in
`rows` and `decoys`.

| column | counts |
|---|---|
| `route`, `sample_id` | which run the row describes |
| `rows` | all rows for this sample |
| `decoys` | decoy rows |
| `targeted` | rows in your reference list |
| `identified` | rows passing your `fdr` threshold |
| `quantified` | rows with an abundance from any tool |
| `identified_and_quantified` | rows that are both |
| `no_signal` | rows with neither an identification nor an abundance |
| `targeted_peptides`, `identified_peptides`, `quantified_peptides` | the same counts as distinct bare peptide sequences, so heavy and light, and different charge states, count once |

`identified` and `quantified` are not limited to your reference list.

### Columns of the quant tables

**Identity**

| column | meaning |
|---|---|
| `sample_id` | from your samplesheet |
| `peptide`, `peptidoform` | the sequence; `peptidoform` carries the modifications |
| `stripped_seq` | the plain sequence, modifications removed |
| `charge` | precursor charge |
| `channel` | `heavy` or `light`. Label-free runs are all `light` |
| `protein` | protein or construct accession |

**Abundance**. One column per tool. A column is empty where that tool
did not run or had no value.

| column | source |
|---|---|
| `abundance_encyclopedia` | EncyclopeDIA: integrated area of the quantitative fragment ions |
| `abundance_openswath` | OpenSWATH: MS2 peak area |
| `abundance_diathem` | diathem: area of the interference-corrected elution profile |
| `abundance_tric` | TRIC, after cross-run alignment. OpenSWATH DIA route only |
| `abundance_primary` | a copy of the one tool `primary_abundance` names |
| `abundance_primary_source` | which tool `abundance_primary` came from |

**Quality and status**

| column | meaning |
|---|---|
| `id_qvalue` | $q$ value, from the rescoring step |
| `id_pep` | posterior error probability (PEP) for each individual identification, from the rescoring step |
| `passes_fdr` | `id_qvalue` <= `fdr` |
| `is_decoy` | a decoy peptide |
| `in_reference_list` | this peptide is in your reference list. When no reference list was supplied, every peptide is regarded as a reference peptide |
| `status` | what happened to this peptide in this sample. See the next table |
| `consistency` | diathem: how well the measurement fits the expected fragment pattern. 0 uninformative, 1 perfect |
| `n_effective_transitions` | diathem: how many fragments actually support the measurement. Near 1 means it rests on a single fragment; exclude below 2 |
| `rt_apex_seconds` | retention time at the peak apex |
| `calibrated` | both channels are present for this peptide, so a light/heavy ratio can be computed |
| `peak_group_rank` | OpenSWATH routes only. OpenSWATH reports several candidate peaks per precursor. `quant.full.tsv` keeps them all; the two filtered tables keep only the best, so each precursor appears once |

**`status` values**

| value | meaning |
|---|---|
| `identified_quantified` | passed your `fdr` threshold, and has an abundance |
| `identified_only` | passed your `fdr` threshold, but no tool gave an abundance |
| `quantified_only` | has an abundance, but did not pass your `fdr` threshold |
| `not_detected` | a tool reported it, but it neither passed the `fdr` threshold nor has an abundance |
| `no_feature` | in your reference list, but no tool reported anything for it |

**Row identity.** `(sample_id, stripped_seq, channel, charge)` is unique in
`quant.reference.tsv`.

---

## Absolute amounts (runs with SIL internal standards): `results/quant/<route>/`

Written when `heavy_label` is set and `standard_amounts` is supplied.

| file | contents |
|---|---|
| `peptide.pairs.confident.tsv` | light/heavy pairs with both channels passing your `fdr` threshold, quantified by the `primary_abundance` tool. One row per pair |
| `peptide.pairs.tsv` | the same pairs without the FDR filter, from every quantification  tool, each in its own rows with a `method` column |
| `peptide.abs_quant.long.tsv` | every peptide the pipeline tried to pair, including the ones it could not, with a `status` saying why |
| `peptide.abs_quant.matrix.tsv` | `abs_pmol_light` from the primary quantification tool, peptide × sample |
| `protein.abs_quant.long.tsv` | amounts per protein, with the number of supporting peptides (`n_peptides`) |
| `abs_quant_qc.tsv` | rows matched to a known amount, for auditing |

`quantify_heavy_light.log` records how many pairs each tool produced.

### The columns that matter

`light_amp` and `heavy_amp` are the two channels' abundances.
`ratio_L_H = light / heavy`, and

```
abs_pmol_light = ratio_L_H × heavy_pmol_known
```

`ratio_H_L` is written alongside because vendor documentation and Skyline
usually quote it the other way round.

The $q$ values come from the rescoring step. `light_qvalue`, `heavy_qvalue`
and `max_qvalue` are the identification $q$ values carried onto the row.
`both_identified` is true when the worse of the two passed your `fdr` threshold.

### `status` in `peptide.abs_quant.long.tsv`

| value | meaning |
|---|---|
| `detected` | both channels measured. The only rows carrying an `abs_pmol_light` |
| `light_missing` | the SIL internal standard was seen, the endogenous peptide was not |
| `heavy_missing` | the endogenous peptide was seen, the SIL internal standard was not |
| `both_missing` | neither channel |
| `skipped_standard` | this peptide belongs to a protein you named in `complete_channels_skip_proteins` (a vendor kit's shared iRT tag, say). It has no light counterpart by design, so its absent light channel is the expected result, not a failed measurement |

`peptide.pairs.tsv` is the `detected` rows, one set per tool.

### Protein amounts are a median, not a sum

Protein-level SIL internal standards assign identical molar amounts to all
constituent peptides. The pipeline takes the median across charge states within
a peptide, then across peptides within a protein. `n_peptides` tells you how
many supported it.

---

## Checks

Checks the pipeline ran on itself. If your run finished, these passed, and the
files are the evidence. If a check stopped the run, its report is in the error
message instead.

| file | route | what it catches |
|---|---|---|
| `qc/<sample>/<sample>.<engine>.search_qc.tsv` | PRM routes | the search identified nothing, or every $q$ value came out 1 |
| `qc/<engine>.quant_qc.tsv` | EncyclopeDIA routes | the quantification pass matched no fragment ions, usually because it ran at a different mass tolerance from the search |
| `openswath_feasibility/<sample>.openswath_feasibility.txt` | OpenSWATH PRM route | spectra too crowded for OpenSWATH to separate targets from decoys. See [PARAMETERS.md §10](PARAMETERS.md#10-acquisition-checks) |

## Other directories

| path | contents |
|---|---|
| `base_quant/<route>/base_quant.long.tsv` | the merged table before the reporting layer annotates it. `quant.full.tsv` is this plus status columns. Intermediate, not a result |
| `base_quant/<route>/base_quant_by_sample/` | the same, one file per sample |
| `diathem_quant/<mode>/` | diathem's own output |
| `encyclopedia_libexport/analysis/analysis.elib` | the EncyclopeDIA quantification database |
| `encyclopedia_libexport/analysis/analysis.elib.peptides.txt` | peptide × sample matrix. Written only when cross-run alignment succeeded |
| `search_dia/<sample>/`, `context_search/<sample>/`, `openswath_workflow/<sample>/`, `context_ms_run/<sample>/`, `pyprophet/`, `pyprophet_prm/<sample>/` | per-sample search and rescoring output, with logs |
| `prepare_library*/`, `complete_channels/`, `harmonize_rt/`, `merge_libraries/` | the libraries the pipeline built, and the checks on them |
| `msconvert/<sample>/` | converted mzML, if you supplied vendor files |
| `pipeline_info/` | Nextflow report, timeline, trace, DAG, log, and the work directory |

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
