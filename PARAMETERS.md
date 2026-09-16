# TEAPOT parameters

Every setting, its default, and when to change it.

Set them in a YAML params file, or as `--flags` on the command line:

```yaml
input:  my_samplesheet.csv
fasta:  my_proteins.fasta
fdr:    0.01
```

```bash
nextflow run /path/to/teapot/main.nf -profile apptainer -params-file my_run.yml
```

This page lists every setting you might want to change. `nextflow.config` holds
a few more that the pipeline needs internally; you can ignore them.

---

## 1. Required

| setting | default | |
|---|---|---|
| `input` | required | your samplesheet CSV |
| `fasta` | required | protein FASTA |
| `outdir` | `results` | where results go |

## 2. Run files

Vendor formats are converted for you. These two settings are passed to
msconvert.

| setting | default | |
|---|---|---|
| `peakpicking` | `--filter "peakPicking vendor msLevel=1-"` | centroids the spectra, using the instrument vendor's own algorithm |
| `otherfilters` | `''` | further msconvert flags, usually more `--filter` arguments |

## 3. The assay library

Supply exactly **one** source. Both engines accept all of these.

| setting | default | |
|---|---|---|
| `elib` | `null` | an EncyclopeDIA chromatogram library. Best quality |
| `dlib` | `null` | an EncyclopeDIA spectrum library |
| `blib` | `null` | a Skyline / BiblioSpec library |
| `library_sheet` | `null` | a CSV of GPF runs; the pipeline builds a chromatogram library from them |
| `use_koina` | `false` | predict a library from your `fasta` |
| `koina_models` | `Prosit_2020_HCD;IM2Deep_CCS;Prosit_2019_iRT` | prediction models |
| `koina_url` | `https://koina.wilhelmlab.org:443/` | |

OpenSWATH also reads three of its own formats directly:

| setting | default | |
|---|---|---|
| `openswath_pqp` | `null` | a ready `.pqp`, used as-is |
| `openswath_tsv` | `null` | a transition list |
| `openswath_traml` | `null` | a TraML |
| `openswath_blib` | `null` | equivalent to `blib` on the OpenSWATH route |

## 4. The background library

Only Context and Context-MS need one, only when your assay library contains
your reference peptides and nothing else. A vendor kit library is exactly
the case.

| setting | default | |
|---|---|---|
| `background_source` | `none` | `none` \| `koina` \| `dlib` \| `blib` \| `elib` \| `gpf` |
| `background_library` | `null` | the file, for `dlib` / `blib` / `elib` |
| `background_fasta` | `null` | required for `koina` |
| `background_min_targets` | `500` | the run stops if the number of background peptides is lower than this threshold |

## 5. Modifications

How a mass shift in your library is matched to a named UniMod modification.

| setting | default | |
|---|---|---|
| `mod_mass_tolerance` | `0.01` | mass tolerance in Da, for matching a mass shift to a UniMod modification |
| `unmapped_mod_max_fraction` | `0.05` | stop if more than this fraction cannot be matched |

## 6. Instrument settings

Both engines default to Orbitrap fragment tolerances. Adjust them to your
instrument resolution. For example:

| your instrument | EncyclopeDIA | OpenSWATH | diathem |
|---|---|---|---|
| ion trap (Thermo Stellar, LTQ) | `-ftol 0.4 -ftolunits AMU -lftol 0.4 -lftolunits AMU -frag HCD` | `-mz_extraction_window 0.4 -mz_extraction_window_unit Th` | `diathem_fragment_tol_th: 0.4` |
| TOF (SCIEX TripleTOF) | `-ftol 20 -ftolunits PPM -frag CID` | `-mz_extraction_window 50 -mz_extraction_window_unit ppm -mz_extraction_window_ms1 20 -mz_extraction_window_ms1_unit ppm` | `diathem_fragment_ppm: 20` |
| Orbitrap | the default | the default | the default |

| setting | default | |
|---|---|---|
| `dia_search_args` | `''` | EncyclopeDIA search, `mode: DIA` |
| `context_search_args` | `''` | EncyclopeDIA search, `mode: PRM` |
| `gpf_search_args` | `''` | the GPF library build by EncyclopeDIA DIA search |
| `libexport_args` | `''` | EncyclopeDIA quantification. Instrument flags are inherited from the search automatically |
| `openswath_search_args` | `''` | extra [OpenSwathWorkflow](https://openms.de/documentation/html/TOPP_OpenSwathWorkflow.html) flags |
| `openswath_assaygen_args` | `''` | extra [OpenSwathAssayGenerator](https://openms.de/documentation/html/TOPP_OpenSwathAssayGenerator.html) flags |
| `openswath_decoygen_args` | `''` | extra [OpenSwathDecoyGenerator](https://openms.de/documentation/html/TOPP_OpenSwathDecoyGenerator.html) flags |
| `java_mem` | `32g` | |

## 7. Retention time

Libraries from different sources arrive on different retention time (RT) scales.
The pipeline puts them on one scale before use, and needs anchor peptides to do it.
See [INPUTS.md](INPUTS.md#rt-anchors) for where the anchors come from.

| setting | default | |
|---|---|---|
| `irt_library` | `null` | your own iRT normalization library (TraML, TSV or PQP) |
| `irt_anchors` | `null` | your own anchor peptides. **Required together with `irt_library`** |
| `irt_top_transitions` | `6` | transitions kept per anchor peptide |
| `rt_harmonize_fallback` | `fail` | `fail` \| `full_window`. `full_window` searches the entire gradient instead of aligning: correct, but much slower |
| `rt_harmonize_min_anchors` | `5` | fewest anchors accepted for a fit |
| `rt_harmonize_min_r2` | `0.98` | lowest fit quality accepted |
| `rt_harmonize_min_anchor_coverage` | `0.8` | warns if the anchors span less of the library than this |
| `rt_coherence_max_range_ratio` | `3.0` | check that the two halves of a merged library did align |
| `rt_coherence_require_iqr_overlap` | `true` | the same check |

## 8. Heavy/light channels and absolute quantification

Leave `heavy_label` empty for label-free runs.

| setting | default | |
|---|---|---|
| `heavy_label` | `[]` | e.g. `[{site: K, unimod: 'UniMod:259'}, {site: R, unimod: 'UniMod:267'}]`. `site` is a residue, or `n_term` / `c_term` |
| `standard_amounts` | `null` | CSV of known pmol of SIL internal standards |
| `min_consistency` | `0.0` | drop measurements whose diathem consistency is below this before computing amounts |
| `complete_channels_skip_proteins` | `[]` | accessions of proteins that exist only as heavy, e.g. `['iRT_Tag']` |
| `complete_channels_ion_types` | `['b','y']` | add `c`, `z` for ETD / EThcD |
| `complete_channels_max_frag_charge` | `2` | |
| `complete_channels_tolerance_ppm` | `20.0` | |
| `complete_channels_coverage_hard_floor` | `0.5` | stop if fewer than this fraction of peptides got a usable twin |
| `complete_channels_coverage_warn_below` | `0.9` | warn below this |

## 9. Confidence estimation

| setting | default | |
|---|---|---|
| `fdr` | `0.01` | the $q$ value threshold |
| `openswath_prm_fdr` | `context` | which confidence estimator the OpenSWATH PRM route uses: `context` (needs a background) or `pyprophet` (no background needed) |
| `add_decoys_to_reference_list` | `true` | if your reference list has no decoys the pipeline generates them automatically |
| `search_min_reference_ids` | `10` | warn if fewer reference peptides pass. Zero passing always stops the run |
| `pyprophet_level` | `ms2` | `ms2` \| `ms1ms2` |
| `pyprophet_learn_args` | XGBoost classifier | add `--pi0_lambda 0 0 0` if the run stops with `pi0 <= 0` |
| `pyprophet_apply_args` | `--classifier=XGBoost` | |
| `pyprophet_export_args` | `''` | |
| `context_engine` | `mprophet` | `percolator` \| `mprophet` |
| `context_input_profile` | `auto` | |
| `context_seed_coefficients` | `encyclopedia` | |
| `context_args` | `''` | extra Context-MS flags |

## 10. Acquisition checks

Scheduled PRM can interleave wide survey scans with the narrow targeted ones.
The pipeline automatically detects this and extracts from the targeted windows
only.

On the OpenSWATH PRM route it also measures how crowded your spectra are, and
stops the run if they are too crowded for OpenSWATH to separate targets from
decoys. The settings below control that check:

| setting | default | |
|---|---|---|
| `openswath_occupancy_tolerance` | `0.4` | in **Th**. Match it to your extraction window |
| `openswath_max_occupancy` | `0.35` | stop above this |
| `openswath_occupancy_window_width` | `null` | which isolation width to measure. `null` picks the narrowest, which is normally right |
| `openswath_feasibility_warn_only` | `false` | report and continue instead of stopping |
| `openswath_auto_matching_window` | `true` | let the pipeline handle interleaved survey windows. Setting `-matching_window_only` yourself in `openswath_search_args` overrides it, either value |
| `skip_openswath_feasibility` | `false` | skip the crowding check and the window detection |

## 11. Quantification

| setting | default | |
|---|---|---|
| `primary_abundance` | `diathem` | which tool fills the headline `abundance_primary` column: `openswath` \| `encyclopedia` \| `diathem` \| `tric` |
| `run_diathem` | `true` | diathem interference-corrected quantification |
| `run_tric` | `false` | cross-run feature alignment. OpenSWATH DIA route only |
| `tric_args` | `''` | |
| `quant_min_runs` | `2` | warn when a route has fewer runs than this |
| `quant_require_min_runs` | `false` | make that warning stop the run. Set both for production |
| `quant_max_zero_ion_fraction` | `0.5` | warn when more than this fraction of peptides matched no fragment ion. All-zero always stops the run |
| `augment_context_with_diathem` | `false` | feed diathem's statistics to Context-MS as extra rescoring features. OpenSWATH PRM route only |
| `openswath_prm_encyclopedia_quant` | `false` | also fill `abundance_encyclopedia` on the OpenSWATH PRM route, for comparing tools. Works with any library source except `openswath_pqp` / `openswath_tsv` / `openswath_traml`, which EncyclopeDIA cannot read |

### diathem

Left at `null`, the two parameters default based on the
samplesheet's `mode`:

* PRM: `diathem_precursor_tol: 0.5`, `diathem_rt_window: 60`
* DIA: `diathem_precursor_tol: 12.5`, `diathem_rt_window: 0`

Override them when the mode does not match the acquisition.

| setting | default | |
|---|---|---|
| `diathem_precursor_tol` | `null` | isolation window half-width in Da |
| `diathem_rt_window` | `null` | RT window half-width in seconds. `0` uses the whole run |
| `diathem_fragment_tol_th` | `null` | fragment tolerance in Th. Suits ion traps; overrides `diathem_fragment_ppm` |
| `diathem_fragment_ppm` | `null` | fragment tolerance in ppm. diathem default: 20, suits Orbitrap |
| `diathem_library_ppm` | `null` | ppm tolerance for matching fragments to the library's peaks |
| `diathem_library` | `null` | library source for expected fragment ratios. Defaults to the assay library |
| `diathem_targets_source` | `reference` | `reference` \| `library`. Quantify the `reference` list, or the whole `library`. OpenSWATH PRM only |
| `diathem_args` | `''` | extra [diathem](https://github.com/statisticalbiotechnology/diathem/pkgs/container/diathem) flags |

## 12. Containers

Every tool runs in a container. All of them are pulled for you; nothing is
built by hand.

| setting | default | |
|---|---|---|
| `encyclopedia_container` | `ghcr.io/shannon225/encyclopedia:6.8.29-teapot` | EncyclopeDIA |
| `openswath_container` | `ghcr.io/openswath/openswath@sha256:509feb438a20252585ff77e05c08a2cd41227ea290e07a959f842a144b6703b8` | OpenMS |
| `pyprophet_container` | `ghcr.io/pyprophet/pyprophet@sha256:fb473fe3222305a94ffd7f9b2900d2d863c182d44108e02daff445e0b7ad08b8` | PyProphet |
| `context_container` | `ghcr.io/shannon225/context-ms:main` | Context-MS |
| `diathem_container` | `ghcr.io/statisticalbiotechnology/diathem:quant` | diathem |
| `python_container` | `quay.io/biocontainers/pandas:2.2.1` | Python |
| `msconvert_container` | `proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses` | Proteowizard msconvert |
| `container_cache_dir` | `<pipeline>/apptainer_cache` | where images are stored. Kept outside `outdir`, so they are pulled once and reused by later runs |
| `container_tmp_dir` | `<pipeline>/apptainer_tmp` | scratch for unpacking images |

## 13. Execution

| setting | default | |
|---|---|---|
| `max_cpus` | detected cores | caps a single-machine run |
| `cluster_options` | `''` | options added to every scheduler job, e.g. `-A <account>` |
| `work_dir` | under `outdir` | where intermediate files go. They can be deleted once the run has finished |
| `publish_mode` | `copy` | `copy` \| `symlink`. How results are placed in `outdir`. `symlink` saves space, but the links break once `work_dir` is deleted |

## 14. CPU, memory and time are not parameters

They are Nextflow *process directives*, set in `conf/base.config` (default) and
`conf/slurm.config` (loaded automatically by the SLURM launcher). If a step runs
out of CPUs, memory or time, raise it there.
