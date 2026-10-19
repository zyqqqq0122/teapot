# TEAPOT input files

What each file you supply actually contains, and what the pipeline does with it.

---

## The samplesheet: `--input`

One row per run. CSV with a header:

```csv
sample_id,file,mode,reference_list,engine
sample_A,/data/sample_A.raw,PRM,/data/reference.csv,encyclopedia
sample_B,/data/sample_B.raw,PRM,/data/reference.csv,encyclopedia
```

| column | what it is |
|---|---|
| `sample_id` | your name for the run, a different one on every row. Letters, digits and `.` `_` `+` `-` only |
| `file` | one acquisition. See run files below |
| `mode` | `PRM` or `DIA` |
| `reference_list` | your reference peptides: the same file on every row of an engine and mode, or blank on all of them. See the reference list below |
| `engine` | `encyclopedia` or `openswath`. Defaults to `encyclopedia` |

## Run files: the `file` column

Either an open format (`.mzML`; `.mzXML` with OpenSWATH only) or a vendor format the pipeline
converts for you inside a container: `.raw` (Thermo), `.wiff` (SCIEX), `.d`
(Bruker/Agilent), `.lcd` (Shimadzu). No local ProteoWizard install.

A `.wiff` is not one file. Name only the `.wiff` in the samplesheet; its
`.wiff.scan` (and `.wiff2.scan`, `.timeseries.data` where present) must sit in
the same directory and are staged automatically. A `.wiff` holding several
samples fans out into several mzML.

## The assay library

Supply exactly one library source:

| you supply | what it is |
|---|---|
| `--blib` or `--dlib` | a spectrum library: peptides with their fragment spectra |
| `--elib` | a chromatogram library, built from GPF runs, so it carries measured retention times (RT) |
| `--library_sheet` | the GPF runs themselves, listed in a CSV of `condition` and `fraction_file` ([template](templates/library_sheet.csv)), as mzML or in a vendor format the pipeline converts. The pipeline searches them, with `--dlib` or with a library predicted from your FASTA (`--use_koina`), to build an `.elib`. Files sharing a `condition` are merged into one run first, as the fractions of one sample |
| `--use_koina` | a switch: the pipeline predicts a spectrum library from your FASTA |

Those all work on either engine. OpenSWATH additionally accepts three formats
of its own:

| you supply | what it is |
|---|---|
| `--openswath_tsv` | an OpenSWATH transition list: one row per fragment ion |
| `--openswath_traml` | a TraML |
| `--openswath_pqp` | OpenSWATH's own library format, with transitions and decoys already in it |

A vendor kit library is a special case, covered
[below](#not-everything-in-your-library-and-your-fasta-is-an-analyte).

### RT anchors

Whichever library source you use, the pipeline puts the library onto a single
RT scale, and needs anchors to do so. They come from the first of these
available:

| anchors from | you supply | when |
|---|---|---|
| `--irt_library` with `--irt_anchors` | two files: the library as a PQP, TraML or transition TSV, and the anchors as the TSV below | whenever you supply them, and they win over the two below |
| the library's own `IrtLibrary` table | nothing | your library is a vendor kit `.blib` |
| the RT column of your `reference_list` | nothing | your library is `--dlib`, `--elib`, `--use_koina` or `--library_sheet` |

An anchors file is tab-separated, one row per anchor peptide:

```tsv
PeptideModSeq	Irt	Standard
LGGNEQVTR	-24.92	1
GAGSSEPVTGLDAK	0.00	1
```

## The background: `background_source`

Needed only by Context and Context-MS, and only when your assay library contains
your reference peptides and nothing else. A vendor kit library is exactly that
case. The background is a library: one you supply, one predicted from the proteins
of your samples, or one built from a few of your runs; [PARAMETERS.md
§4](PARAMETERS.md#4-the-background-library) explains how to choose. What you supply:

| you supply | what it is |
|---|---|
| `background_library` | a background library you already have: an `.elib`, `.dlib` or `.blib` |
| `background_fasta` | a protein FASTA of the proteins potentially present in your samples. With `background_peptides`, your species' whole proteome works as well. Not the vendor kit's FASTA, which holds only the standards |
| `background_peptides` | the peptides to keep of a predicted background: the tables of a PeptideAtlas build, or a list of peptides. See below |
| `background_sheet` | a few runs of samples like yours to build the background from, listed as in `library_sheet` ([template](templates/background_sheet.csv)), each under a `condition` of its own. Your PRM runs can be used as they are, as mzML or in a vendor format the pipeline converts |
| `background_search_library` | a spectrum library to search the `background_sheet` runs with, instead of a prediction from `background_fasta`: an `.elib`, `.dlib` or `.blib` |

### Peptides to keep: `background_peptides`

Either of two files. The pipeline predicts the peptides of the proteins in
`background_fasta`, then keeps only those this file names, so it always comes
together with `background_fasta`.

**The tables of a PeptideAtlas build.** [PeptideAtlas](https://peptideatlas.org/builds/)
counts how often each peptide was observed in the public data of one sample type.
Give the build file as it is; unpacking it is not needed. The pipeline reads `peptide.tsv`
and `peptide_instance.tsv`, and keeps the `background_peptides_top` peptides
observed most often.

**A list of peptides**, any you choose: a text file with one peptide per line (a
header line such as `peptide` is skipped; in a CSV or TSV the first column is read):

```
peptide
LGGNEQVTR
GAGSSEPVTGLDAK
```

or a FASTA of peptides:

```
>pep1
LGGNEQVTR
>pep2
GAGSSEPVTGLDAK
```

Modifications and flanking residues are ignored. Every listed peptide the
prediction holds is kept.

## The FASTA: `--fasta`

Used to assign peptides to proteins, to generate decoys, and as the source for
Koina prediction. A background has its own FASTA, `background_fasta`. On a
label-free or DIA run, supply the proteins potentially
present in your samples. If your SIL internal standards came as a kit, supply
the vendor's own FASTA; see
[below](#not-everything-in-your-library-and-your-fasta-is-an-analyte).

## The reference list: `reference_list`

Your reference peptides: what you told the mass spectrometer to look for. The
pipeline reads it to know what *should* have been found, and the list you gave
the instrument usually works unedited (A CSV or TSV whose first seven columns
are these, in this order):

```csv
Compound,Formula,Adduct,m/z,z,RT Time (min),Window (min)
LGGNEQVTR,,(no adduct),487.2567,2,23.41,5
GAGSSEPVTGLDAK,,(no adduct),644.8232,2,28.73,5
```

`Formula` and `Adduct` may be empty. An eighth column, `isDecoy` (`true` or `false`),
is optional: without decoys the pipeline makes its own. Further columns are ignored.

Name the reference list in the samplesheet, or leave the column blank:

| left blank | then |
|---|---|
| your library is a vendor kit `.blib` | the list is taken from the library |
| PRM runs with any other library | the pipeline cannot start |
| DIA runs with any other library | the results list the peptides your runs detected |

## Known amounts: `--standard_amounts`

Needed only for absolute quantification. A CSV keyed either **per protein**:

```csv
Uniprot,Gene,Amount per well [pmol]
P04217,A1BG,1.327
P01023,A2M,1.971
```

or **per peptide**:

```csv
Sequence,Charge,pmol
LGGNEQVTR,2,0.500
GAGSSEPVTGLDAK,2,0.500
```

Column names are matched case-insensitively from a set of aliases, so a vendor's
own certificate of analysis usually works unedited. Extra columns are ignored:

| meaning | accepted names |
|---|---|
| amount | `pmol`, `amount`, `heavy_pmol`, `known_pmol`, `Amount per well [pmol]` |
| peptide | `Sequence`, `Peptide`, `PeptideSequence`, `Compound`, `Peptidoform` |
| protein | `Uniprot`, `Protein`, `Protein_ID`, `Accession` |
| charge | `Charge`, `z`, `PrecursorCharge` (optional) |

The amount is always on the **heavy** side:
`abs_pmol_light = (light / heavy) × known heavy pmol`.

A protein-keyed file is matched by accession token, and by the part before the
first `_`, so a library accession like `QR3140513_A1BG` matches `QR3140513` in
the amounts file.

---

## Not everything in your library and your FASTA is an analyte

A vendor SIL kit supplies heavy-labelled standards: recombinant fragments of
the proteins you want to quantify, each covering a few of their peptides. Every
standard carries the same tag, fused to it, which the kit calls the iRT tag:
its peptides give the run its RT scale in place of a separate iRT kit.

The library, FASTA and amounts file that come with the kit describe those
standard proteins:

| the kit's file | what it holds |
|---|---|
| assay library | your SIL internal standards and the iRT tag's peptides, **no background** |
| FASTA | one entry per artificial protein, iRT tag included |
| known amounts | the pmol of each artificial protein |

Your reference peptides exist in both channels: heavy (SIL internal standards)
and light (endogenous peptides). When the library only includes one channel,
the pipeline synthesizes the missing twin so both are analyzed. The iRT tag's
peptides exist **only as heavy** since the iRT tag is a synthetic sequence
with no endogenous counterpart. Synthesizing a light twin would force the
pipeline to extract a precursor the instrument never isolated, merely
collecting whatever randomly co-elutes in its window.

Name the iRT tag and the pipeline skips its peptides: it gives them no light twin,
and leaves them out of the results (only `quant.full.tsv` lists them):

```yaml
complete_channels_skip_proteins:
  - 'iRT_Tag'          # its accession in the FASTA
```
