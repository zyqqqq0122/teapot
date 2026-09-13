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
| `sample_id` | your name for the run |
| `file` | one acquisition. See run files below |
| `mode` | `PRM` or `DIA` |
| `reference_list` | the reference peptides for this run. Can be blank |
| `engine` | `encyclopedia` or `openswath`. Defaults to `encyclopedia` |

## Run files: the `file` column

Either an open format (`.mzML`, `.mzXML`) or a vendor format the pipeline
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
| `--library_sheet` | the GPF runs themselves, listed in a CSV of `condition` and `fraction_file` ([template](templates/library_sheet.csv)). The pipeline searches them to build an `.elib` |
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

## The background: `--background_library` or `--background_fasta`

Needed only by Context and Context-MS, and only when your assay library contains
your reference peptides and nothing else. A vendor kit library is exactly that
case. Supply either:

| source | what it is |
|---|---|
| `--background_library` | an `.elib`, `.blib` or `.dlib`, matching `background_source` |
| `--background_fasta` | a FASTA of the proteins potentially present in your samples |

## The FASTA: `--fasta`

Used to assign peptides to proteins, to generate decoys, and as the source for
Koina prediction. On a label-free or DIA run, supply the proteins potentially
present in your samples. If your SIL internal standards came as a kit, supply
the vendor's own FASTA; see
[below](#not-everything-in-your-library-and-your-fasta-is-an-analyte).

## The reference list: `reference_list`

Your reference peptides: what you told the mass spectrometer to look for. The
pipeline reads it to know what *should* have been found, and the list you gave
the instrument usually works unedited:

```csv
Compound,m/z,z,RT Time (min)
LGGNEQVTR,487.2567,2,23.41
GAGSSEPVTGLDAK,644.8232,2,28.73
```

CSV or TSV, detected from the header. Column names are matched
case-insensitively from a set of aliases:

| meaning | accepted names |
|---|---|
| peptide | `Compound`, `Sequence`, `Peptide`, `PeptideSequence`, `ModifiedPeptideSequence` |
| RT | `RT`, `RT Time (min)`, `Retention Time`, `NormalizedRetentionTime` |
| charge | `z`, `Charge`, `PrecursorCharge` |
| *m/z* | `m/z`, `mz`, `PrecursorMz` |
| decoy flag | `isDecoy`, `decoy`, `is_decoy`, with `true`, `1`, `yes` or `decoy` meaning a decoy |

An RT column is optional. On the OpenSWATH PRM route it enables the pipeline to
derive RT anchors itself, without extra files, unless you supply a vendor `.blib`,
which uses its own `IrtLibrary` instead.

The decoy column is also optional. $q$ value and posterior error probability (PEP)
estimates come from target-decoy competition. If your list has no decoys, the
pipeline detects that and generates them for you. OpenSWATH takes its decoys from
the `.pqp` instead.

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

A vendor SIL kit supplies artificial, heavy-labelled standard proteins covering
the peptides you want to quantify. The standard proteins share one extra stretch
of sequence, an **iRT tag**, whose peptides give the run its RT scale in place of
a separate iRT kit.

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

Name the iRT tag and the pipeline skips its peptides:

```yaml
complete_channels_skip_proteins:
  - 'iRT_Tag'          # its accession in the FASTA
```
