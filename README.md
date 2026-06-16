# TEAPOT

**T**argeted **E**valuation & **A**nalysis **P**ipeline for **P**roteomics **O**n **T**andem mass spectrometry

---

## Overview

TEAPOT is a pipeline for targeted evaluation and analysis of proteomics data derived from tandem mass spectrometry (MS/MS) experiments. It is designed to streamline the processing, quality assessment, and interpretation of targeted proteomics workflows such as SRM/MRM, PRM, and DIA.

## Features

- Targeted peptide/protein quantification from MS/MS data
- Quality control and evaluation metrics
- Flexible input support for common data formats
- Reproducible analysis pipeline

## Requirements

- Python 3.8+

## Tools to consider 
- OpenSWATH 
- EncyclopeDIA 
- Skyline

## Installation

```bash
git clone https://github.com/thanadol-git/teapot.git
cd teapot
pip install -r requirements.txt
```

## Usage

```bash
python teapot.py --input <input_file> --output <output_dir>
```

## License

MIT
