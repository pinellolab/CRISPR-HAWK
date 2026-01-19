# Instructions to reproduce the results presented in "CRISPR-HAWK: Haplotype- and Variant-Aware Guide Design Toolkit for CRISPR-Cas"

This folder provides a fully reproducible workflow for regenerating the analyses and figures presented in the manuscript. All required input datasets, pre-computed CRISPR-HAWK outputs, and supporting scripts are included or referenced here.

The workflow illustrates how haplotype- and variant-aware CRISPR guide design changes predicted on-target and off-target activity across human populations. In particular, the analyses demonstrate how genetic diversity, both single variants and phased haplotypes, can modify target-site sequences and, in turn, influence CRISPR activity predictions.

## Downloading the Data

All datasets and pre-computed outputs needed to reproduce the manuscript results are archived on Zenodo:

[https://zenodo.org/records/18070463](https://zenodo.org/records/18070463)

Download the archive and extract it directly into this directory. After unpacking, the folder structure will appear as:

```
crisprhawk-data
├── data
│   ├── annotations
│   ├── genomes
│   ├── regions
│   └── variants
├── results
└── results_nsamples
```

### Folder Contents

#### `data/annotations/`

Genomic feature annotations used to characterize candidate guide locations, including:

* **COSMIC**: catalog of somatic cancer mutations
* **DHS**: DNase I hypersensitive sites
* **ENCODE**: experimentally annotated regulatory elements
* **GENCODE**: protein-coding gene annotations

These files are provided in indexed BED format to enable efficient querying.

#### `data/genomes/hg38/`

Reference genome FASTA files and index files (`.fa` + `.fai`) for the hg38 chromosomes analyzed in the study.

#### `data/regions/`

BED files defining the genomic loci supplied as input to CRISPR-HAWK:

* `cas9.bed` — regions targeted by **SpCas9** (NGG PAM; 20-nt guides)
* `cpf1.bed` — regions targeted by **Cas12a/Cpf1** (TTTV PAM; 23-nt guides)

These include both **therapeutic CRISPR targets** and **benchmark regions** referenced in the manuscript.

#### `data/variants/`

Population-scale human genetic variation datasets:

* **1000 Genomes Project (1000G)**
* **HGDP**
* **gnomAD v4.1 genomes**
  (pre-processed using `crisprhawk convert-gnomad-vcfs`)

CRISPR-HAWK uses these datasets to generate variant- and haplotype-resolved target sequences.

#### `results/`

Pre-computed CRISPR-HAWK outputs for all combinations of:

* Cas nucleases
* variant datasets
* genomic regions

They can be used directly to regenerate the manuscript figures.

## Re-Running the Full Analysis

To recompute all results from scratch:

1. Ensure **CRISPR-HAWK** is installed and available in your environment.
2. Run within `crisprhawk-data`:

```bash
bash reproduce-results.sh
```

where `crisprhawk-data` is the folder downloaded from Zenodo.

This script automatically:

* launches all CRISPR-HAWK analyses
* for **SpCas9 and Cas12a/Cpf1**
* across **all therapeutic and benchmark regions**
* using **all population variant datasets**

Outputs are organized by nuclease and dataset for convenience.

## Computing the Number of Supporting Individuals

Several manuscript figures report the **number of individuals** carrying each variant- or haplotype-matched CRISPR target.

To reproduce these statistics:

```bash
python compute_samples_number.py
```

Run this script **inside the `crisprhawk-data` directory**.

This script:

* parses CRISPR-HAWK haplotype tables
* determines how many individuals support each alternative gRNAs
* appends this information directly into the result tables

No manual intervention is required.

## Generating the Figures

Once CRISPR-HAWK analyses are complete and supporting-sample counts have been added you can regenerate all manuscript visualizations using the provided reproducibility notebook (`generate-figures.ipynb`).

The notebook:

* loads the annotated CRISPR-HAWK outputs
* reproduces each figure panel step-by-step
* documents the corresponding data-processing logic
