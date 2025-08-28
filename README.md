<p align="center">
    <img src="assets/readme/logo.jpeg", alt="log.jpeg", height=200, width=200>
</p>

# CRISPR-HAWK

CRISPR-HAWK is a comprehensive and scalable tool for designing guide RNAs (gRNAs) in CRISPR-Cas systems. Available as an offline tool with a user-friendly command-line interface, CRISPR-HAWK integrates large-scale human genetic variation datasets—such as the 1000 Genomes Project, the Human Genome Diversity Project (HGDP), and gnomAD—with orthogonal genomic annotations to systematically prioritize gRNAs targeting regions of interest. CRISPR-HAWK is Cas system-independent, supporting a wide range of nucleases including Cas9, SaCas9, Cpf1 (Cas12a), and others. It offers users full flexibility to define custom PAM sequences and guide lengths, enabling compatibility with emerging CRISPR technologies and tailored experimental requirements. The tool accounts for both single-nucleotide variants (SNVs) and small insertions and deletions (indels), and it is capable of handling individual- and population-specific haplotypes. This makes CRISPR-HAWK particularly suitable for both personalized and population-wide gRNA design. CRISPR-HAWK automates the entire workflow—from variant-aware preprocessing to gRNA discovery—delivering comprehensive outputs including ranked tables, annotated sequences, and high-quality figures. Its modular design ensures easy integration with existing pipelines and tools, such as [CRISPRme](https://github.com/pinellolab/CRISPRme) or [CRISPRitz](https://github.com/pinellolab/CRISPRitz), for subsequent off-target prediction and analysis of prioritized gRNAs.

## Table of Contents

0 [System Requirements](#0-system-requirements)
<br>1 [Installation](#1-installation)
<br>&nbsp;&nbsp;1.1 [Install CRISPR-HAWK from Mamba/Conda](#11-install-crispr-hawk-from-mambaconda)
<br>&nbsp;&nbsp;1.2 [Install CRISPR-HAWK from Docker](#12-install-crispr-hawk-from-docker)
<br>&nbsp;&nbsp;1.3 [Install CRISPR-HAWK from PyPI](#13-install-crispr-hawk-from-pypi)
<br>&nbsp;&nbsp;1.4 [Install CRISPR-HAWK from Source Code](#14-install-crispr-hawk-from-source-code)

## 0 System Requirements

To ensure optimal performance, CRISPR-HAWK requires the following system specifications:

- **Operating System**:
<br>macOS or any modern Linux distribution (e.g., Ubuntu, CentOS)

- **Minimum RAM**:
<br>16 GB — sufficient for standard use cases and small to medium-sized datasets

- **Recommended RAM for Large-Scale Analyses**:
<br>32 GB or more — recommended for memory-intensive tasks such as:

    - Scanning regions larger than 1 Mb

    - Processing large-scale variant datasets (e.g., gnomAD data)

⚠️ Note: For optimal performance and stability, especially when dealing with large-scale variant datasets, ensure that your system meets or exceeds the recommended specifications.

## 1 Installation

This section provides step-by-step instructions to install CRISPR-HAWK. Choose the method that best suits your environment and preferences:

- **[Install CRISPR-HAWK from Mamba/conda](#11-install-crispr-hawk-from-mambaconda)** (recommended)
<br>Best for users seeking an isolated and reproducible environment with minimal manual dependency handling.

- **[Install CRISPR-HAWK from Docker](#12-install-crispr-hawk-from-docker)**
<br>Ideal for users who prefer containerized deployments or want to avoid configuring the environment manually.

- **[Install CRISPR-HAWK from PyPI](#13-install-crispr-hawk-from-pypi)**
<br>Quick option for Python users already working within a virtual environment. May require manual handling of some dependencies.

- **[Install CRISPR-HAWK from source code](#14-install-crispr-hawk-from-source-code)**
<br>Suitable for developers or contributors who want full control over the codebase or plan to customize CRISPR-HAWK.

⚠️ **Note:** We recommend using the Mamba/Conda or Docker installation methods for most users, as they ensure the highest compatibility and stability across systems.

### 1.1 Install CRISPR-HAWK from Mamba/Conda

TBA

### 1.2 Install CRISPR-HAWK from Docker

TBA

### 1.3 Install CRISPR-HAWK from PyPI

TBA

### 1.4 Install CRISPR-HAWK from Source Code

Installing CRISPR-HAWK from source is ideal for developers, contributors, or users who wish to inspect or customize the codebase.

This method assumes you already have **Python 3.8** installed and accessible from your system’s environment.

✅ **Prerequisites**

- Python **3.8** (strictly required)

- `git`

- A virtual environment (optional but recommended)

🧱 **Installation Steps**

**1. Clone the Repository**
```bash
git clone https://github.com/ManuelTgn/CRISPR-HAWK.git
cd CRISPR-HAWK
```

**2. (Optional) Create and Activate a Virtual Environment**
```bash
python3.8 -m venv crisprhawk-env
source crisprhawk-env/bin/activate
```

**3. Install CRISPR-HAWK and Its Dependencies**
```bash
pip install .  # regular installation
pip install -e .  # development-mode installation
```

The `.` tells `pip` to install the current directory as a package, including all dependencies specified in `setup.py` or `pyproject.toml`.

🚀 **Quick Test**

Once installation is complete, verify that the command-line interface is working:
```bash
crisprhawk -h
```

If the help message is displayed correctly, CRISPR-HAWK is successfully installed and callable from any directory in your system.

## 2 Usage

This section describes CRISPR-HAWK's fubctionalities

### 2.1 Search

The `crisprhawk search` command is the core functionality of CRISPR-HAWK, designed to identify and annotate candidate gRNAs in both reference and variant genomes.
It integrates variant-aware search, functional annotation, and predictive scoring to help you prioritize the most robust and context-aware guides for CRISPR editing.

The search includes:

* Support for any Cas system (Cas9, Cpf1, SaCas9, etc.)
* Compatibility with custom PAM sequences and guide lengths
* Consideration of individual or population-level variants (SNVs/indels) from VCF files
* Scoring using **Azimuth**, **RS3**, **CFDon**, and **DeepCpf1**
* Functional and gene annotation using user-specified BED files <!--* Optional estimation and reporting of **off-targets**-->
* Output in detailed and structured reports (TSV, haplotype tables)

Usage:
```bash
crisprhawk search -f <fasta> -r <bedfile> -v <vcf> -p <pam> -g <guide-length> -o <output-dir>
```

---

Arguments and Options:

#### Required

| Option                       | Description                                                                                                                                    |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `-f`, `--fasta <FASTA-FILE>` | Path to the **reference genome** in FASTA format. Guides will be searched on this genome.                                                      |
| `-r`, `--regions <BED-FILE>` | BED file specifying the **target genomic regions** where gRNAs will be searched (e.g., promoters, exons, or custom loci).                      |
| `-v`, `--vcf <VCF-DIR>`      | *(Optional but recommended)* Folder containing **VCF files** for the variant-aware design. If omitted, only the reference genome will be used. |
| `-p`, `--pam <PAM>`          | PAM sequence used to define valid gRNA targets (e.g., `NGG` for SpCas9, `TTTV` for Cpf1).                                                      |
| `-g`, `--guide-len <LENGTH>` | Length of the guide RNA, excluding the PAM (typically 20 for SpCas9).                                                                          |
| `-o`, `--outdir <DIR>`       | Output directory to store results, reports, and intermediate files. If not specified, defaults to the current working directory.               |

#### Optional

| Option                                         | Description                                                                                                                                                                       |
| ---------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-i`, `--fasta-idx <FAI>`                      | Optional FASTA index file (`.fai`). If not provided, it will be computed automatically.                                                                                           |
| `--right`                                      | By default, guides are extracted **upstream** of the PAM. Use this flag to extract them **downstream** (right side), useful for Cpf1 or other reverse-PAM systems.                |
| `--no-filter`                                  | By default, only VCF variants with `FILTER == PASS` are considered. Use this flag to include **all variants**, regardless of filter status.                                       |
| `--annotation <BED1 BED2 ...>`                 | Provide one or more **BED files** with custom genomic features (e.g., enhancers, DHS, regulatory elements). Must follow BED format and have at least 4 columns (name in the 4th). |
| `--annotation-colnames <name1 name2 ...>`      | Custom column names for the annotations from the `--annotation` files. Must match the number and order of the provided BEDs.                                                      |
| `--gene-annotation <GENE-BED1 GENE-BED2 ...>`  | One or more **gene annotation BED files** (9-column format). Must follow GENCODE-style structure, with gene name in the 9th column and gene feature (e.g., exon, UTR) in the 7th. |
| `--gene-annotation-colnames <name1 name2 ...>` | Custom column names for gene annotations, matching the files in `--gene-annotation`.                                                                                              |
| `--haplotype-table`                            | When enabled, a TSV file reporting haplotype-aware variants and their associated guide matches will be produced.                                                                  |
| `--estimate-offtargets`                        | Enables **off-target site estimation** for each candidate guide (uses efficient genome-wide search).                                                                              |
| `--write-offtargets-report`                    | When enabled, generates a detailed **off-target report** for each guide. Disabled by default for performance.                                                                     |
| `-t`, `--threads <INT>`                        | Number of threads to use for parallel processing. Use `-t 0` to utilize **all available cores**. *(Default: 1)*                                                                   |
| `--verbosity <LEVEL>`                          | Controls the verbosity of logs. Options: `0` = Silent, `1` = Normal, `2` = Verbose, `3` = Debug. *(Default: 1)*                                                                   |
| `--debug`                                      | Enables **debug mode**, showing full stack traces and internal logs for troubleshooting.                                                                                          |

---

Example:

```bash
crisprhawk search \
  -f hg38.fa \
  -r targets.bed \
  -v gnomAD_chr22/ \
  -p NGG \
  -g 20 \
  -o results/ \
  --annotation regulatory_regions.bed \
  --gene-annotation gencode.bed \
  --haplotype-table \
  -t 8
```

This will:

* Search 20 bp gRNAs with an NGG PAM in `targets.bed`
* Consider population variants from VCFs in `gnomAD_chr22/`
* Annotate guides using custom regulatory and gene annotations <!--* Estimate and report off-targets-->
* Run in parallel using 8 threads






