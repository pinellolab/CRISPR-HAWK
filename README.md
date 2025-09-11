<p align="center">
    <img src="assets/readme/logo.jpeg", alt="log.jpeg", height=200, width=200>
</p>

# CRISPR-HAWK

CRISPR-HAWK is a comprehensive and scalable tool for designing guide RNAs (gRNAs) and assessing genetic variants impact on on-target sites in CRISPR-Cas systems. Available as an offline tool with a user-friendly command-line interface, CRISPR-HAWK integrates large-scale human genetic variation datasets—such as the 1000 Genomes Project, the Human Genome Diversity Project (HGDP), and gnomAD—with orthogonal genomic annotations to systematically prioritize gRNAs targeting regions of interest. CRISPR-HAWK is Cas system-independent, supporting a wide range of nucleases including Cas9, SaCas9, Cpf1 (Cas12a), and others. It offers users full flexibility to define custom PAM sequences and guide lengths, enabling compatibility with emerging CRISPR technologies and tailored experimental requirements. The tool accounts for both single-nucleotide variants (SNVs) and small insertions and deletions (indels), and it is capable of handling individual- and population-specific haplotypes. This makes CRISPR-HAWK particularly suitable for both personalized and population-wide gRNA design. CRISPR-HAWK automates the entire workflow—from variant-aware preprocessing to gRNA discovery—delivering comprehensive outputs including ranked tables, annotated sequences, and high-quality figures. Its modular design ensures easy integration with existing pipelines and tools, such as [CRISPRme](https://github.com/pinellolab/CRISPRme) or [CRISPRitz](https://github.com/pinellolab/CRISPRitz), for subsequent off-target prediction and analysis of prioritized gRNAs.

## Table of Contents

0 [System Requirements](#0-system-requirements)
<br>1 [Installation](#1-installation)
<br>&nbsp;&nbsp;1.1 [Install CRISPR-HAWK from Mamba/Conda](#11-install-crispr-hawk-from-mambaconda)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.1.1 [Install Conda or Mamba](#111-install-conda-or-mamba)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.1.2 [Install CRISPR-HAWK](#112-install-crispr-hawk)
<br>&nbsp;&nbsp;1.2 [Install CRISPR-HAWK from Docker](#12-install-crispr-hawk-from-docker)
<br>&nbsp;&nbsp;1.3 [Install CRISPR-HAWK from PyPI](#13-install-crispr-hawk-from-pypi)
<br>&nbsp;&nbsp;1.4 [Install CRISPR-HAWK from Source Code](#14-install-crispr-hawk-from-source-code)
<br>&nbsp;&nbsp;1.5 [Install External Software Dependencies](#15-install-external-software-dependencies)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.5.1 [Install CRISPRitz (for Off-target Estimation)](#151-install-crispritz-for-off-target-estimation)
<br>2 [Usage](#2-usage)
<br>&nbsp;&nbsp;2.1 [General Syntax](#21-general-syntax)
<br>&nbsp;&nbsp;2.2 [Search](#22-search)
<br>&nbsp;&nbsp;2.3 [Convert gnomAD VCF](#23-convert-gnomad-vcf)
<br>&nbsp;&nbsp;2.4 [Prepare Data for CRISPRme](#24-prepare-data-for-crisprme)
<br>3 [Test](#3-test)
<br>4 [Citation](#4-citation)
<br>5 [Contacts](#5-contacts)
<br>6 [License](#6-license)

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

> ⚠️ **Note**: For optimal performance and stability, especially when dealing with large-scale variant datasets, ensure that your system meets or exceeds the recommended specifications.

## 1 Installation

This section provides step-by-step instructions to install CRISPR-HAWK and external dependencies. Choose the method that best suits your environment and preferences:

- **[Install CRISPR-HAWK from Mamba/conda](#11-install-crispr-hawk-from-mambaconda)** (recommended)
<br>Best for users seeking an isolated and reproducible environment with minimal manual dependency handling.

- **[Install CRISPR-HAWK from Docker](#12-install-crispr-hawk-from-docker)**
<br>Ideal for users who prefer containerized deployments or want to avoid configuring the environment manually.

- **[Install CRISPR-HAWK from PyPI](#13-install-crispr-hawk-from-pypi)**
<br>Quick option for Python users already working within a virtual environment. May require manual handling of some dependencies.

- **[Install CRISPR-HAWK from source code](#14-install-crispr-hawk-from-source-code)**
<br>Suitable for developers or contributors who want full control over the codebase or plan to customize CRISPR-HAWK.

> ⚠️ **Note:** We recommend using the Mamba/Conda or Docker installation methods for most users, as they ensure the highest compatibility and stability across systems.

### 1.1 Install CRISPR-HAWK from Mamba/Conda

#### 1.1.1 Install Conda or Mamba

Before installing CRISPR-HAWK, ensure that either **Conda** or **Mamba** is installed on your system. Based on recommendations from the Bioconda community and performance testing during CRISPR-HAWK development, we **recommend using [Mamba](https://mamba.readthedocs.io/en/latest/index.html)** over Conda. Mamba is a fast, efficient drop-in replacement for Conda, built with a high-performance dependency solver in C++.

**Installation Steps**

**1. Install Conda or Mamba**

* To install **Conda**, follow the official instructions:
<br>[Conda Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

* To install **Mamba**, follow the official instructions:
<br>[Mamba Installation Guide](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

**2. Configure Bioconda Channels**
Once Mamba (or Conda) is installed, configure your environment with the appropriate channels used by CRISPR-HAWK:
```bash
mamba config --add channels bioconda
mamba config --add channels defaults
mamba config --add channels conda-forge
mamba config --set channel_priority strict
```

> 💡 **Tip**: If you are using Conda instead of Mamba, simply replace `mamba` with `conda` in the commands above.

By completing these steps your system will be correctly configured to install CRISPR-HAWK and all required dependencies via Bioconda.


#### 1.1.2 Install CRISPR-HAWK

TBA

### 1.2 Install CRISPR-HAWK from Docker

TBA

### 1.3 Install CRISPR-HAWK from PyPI

TBA

### 1.4 Install CRISPR-HAWK from Source Code

Installing CRISPR-HAWK from source is ideal for developers, contributors, or users who wish to inspect or customize the codebase.

This method assumes you already have **Python 3.8** installed and accessible from your system’s environment.

**Prerequisites**

- Python **3.8** (strictly required)

- `git`

- A virtual environment (optional but recommended)

**Installation Steps**

**1. Clone the Repository**
```bash
git clone https://github.com/ManuelTgn/CRISPR-HAWK.git
cd CRISPR-HAWK
```

**2. (Optional) Create and Activate a Virtual Environment**
```bash
mamba create -n crisprhawk-env python=3.8 -y
mamba activate crisprhawk-env
```

**3. Install CRISPR-HAWK and Its Dependencies**
```bash
pip install .  # regular installation
pip install -e .  # development-mode installation
```

The `.` tells `pip` to install the current directory as a package, including all dependencies specified in `setup.py` or `pyproject.toml`.

**Quick Test**

Once installation is complete, verify that the command-line interface is working:
```bash
crisprhawk -h
```

If the help message is displayed correctly, CRISPR-HAWK is successfully installed and callable from any directory in your system.

### 1.5 Install External Software Dependencies

CRISPR-HAWK relies on a few external tools for certain optional features, such as genome-wide **off-target nomination**. These dependencies are **not bundled** with the core CRISPR-HAWK installation and must be installed separately if you wish to enable advanced features.

> ⚠️ **Note**: External tools are currently only supported on Linux-based systems. Windows and macOS users can still run the core pipeline (variant-aware gRNA search, scoring, annotation), but **off-target estimation** will not be available.

> 💡 **Tip**: Installing these tools requires Mamba or Conda to be available on your system. If you haven't installed Mamba yet, refer to the instructions in [Section 1.1.1](#111-install-conda-or-mamba).

#### 1.5.1 Install CRISPRitz (for Off-target Estimation)

[CRISPRitz](https://github.com/pinellolab/CRISPRitz) is an efficient tool for nominating CRISPR-Cas off-target sites across large genomes accounting for mismatches and DNA/RNA bulges. CRISPR-HAWK uses it to enable **fast, high-throughput off-target estimation** in the reference genome for the identified candidate gRNAs.

> ⚠️ **Note**: CRISPRitz is required only if you plan to use the `--estimate-offtargets` feature in the `crisprhawk search` command.

**Installation Steps**:

1. Create a dedicated environment for CRISPRitz:
```bash
mamba create -n crispritz-crisprhawk python=3.8 -y
```

2. Activate the environment::
```bash
mamba activate crispritz-crisprhawk
```

3. Install CRISPRitz (latest version):
```bash
mamba install crispritz=2.6.6 -y
```

> 📝 **Why use a separate environment?**
<br>This prevents potential conflicts between CRISPR-HAWK's and CRISPRitz’s dependencies, and ensures better reproducibility.


**Test your installation**

You can confirm that CRISPRitz is correctly installed and functional by running:
```bash
mamba run -n crispritz-crisprhawk crispritz.py --help
```

If everything is working, the CRISPRitz help menu should appear, displaying available options and usage instructions.

## 2 Usage

CRISPR-HAWK provides multiple functionalities designed to support variant- and haplotype-aware CRISPR guide design, gRNA efficiency evaluation, and integration with downstream analysis pipelines. Each command serves a distinct role in the workflow.

### 2.1 General Syntax
```bash
crisprhawk <command> [options]
```

To view available commands:
```bash
crisprhawk --help
```

To check version:
```bash
crisprhawk --version
```

### 2.2 Search

The `crisprhawk search` command is the core functionality of CRISPR-HAWK, designed to identify and annotate candidate gRNAs in both reference and variant genomes.
It integrates variant-aware search, functional annotation, and predictive scoring to help you prioritize the most robust and context-aware guides for CRISPR editing.

The search includes:

* Support for any Cas system (Cas9, Cpf1, SaCas9, etc.)
* Compatibility with custom PAM sequences and guide lengths
* Variant-aware design from individual or population-level VCF files (SNVs and indels)
* Scoring using **Azimuth**, **RS3**, **CFDon**, **Elevation-on**, and **DeepCpf1**
* Functional and gene annotation using user-specified BED files 
* Optional estimation and reporting of **off-targets**
* Output in detailed and structured reports (TSV, haplotype tables, off-target tables)

Usage:
```bash
crisprhawk search -f <fasta-dir> -r <bedfile> -v <vcf-dir> -p <pam> -g <guide-length> -o <output-dir>
```

> ⚠️ **Note**: All FASTA files in `<fasta-dir>` must be one per chromosome (e.g., chr1.fa, chr2.fa, etc.).

---

#### Required Arguments

| Option                       | Description                                                                                                                                    |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `-f`, `--fasta <FASTA-DIR>` | Directory containing chromosome-separated FASTA files for the reference genome. All files will be used as reference input.                                                      |
| `-r`, `--regions <BED-FILE>` | BED file defining the target regions where candidate gRNAs will be searched (e.g., promoters, exons, enhancers).                      |
| `-v`, `--vcf <VCF-DIR>`      | *(Optional but recommended)* Directory containing per-chromosome VCF files for variant-aware guide design. If omitted, the tool performs reference-only analysis. |
| `-p`, `--pam <PAM>`          | PAM sequence used to define valid gRNA targets (e.g., `NGG` for SpCas9, `TTTV` for Cpf1).                                                      |
| `-g`, `--guide-len <LENGTH>` | Length of the guide RNA (excluding PAM), e.g., 20 for SpCas9.                                                                         |

#### Optional Arguments

| Option                                         | Description                                                                                                                                                                       |
| ---------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-i`, `--fasta-idx <FAI>`                      | Optional FASTA index file (`.fai`). If not provided, it will be computed automatically.                                                                                           |
| `--right`                                      | By default, guides are extracted **upstream** of the PAM. Use this flag to extract them **downstream** (right side), useful for Cpf1 or other reverse-PAM systems.                |
| `--no-filter`                                  | By default, only VCF variants with `FILTER == PASS` are considered. Use this flag to include **all variants**, regardless of filter status.                                       |
| `--annotation <BED1 BED2 ...>`                 | Provide one or more BED files with **custom genomic** features (e.g., enhancers, DHS, regulatory elements). Must include a 4th column with annotation name. |
| `--annotation-colnames <name1 name2 ...>`      | Custom column names for the annotations from the `--annotation` files. Must match the number and order of BED files.                                                      |
| `--gene-annotation <GENE-BED1 GENE-BED2 ...>`  | One or more **gene annotation BED files** (9-column format). Must follow GENCODE-style structure, with gene name in the 9th column and gene feature (e.g., exon, UTR) in the 7th. |
| `--gene-annotation-colnames <name1 name2 ...>` | Custom column names for gene annotations, matching the files in `--gene-annotation`.                                                                                              |
| `--haplotype-table`                            | When enabled, a TSV file reporting haplotype-aware variants and their associated guide matches will be produced.                                                                  |
| `-t`, `--threads <INT>`                        | Number of threads to use for parallel processing. Use `-t 0` to utilize **all available cores**. *(Default: 1)*                                                                   |
| `-o`, `--outdir <DIR>`       | Output directory to store results, reports, and intermediate files. If not specified, defaults to the current working directory.               |
| `--verbosity <LEVEL>`                          | Controls the verbosity of logs. Options: `0` = Silent, `1` = Normal, `2` = Verbose, `3` = Debug. *(Default: 1)*                                                                   |
| `--debug`                                      | Enables **debug mode**, showing full stack traces and internal logs for troubleshooting.                                                                                          |

---

Example:

```bash
crisprhawk search \
  -f hg38.fa \
  -r targets.bed \
  -v 1000G/ \
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
* Consider population variants from VCFs in `1000G/`
* Annotate guides using custom regulatory and gene annotations 
* Run in parallel using 8 threads

#### Off-targets Estimation (Optional)

> 🐧 **Note**: Off-target estimation is **only available on Linux-based operating systems** and assumes you have successfully followed the installation instructions in the [Install CRISPRitz](#151-install-crispritz-for-off-target-estimation) section.

CRISPR-HAWK supports genome-wide off-target nomination through integration with [CRISPRitz](https://github.com/pinellolab/CRISPRitz), enabling the identification of potential unintended gRNA binding sites across the reference genome.

> ⚠️ **Note**: Off-target nomination in CRISPR-HAWK is limited to the reference genome only, for performance and scalability reasons. If you need to estimate off-targets while accounting for genetic variants (e.g., SNVs, indels, population haplotypes), we recommend using [CRISPRme](https://github.com/pinellolab/CRISPRme) — a specialized, variant-aware off-target analysis tool.

When enabled, the off-target module allows:

* Comprehensive search of potential off-targets in the reference genome

* Support for mismatches and DNA/RNA bulges

* Output of structured off-target reports for guide prioritization

*Enabling Off-Target Search*

To run off-target estimation, you must enable `--estimate-offtargets` and provide a pre-built genome index compatible with CRISPRitz.

| Option                      | Description                                                                                          |
| --------------------------- | ---------------------------------------------------------------------------------------------------- |
| `--estimate-offtargets`     | Activates the off-target search pipeline for all candidate guides using CRISPRitz.                   |
| `--crispritz-index <DIR>`   | Directory containing the CRISPRitz genome index. Must match the FASTA files provided with `--fasta`. |
| `--mm <INT>`                | Max number of mismatches allowed in off-target search. *(Default: 4)*                                  |
| `--bdna <INT>`              | Max number of DNA bulges allowed. *(Default: 0)*                                                       |
| `--brna <INT>`              | Max number of RNA bulges allowed. *(Default: 0)*                                                       |

Example:
```bash
crisprhawk search \
  -f genome_fasta/ \
  -r targets.bed \
  -v 1000G/ \
  -p NGG \
  -g 20 \
  -o results/ \
  --estimate-offtargets \
  --crispritz-index genome_library/NGG_2_hg38/ \
  --mm 6 \
  --bdna 1 \
  --brna 1 \
  -t 8
```

This will:

* Design guides from `targets.bed`, using 1000G variants

* Estimate genome-wide off-target sites (up to 6 mismatches, 1 bulge each)

* Save off-targets in a detailed TSV file per guide

* Use 8 threads for faster processing


**Build a CRISPRitz Genome Index**

To use CRISPRitz for off-target estimation, you must first generate an indexed version of your reference genome. The `index-genome` command in CRISPRitz precomputes a genome-wide searchable index for a specific PAM and guide configuration — similar to building a BWA index. This index is essential for fast and efficient off-target search, especially when allowing bulges (RNA or DNA) and scanning across large genomes or thousands of guide RNAs. 

> 💡 **Tip**: CRISPRitz indexes can be reused multiple times and do not need to be generated before each run

*Required Inputs*

1. **Output name of the genome index**
  <br>This is the name of the directory that will store the generated index (e.g., `hg38_index/`).

2. **Directory containing genome FASTA files**
  <br>The genome must be split by chromosome — i.e., one file per chromosome (e.g., `chr1.fa`, `chr2.fa`, ...).

3. **PAM configuration file**
  <br>A text file containing:

    * The full gRNA+PAM string with `N`s representing the guide length
    <br>(e.g., `NNNNNNNNNNNNNNNNNNNNGG`)

    * A space-separated number indicating the PAM length
    <br>(e.g., `3` for NGG where the PAM is 3 bp)

Example file content (`pamNGG.txt`):
```
NNNNNNNNNNNNNNNNNNNNGG 3
```

4. Maximum number of bulges to index (`-bMax`)
  <br>This sets the maximum number of **DNA and RNA bulges** that can be used in future off-target searches using this index.

5. Number of threads (`-th`, optional)
Number of threads for parallelization during indexing.

*Output*

A new directory (named after the first argument, e.g., `hg38_index/`) containing:

* Indexed `.bin` files (one per chromosome)

* These include all candidate target sites for the selected PAM, including padding required to enable bulge-aware search.

*Example Usage*
```bash
crispritz.py index-genome hg19_ref hg19_ref/ pam/pamNGG.txt -bMax 2 -th 4
```

| Parameter        | Description                                                                   |
| ---------------- | ----------------------------------------------------------------------------- |
| `hg19_ref`       | Name of the output directory that will store the indexed genome.              |
| `hg19_ref/`      | Folder containing FASTA files (one per chromosome).                           |
| `pam/pamNGG.txt` | Path to the PAM specification file (see format above).                        |
| `-bMax 2`        | Index will support searches with up to **2 RNA bulges** and **2 DNA bulges**. |
| `-th 4`          | Use 4 threads during index creation (optional, increases performance).        |


>💡 **Tip**: Need more help generating a CRISPRitz index?
<br>Refer to the [CRISPRitz documentation on indexing](https://github.com/pinellolab/CRISPRitz#crispritz-installation-and-usage) for complete details, tips, and supported genome formats.

### 2.3 Convert gnomAD VCF

The `crisprhawk convert-gnomad-vcf` command is a utility designed to convert gnomAD VCF files (version ≥ 3.1) into a format compatible with CRISPR-HAWK. This step is essential when working with large-scale population datasets (e.g., gnomAD v3.1, v4.1), ensuring proper variant normalization, filtering, and indexing.

This module:

* Supports both `.vcf.bgz` and `.vcf.gz` formats

* Extracts and preserves **sample-level** genotypes for population-aware haplotype modeling

* Can optionally handle **joint allele frequency files** from gnomAD v4.1 (joint genome/exome releases)

* Supports parallel processing for large-scale datasets

Usage:
```bash
crisprhawk convert-gnomad-vcf -d <vcf-dir> -o <output-dir>
```

> ⚠️ **Note**: All `.vcf.bgz` (or `.vcf.gz`) files in `<vcf-dir>` will be processed. Ensure tabix indices (`.tbi`) are present in the same folder.

#### Required Arguments
| Option                  | Description                                                                                               |
| ----------------------- | --------------------------------------------------------------------------------------------------------- |
| `-d`, `--vcf-dir <DIR>` | Directory containing per-chromosome **gnomAD VCF files** (compressed as `.vcf.bgz` or `.vcf.gz`).         |
| `-o`, `--outdir <DIR>`  | Output directory for the converted VCF files. If not provided, defaults to the current working directory. |


#### Optional Arguments
| Option                  | Description                                                                                                           |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `--joint`               | Enable this flag when using **gnomAD v4.1 joint** genome/exome files. Ensures correct parsing of allele frequencies.  |
| `--keep`                | Include **all variants**, regardless of the VCF `FILTER` field. By default, only `FILTER=PASS` variants are retained. |
| `--suffix <SUFFIX>`     | Optional suffix to append to the output VCF filenames (e.g., `_converted`). Useful for traceability.                  |
| `-t`, `--threads <INT>` | Number of threads to use. Set `-t 0` to use all available CPU cores. *(Default: 1)*                                   |
| `--verbosity <LEVEL>`   | Logging verbosity. Options: `0` = Silent, `1` = Normal, `2` = Verbose, `3` = Debug *(Default: 1)*                     |
| `--debug`               | Enable debug mode and print full error tracebacks.                                                                    |


**Example**
```bash
crisprhawk convert-gnomad-vcf \
  -d gnomad_v4.1/ \
  -o converted_vcfs/ \
  --suffix _crisprhawk \
  -t 4
```

This command will:

* Convert all `.vcf.bgz` files in `gnomad_v4.1/`

* Retain only variants with `FILTER=PASS`

* Append `_crisprhawk` to each output filename

* Run using 4 threads

> ⚠️ **Note**: Make sure your input files are correctly indexed (`.tbi`) and are chromosome-specific. CRISPR-HAWK expects one VCF per chromosome.

### 2.4 Prepare Data for CRISPRme

The `crisprhawk prepare-data-crisprme` command transforms a CRISPR-HAWK report into guide files compatible with **[CRISPRme](https://github.com/pinellolab/CRISPRme)**, enabling downstream variant- and haplotype-aware off-target prediction using CRISPRme's framework.

This module:

* Extracts gRNA sequences from a CRISPR-HAWK report

* Creates per-guide FASTA files required by CRISPRme

* Optionally generates a **PAM specification file** compatible with the selected CRISPR system

Usage:
```bash
crisprhawk prepare-data-crisprme --report <crisprhawk-report> -o <output-dir>
```

#### Required Arguments

| Option                         | Description                                                                                      |
| ------------------------------ | ------------------------------------------------------------------------------------------------ |
| `--report <CRISPRHAWK-REPORT>` | Path to a CRISPR-HAWK report file containing guide sequences.        |
| `-o`, `--outdir <DIR>`         | Output directory for the CRISPRme-compatible guide files. Defaults to current working directory. |

#### Optional Arguments

| Option              | Description                                                                            |
| ------------------- | -------------------------------------------------------------------------------------- |
| `--create-pam-file` | If set, also generates a **PAM file** in CRISPRme format in the same output directory. |
| `--debug`           | Enable debug mode and print full error tracebacks.                                     |

**Output Structure**

The command will generate:

* One FASTA file per gRNA (named using the guide label or sequence)

* Optionally, a `pam.txt` file for CRISPRme if `--create-pam-file` is specified

These files can be directly used as input to CRISPRme's `--guide` and `--pam` options.

**Example**

```bash
crisprhawk prepare-data-crisprme \
  --report results/crisprhawk_report.tsv \
  -o crisprme_inputs/ \
  --create-pam-file
```

This command will:

* Parse all guides listed in the `crisprhawk_report.tsv` report

* Write per-guide FASTA files to the `crisprme_inputs/` folder

* Generate a PAM file (`pam.txt`) compatible with CRISPRme

> 💡 **Tip**: This is especially useful when transitioning from **on-target selection** in CRISPR-HAWK to **off-target analysis** in CRISPRme.

## 3 Test

## 4 Citation

If you use CRISPR-HAWK in your research, please cite:

```bibtex
@software{crisprhawk2025,
  title = {CRISPR-HAWK: Haplotype- and vAriant-aWare guide design toolKit},
  author = {Manuel Tognon},
  year = {2025},
  url = {https://github.com/ManuelTgn/CRISPR-HAWK}
}
```

## 5 Contacts

* Manuel Tognon
  <br>manuel.tognon@univr.it

* Rosalba Giugno
  <br>rosalba.giugno@univr.it

* Luca Pinello
  <br>lpinello@mgh.harvard.edu

## 6 License

CRISPR-HAWK is licensed under the AGPL-3.0 license, which permits its use for academic research purposes only.

For any commercial or for-profit use, please contact the authors.


