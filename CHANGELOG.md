# Changelog

All notable changes to this project will be documented in this file.
This is the first release of **CRISPR‑HAWK** (v0.1.0) — initial fully functional version.

---

## v0.1.0 — 2025-09-16

### Added

* Command‑line interface `crisprhawk` with core functionality for variant‑ and haplotype‑aware gRNA (guide RNA) design.
* Support for multiple Cas systems: Cas9, SaCas9, Cpf1 (Cas12a), etc., with custom PAM sequences and guide lengths.
* Integration of large‑scale human variation datasets (e.g. 1000 Genomes, HGDP, gnomAD) for variant‑aware guide design.
* Ability to handle both single nucleotide variants (SNVs) and small insertions/deletions (indels).
* Gene and feature annotation via user‑supplied BED files (custom genomic features, regulatory regions, etc.).
* Haplotype‑aware processing and output of matches per haplotype.
* Output formats include ranked tables, annotated sequences, reports.

### Features / Options

* `--right`: extract guides downstream of the PAM (useful for reverse‑PAM systems)
* `--no-filter`: include variants regardless of filter status (not just `FILTER == PASS`)
* `--annotation` / `--annotation-colnames`: custom annotation files and column naming
* `--gene-annotation` / `--gene-annotation-colnames`: support for GENCODE‑style gene annotation (9‑column format)
* `--haplotype-table`: option to produce TSV files reporting haplotype‑aware variants and associated guide matches
* `--threads` / `-t`: parallelism, ability to use all cores via `-t 0`
* `--outdir`: configurable output directory for all results and reports
* `--verbosity`: logging levels (0 = Silent, up to 3 = Debug)
* `--debug`: full stack traces and verbose logs for troubleshooting

### Enhancements

* Added `--compute-elevation-score`: to compute Elevation & Elevation-On scores for guide efficiency (when combined guide + PAM length = 23 bp and guide is downstream of PAM)
* Added `--graphical-reports`: generate figures (pie chart of guide‐type distribution; delta plots showing effect of genetic diversity on guide efficiency and on‑target activity)

### Fixes / Internal Improvements

* Project structure with source code under `src/`, tests under `tests/`
* Packaging files: `setup.py`, `setup.cfg`, `pyproject.toml`, `requirements.txt` included for installation and dependency management
* Documentation: README covers installation, usage, commands; license file included (AGPL‑3.0 - academic use only)

### Known Limitations

* Off‑target estimation currently only works against the reference genome (not variant or haplotype aware)
* External dependencies (e.g. CRISPRitz for off‑target) must be installed separately
* Performance: large genome‑scale operations may need high memory / compute resources

---

If you want, I can prepare a version of this CHANGELOG with suggested dates, or polish it for your release on PyPI / GitHub.
