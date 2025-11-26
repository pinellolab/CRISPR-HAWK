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

## v0.1.1 – 2025-11-13

### Added

* Support for off-target annotation via CLI and core logic: users can now supply BED files + custom column names for enriched off-target reporting. (PR #3)  
* Workflow reorganised into discrete modules (`annotation.py`, `scoring.py`, `search_offtargets.py`), with the main driver updated for better clarity and maintainability. (PR #3)  
* Enhanced CLI argument parsing/validation in `crisprhawk_argparse.py`, including new flag `--candidate-guides` for detailed guide analysis and graphical reporting. (PR #3)  
* Introduced `CandidateGuide` class and associated logic (`candidate_guides.py`). Added new reporting features (dot-plot, scatter-plot) for variant‐effect analysis of candidate guides. (PR #3)  
* Unified `CrisprHawkSearchInputArgs` object for major modules to streamline parameter passing and reduce redundancy. (PR #3)  
* Removed requirement for manual FASTA index file (FAI) — index handling is now automatic within the Fasta class. (PR #3)  
* Improved handling of IUPAC nucleotide codes in micro-homology scoring and guide‐search routines (supports ambiguous bases better). (PR #3)  

### Fixed / Improved

* Added parallel computation for DeepCpf1 scoring (chunking guide sequences + `ProcessPoolExecutor`) to improve performance on large guide sets. (PR #4)  
* Refactored candidate‐guide handling and graphical reports: updated coordinate parsing to require strand information; improved helper functions for guide-ID assignment, delta calculations and colour palette generation. (PR #4)  
* Fixed IUPAC encoding logic in `haplotypes.py` (use of IUPACTABLE for reference base conversion) and cleaned up variant addition logic. (PR #4)  
* Fixed bugs in variant annotation and guide mapping logic: guide retrieval and annotation functions now use position‐maps and variant normalisation; handling of multiallelic variants improved. (PR #4)  
* Improved formatting for allele frequencies (scientific notation) and pie-chart representation in graphical reports; corrected parsing of AF field when no additional info present. (PR #4)  
* Refactored CLI argument groups: made `--outdir` optional in search & converter subcommands; updated README logo file. (PR #4)  
* Adjusted test assertions to reflect changes in allele frequency formatting and scoring logic. (PR #4)  

### Breaking Changes / Migration Notes

* The manual provision of the FASTA index (`.fai`) is no longer required – indexing is handled internally.  
* Some function signatures changed due to the move to the unified `CrisprHawkSearchInputArgs` object; if you extend CRISPR-HAWK modules or have custom downstream code, you may need to update accordingly.  
* Existing workflows remain supported and functional, but we recommend verifying your pipelines after update.

### Version Bump
- Version updated to **v0.1.1**.

---

## v0.1.2 – 2025-11-26

### Fixed / Improved

* Corrected the formatting and ordering of report columns in off-target estimation output for each reported guide. This fix ensures consistent column alignment across guides, resolves mismatches between header and values, and improves compatibility with downstream parsing/analysis tools (PR #5).

### Breaking Changes / Migration Notes

* Existing workflows remain supported and functional, but we recommend verifying your pipelines after update.

### Version Bump
- Version updated to **v0.1.2**.


