"""
CRISPR-HAWK {version}

Copyright (C) 2025 Manuel Tognon <manu.tognon@gmail.com> <manuel.tognon@univr.it> <mtognon@mgh.harvard.edu>

CRISPR-HAWK: Haplotype- and vAriant-aWare guide design toolKit

CRISPR-HAWK is a tool for haplotype- and variant-aware guide RNAs design (support all CRISPR systems), gRNA
efficiency assessment (support for Cas9 and Cpf1 systems), and analysis of genetic diversity impact on
on-targets specificity.

Usage:
    crisprhawk search -f <fasta> -r <bedfile> -v <vcf> -p <pam> -g <guide-length> -o <output-dir>
    crisprhawk convert-gnomad-vcf -d <vcf-dir> -o <output-dir>
    crisprhawk prepare-data-crisprme --report <crisprhawk-report> -o <output-dir>

Run 'crisprhawk -h/--help' to display the complete help
"""

from .crisprhawk_argparse import (
    CrisprHawkArgumentParser,
    CrisprHawkSearchInputArgs,
    CrisprHawkConverterInputArgs,
    CrisprHawkPrepareDataInputArgs,
)
from .crisprhawk import (
    crisprhawk_search,
    crisprhawk_converter,
    crisprhawk_prepare_data_crisprme,
)
from .exception_handlers import sigint_handler
from .crisprhawk_version import __version__
from .utils import prepare_package, TOOLNAME

from argparse import _SubParsersAction
from time import time

import sys
import os

# crisprhawk commands
SEARCH = "search"
CONVERTGNOMADVCF = "convert-gnomad-vcf"
PREPAREDATACRISPRME = "prepare-data-crisprme"
COMMANDS = [SEARCH, CONVERTGNOMADVCF, PREPAREDATACRISPRME]


def create_parser_crisprhawk() -> CrisprHawkArgumentParser:
    # force displaying docstring at each usage display and force
    # the default help to not being shown
    parser = CrisprHawkArgumentParser(usage=__doc__, add_help=False)  # type: ignore
    group = parser.add_argument_group("Options")  # arguments group
    # add help and version arguments
    group.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit"
    )
    group.add_argument(
        "--version",
        action="version",
        help=f"Show {TOOLNAME} version and exit",
        version=__version__,
    )
    # create subparsers for different functionalities
    subparsers = parser.add_subparsers(
        dest="command",
        title="Available commands",
        metavar="",  # needed for help formatting (avoid <command to be displayed>)
        description=None,
    )
    # crisprhawk search command
    parser_search = create_search_parser(subparsers)
    # crisprhawk graphical-reports command
    # crisprhawk convert-gnomad-vcf command
    parser_converter = create_converter_parser(subparsers)
    # crisprhawk prepare-data-crisprme command
    parser_prepare = create_parser_prepare_data(subparsers)
    return parser


def create_search_parser(subparser: _SubParsersAction) -> _SubParsersAction:
    parser_search = subparser.add_parser(
        SEARCH,
        usage="CRISPR-HAWK search {version}\n\nUsage:\n"
        "\tcrisprhawk search -f <fasta> -r <bedfile> -v <vcf> -p <pam> -g "
        "<guide-length> -o <output-dir>\n\n",
        description="Automated end-to-end search pipeline that processes raw input "
        "data through gRNA identification, scoring, and annotation of results",
        help="perform a comprehensive gRNA search across the reference genome "
        "and optionally variant-aware genomes. Includes Azimuth and RS3 (for "
        "Cas9 systems), and DeepCpf1 (for Cpf1 systems) scores, CFDon score (for "
        "Cas systems) to evaluate genetic diversity impact on on-targets, and "
        "automated gRNA annotation",
    )
    parser_search.add_argument(
        "-f",
        "--fasta",
        type=str,
        metavar="FASTA-FILE",
        dest="fasta",
        required=True,
        help="Reference genome in FASTA format used for guide search",
    )
    parser_search.add_argument(
        "-i",
        "--fasta-idx",
        type=str,
        metavar="FASTA-IDX",
        dest="fasta_idx",
        nargs="?",
        default="",
        help="Optional FASTA index file (FAI) for the input reference (default: "
        "compute FAI)",
    )
    parser_search.add_argument(
        "-r",
        "--regions",
        type=str,
        metavar="GENOMIC-REGIONS-BED",
        dest="bedfile",
        required=True,
        help="BED file specifying genomic regions where guides will be searched",
    )
    parser_search.add_argument(
        "-v",
        "--vcf",
        type=str,
        metavar="VCF-DIR",
        dest="vcf",
        nargs="?",
        default="",
        help="Optional folder storing VCF files to consider in the guide design. "
        "(default: no variant-aware analysis)",
    )
    parser_search.add_argument(
        "-p",
        "--pam",
        type=str,
        metavar="PAM",
        dest="pam",
        required=True,
        help="PAM sequence used to identify candidate guides (e.g., NGG, NAG, " "etc.)",
    )
    parser_search.add_argument(
        "-g",
        "--guide-len",
        type=int,
        metavar="GUIDE-LENGTH",
        dest="guidelen",
        required=True,
        help="Length of the guide (excluding the PAM)",
    )
    parser_search.add_argument(
        "--right",
        action="store_true",
        dest="right",
        default=False,
        help="If set, guides are extracted downstream (right side) of the PAM "
        "site. (default: guides are extracted upstream (left side))",
    )
    parser_search.add_argument(
        "-o",
        "--outdir",
        type=str,
        metavar="OUTDIR",
        dest="outdir",
        nargs="?",
        default=os.getcwd(),
        help="Output directory where reports and results will be saved. "
        "(default: current working directory)",
    )
    parser_search.add_argument(
        "--no-filter",
        action="store_true",
        dest="no_filter",
        default=False,
        help="If set, all variants in the input VCF file will be considered "
        "regardless of FILTER status (default: only variants with FILTER == "
        "'PASS' are used)",
    )
    parser_search.add_argument(
        "--annotation",
        type=str,
        metavar="ANNOTATION-BED",
        dest="annotations",
        nargs="*",
        default=[],
        help="One or more BED files specifying genomic regions used to annotate "
        "guide candidates. Each file should follow the standard BED format "
        "(at least: chrom, start, end), and should include additional annotation "
        "on the 4th column (default: no annotation)",
    )
    parser_search.add_argument(
        "--annotation-colnames",
        type=str,
        metavar="ANNOTATION-COLNAMES",
        dest="annotation_colnames",
        nargs="*",
        default=[],
        help="List of custom column names to use in the final report. Each name "
        "corresponds to one of the input BED files provided with '--annotation'. "
        "Must match the number and order of files in '--annotation' (default: "
        "annotation columns are named 'annotation_<i>')",
    )
    parser_search.add_argument(
        "--gene-annotation",
        type=str,
        metavar="GENE-ANNOTATION-BED",
        dest="gene_annotations",
        nargs="*",
        default=[],
        help="One or more BED files specifying gene regions used to annotate guide "
        "candidates. The file should follow standard BED format (chrom, start, "
        "end) and should include 9 columns. The 7th column should indicate the "
        "gencode feature (e.g., start_codon, exon, etc.). The 9th column should "
        "be a semicolon-separated list with the gene name identified by "
        "gene_name (e.g., gene_id=ENSG00000281518;gene_name=FOXO6;...;) "
        "(default: no gene annotation)",
    )
    parser_search.add_argument(
        "--gene-annotation-colnames",
        type=str,
        metavar="GENE-ANNOTATION-COLNAMES",
        dest="gene_annotation_colnames",
        nargs="*",
        default=[],
        help="Custom column names to assign to the gene annotation fields in the "
        "final report. These should correspond to the columns present in the BED "
        "file provided via '--gene-annotation' (default: column names assigned "
        "as 'gene_annotation_<i>')",
    )
    parser_search.add_argument(
        "--haplotype-table",
        action="store_true",
        dest="haplotype_table",
        default=False,
        help="When enabled, the haplotype table is returned in the output folder "
        "as TSV file (default: disabled)",
    )
    parser_search.add_argument(
        "--estimate-offtargets",
        action="store_true",
        dest="estimate_offtargets",
        default=False,
        help="estimate potential off-target sites on reference genome for each "
        "guide RNA candidate using CRISPRitz. This feature is only supported on "
        "Linux-based systems (default: disabled)",
    )
    parser_search.add_argument(
        "-t",
        "--threads",
        type=int,
        metavar="THREADS",
        dest="threads",
        required=False,
        default=1,
        help="Number of threads. Use 0 for using all available cores (default: 1)",
    )
    parser_search.add_argument(
        "--verbosity",
        type=int,
        metavar="VERBOSITY",
        dest="verbosity",
        nargs="?",
        default=1,  # minimal output
        help="Verbosity level of output messages: 0 = Silent, 1 = Normal, 2 = "
        "Verbose, 3 = Debug (default: 1)",
    )
    parser_search.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Enter debug mode and trace the full error stack",
    )
    return parser_search


def create_converter_parser(subparser: _SubParsersAction) -> _SubParsersAction:
    parser_converter = subparser.add_parser(
        CONVERTGNOMADVCF,
        usage="CRISPR-HAWK convert-gnomad-vcf {version}\n\nUsage:\n"
        "\tcrisprhawk convert-gnomad-vcf -d <vcf-dir> -o <output-dir>\n\n",
        description="Convert gnomAD VCF files (version ≥ 3.1) into a format "
        f"compatible with {TOOLNAME}. This utility preprocesses gnomAD VCFs to "
        "ensure both structural and content compatibility, and incorporates "
        "sample-level information to enable population-aware variant representation.",
        help=f"convert gnomAD VCFs (v3.1 or newer) into {TOOLNAME}-compatible "
        "format",
    )
    parser_converter.add_argument(
        "-d",
        "--vcf-dir",
        type=str,
        dest="gnomad_vcf_dir",
        metavar="GNOMAD-VCF-DIR",
        required=True,
        help="path to the directory containing gnomAD VCF files (with .vcf.bgz "
        "or vcf.gz extension). All .vcf.bgz files in the directory will be automatically "
        "processed",
    )
    parser_converter.add_argument(
        "-o",
        "--outdir",
        type=str,
        metavar="OUTDIR",
        dest="outdir",
        nargs="?",
        default=os.getcwd(),
        help="Output directory where converted VCF files will be saved "
        "(default: current working directory)",
    )
    parser_converter.add_argument(
        "--joint",
        action="store_true",
        dest="joint",
        help="Set this flag if the input VCFs contain joint allele frequencies, "
        "as in gnomAD v4.1 joint exomes/genomes releases (default: disabled)",
    )
    parser_converter.add_argument(
        "--keep",
        action="store_true",
        dest="keep",
        help="Retain all variants regardless of their FILTER status. "
        "By default, only variants with FILTER=PASS are included (default: "
        "disabled)",
    )
    parser_converter.add_argument(
        "--suffix",
        type=str,
        dest="suffix",
        required=False,
        default="",
        help="Optional suffix to append to the names of the converted VCF files. "
        "Useful for distinguishing output files (default: no suffix)",
    )
    parser_converter.add_argument(
        "-t",
        "--threads",
        type=int,
        metavar="THREADS",
        dest="threads",
        required=False,
        default=1,
        help="Number of threads. Use 0 for using all available cores (default: 1)",
    )
    parser_converter.add_argument(
        "--verbosity",
        type=int,
        metavar="VERBOSITY",
        dest="verbosity",
        nargs="?",
        default=1,  # minimal output
        help="Verbosity level of output messages: 0 = Silent, 1 = Normal, 2 = "
        "Verbose, 3 = Debug (default: 1)",
    )
    parser_converter.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Enter debug mode and trace the full error stack",
    )
    return parser_converter


def create_parser_prepare_data(subparser: _SubParsersAction) -> _SubParsersAction:
    parser_prepare = subparser.add_parser(
        PREPAREDATACRISPRME,
        usage="CRISPR-HAWK prepare-data-crisprme {version}\n\nUsage:\n"
        "\tcrisprhawk prepare-data-crisprme --report <crisprhawk-report> -o "
        "<output-dir>\n\n",
        description="Generate guide files from a CRISPR-HAWK report for "
        "downstream analysis with CRISPRme. For each guide listed in the report, "
        "this utility creates a guide file compatible with CRISPRme, enabling "
        "variant- and haplotype-aware off-target prediction",
        help=f"generate CRISPRme-compatible guide files from a CRISPR-HAWK report",
    )
    parser_prepare.add_argument(
        "--report",
        type=str,
        dest="report",
        metavar="CRISPRHAWK-REPORT",
        required=True,
        help="path to the CRISPR-HAWK report file containing guide sequences",
    )
    parser_prepare.add_argument(
        "--create-pam-file",
        action="store_true",
        dest="create_pam",
        help="If set, a PAM file suitable for CRISPRme will also be generated "
        "in the same output directory (default: disabled)",
    )
    parser_prepare.add_argument(
        "-o",
        "--outdir",
        type=str,
        metavar="OUTDIR",
        dest="outdir",
        nargs="?",
        default=os.getcwd(),
        help="Directory where the guide (and optional PAM) files will be saved "
        "(default: current working directory)",
    )
    parser_prepare.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Enter debug mode and trace the full error stack",
    )
    return parser_prepare


def main():
    prepare_package()  # check if models and data are available and uncompressed
    start = time()  # track elapsed time
    try:
        parser = create_parser_crisprhawk()  # parse input argument using custom parser
        if not sys.argv[1:]:  # no input args -> print help and exit
            parser.error_noargs()
        args = parser.parse_args(sys.argv[1:])  # parse input args
        if args.command == SEARCH:  # search command
            crisprhawk_search(CrisprHawkSearchInputArgs(args, parser))
        elif args.command == CONVERTGNOMADVCF:  # convert-gnoamd-vcf command
            crisprhawk_converter(CrisprHawkConverterInputArgs(args, parser))
        elif args.command == PREPAREDATACRISPRME:  # prepare-data-crisprme command
            crisprhawk_prepare_data_crisprme(
                CrisprHawkPrepareDataInputArgs(args, parser)
            )
    except KeyboardInterrupt as e:
        sigint_handler()  # catch SIGINT and exit gracefully
    sys.stdout.write(f"{TOOLNAME} - Elapsed time {(time() - start):.2f}s\n")


# --------------------------------> ENTRY POINT <--------------------------------
if __name__ == "__main__":
    main()
