"""
CRISPR-HAWK {version}

Copyright (C) 2025 Manuel Tognon <manu.tognon@gmail.com> <manuel.tognon@univr.it> <mtognon@mgh.harvard.edu>

CRISPR-HAWK: Haplotype- and vAriant-aWare guide design toolKit

CRISPR-HAWK is a tool for haplotype- and variant-aware guide RNAs design (support all CRISPR systems), gRNA
efficiency assessment (support for Cas9 and Cpf1 systems), and analysis of genetic diversity impact on
on-targets specificity.

Usage:
    crisprhawk search -f <fasta-dir> -r <bedfile> -v <vcf-dir> -p <pam> -g <guide-length> -o <output-dir>
    crisprhawk convert-gnomad-vcf -d <vcf-dir> -o <output-dir>
    crisprhawk prepare-data-crisprme --report <crisprhawk-report> -o <output-dir>
    crisprhawk crispritz-config --env <crispritz-env>

Run 'crisprhawk -h/--help' to display the complete help
"""

from .crisprhawk_argparse import (
    CrisprHawkArgumentParser,
    CrisprHawkSearchInputArgs,
    CrisprHawkConverterInputArgs,
    CrisprHawkPrepareDataInputArgs,
    CrisprHawkCrispritzConfigInputArgs,
)
from .crisprhawk import (
    crisprhawk_search,
    crisprhawk_converter,
    crisprhawk_prepare_data_crisprme,
    crisprhawk_crispritz_config,
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
CRISPRPITZCONFIG = "crispritz-config"
COMMANDS = [SEARCH, CONVERTGNOMADVCF, PREPAREDATACRISPRME, CRISPRPITZCONFIG]


def create_parser_crisprhawk() -> CrisprHawkArgumentParser:
    """Creates and configures the main argument parser for the CRISPR-HAWK CLI.

    This function sets up the command-line interface, including all available
    commands and their arguments, for the CRISPR-HAWK toolkit.

    Returns:
        CrisprHawkArgumentParser: The configured argument parser for CRISPR-HAWK.
    """
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
    # crisprhawk convert-gnomad-vcf command
    parser_converter = create_converter_parser(subparsers)
    # crisprhawk prepare-data-crisprme command
    parser_prepare = create_parser_prepare_data(subparsers)
    # crisprhawk crispritz-config command
    parser_crispritz_config = create_crispritz_config_parser(subparsers)
    return parser


def create_search_parser(subparser: _SubParsersAction) -> _SubParsersAction:
    """Creates the argument parser for the CRISPR-HAWK search command.

    This function defines and configures all arguments and options available for
    the search functionality of CRISPR-HAWK.

    Args:
        subparser (_SubParsersAction): The subparsers object to which the search
            parser will be added.

    Returns:
        _SubParsersAction: The configured search command parser.
    """
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
        add_help=False,
    )
    general_group = parser_search.add_argument_group("General options")
    general_group.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    required_group = parser_search.add_argument_group("Options")
    required_group.add_argument(
        "-f",
        "--fasta",
        type=str,
        metavar="FASTA-FILE",
        dest="fasta",
        required=True,
        help="folder containing genome FASTA files for guide search. Each "
        "chromosome must be in a separate FASTA file (e.g., chr1.fa, chr2.fa). "
        "All files in the folder will be used as the reference genome",
    )
    required_group.add_argument(
        "-r",
        "--regions",
        type=str,
        metavar="GENOMIC-REGIONS-BED",
        dest="bedfile",
        required=True,
        help="BED file specifying genomic regions where guides will be searched",
    )
    required_group.add_argument(
        "-p",
        "--pam",
        type=str,
        metavar="PAM",
        dest="pam",
        required=True,
        help="PAM sequence used to identify candidate guides (e.g., NGG, NAG, " "etc.)",
    )
    required_group.add_argument(
        "-g",
        "--guide-len",
        type=int,
        metavar="GUIDE-LENGTH",
        dest="guidelen",
        required=True,
        help="length of the guide (excluding the PAM)",
    )
    required_group.add_argument(
        "-o",
        "--outdir",
        type=str,
        metavar="OUTDIR",
        dest="outdir",
        nargs="?",
        default=os.getcwd(),
        help="output directory where reports and results will be saved. "
        "(default: current working directory)",
    )
    optional_group = parser_search.add_argument_group("Optional arguments")
    optional_group.add_argument(
        "-i",
        "--fasta-idx",
        type=str,
        metavar="FASTA-IDX",
        dest="fasta_idx",
        nargs="?",
        default="",
        help="optional FASTA index file (FAI) for the input reference (default: "
        "compute FAI)",
    )
    optional_group.add_argument(
        "-v",
        "--vcf",
        type=str,
        metavar="VCF-DIR",
        dest="vcf",
        nargs="?",
        default="",
        help="optional folder storing VCF files to consider in the guide design. "
        "(default: no variant-aware analysis)",
    )
    optional_group.add_argument(
        "--right",
        action="store_true",
        dest="right",
        default=False,
        help="if set, guides are extracted downstream (right side) of the PAM "
        "site. (default: guides are extracted upstream (left side))",
    )
    optional_group.add_argument(
        "--no-filter",
        action="store_true",
        dest="no_filter",
        default=False,
        help="if set, all variants in the input VCF file will be considered "
        "regardless of FILTER status (default: only variants with FILTER == "
        "'PASS' are used)",
    )
    optional_group.add_argument(
        "--annotation",
        type=str,
        metavar="ANNOTATION-BED",
        dest="annotations",
        nargs="*",
        default=[],
        help="one or more BED files specifying genomic regions used to annotate "
        "guide candidates. Each file should follow the standard BED format "
        "(at least: chrom, start, end), and should include additional annotation "
        "on the 4th column (default: no annotation)",
    )
    optional_group.add_argument(
        "--annotation-colnames",
        type=str,
        metavar="ANNOTATION-COLNAMES",
        dest="annotation_colnames",
        nargs="*",
        default=[],
        help="list of custom column names to use in the final report. Each name "
        "corresponds to one of the input BED files provided with '--annotation'. "
        "Must match the number and order of files in '--annotation' (default: "
        "annotation columns are named 'annotation_<i>')",
    )
    optional_group.add_argument(
        "--gene-annotation",
        type=str,
        metavar="GENE-ANNOTATION-BED",
        dest="gene_annotations",
        nargs="*",
        default=[],
        help="one or more BED files specifying gene regions used to annotate guide "
        "candidates. The file should follow standard BED format (chrom, start, "
        "end) and should include 9 columns. The 7th column should indicate the "
        "gencode feature (e.g., start_codon, exon, etc.). The 9th column should "
        "be a semicolon-separated list with the gene name identified by "
        "gene_name (e.g., gene_id=ENSG00000281518;gene_name=FOXO6;...;) "
        "(default: no gene annotation)",
    )
    optional_group.add_argument(
        "--gene-annotation-colnames",
        type=str,
        metavar="GENE-ANNOTATION-COLNAMES",
        dest="gene_annotation_colnames",
        nargs="*",
        default=[],
        help="custom column names to assign to the gene annotation fields in the "
        "final report. These should correspond to the columns present in the BED "
        "file provided via '--gene-annotation' (default: column names assigned "
        "as 'gene_annotation_<i>')",
    )
    optional_group.add_argument(
        "--compute-elevation-score",
        action="store_true",
        dest="compute_elevation",
        default=False,
        help="compute Elevation and Elevation-on scores to evaluate guide "
        "efficiency. This requires that the combined length of the guide and "
        "PAM is exactly 23 bp, and that the guide sequence is located upstream "
        "of the PAM (default: disabled)",
    )
    optional_group.add_argument(
        "--haplotype-table",
        action="store_true",
        dest="haplotype_table",
        default=False,
        help="when enabled, the haplotype table is returned in the output folder "
        "as TSV file (default: disabled)",
    )
    optional_group.add_argument(
        "--estimate-offtargets",
        action="store_true",
        dest="estimate_offtargets",
        default=False,
        help="estimate potential off-target sites on reference genome for each "
        "guide RNA candidate using CRISPRitz. This feature is only supported on "
        "Linux-based systems (default: disabled)",
    )
    optional_group.add_argument(
        "--crispritz-index",
        type=str,
        dest="crispritz_index",
        required=False,
        default="",
        help="path to the genome index directory generated by CRISPRitz. "
        "Required only when using --estimate-offtargets (default: ignored)",
    )
    optional_group.add_argument(
        "--mm",
        type=int,
        dest="mm",
        required=False,
        default=4,
        help="maximum number of mismatches to consider during off-target "
        "estimation. Only used if --estimate-offtargets is enabled (default: 4)",
    )
    optional_group.add_argument(
        "--bdna",
        type=int,
        dest="bdna",
        required=False,
        default=0,
        help="maximum number of DNA bulges to consider during off-target "
        "estimation. Only used if --estimate-offtargets is enabled (default: 0)",
    )
    optional_group.add_argument(
        "--brna",
        type=int,
        dest="brna",
        required=False,
        default=0,
        help="maximum number of RNA bulges to consider during off-target "
        "estimation. Only used if --estimate-offtargets is enabled (default: 0)",
    )
    optional_group.add_argument(
        "--offtargets-annotation",
        type=str,
        metavar="ANNOTATION-BED",
        dest="offtargets_annotations",
        nargs="*",
        default=[],
        help="one or more BED files specifying genomic regions used to annotate "
        "estimated offtargets. Each file should follow the standard BED format "
        "(at least: chrom, start, end), and should include additional annotation "
        "on the 4th column (default: no annotation)",
    )
    optional_group.add_argument(
        "--offtargets-annotation-colnames",
        type=str,
        metavar="ANNOTATION-COLNAMES",
        dest="offtargets_annotation_colnames",
        nargs="*",
        default=[],
        help="list of custom column names to use in the offtargets report. Each "
        "name corresponds to one of the input BED files provided with "
        "'--offtargets-annotation'. Must match the number and order of files in "
        "'--offtargets-annotation' (default: annotation columns are named "
        "'annotation_<i>')",
    )
    optional_group.add_argument(
        "--graphical-reports",
        action="store_true",
        dest="graphical_reports",
        default=False,
        help="generate graphical reports to summarize findings. Includes a pie "
        "chart showing the distribution of guide types (reference, spacer+PAM "
        "alternative, spacer alternative, PAM alternative) and delta plots "
        "illustrating the impact of genetic diversity on guide efficiency and "
        "on-target activity (default: disabled)",
    )
    optional_group.add_argument(
        "-t",
        "--threads",
        type=int,
        metavar="THREADS",
        dest="threads",
        required=False,
        default=1,
        help="number of threads. Use 0 for using all available cores (default: 1)",
    )
    optional_group.add_argument(
        "--verbosity",
        type=int,
        metavar="VERBOSITY",
        dest="verbosity",
        nargs="?",
        default=1,  # minimal output
        help="verbosity level of output messages: 0 = Silent, 1 = Normal, 2 = "
        "Verbose, 3 = Debug (default: 1)",
    )
    optional_group.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="enter debug mode and trace the full error stack",
    )
    return parser_search


def create_converter_parser(subparser: _SubParsersAction) -> _SubParsersAction:
    """Creates the argument parser for the CRISPR-HAWK convert-gnomad-vcf command.

    This function defines and configures all arguments and options available for
    converting gnomAD VCF files to a CRISPR-HAWK compatible format.

    Args:
        subparser (_SubParsersAction): The subparsers object to which the converter
            parser will be added.

    Returns:
        _SubParsersAction: The configured convert-gnomad-vcf command parser.
    """
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
        add_help=False,
    )
    general_group = parser_converter.add_argument_group("General options")
    general_group.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
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
    """Creates the argument parser for the CRISPR-HAWK prepare-data-crisprme command.

    This function defines and configures all arguments and options available for
    generating CRISPRme-compatible guide files from a CRISPR-HAWK report.

    Args:
        subparser (_SubParsersAction): The subparsers object to which the
            prepare-data-crisprme parser will be added.

    Returns:
        _SubParsersAction: The configured prepare-data-crisprme command parser.
    """
    parser_prepare = subparser.add_parser(
        PREPAREDATACRISPRME,
        usage="CRISPR-HAWK prepare-data-crisprme {version}\n\nUsage:\n"
        "\tcrisprhawk prepare-data-crisprme --report <crisprhawk-report> -o "
        "<output-dir>\n\n",
        description="Generate guide files from a CRISPR-HAWK report for "
        "downstream analysis with CRISPRme. For each guide listed in the report, "
        "this utility creates a guide file compatible with CRISPRme, enabling "
        "variant- and haplotype-aware off-target prediction",
        help="generate CRISPRme-compatible guide files from a CRISPR-HAWK report",
        add_help=False,
    )
    general_group = parser_prepare.add_argument_group("General options")
    general_group.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
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


def create_crispritz_config_parser(subparser: _SubParsersAction) -> _SubParsersAction:
    """Creates the argument parser for the CRISPR-HAWK crispritz-config command.

    This function defines and configures all arguments and options for managing
    CRISPRitz integration settings, including environment and target directory
    configuration.

    Args:
        subparser (_SubParsersAction): The subparsers object to which the
            crispritz-config parser will be added.

    Returns:
        _SubParsersAction: The configured crispritz-config command parser.
    """
    parser_crispritz_config = subparser.add_parser(
        CRISPRPITZCONFIG,
        usage="CRISPR-HAWK crispritz-config {version}\n\nUsage:\n"
        "\tcrisprhawk crispritz-config --env <crispritz-env-name> --targets-dir "
        "<crispritz-targets-dir>\n\n",
        description="Configure CRISPRitz integration settings including "
        "environment name and output directories. This command manages the "
        "configuration for CRISPRitz integration, allowing the user to specify "
        "the conda/mamba environment where CRISPRitz is installed and customize "
        "where target files will be stored. Configuration is automatically saved "
        "and persists across sessions. The configuration is stored in a JSON "
        "file and can be modified at any time. Use --show to display current "
        "settings or --reset to restore defaults.",
        help="configure CRISPRitz environment and target storage settings",
        add_help=False,
    )
    general_group = parser_crispritz_config.add_argument_group("General options")
    general_group.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    # main configuration arguments
    config_group = parser_crispritz_config.add_argument_group("Configuration settings")
    config_group.add_argument(
        "--env",
        type=str,
        metavar="CRISPRITZ-ENV-NAME",
        dest="env_name",
        required=False,
        default=None,
        help="name of the conda/mamba environment where CRISPRitz is installed. "
        "This environment will be activated when running CRISPRitz commands. "
        "Must be a valid environment name accessible via 'conda activate <name>' "
        "or 'mamba activate <name>'. (default: 'crispritz')",
    )
    config_group.add_argument(
        "--targets-dir",
        type=str,
        metavar="CRISPRITZ-TARGETS-DIR",
        dest="targets_dir",
        required=False,
        default=None,
        help="Directory path where CRISPRitz target files will be stored. "
        "Can be absolute (/path/to/targets) or relative (./targets). "
        "Directory will be created if it doesn't exist. (default: "
        "'.crispritz_targets' - hidden folder in search output directory)",
    )
    # action arguments
    action_group = parser_crispritz_config.add_argument_group("Configuration actions")
    action_group.add_argument(
        "--show",
        dest="show",
        action="store_true",
        help="Display current CRISPRitz configuration settings without making "
        "changes. Shows environment name, targets directory, and configuration "
        "file location.",
    )
    action_group.add_argument(
        "--reset",
        dest="reset",
        action="store_true",
        help="Reset CRISPRitz configuration to default values. This will restore "
        "environment name to 'crispritz' and targets directory to "
        "'.crispritz_targets'. Use with caution as this will overwrite current "
        "settings",
    )
    action_group.add_argument(
        "--validate",
        dest="validate",
        action="store_true",
        help="Validate current configuration settings. Checks if the specified "
        "environment exists and if the targets directory is accessible. Reports "
        "any configuration issues found",
    )
    return parser_crispritz_config


def main():
    """Entry point for the CRISPR-HAWK command-line interface.

    This function parses command-line arguments and dispatches execution to the
    appropriate CRISPR-HAWK command handler.
    """
    start = time()  # track elapsed time
    try:
        parser = create_parser_crisprhawk()  # parse input argument using custom parser
        if not sys.argv[1:]:  # no input args -> print help and exit
            parser.error_noargs()
        args = parser.parse_args(sys.argv[1:])  # parse input args
        prepare_package()  # check if models and data are available and uncompressed
        if args.command == SEARCH:  # search command
            crisprhawk_search(CrisprHawkSearchInputArgs(args, parser))
        elif args.command == CONVERTGNOMADVCF:  # convert-gnoamd-vcf command
            crisprhawk_converter(CrisprHawkConverterInputArgs(args, parser))
        elif args.command == PREPAREDATACRISPRME:  # prepare-data-crisprme command
            crisprhawk_prepare_data_crisprme(
                CrisprHawkPrepareDataInputArgs(args, parser)
            )
        elif args.command == CRISPRPITZCONFIG:  # crispritz-configure command
            crisprhawk_crispritz_config(
                CrisprHawkCrispritzConfigInputArgs(args, parser)
            )
    except KeyboardInterrupt:
        sigint_handler()  # catch SIGINT and exit gracefully
    sys.stdout.write(f"{TOOLNAME} - Elapsed time {(time() - start):.2f}s\n")


# --------------------------------> ENTRY POINT <--------------------------------
if __name__ == "__main__":
    main()
