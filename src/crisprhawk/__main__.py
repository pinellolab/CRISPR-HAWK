"""
CRISPR-HAWK {version}

Copyright (C) 2025 Manuel Tognon <manu.tognon@gmail.com> <manuel.tognon@univr.it> <mtognon@mgh.harvard.edu>

CRISPR-HAWK: Haplotype- and vAriant-aWare guide design toolKit

Usage:
    crisprhawk -f <fasta> -r <bedfile> -v <vcf> -p <pam> -g <guide-length> -o <output-dir>

Run 'crisprhawk -h/--help' to display the complete help
"""

from .crisprhawk_argparse import CrisprHawkArgumentParser, CrisprHawkInputArgs
from .exception_handlers import sigint_handler
from .crisprhawk_version import __version__
from .crisprhawk import crisprhawk
from .utils import TOOLNAME

from time import time

import sys
import os


def create_parser_crisprhawk() -> CrisprHawkArgumentParser:
    # force displaying docstring at each usage display and force
    # the default help to not being shown
    parser = CrisprHawkArgumentParser(usage=__doc__, add_help=False)  # type: ignore
    group = parser.add_argument_group("Options")  # arguments group
    # input arguments
    group.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit"
    )
    group.add_argument(
        "--version",
        action="version",
        help=f"Show {TOOLNAME} version and exit",
        version=__version__,
    )
    group.add_argument(
        "-f",
        "--fasta",
        type=str,
        metavar="FASTA-FILE",
        dest="fasta",
        required=True,
        help="Reference genome in FASTA format used for guide search",
    )
    group.add_argument(
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
    group.add_argument(
        "-r",
        "--regions",
        type=str,
        metavar="GENOMIC-REGIONS-BED",
        dest="bedfile",
        required=True,
        help="BED file specifying genomic regions where guides will be searched",
    )
    group.add_argument(
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
    group.add_argument(
        "-p",
        "--pam",
        type=str,
        metavar="PAM",
        dest="pam",
        required=True,
        help="PAM sequence used to identify candidate guides (e.g., NGG, NAG, " "etc.)",
    )
    group.add_argument(
        "-g",
        "--guide-len",
        type=int,
        metavar="GUIDE-LENGTH",
        dest="guidelen",
        required=True,
        help="Length of the guide (excluding the PAM)",
    )
    group.add_argument(
        "--right",
        action="store_true",
        dest="right",
        default=False,
        help="If set, guides are extracted downstream (right side) of the PAM "
        "site. (default: guides are extracted upstream (left side))",
    )
    group.add_argument(
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
    group.add_argument(
        "--no-filter",
        action="store_true",
        dest="no_filter",
        default=False,
        help="If set, all variants in the input VCF file will be considered "
        "regardless of FILTER status (default: only variants with FILTER == "
        "'PASS' are used)",
    )
    group.add_argument(
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
    group.add_argument(
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
    group.add_argument(
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
    group.add_argument(
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
    group.add_argument(
        "--haplotype-table",
        action="store_true",
        dest="haplotype_table",
        default=False,
        help="When enabled, the haplotype table is returned in the output folder "
        "as TSV file (default: disabled)",
    )
    # group.add_argument(
    #     "--estimate-offtargets",
    #     action="store_true",
    #     dest="estimate_offtargets",
    #     default=False,
    #     help="When enabled, the off-targets are estimated for each guide RNA "
    #     "candidate ()",
    # )
    # group.add_argument(
    #     "--write-offtargets-report",
    #     action="store_true",
    #     dest="write_offtargets_report",
    #     default=False,
    #     help="When enabled, write a report for all off-targets found for each "
    #     "guide RNA candidate. By default, off-targets are not reported",
    # )
    group.add_argument(
        "-t",
        "--threads",
        type=int,
        metavar="THREADS",
        dest="threads",
        nargs="?",
        default=1,
        help="Number of threads. Use 0 for using all available cores (default: 1)",
    )
    group.add_argument(
        "--verbosity",
        type=int,
        metavar="VERBOSITY",
        dest="verbosity",
        nargs="?",
        default=1,  # minimal output
        help="Verbosity level of output messages: 0 = Silent, 1 = Normal, 2 = "
        "Verbose, 3 = Debug (default: 1)",
    )
    group.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Enter debug mode and trace the full error stack",
    )
    return parser


def main():
    start = time()  # track elapsed time
    try:
        parser = create_parser_crisprhawk()  # parse input argument using custom parser
        if not sys.argv[1:]:  # no input args -> print help and exit
            parser.error_noargs()
        crisprhawk(CrisprHawkInputArgs(parser.parse_args(sys.argv[1:]), parser))  # type: ignore
    except KeyboardInterrupt as e:
        sigint_handler()  # catch SIGINT and exit gracefully
    sys.stdout.write(f"{TOOLNAME} - Elapsed time {(time() - start):.2f}s\n")


# --------------------------------> ENTRY POINT <--------------------------------
if __name__ == "__main__":
    main()
