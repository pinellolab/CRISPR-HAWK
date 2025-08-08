"""
CRISPR-HAWK {version}

Copyright (C) 2025 Manuel Tognon <manu.tognon@gmail.com> <manuel.tognon@univr.it> <mtognon@mgh.harvard.edu>

Usage:
    crisprhawk -f <fasta> -r <bedfile> -v <vcf>

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
    parser = CrisprHawkArgumentParser(usage=__doc__, add_help=False)
    group = parser.add_argument_group("Options")  # arguments group
    # input arguments
    group.add_argument("-h", "--help", action="help", help="Show this message and exit")
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
        help="Input FASTA file",
    )
    group.add_argument(
        "-i",
        "--fasta-idx",
        type=str,
        metavar="FASTA-IDX",
        dest="fasta_idx",
        nargs="?",
        default="",
        help="Fasta index (FAI), indexing the input fastafile",
    )
    group.add_argument(
        "-r",
        "--regions",
        type=str,
        metavar="GENOMIC-REGIONS-BED",
        dest="bedfile",
        required=True,
        help="BED file specifying genomic regions to search for guide RNAs",
    )
    group.add_argument(
        "-v",
        "--vcf",
        type=str,
        metavar="VCF-FILE",
        dest="vcf",
        nargs="?",
        default="",
        help="VCF file containing genetic variants to account for during guide RNA search",
    )
    group.add_argument(
        "-p",
        "--pam",
        type=str,
        metavar="PAM",
        dest="pam",
        required=True,
        help="PAM sequence used during candidate guide nomination",
    )
    group.add_argument(
        "-g",
        "--guide-len",
        type=int,
        metavar="GUIDE-LENGTH",
        dest="guidelen",
        required=True,
        help="Length of guides to nominate",
    )
    group.add_argument(
        "--right",
        action="store_true",
        dest="right",
        default=False,
        help="If selected, extract guides downstream PAM matching positions",
    )
    group.add_argument(
        "-o",
        "--outdir",
        type=str,
        metavar="OUTDIR",
        dest="outdir",
        nargs="?",
        default=os.getcwd(),
        help="Path to output directory, by default the reports are stored in the "
        "current working directory",
    )
    group.add_argument(
        "--no-filter",
        action="store_true",
        dest="no_filter",
        default=False,
        help="When enabled, all variants in the input VCF files will be "
        "considered. By default, only variants marked with 'PASS' in the FILTER "
        "field are processed, while others are skipped",
    )
    group.add_argument(
        "--functional-annotation",
        type=str,
        metavar="FUNC-ANN-BED",
        dest="functional_annotation",
        nargs="?",
        default="",
        help="BED file specifying functional genomic regions to functionally "
        "annotate guide candidates"
    )
    group.add_argument(
        "--gene-annotation",
        type=str,
        metavar="GENE-ANN-BED",
        dest="gene_annotation",
        nargs="?",
        default="",
        help="BED file specifying gene regions to annotate guide candidates"
    )
    group.add_argument(
        "--haplotype-table",
        action="store_true",
        dest="haplotype_table",
        default=False,
        help="When enabled, the haplotype table is returned in the output folder "
        "as TSV file. By default, the haplotype table is not returned",
    )
    group.add_argument(
        "--estimate-offtargets",
        action="store_true",
        dest="estimate_offtargets",
        default=False,
        help="When enabled, the off-targets are estimated for each guide RNA "
        "candidate. By default, off-targets are not estimated",
    )
    group.add_argument(
        "--write-offtargets-report",
        action="store_true",
        dest="write_offtargets_report",
        default=False,
        help="When enabled, write a report for all off-targets found for each "
        "guide RNA candidate. By default, off-targets are not reported"
    )
    group.add_argument(
        "--verbosity",
        type=int,
        metavar="VERBOSITY",
        dest="verbosity",
        nargs="?",
        default=1,  # minimal output
        help="Sets the level of detail in the output messages. Available levels "
        "are: 0 (Silent), 1 (Default), 2 (Verbose), 3 (Debug)",
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
