"""
"""

from crisprhawk_argparse import CrisprHawkArgumentParser
from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkError
from bedfile import Bed
from fasta import Fasta

from argparse import Namespace

import os


def enrichment(
    fasta: str,
    bedfile: str,
    fasta_idx: str,
    debug: bool,
    parser: CrisprHawkArgumentParser,
):
    # check fasta file existance and content
    if not os.path.isfile(fasta):
        parser.error(f"Unable to find {fasta}")
    if os.stat(fasta).st_size <= 0:
        parser.error(f"The input FASTA file {fasta} seems empty")
    if not os.path.isfile(bedfile):
        parser.error(f"Unable to find {bedfile}")
    if os.stat(bedfile).st_size <= 0:
        parser.error(f"The input BED file {bedfile} seems empty")
    # create fasta and bedfile objects to represent the input fasta and regions
    fasta = Fasta(fasta, debug, fasta_idx)
    bedfile = Bed(bedfile, debug)
    # extract regions from


def crisprhawk(args: Namespace, parser: CrisprHawkArgumentParser) -> None:
    # sequence enrichment -> add genetic variants to input sequences
    enrichment(args.fasta, args.bedfile, args.fasta_idx, args.debug, parser)
    exception_handler(CrisprHawkError, "error raised", os.EX_USAGE, False)
