"""
"""

from crisprhawk_argparse import CrisprHawkInputArgs
from regions import construct_regions


def crisprhawk(args: CrisprHawkInputArgs) -> None:
    # extract genomic regions defined in input bed file
    regions = construct_regions(args.fasta, args.bedfile, args.fasta_idx, args.guidelen, args.verbosity, args.debug)