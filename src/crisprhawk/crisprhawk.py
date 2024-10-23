"""
"""

from crisprhawk_argparse import CrisprHawkArgumentParser
from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkError
from search_guides import search
from bedfile import Bed, Region
from sequences import Fasta, PAM

from typing import List
from argparse import Namespace

import os


def enrichment(
    fasta: str,
    bedfile: str,
    fasta_idx: str,
    vcf: str, 
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
    bedfile = Bed(bedfile, debug)  # read regions from input bedfile 
    regions = bedfile.extract(Fasta(fasta, debug, fasta_idx))  # extract regions
    if not vcf:  # no variants in input, skip enrichment and go to encoding
        return regions
    
    # TODO: enrichment
    print("enrichment") 

def encoding(regions: List[Region]) -> List[Region]:
    for region in regions:
        region.encode()  # encode input regions as sequences of bits 
    return regions


def guide_search(pamseq: str, regions: List[Region], guidelen: int, right: bool, debug: bool):
    pam = PAM(pamseq, debug)
    for region in regions:
        search(pam, region, guidelen, right, debug)
   


def crisprhawk(args: Namespace, parser: CrisprHawkArgumentParser) -> None:
    # sequence enrichment -> add genetic variants to input sequences if input 
    # vcf is given, return reference sequences otherwise
    regions = enrichment(args.fasta, args.bedfile, args.fasta_idx, args.vcf, args.debug, parser)
    # encode sequences in bit for efficient candidate guide search
    regions = encoding(regions)
    # search guides in the input regions
    guide_search(args.pam, regions, args.guidelen, args.right, args.debug)
