"""
"""

from crisprhawk_argparse import CrisprHawkArgumentParser
from search_guides import search
from bedfile import Bed, Region
from sequences import Fasta, PAM
from reports import report_guides

from typing import List, Dict, Tuple
from argparse import Namespace

import os


def enrichment(
    fasta: str,
    bedfile: str,
    fasta_idx: str,
    vcf: str,
    guidelen: int,
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
    bedfile = Bed(bedfile, guidelen, debug)  # read regions from input bedfile
    regions = bedfile.extract(Fasta(fasta, debug, fasta_idx))  # extract regions
    if not vcf:  # no variants in input, skip enrichment and go to encoding
        return regions
    # TODO: enrichment
    print("enrichment")


def encoding(regions: List[Region]) -> List[Region]:
    for region in regions:
        region.encode()  # encode input regions as sequences of bits
    return regions


def guide_search(
    pam: PAM, regions: List[Region], guidelen: int, right: bool, debug: bool
) -> Dict[
    Region, Tuple[Tuple[List[int], List[int]], Tuple[List[List[str]], List[List[str]]]]
]:
    return {region: search(pam, region, guidelen, right, debug) for region in regions}


def crisprhawk(args: Namespace, parser: CrisprHawkArgumentParser) -> None:
    # sequence enrichment -> add genetic variants to input sequences if input
    # vcf is given, return reference sequences otherwise
    regions = enrichment(
        args.fasta,
        args.bedfile,
        args.fasta_idx,
        args.vcf,
        args.guidelen,
        args.debug,
        parser,
    )
    # encode sequences in bit for efficient candidate guide search
    regions = encoding(regions)
    # search guides in the input regions
    pam = PAM(args.pam, args.debug)  # initialize pam
    guides = guide_search(pam, regions, args.guidelen, args.right, args.debug)
    # report guides in output directory
    for region, (positions, guides) in guides.items():
        report_guides(
            args.outdir,
            region,
            guides,
            positions,
            pam,
            args.right,
            args.guidelen,
            args.debug,
        )
