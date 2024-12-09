"""
"""

from crisprhawk_argparse import CrisprHawkInputArgs
from utils import print_verbosity, warning, VERBOSITYLVL
from search_guides import search
from bedfile import Bed, Region, RegionList
from sequences import Fasta, PAM
from reports import report_guides
from enrichment import enricher
from variants import VariantRecord
from haplotypes import track_haplotypes

from typing import List, Dict, Tuple


def enrichment(
    fasta: str,
    bedfile: str,
    fasta_idx: str,
    vcfs: List[str],
    guidelen: int,
    no_filter: bool,
    verbosity: int,
    debug: bool,
) -> Tuple[RegionList, Dict[Region, Dict[int, VariantRecord]], bool]:
    print_verbosity(f"Parsing input BED file {bedfile}", verbosity, VERBOSITYLVL[2])
    bed = Bed(bedfile, guidelen, debug)  # read regions from input bedfile
    print_verbosity(f"Parsed regions number: {len(bed)}", verbosity, VERBOSITYLVL[3])
    print_verbosity(f"Extracting regions from {fasta}", verbosity, VERBOSITYLVL[2])
    regions = bed.extract(Fasta(fasta, verbosity, debug, fasta_idx))  # extract regions
    print_verbosity(
        f"Extracted regions:\n{regions.format(pad=guidelen)}",
        verbosity,
        VERBOSITYLVL[3],
    )
    if not vcfs:  # no variants in input, skip enrichment and go to encoding
        warning("Skipping enrichment (no input VCF)", verbosity)
        return regions
    # TODO: enrichment - haplotype tracking + indels
    return enricher(regions, vcfs, guidelen, no_filter, verbosity, debug)


def encoding(regions: RegionList) -> RegionList:
    for region in regions:
        region.encode()  # encode input regions as sequences of bits
    return regions


def guide_search(
    pam: PAM, regions: List[Region], guidelen: int, right: bool, debug: bool
) -> Dict[
    Region, Tuple[Tuple[List[int], List[int]], Tuple[List[List[str]], List[List[str]]]]
]:
    return {region: search(pam, region, guidelen, right, debug) for region in regions}


def crisprhawk(args: CrisprHawkInputArgs) -> None:
    # sequence enrichment -> add genetic variants to input sequences if input
    # vcf is given, return reference sequences otherwise
    regions, variants_maps, phased = enrichment(
        args.fasta,
        args.bedfile,
        args.fasta_idx,
        args.vcfs,
        args.guidelen,
        args.no_filter,
        args.verbosity,
        args.debug,
    )
    # encode sequences in bit for efficient candidate guide search
    regions = encoding(regions)
    # search guides in the input regions
    pam = PAM(args.pam, args.debug)  # initialize pam
    guides = guide_search(pam, regions, args.guidelen, args.right, args.debug)
    # track haplotypes on guides
    for r in guides:  # update current region guides data
        guides[r] = track_haplotypes(guides[r], variants_maps[r], args.guidelen, phased)
    # report guides in output directory
    for region, (positions, guides, samples) in guides.items():
        report_guides(
            args.outdir,
            region,
            guides,
            positions,
            samples,
            pam,
            args.right,
            args.guidelen,
            args.debug,
        )
