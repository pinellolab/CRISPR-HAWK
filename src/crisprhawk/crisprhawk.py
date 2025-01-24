"""
"""

from crisprhawk_argparse import CrisprHawkInputArgs
from utils import print_verbosity, warning, VERBOSITYLVL
from search_guides import search
from bedfile import Bed, Region, IndelRegion, RegionList
from sequences import Fasta, PAM
from reports import report_guides
from enrichment import enricher
from variants import VariantRecord
from haplotypes import track_haplotypes
from guide import Guide

from typing import List, Dict, Tuple
from time import time


def construct_regions(bedfile: str, fastafile: str, fasta_idx: str, guidelen: int, verbosity: int, debug: bool) -> RegionList:
    # read input bed and construct Region object for each genomic region
    print_verbosity(f"Parsing input BED file {bedfile}", verbosity, VERBOSITYLVL[2])
    start = time()  # track processing time
    bed = Bed(bedfile, guidelen, debug)  # guidelen used to pad regions
    print_verbosity(f"Parsed regions number: {len(bed)}", verbosity, VERBOSITYLVL[3])
    # extract genomic sequences foe each input region
    print_verbosity(f"Extracting regions from {fastafile}", verbosity, VERBOSITYLVL[2])
    regions = bed.extract(Fasta(fastafile, verbosity, debug, fasta_idx)) 
    print_verbosity(
        f"Extracted regions:\n{regions.format(pad=guidelen, string=True)}",
        verbosity,
        VERBOSITYLVL[3],
    )
    print_verbosity(f"Regions extracted in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3])
    return regions  # regions with associated genomic sequence (see Region)


def enrichment(
    fasta: str,
    bedfile: str,
    fasta_idx: str,
    vcfs: List[str],
    guidelen: int,
    pamlen: int,
    onlyref: bool,
    no_filter: bool,
    verbosity: int,
    debug: bool,
) -> Tuple[RegionList, Dict[Region, Dict[int, VariantRecord]], bool]:
    # read and parse regions from input bedfile
    print_verbosity(f"Parsing input BED file {bedfile}", verbosity, VERBOSITYLVL[2])
    bed = Bed(bedfile, guidelen, debug)  
    print_verbosity(f"Parsed regions number: {len(bed)}", verbosity, VERBOSITYLVL[3])
    # extract regions sequences for subsequent steps
    print_verbosity(f"Extracting regions from {fasta}", verbosity, VERBOSITYLVL[2])
    regions = bed.extract(Fasta(fasta, verbosity, debug, fasta_idx))  # extract regions
    print_verbosity(
        f"Extracted regions:\n{regions.format(pad=guidelen, string=True)}",
        verbosity,
        VERBOSITYLVL[3],
    )
    if onlyref:  # no variants in input, skip enrichment and go to encoding
        warning("Skipping enrichment (no input VCF)", verbosity)
        return regions, None, None
    # variants in input -> proceed with genome enrichment
    return enricher(regions, vcfs, guidelen, pamlen, no_filter, verbosity, debug)


def encoding(regions: RegionList, verbosity: int) -> RegionList:
    print_verbosity(f"Encoding {len(regions)} regions in bits", verbosity, VERBOSITYLVL[2])
    start = time()  # encoding start time
    for region in regions:
        region.encode()  # encode input regions as sequences of bits
    print_verbosity(f"Encoding completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3])
    return regions


def guide_search(
    pam: PAM, regions: List[Region], guidelen: int, right: bool, verbosity: int, debug: bool
) -> Dict[Region, List[Guide]]:
    print_verbosity(f"Begin guide candidates search on {len(regions)} regions", verbosity, VERBOSITYLVL[3])
    return {region: search(pam, region, guidelen, right, verbosity, debug) for region in regions}

def initialize_pam(pam: str, verbosity: int, debug: bool) -> PAM:
    # construct PAM object 
    print_verbosity(f"Creating PAM object for PAM {pam}", verbosity, VERBOSITYLVL[2])
    return PAM(pam, debug)

def crisprhawk_ref(regions: RegionList, pam: PAM, guidelen: int, right: bool, outdir: str, verbosity: int, debug: bool) -> Dict[Region, List[Guide]]:
    # encode sequences in bit for efficient candidate guide search
    regions = encoding(regions, verbosity)  
    return guide_search(pam, regions, guidelen, right, verbosity, debug)  # search guides
    

def crisprhawk(args: CrisprHawkInputArgs) -> None:
    # extract genomic regions defined in input bed 
    regions = construct_regions(args.bedfile, args.fasta, args.fasta_idx, args.guidelen, args.verbosity, args.debug)
    pam = initialize_pam(args.pam, args.verbosity, args.debug)  # initialize pam object
    onlyref = not args.vcfs  # establish whether variants have been given
    if onlyref:  # reference-only data in input, follow the reference only pipeline
        guides = crisprhawk_ref(regions, pam, args.guidelen, args.right, args.outdir, args.verbosity, args.debug)
    else:
        pass
    # report found guides in output directory
    report_guides(guides, args.guidelen, pam, args.outdir, args.right, args.debug)

    # # sequence enrichment -> add genetic variants to input sequences if input
    # # vcf is given, return reference sequences otherwise
    # regions, variants_maps, phased = enrichment(
    #     args.fasta,
    #     args.bedfile,
    #     args.fasta_idx,
    #     args.vcfs,
    #     args.guidelen,
    #     len(args.pam),
    #     onlyref,
    #     args.no_filter,
    #     args.verbosity,
    #     args.debug,
    # )
    
    # regions = encoding(regions)
    # # search guides in the input regions
    # pam = PAM(args.pam, args.debug)  # initialize pam
    # guides = guide_search(pam, regions, args.guidelen, args.right, args.debug)
    # if not onlyref: # track haplotypes on guides
    #     for r in guides:  # update current region guides data
    #         guides[r] = track_haplotypes(
    #             guides[r], variants_maps[r], args.guidelen, len(pam), phased, r, args.debug
    #         )
    # # report guides in output directory
    # report_guides(guides, args.outdir, pam, args.right, args.guidelen, args.debug)
