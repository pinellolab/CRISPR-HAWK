"""
"""

from bedfile import Bed, Region, IndelRegion, RegionList
from utils import print_verbosity, warning, VERBOSITYLVL
from crisprhawk_argparse import CrisprHawkInputArgs
from enrichment import enricher, load_vcfs
from haplotypes import track_haplotypes
from variant_map import VariantMap
from reports import report_guides
from search_guides import search
from sequences import Fasta, PAM
from guide import Guide

from typing import List, Dict, Tuple, Union
from time import time


def construct_regions(
    bedfile: str,
    fastafile: str,
    fasta_idx: str,
    guidelen: int,
    verbosity: int,
    debug: bool,
) -> RegionList:
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
    print_verbosity(
        f"Regions extracted in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return regions  # regions with associated genomic sequence (see Region)


def enrichment(
    regions: RegionList,
    vcflist: List[str],
    guidelen: int,
    pamlen: int,
    no_filter: bool,
    verbosity: int,
    debug: bool,
) -> Dict[Union[Region, IndelRegion], VariantMap]:
    vcfs = load_vcfs(vcflist, verbosity, debug)  # load input vcfs
    return enricher(regions, vcfs, guidelen, pamlen, no_filter, verbosity, debug)


def encoding(regions: RegionList, verbosity: int) -> RegionList:
    print_verbosity(
        f"Encoding {len(regions)} regions in bits", verbosity, VERBOSITYLVL[2]
    )
    start = time()  # encoding start time
    for region in regions:
        region.encode()  # encode input regions as sequences of bits
    print_verbosity(
        f"Encoding completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return regions


def guide_search(
    pam: PAM,
    regions: List[Region],
    guidelen: int,
    right: bool,
    verbosity: int,
    debug: bool,
) -> Dict[Union[Region, IndelRegion], List[Guide]]:
    print_verbosity(
        f"Begin guide candidates search on {len(regions)} regions",
        verbosity,
        VERBOSITYLVL[3],
    )
    return {
        region: search(pam, region, guidelen, right, verbosity, debug)
        for region in regions
    }


def initialize_pam(pam: str, verbosity: int, debug: bool) -> PAM:
    # construct PAM object
    print_verbosity(f"Creating PAM object for PAM {pam}", verbosity, VERBOSITYLVL[2])
    return PAM(pam, debug)


def crisprhawk_ref(
    regions: RegionList,
    pam: PAM,
    guidelen: int,
    right: bool,
    verbosity: int,
    debug: bool,
) -> Dict[Region, List[Guide]]:
    # encode sequences in bit for efficient candidate guide search
    regions = encoding(regions, verbosity)
    # search guide candidates within input regions
    return guide_search(
        pam, regions, guidelen, right, verbosity, debug
    )  # search guides


def crisprhawk_alt(
    regions: RegionList,
    vcflist: List[str],
    guidelen: int,
    pam: PAM,
    right: bool,
    no_filter: bool,
    verbosity: int,
    debug: bool,
):
    # enrich input regions with SNVs and indels and assess whether input VCFs are phased
    variant_maps = enrichment(
        regions, vcflist, guidelen, len(pam), no_filter, verbosity, debug
    )
    # encode sequences in bit for efficient candidate guide search
    regions = encoding(regions, verbosity)
    # search guide candidates within input regions
    guides = guide_search(pam, regions, guidelen, right, verbosity, debug)
    # track haplotypes on guides
    for region in regions:
        guides[region] = track_haplotypes(region, guides[region], variant_maps[region])
    return guides


def crisprhawk(args: CrisprHawkInputArgs) -> None:
    # extract genomic regions defined in input bed
    regions = construct_regions(
        args.bedfile,
        args.fasta,
        args.fasta_idx,
        args.guidelen,
        args.verbosity,
        args.debug,
    )
    pam = initialize_pam(args.pam, args.verbosity, args.debug)  # initialize pam object
    onlyref = not args.vcfs  # establish whether variants have been given
    if onlyref:  # reference-only data in input, follow the reference only pipeline
        warning("Skipping enrichment (no input VCF)", args.verbosity)
        guides = crisprhawk_ref(
            regions, pam, args.guidelen, args.right, args.verbosity, args.debug
        )
    else:
        guides = crisprhawk_alt(
            regions,
            args.vcfs,
            args.guidelen,
            pam,
            args.right,
            args.no_filter,
            args.verbosity,
            args.debug,
        )
    # report found guides in output directory
    report_guides(guides, args.guidelen, pam, args.outdir, args.right, args.debug)

