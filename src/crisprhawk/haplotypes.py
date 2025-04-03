""" """

from exception_handlers import exception_handler
from utils import print_verbosity, VERBOSITYLVL

from hapsolver import (
    VCF,
    VariantRecord,
    Region,
    RegionList,
    HaplotypeGraph,
    Coordinate,
    Haplotype,
    Sequence,
    flatten_list,
)
from typing import List, Dict
from time import time

import os


def read_vcf(vcflist: List[str], verbosity: int, debug: bool) -> Dict[str, VCF]:
    # load vcf files and map each vcf to its contig (assume on vcf per contig)
    print_verbosity("Loading VCF files", verbosity, VERBOSITYLVL[3])
    start = time()  # track vcf parsing time
    try:  # create vcf dictionary
        vcfs = {vcf.contig: vcf for vcf in [VCF(f) for f in vcflist]}
    except FileNotFoundError as e:
        exception_handler(
            Exception, f"Failed parsing VCF files", os.EX_DATAERR, debug, e
        )
    print_verbosity(
        f"Loaded {len(vcfs)} VCFs in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return vcfs


def fetch_variants(
    vcfs: Dict[str, VCF], regions: RegionList, verbosity: int, debug: bool
) -> Dict[Region, List[VariantRecord]]:
    # recover variants mapped on the query region
    print_verbosity("Fetching variants", verbosity, VERBOSITYLVL[3])
    start = time()  # track variants fetching time
    try:  # fecth variants in each region
        variants = {
            region: flatten_list(
                [v.split() for v in vcfs[region.contig].fetch(region.coordinates)]
            )
            for region in regions
        }
    except (ValueError, IndexError, Exception) as e:
        exception_handler(
            Exception, "Failed fecthing variants", os.EX_DATAERR, debug, e
        )
    print_verbosity(
        f"Fetched {sum(len(v) for v in variants.values())} variants in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return variants


def construct_hapgraph(
    variants: List[VariantRecord],
    coordinate: Coordinate,
    chromcopy: int,
    verbosity: int,
    debug: bool,
) -> HaplotypeGraph:
    # construct haplotype graph from input variants
    print_verbosity(
        f"Constructing haplotype graph for region {coordinate}",
        verbosity,
        VERBOSITYLVL[3],
    )
    start = time()  # track haplotype graph construction time
    try:  # create haplotype graph
        hapgraph = HaplotypeGraph(variants, chromcopy)
        hapgraph.compute_graph()  # populate haplotype graph
    except ValueError as e:
        exception_handler(
            Exception,
            f"Failed constructing haplotype graph for region {coordinate}",
            os.EX_DATAERR,
            debug,
            e,
        )
    print_verbosity(
        f"Haplotype graph constructed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return hapgraph


def initialize_haplotype(
    sequence: str, coord: Coordinate, phased: bool, chromcopy: int
) -> Haplotype:
    # initialize haplotype object
    return Haplotype(Sequence(sequence), coord, phased, chromcopy)


def retrieve_haplotypes(
    hapgraph: HaplotypeGraph,
    refseq: str,
    coordinate: Coordinate,
    phased: bool,
    chromcopy: int,
    verbosity: int,
    debug: bool,
) -> List[Haplotype]:
    # retrieve haplotypes from haplotype graph
    print_verbosity(
        f"Retrieving haplotypes for region {coordinate}", verbosity, VERBOSITYLVL[3]
    )
    start = time()  # track haplotype retrieval time
    haplotypes = []  # list of haplotypes
    try:  # retrieve haplotypes
        for variants, samples in hapgraph.retrieve_haplotypes():
            hap = initialize_haplotype(refseq, coordinate, phased, chromcopy)
            hap.add_variants(variants, samples)  # add variants to haplotype
            haplotypes.append(hap)  # add haplotype to haplotypes list
    except ValueError as e:
        exception_handler(
            Exception, "Failed retrieving haplotypes", os.EX_DATAERR, debug, e
        )
    print_verbosity(
        f"Haplotypes retrieved in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return haplotypes


def reconstruct_haplotypes(
    vcflist: List[str], regions: RegionList, verbosity: int, debug: bool
) -> Dict[Region, List[Haplotype]]:
    # read input vcf files and fetch variants in each region
    print_verbosity("Reconstructing haplotypes", verbosity, VERBOSITYLVL[1])
    start = time()  # track haplotypes reconstruction time
    vcfs = read_vcf(vcflist, verbosity, debug)
    variants = fetch_variants(vcfs, regions, verbosity, debug)
    # reconstruct haplotypes for each region
    phased = vcfs[regions[0].contig].phased  # assess VCF phasing
    chromcopies = [0, 1] if phased else [0]
    # initialize haplotypes list with reference sequence haplotype
    haplotypes = {r: [] for r in regions}
    for region in regions:
        # region haplotypes initialized with reference haplotype
        region_haps = [
            initialize_haplotype(
                region.sequence.sequence, region.coordinates, phased, 0
            )
        ]
        for chromcopy in chromcopies:  # iterate over chromosome copies
            hapgraph = construct_hapgraph(
                variants[region], region.coordinates, chromcopy, verbosity, debug
            )
            region_haps.extend(
                retrieve_haplotypes(
                    hapgraph,
                    region.sequence.sequence,
                    region.coordinates,
                    phased,
                    chromcopy,
                    verbosity,
                    debug,
                )
            )
        haplotypes[region] = region_haps
    print_verbosity(
        f"Haplotypes reconstructed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )
    return haplotypes


def reconstruct_haplotypes_ref(
    regions: RegionList, verbosity: int, debug: bool
) -> Dict[Region, List[Haplotype]]:
    # initialize haplotypes list with reference sequence haplotype
    print_verbosity("Reconstructing haplotypes (REF only)", verbosity, VERBOSITYLVL[1])
    start = time()  # track haplotype reconstruction time
    haplotypes = {
        r: [initialize_haplotype(r.sequence.sequence, r.coordinates, False, 0)]
        for r in regions
    }
    print_verbosity(
        f"Haplotypes reconstructed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )
    return haplotypes
