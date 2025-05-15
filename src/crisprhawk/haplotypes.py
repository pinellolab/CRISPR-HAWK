""" """

from exception_handlers import exception_handler
from utils import print_verbosity, flatten_list, VERBOSITYLVL
from region import Region, RegionList
from variant import VCF, VariantRecord, VTYPES
from coordinate import Coordinate
from sequence import Sequence
from haplotype import Haplotype


from hapsolver import (
    HaplotypeGraph,
)
from typing import List, Dict, Tuple
from time import time

import os




    
    

    
        


def read_vcf(vcflist: List[str], verbosity: int, debug: bool) -> Dict[str, VCF]:
    """Load VCF files and map each VCF to its contig.

    Reads a list of VCF file paths, loads each as a VCF object, and returns a 
    dictionary mapping contig names to their corresponding VCF objects. Assumes 
    one VCF per contig.

    Args:
        vcflist: List of VCF file paths.
        verbosity: The verbosity level for logging.
        debug: Whether to enable debug mode for exception handling.

    Returns:
        A dictionary mapping contig names to VCF objects.

    Raises:
        Exception: If parsing any VCF file fails.
    """
    # load vcf files and map each vcf to its contig (assume on vcf per contig)
    print_verbosity("Loading VCF files", verbosity, VERBOSITYLVL[3])
    start = time()  # track vcf parsing time
    try:  # create vcf dictionary
        vcfs = {vcf.contig: vcf for vcf in [VCF(f, verbosity, debug) for f in vcflist]}
    except FileNotFoundError as e:
        exception_handler(
            Exception, "Failed parsing VCF files", os.EX_DATAERR, debug, e
        )
    print_verbosity(
        f"Loaded {len(vcfs)} VCFs in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return vcfs


def fetch_variants(
    vcfs: Dict[str, VCF], regions: RegionList, verbosity: int, debug: bool
) -> Dict[Region, List[VariantRecord]]:
    """Fetch variants mapped to each query region.

    Iterates over the provided regions and retrieves variant records from the 
    corresponding VCFs, returning a dictionary mapping each region to its list 
    of variants.

    Args:
        vcfs: Dictionary mapping contig names to VCF objects.
        regions: List of regions to fetch variants for.
        verbosity: The verbosity level for logging.
        debug: Whether to enable debug mode for exception handling.

    Returns:
        A dictionary mapping each region to a list of VariantRecord objects.

    Raises:
        Exception: If fetching variants fails.
    """
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
    except Exception as e:
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

def compute_haplotypes_phased(variants: List[VariantRecord], samples: List[str]):
    # initialize sample-variant map for both copies
    sample_variants = {s: ([], []) for s in samples}
    for variant in variants:
        assert len(variant.samples) == 1
        for chromcopy in [0, 1]:  # iterate over chromosome copies
            for s in variant.samples[0][chromcopy]:  # add variant to sample-variant map
                sample_variants[s][chromcopy].append(variant)
    return sample_variants

def resolve_haplotypes(haplotypes: Dict[str, Tuple[List[VariantRecord], List[VariantRecord]]], refseq: str, region: Region, phased: bool):
    for sample in haplotypes:
        for i in [0,1]:
            h = Haplotype(refseq, region.coordinates, phased, i)


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

    for region in regions:
        samples = vcfs[region.contig].samples
        if phased:
            haps =  compute_haplotypes_phased(variants[region], samples)
            resolve_haplotypes(haps, region.sequence.sequence, region, phased)
            

    exit()

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
