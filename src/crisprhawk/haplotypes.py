""" """

from exception_handlers import exception_handler
from utils import print_verbosity, flatten_list, VERBOSITYLVL
from region import Region, RegionList
from variant import VCF, VariantRecord, VTYPES
from coordinate import Coordinate
from sequence import Sequence
from haplotype import Haplotype

from typing import List, Dict, Tuple, Optional
from collections import defaultdict
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



def initialize_haplotypes(regions: RegionList, debug: bool) -> Dict[Region, List[Haplotype]]:
    # initialize haplotype object with REF haplotype
    return {r: [Haplotype(Sequence(r.sequence.sequence, debug), r.coordinates, False, 0, debug)] for r in regions}


def compute_haplotypes_phased(variants: List[VariantRecord], samples: List[str]):
    # initialize sample-variant map for both copies
    sample_variants = {s: ([], []) for s in samples}
    for variant in variants:
        assert len(variant.samples) == 1
        for chromcopy in [0, 1]:  # iterate over chromosome copies
            for s in variant.samples[0][chromcopy]:  # type: ignore 
                # add variant to sample-variant map
                sample_variants[s][chromcopy].append(variant)
    return sample_variants

def ishomozygous(haplotypes: List[Haplotype]) -> bool:
    return len({h.sequence.sequence for h in haplotypes}) == 1

def _solve_haplotypes(sequence: str, coordinates: Coordinate, phased: bool, variants: Tuple[List[VariantRecord], List[VariantRecord]], sample: str, debug: bool) -> List[Haplotype]:
    # solve haplotypes for diploid samples
    h0 = Haplotype(Sequence(sequence, debug), coordinates, phased, 0, debug)  # first copy
    h0.add_variants(variants[0], sample)  # add variants to haplotype
    h1 = Haplotype(Sequence(sequence, debug), coordinates, phased, 1, debug)  # second copy
    h1.add_variants(variants[1], sample)  # add variants to haplotype
    # check for homozygous haplotypes
    if ishomozygous([h0, h1]):
        h0.homozygous_samples()
        return [h0]  # return only one haplotype (homozygous sample)
    return [h0, h1]  # return both haplotypes

def _collapse_haplotypes(sequence: str, haplotypes: List[Haplotype], debug: bool) -> Haplotype:
    hap = Haplotype(Sequence(sequence, debug), haplotypes[0].coordinates, haplotypes[0].phased, 0, debug)
    hap.set_samples(",".join([h.samples for h in haplotypes]))
    hap.set_variants(",".join([h.variants for h in haplotypes]))
    return hap
    

def collapse_haplotypes(haplotypes: List[Haplotype], debug: bool) -> List[Haplotype]:
    haplotypes_dict = [(h.sequence.sequence, h) for h in haplotypes]
    haplotypes_collapsed = defaultdict(list)
    for seq, hap in haplotypes_dict:
        haplotypes_collapsed[seq].append(hap)
    return [_collapse_haplotypes(seq, haplist, debug) for seq, haplist in haplotypes_collapsed.items()]

def solve_haplotypes(sample_variants: Dict[str, Tuple[List[VariantRecord], List[VariantRecord]]], hapseqs: List[Haplotype], refseq: str, coordinates: Coordinate, phased: bool, debug: bool) -> List[Haplotype]:
    # solve haplotypes for each sample (assumes diploid samples)
    for sample, variants in sample_variants.items():
        hapseqs += _solve_haplotypes(refseq, coordinates, phased, variants, sample, debug)
    return hapseqs

def add_variants(vcflist: List[str], regions: RegionList, haplotypes: Dict[Region, List[Haplotype]], verbosity: int, debug: bool) -> Dict[Region, List[Haplotype]]:
    # read VCF files and extract variants located within the region
    vcfs = read_vcf(vcflist, verbosity, debug)  
    variants = fetch_variants(vcfs, regions, verbosity, debug)  
    phased = vcfs[regions[0].contig].phased  # assess VCF phasing
    for region in regions:  # reconstruct haplotypes for each region
        if phased:  # phased VCFs
            # recover variants combinations on each chromosome copy
            samples_variants = compute_haplotypes_phased(variants[region], vcfs[regions[0].contig].samples)
            # solve haplotypes for each sample
            haplotypes[region] = solve_haplotypes(samples_variants, haplotypes[region], region.sequence.sequence, region.coordinates, phased, debug)
    return haplotypes

def reconstruct_haplotypes(
    vcflist: List[str], regions: RegionList, verbosity: int, debug: bool
) -> Dict[Region, List[Haplotype]]:
    # read input vcf files and fetch variants in each region
    print_verbosity("Reconstructing haplotypes", verbosity, VERBOSITYLVL[1])
    start = time()  # track haplotypes reconstruction time
    # initialize haplotypes list with reference sequence haplotype
    haplotypes = initialize_haplotypes(regions, debug)
    if vcflist:  # add variants to regions and solve haplotypes
        haplotypes = add_variants(vcflist, regions, haplotypes, verbosity, debug)
    print_verbosity(
        f"Haplotypes reconstructed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )
    return haplotypes
