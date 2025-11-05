"""
This module provides functions for reconstructing haplotypes from VCF files and
genomic regions.

It supports both phased and unphased variant data, enabling the creation, manipulation,
and output of haplotype information for downstream analysis.
"""

from .crisprhawk_argparse import CrisprHawkSearchInputArgs
from .crisprhawk_error import CrisprHawkHaplotypeError
from .exception_handlers import exception_handler
from .utils import print_verbosity, flatten_list, VERBOSITYLVL
from .region import Region, RegionList
from .variant import VCF, VariantRecord, VTYPES
from .coordinate import Coordinate
from .sequence import Sequence
from .haplotype import Haplotype, HaplotypeIndel

from itertools import product
from typing import List, Dict, Tuple, Set
from collections import defaultdict
from time import time

import random
import string
import os

HAPTABCNAMES = ["id", "haplotype", "variants", "samples"]


def read_vcf(vcflist: List[str], verbosity: int, debug: bool) -> Dict[str, VCF]:
    """
    Loads VCF files and maps each VCF to its contig.
    Returns a dictionary mapping contig names to VCF objects for further processing.

    Args:
        vcflist (List[str]): List of VCF file paths.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode.

    Returns:
        Dict[str, VCF]: Dictionary mapping contig names to VCF objects.

    Raises:
        Exception: If any VCF file cannot be parsed.
    """
    # load vcf files and map each vcf to its contig (assume on vcf per contig)
    print_verbosity("Loading VCF files", verbosity, VERBOSITYLVL[3])
    start = time()  # track vcf parsing time
    try:  # create vcf dictionary
        vcfs = {vcf.contig: vcf for vcf in [VCF(f, verbosity, debug) for f in vcflist]}
    except FileNotFoundError as e:
        exception_handler(
            Exception, "Failed parsing VCF files", os.EX_DATAERR, debug, e  # type: ignore
        )
    print_verbosity(
        f"Loaded {len(vcfs)} VCFs in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return vcfs


def fetch_variants(
    vcfs: Dict[str, VCF], regions: RegionList, verbosity: int, debug: bool
) -> Dict[Region, List[VariantRecord]]:
    """
    Retrieves variants mapped to each query region from the provided VCF files.
    Returns a dictionary mapping each region to a list of variant records found
    within that region.

    Args:
        vcfs (Dict[str, VCF]): Dictionary mapping contig names to VCF objects.
        regions (RegionList): List of regions to query for variants.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode.

    Returns:
        Dict[Region, List[VariantRecord]]: Dictionary mapping regions to lists of
            variant records.

    Raises:
        Exception: If fetching variants fails for any region.
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
            Exception, "Failed fecthing variants", os.EX_DATAERR, debug, e  # type: ignore
        )
    print_verbosity(
        f"Fetched {sum(len(v) for v in variants.values())} variants in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return variants


def initialize_haplotypes(
    regions: RegionList, debug: bool
) -> Dict[Region, List[Haplotype]]:
    """
    Initializes haplotype objects for each region using the reference sequence.
    Returns a dictionary mapping each region to a list containing its reference
    haplotype.

    Args:
        regions (RegionList): List of regions for which to initialize haplotypes.
        debug (bool): Flag to enable debug mode.

    Returns:
        Dict[Region, List[Haplotype]]: Dictionary mapping regions to lists of
            haplotypes.
    """
    # initialize haplotype object with REF haplotype
    return {
        r: [
            Haplotype(
                Sequence(r.sequence.sequence, debug), r.coordinates, False, 0, debug
            )
        ]
        for r in regions
    }


def compute_haplotypes_phased(
    variants: List[VariantRecord], samples: List[str]
) -> Dict[str, Tuple[List[VariantRecord], List[VariantRecord]]]:
    """
    Assigns variants to each sample for both chromosome copies in a phased manner.
    Returns a dictionary mapping each sample to a tuple of variant lists for each
    chromosome copy.

    Args:
        variants (List[VariantRecord]): List of variant records to assign.
        samples (List[str]): List of sample names.

    Returns:
        Dict[str, Tuple[List[VariantRecord], List[VariantRecord]]]: Dictionary
            mapping sample names to tuples of variant lists for each chromosome
            copy.
    """
    # initialize sample-variant map for both copies
    sample_variants = {s: ([], []) for s in samples}
    for variant in variants:
        assert len(variant.samples) == 1
        for chromcopy in [0, 1]:  # iterate over chromosome copies
            for s in variant.samples[0][chromcopy]:  # type: ignore
                # add variant to sample-variant map
                sample_variants[s][chromcopy].append(variant)
    # remove samples without variants
    sample_variants = {s: v for s, v in sample_variants.items() if v[0] or v[1]}
    return sample_variants


def compute_haplotypes_unphased(
    variants: List[VariantRecord], samples: List[str]
) -> Dict[str, List[VariantRecord]]:
    """
    Assigns unphased SNV variants to each sample without assuming diploid copies.
    Returns a dictionary mapping each sample to a list of variant records.

    Args:
        variants (List[VariantRecord]): List of variant records to assign.
        samples (List[str]): List of sample names.

    Returns:
        Dict[str, List[VariantRecord]]: Dictionary mapping sample names to lists
            of variant records.
    """
    variants = [v for v in variants if v.vtype[0] == VTYPES[0]]
    # initialize sample-variant map
    sample_variants = {s: [] for s in samples}  # do not assume diploid copies
    for variant in variants:
        assert len(variant.samples) == 1
        # always used first copy in unphased
        for s in variant.samples[0][0]:  # type: ignore
            # add variant to sample-variant map
            sample_variants[s].append(variant)
    # remove samples without variants
    sample_variants = {s: v for s, v in sample_variants.items() if v}
    return sample_variants


def compute_indel_haplotypes_unphased(
    variants: List[VariantRecord], samples: List[str]
) -> Dict[str, List[VariantRecord]]:
    """
    Assigns unphased indel variants to each sample without assuming diploid copies.
    Returns a dictionary mapping each sample to a list of indel variant records.

    Args:
        variants (List[VariantRecord]): List of indel variant records to assign.
        samples (List[str]): List of sample names.

    Returns:
        Dict[str, List[VariantRecord]]: Dictionary mapping sample names to lists
            of indel variant records.
    """
    # initialize sample-variant map
    sample_variants = {s: [] for s in samples}  # do not assume diploid copies
    for variant in variants:
        assert len(variant.samples) == 1
        # always used first copy in unphased
        for s in variant.samples[0][0]:  # type: ignore
            if s not in samples:
                continue
            # add variant to sample-variant map
            sample_variants[s].append(variant)
    # remove samples without variants
    sample_variants = {s: v for s, v in sample_variants.items() if v}
    return sample_variants


def ishomozygous(haplotypes: List[Haplotype]) -> bool:
    """
    Determines if all haplotypes in the list have identical sequences.
    Returns True if all haplotypes share the same sequence, otherwise False.

    Args:
        haplotypes (List[Haplotype]): List of haplotype objects to compare.

    Returns:
        bool: True if all haplotypes have the same sequence, False otherwise.
    """
    return len({h.sequence.sequence for h in haplotypes}) == 1


def _collapse_haplotypes(
    sequence: str, haplotypes: List[Haplotype], debug: bool
) -> Haplotype:
    """
    Collapses a list of haplotypes with the same sequence into a single haplotype
    object. Returns a new haplotype representing the merged information from all
    input haplotypes.

    Args:
        sequence (str): The nucleotide sequence for the collapsed haplotype.
        haplotypes (List[Haplotype]): List of haplotype objects to collapse.
        debug (bool): Flag to enable debug mode.

    Returns:
        Haplotype: A new haplotype object representing the collapsed haplotypes.
    """
    hap = Haplotype(
        Sequence(sequence, debug, allow_lower_case=True),
        haplotypes[0].coordinates,
        haplotypes[0].phased,
        0,
        debug,
    )
    samples = "REF" if sequence.isupper() else ",".join({h.samples for h in haplotypes})
    hap.samples = samples
    variants = (
        "NA"
        if sequence.isupper()
        else ",".join(sorted({h.variants for h in haplotypes}))
    )
    hap.variants = variants
    hap.set_afs(haplotypes[0].afs)
    hap.set_posmap(
        haplotypes[0].posmap, haplotypes[0].posmap_rev
    )  # same posmap for all collapsed haplotypes
    hap.set_variant_alleles(haplotypes[0].variant_alleles)  # same variant alleles
    return hap


def collapse_haplotypes(haplotypes: List[Haplotype], debug: bool) -> List[Haplotype]:
    """
    Collapses haplotypes with identical sequences into single representative
    haplotype objects. Returns a list of unique haplotypes, each representing all
    input haplotypes with the same sequence.

    Args:
        haplotypes (List[Haplotype]): List of haplotype objects to collapse.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Haplotype]: List of collapsed, unique haplotype objects.
    """
    haplotypes_dict = [(h.sequence.sequence, h) for h in haplotypes]
    haplotypes_collapsed = defaultdict(list)
    for seq, hap in haplotypes_dict:
        haplotypes_collapsed[seq].append(hap)
    return [
        _collapse_haplotypes(seq, haplist, debug)
        for seq, haplist in haplotypes_collapsed.items()
    ]


def _solve_haplotypes_phased(
    sequence: str,
    coordinates: Coordinate,
    phased: bool,
    variants: Tuple[List[VariantRecord], List[VariantRecord]],
    sample: str,
    debug: bool,
) -> List[Haplotype]:
    """
    Solves haplotypes for diploid samples using phased variant information.
    Returns a list containing one or two haplotype objects depending on zygosity.

    Args:
        sequence (str): The reference sequence for the haplotype.
        coordinates (Coordinate): Genomic coordinates for the haplotype.
        phased (bool): Indicates if the variants are phased.
        variants (Tuple[List[VariantRecord], List[VariantRecord]]): Tuple of
            variant lists for each chromosome copy.
        sample (str): Sample name.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Haplotype]: List containing one or two haplotype objects.
    """
    # solve haplotypes for diploid samples
    h0 = Haplotype(
        Sequence(sequence, debug), coordinates, phased, 0, debug
    )  # first copy
    h0.add_variants_phased(variants[0], sample)  # add variants to haplotype
    h1 = Haplotype(
        Sequence(sequence, debug), coordinates, phased, 1, debug
    )  # second copy
    h1.add_variants_phased(variants[1], sample)  # add variants to haplotype
    # check for homozygous haplotypes
    if ishomozygous([h0, h1]):
        h0.homozygous_samples()
        return [h0]  # return only one haplotype (homozygous sample)
    return [h0, h1]  # return both haplotypes


def solve_haplotypes_phased(
    sample_variants: Dict[str, Tuple[List[VariantRecord], List[VariantRecord]]],
    hapseqs: List[Haplotype],
    refseq: str,
    coordinates: Coordinate,
    phased: bool,
    debug: bool,
) -> List[Haplotype]:
    """
    Solves haplotypes for each sample using phased variant information in diploid
    samples. Returns a list of unique haplotype objects after collapsing those with
    identical sequences.

    Args:
        sample_variants (Dict[str, Tuple[List[VariantRecord], List[VariantRecord]]]):
            Dictionary mapping sample names to tuples of variant lists for each
            chromosome copy.
        hapseqs (List[Haplotype]): List of initial haplotype objects.
        refseq (str): Reference sequence for the haplotypes.
        coordinates (Coordinate): Genomic coordinates for the haplotypes.
        phased (bool): Indicates if the variants are phased.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Haplotype]: List of unique haplotype objects.
    """
    # solve haplotypes for each sample (assumes diploid samples)
    for sample, variants in sample_variants.items():
        hapseqs += _solve_haplotypes_phased(
            refseq, coordinates, phased, variants, sample, debug
        )
    return collapse_haplotypes(hapseqs, debug)  # collapse haplotypes with same sequence


def _split_id(vid: str) -> str:
    """
    Splits a variant ID string into its chromosome, position, and reference allele
    components. Returns a string combining the chromosome, position, and reference
    allele for use in identifying multiallelic sites.

    Args:
        vid (str): Variant ID string in the format 'chrom-pos-ref/alt'.

    Returns:
        str: String in the format 'chrom-pos-ref' representing the variant's origin.
    """
    chrom, pos, ref_alt = vid.split("-")  # retrieve variant ids fields
    ref, alt = ref_alt.split("/")
    return f"{chrom}-{pos}-{ref}"  # to recover origin of multiallelic sites


def generate_variants_combinations(
    variants: List[VariantRecord], sample: str
) -> List[List[VariantRecord]]:
    """
    Generates all possible combinations of variants for a given sample, accounting
    for multiallelic sites. Returns a list of variant combinations, where each
    combination represents a possible haplotype configuration.

    Args:
        variants (List[VariantRecord]): List of variant records to combine.
        sample (str): Sample name for which to generate variant combinations.

    Returns:
        List[List[VariantRecord]]: List of variant combinations, each as a list
            of variant records.
    """
    variant_groups = {variant.id[0]: [variant] for variant in variants}
    return [list(comb) for comb in product(*variant_groups.values())]


def _solve_haplotypes_unphased(
    sequence: str,
    coordinates: Coordinate,
    phased: bool,
    variants: List[VariantRecord],
    sample: str,
    debug: bool,
) -> List[Haplotype]:
    """
    Solves haplotypes for a sample using unphased variant information.
    Returns a list of haplotype objects representing all possible variant
    combinations for the sample.

    Args:
        sequence (str): The reference sequence for the haplotype.
        coordinates (Coordinate): Genomic coordinates for the haplotype.
        phased (bool): Indicates if the variants are phased.
        variants (List[VariantRecord]): List of variant records to assign.
        sample (str): Sample name.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Haplotype]: List of haplotype objects representing possible variant
            combinations.
    """
    variants_combinations = generate_variants_combinations(variants, sample)
    haps = []  # haplotype list
    for variant_combination in variants_combinations:
        variant_combination = [v for v in variant_combination if v is not None]
        h = Haplotype(Sequence(sequence, debug), coordinates, phased, 0, debug)
        # add variants to sample haplotypes
        h.add_variants_unphased(variant_combination, sample)  
        haps.append(h)  # insert haplotype to haplotypes list
    return [haps[0]] if ishomozygous(haps) else haps


def solve_haplotypes_unphased(
    sample_variants: Dict[str, List[VariantRecord]],
    hapseqs: List[Haplotype],
    refseq: str,
    coordinates: Coordinate,
    phased: bool,
    debug: bool,
) -> List[Haplotype]:
    """
    Solves haplotypes for each sample using unphased variant information.
    Returns a list of unique haplotype objects after collapsing those with
    identical sequences.

    Args:
        sample_variants (Dict[str, List[VariantRecord]]): Dictionary mapping
            sample names to lists of variant records.
        hapseqs (List[Haplotype]): List of initial haplotype objects.
        refseq (str): Reference sequence for the haplotypes.
        coordinates (Coordinate): Genomic coordinates for the haplotypes.
        phased (bool): Indicates if the variants are phased.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Haplotype]: List of unique haplotype objects.
    """
    # solve haplotypes for each sample
    for sample, variants in sample_variants.items():
        hapseqs += _solve_haplotypes_unphased(
            refseq, coordinates, phased, variants, sample, debug
        )
    return collapse_haplotypes(hapseqs, debug)


def classify_variants(
    variants: List[VariantRecord],
) -> Tuple[List[VariantRecord], List[VariantRecord]]:
    """
    Splits a list of variant records into SNVs and indels based on their type.
    Returns a tuple containing two lists: one for SNVs and one for indels.

    Args:
        variants (List[VariantRecord]): List of variant records to classify.

    Returns:
        Tuple[List[VariantRecord], List[VariantRecord]]: Tuple containing a list
            of SNVs and a list of indels.
    """
    # split variants according to their type (required for unphased variants processing)
    snvs, indels = [], []
    for v in variants:
        (snvs if v.vtype[0] == VTYPES[0] else indels).append(v)
    assert len(snvs) + len(indels) == len(variants)
    return snvs, indels


def compute_snvs_haplotype_unphased(
    snvs: List[VariantRecord],
    samples: List[str],
    refseq: str,
    coordinates: Coordinate,
    phased: bool,
    debug: bool,
) -> List[Haplotype]:
    """
    Computes and solves haplotypes using only unphased SNV variants for the given
    samples. Returns a list of unique haplotype objects representing all possible
    SNV configurations.

    Args:
        snvs (List[VariantRecord]): List of SNV variant records to assign.
        samples (List[str]): List of sample names.
        refseq (str): Reference sequence for the haplotypes.
        coordinates (Coordinate): Genomic coordinates for the haplotypes.
        phased (bool): Indicates if the variants are phased.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Haplotype]: List of unique haplotype objects.
    """
    # compute and solve snvs-only haplotype (only unphased)
    samples_variants_snvs = compute_haplotypes_unphased(snvs, samples)
    return solve_haplotypes_unphased(
        samples_variants_snvs,
        [],  # start with empty list instead of reference
        refseq,
        coordinates,
        phased,
        debug,
    )


def create_indel_window(
    indel: VariantRecord, region: Region
) -> Tuple[str, Coordinate, int, int, int]:
    """
    Creates a window around an indel variant within a genomic region, including
    additional bases upstream and downstream. Returns the sequence of the window,
    its coordinates, window start and stop positions, and the indel length.

    Args:
        indel (VariantRecord): The indel variant for which to create the window.
        region (Region): The genomic region containing the indel.

    Returns:
        Tuple[str, Coordinate, int, int, int]: A tuple containing the window sequence,
            its coordinates, window start, window stop, and indel length.
    """
    # compute indel window; additional 50 bp upstream and downstream
    indel_length = len(indel.alt[0]) - len(indel.ref) if indel.alt else 0
    # compute window boundaries
    window_start = max(region.start, indel.position - 100)
    window_stop = min(region.stop, indel.position - (indel_length) + 100)
    indel_coords = Coordinate(region.contig, window_start, window_stop, 0)
    # extract sequence for the indel window
    startrel = window_start - region.start  # relative window start position
    stoprel = window_stop - region.start + 1  # relative window stop position
    indel_sequence = region.sequence.sequence[startrel:stoprel]
    return indel_sequence, indel_coords, window_start, window_stop, indel_length


def find_overlapping_snvs(
    indel_start: int, indel_stop: int, snvs: List[VariantRecord]
) -> List[VariantRecord]:
    """
    Identifies SNV variants that overlap with a specified indel window.
    Returns a list of SNV variant records whose positions fall within the indel
    window.

    Args:
        indel_start (int): Start position of the indel window.
        indel_stop (int): Stop position of the indel window.
        snvs (List[VariantRecord]): List of SNV variant records to check.

    Returns:
        List[VariantRecord]: List of SNV variant records overlapping the indel window.
    """
    # find snvs overlapping the current indel window
    return [snv for snv in snvs if indel_start <= snv.position <= indel_stop]


def retrieve_indel_samples(indel: VariantRecord) -> Set[str]:
    """
    Retrieves the set of sample names that carry the given indel variant.
    Returns a set of sample names if present, otherwise returns an empty set.

    Args:
        indel (VariantRecord): The indel variant from which to extract sample names.

    Returns:
        Set[str]: Set of sample names carrying the indel variant.
    """
    return set(indel.samples[0][0]) if indel and len(indel.samples) > 0 else set()


def set_haplotypes_samples(
    indel_haplotypes: List[Haplotype], indel_samples: Set[str]
) -> List[Haplotype]:
    """
    Updates the samples attribute of each haplotype to reflect only those samples
    that carry the indel. Returns the list of haplotypes with updated sample
    information for indel carriers.

    Args:
        indel_haplotypes (List[Haplotype]): List of haplotype objects to update.
        indel_samples (Set[str]): Set of sample names carrying the indel.

    Returns:
        List[Haplotype]: List of haplotype objects with updated sample information.
    """
    # set samples for each haplotype to reflect indel carriers
    for hap in indel_haplotypes:
        if hap.samples != "REF":  # alternative
            # ensure haplotype samples are subset of indel carriers
            hap_samples = set(hap.samples.split(",")) if hap.samples else set()
            final_samples = hap_samples.intersection(indel_samples)
            if final_samples:
                hap.samples = ",".join(sorted(final_samples))
    return indel_haplotypes


def create_indels_haplotype_unphased(
    indel: VariantRecord,
    snvs: List[VariantRecord],
    region: Region,
    phased: bool,
    debug: bool,
) -> List[Haplotype]:
    """
    Computes and solves haplotypes for samples carrying a given indel using unphased
    variant information. Returns a list of haplotype objects representing all possible
    indel and overlapping SNV configurations for indel carriers.

    Args:
        indel (VariantRecord): The indel variant for which to create haplotypes.
        snvs (List[VariantRecord]): List of SNV variant records in the region.
        region (Region): The genomic region containing the indel.
        phased (bool): Indicates if the variants are phased.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Haplotype]: List of haplotype objects for samples carrying the indel.
    """
    # compute indel window and find overlapping snvs
    indel_sequence, indel_coords, indel_start, indel_stop, indel_length = (
        create_indel_window(indel, region)
    )
    snvs_overlapping = find_overlapping_snvs(indel_start, indel_stop, snvs)
    indel_samples = retrieve_indel_samples(indel)  # recover indel's samples
    indel_haplotypes = []  # initialize variable
    # create sample-variant mapping only for samples carrying the indel
    if indel_samples:  # compute and solve indel haplotypes
        indel_variants = [indel] + snvs_overlapping
        samples_variants_indel = compute_indel_haplotypes_unphased(
            indel_variants, list(indel_samples)
        )
        indel_haplotypes = solve_haplotypes_unphased(
            samples_variants_indel,
            [],  # start with empty list
            indel_sequence,
            indel_coords,
            phased,
            debug,
        )
        # set samples for each indel haplotype to reflect indel carriers
        indel_haplotypes = set_haplotypes_samples(indel_haplotypes, indel_samples)
    return indel_haplotypes


def add_variants_unphased(
    haplotypes: List[Haplotype],
    region: Region,
    vcfs: Dict[str, VCF],
    variants: List[VariantRecord],
    phased: bool,
    debug: bool,
) -> List[Haplotype]:
    """
    Adds unphased SNV and indel variants to the list of haplotypes for a given region.
    Returns an updated list of haplotype objects representing all possible variant
    configurations for the region.

    Args:
        haplotypes (List[Haplotype]): List of initial haplotype objects.
        region (Region): The genomic region for which to add variants.
        vcfs (Dict[str, VCF]): Dictionary mapping contig names to VCF objects.
        variants (List[VariantRecord]): List of variant records to add.
        phased (bool): Indicates if the variants are phased.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Haplotype]: Updated list of haplotype objects with added variants.
    """
    snvs, indels = classify_variants(variants)  # split variants according to their type
    if snvs:  # compute snvs-only haplotypes
        haplotypes.extend(
            compute_snvs_haplotype_unphased(
                snvs,
                vcfs[region.contig].samples,
                region.sequence.sequence,
                region.coordinates,
                phased,
                debug,
            )
        )
    for indel in indels:  # create haplotype for each individual indel
        if region.coordinates.startp <= indel.position < region.coordinates.stopp:
            haplotypes.extend(
                create_indels_haplotype_unphased(indel, snvs, region, phased, debug)
            )
    return haplotypes


def add_variants_phased(
    haplotypes: List[Haplotype],
    region: Region,
    vcfs: Dict[str, VCF],
    variants: List[VariantRecord],
    phased: bool,
    debug: bool,
) -> List[Haplotype]:
    """
    Adds phased SNV and indel variants to the list of haplotypes for a given region.
    Returns an updated list of haplotype objects representing all possible phased
    variant configurations for the region.

    Args:
        haplotypes (List[Haplotype]): List of initial haplotype objects.
        region (Region): The genomic region for which to add variants.
        vcfs (Dict[str, VCF]): Dictionary mapping contig names to VCF objects.
        variants (List[VariantRecord]): List of variant records to add.
        phased (bool): Indicates if the variants are phased.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Haplotype]: Updated list of haplotype objects with added phased variants.
    """
    # recover variants combinations on each chromosome copy
    samples_variants = compute_haplotypes_phased(variants, vcfs[region.contig].samples)
    # solve haplotypes for each sample
    return solve_haplotypes_phased(
        samples_variants,
        haplotypes,
        region.sequence.sequence,
        region.coordinates,
        phased,
        debug,
    )


def add_variants(
    vcflist: List[str],
    regions: RegionList,
    haplotypes: Dict[Region, List[Haplotype]],
    verbosity: int,
    debug: bool,
) -> Tuple[Dict[Region, List[Haplotype]], bool]:
    """
    Adds SNV and indel variants from VCF files to haplotypes for each region.
    Returns an updated dictionary of haplotypes for each region and a flag indicating
    whether the VCFs are phased.

    Args:
        vcflist (List[str]): List of VCF file paths.
        regions (RegionList): List of regions to process.
        haplotypes (Dict[Region, List[Haplotype]]): Dictionary mapping regions to
            lists of haplotype objects.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode.

    Returns:
        Tuple[Dict[Region, List[Haplotype]], bool]: Updated haplotypes dictionary
            and a boolean indicating if the VCFs are phased.
    """
    # read VCF files and extract variants located within the region
    vcfs = read_vcf(vcflist, verbosity, debug)
    variants = fetch_variants(vcfs, regions, verbosity, debug)
    phased = vcfs[regions[0].contig].phased  # assess VCF phasing
    for region in regions:  # reconstruct haplotypes for each region
        if phased:  # phased VCFs
            haplotypes[region] = add_variants_phased(
                haplotypes[region], region, vcfs, variants[region], phased, debug
            )
        else:  # unphased VCFs
            haplotypes[region] = add_variants_unphased(
                haplotypes[region], region, vcfs, variants[region], phased, debug
            )
    return haplotypes, phased


def generate_haplotype_ids(
    haplotypes: Dict[Region, List[Haplotype]],
) -> Dict[Region, List[Haplotype]]:
    """
    Generates and assigns unique random IDs to each haplotype in all regions.
    Returns the updated haplotypes dictionary with assigned IDs.

    Args:
        haplotypes (Dict[Region, List[Haplotype]]): Dictionary mapping regions
            to lists of haplotype objects.

    Returns:
        Dict[Region, List[Haplotype]]: Updated dictionary with haplotype IDs assigned.
    """
    chars = string.ascii_letters + string.digits  # generate random characters
    for region, haps in haplotypes.items():
        ids = set()
        while len(ids) < len(haps):
            ids.add("hap_" + "".join(random.choices(chars, k=8)))  # generate random ID
        ids = list(ids)
        for i, hap in enumerate(haps):
            hap.id = ids[i]  # set haplotype ID
    return haplotypes


def haplotypes_table(
    haplotypes: Dict[Region, List[Haplotype]], outdir: str, verbosity: int, debug: bool
) -> None:
    """
    Writes a table of haplotypes for each region to a TSV file in the specified
    output directory. Each table includes haplotype ID, sequence, variants, and
    samples for all haplotypes in the region.

    Args:
        haplotypes (Dict[Region, List[Haplotype]]): Dictionary mapping regions to lists of haplotype objects.
        outdir (str): Output directory where the haplotype tables will be written.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode.

    Returns:
        None
    """
    # write haplotypes table to file
    print_verbosity("Writing haplotypes table", verbosity, VERBOSITYLVL[1])
    start = time()  # track haplotypes table writing time
    for region, haps in haplotypes.items():
        haptable_fname = os.path.join(
            outdir, f"haplotypes_table_{region.contig}_{region.start}_{region.stop}.tsv"
        )
        try:
            with open(haptable_fname, "w") as outfile:
                outfile.write("\t".join(HAPTABCNAMES) + "\n")
                for hap in haps:
                    outfile.write(
                        f"{hap.id}\t{hap.sequence}\t{hap.variants}\t{hap.samples}\n"
                    )
        except OSError as e:
            exception_handler(
                CrisprHawkHaplotypeError,
                f"Failed writing haplotype table for region {region}",
                os.EX_IOERR,
                debug,
                e,
            )
    print_verbosity(
        f"Haplotypes table written in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )


def reconstruct_haplotypes(
    regions: RegionList, args: CrisprHawkSearchInputArgs
) -> Tuple[Dict[Region, List[Haplotype]], bool, bool]:
    """Reconstructs haplotypes for each genomic region using VCF data and input
    arguments.

    This function initializes reference haplotypes, adds variants from VCF files
    if provided, generates unique haplotype IDs, and optionally writes haplotype
    tables to disk.

    Args:
        regions: List of Region objects for which to reconstruct haplotypes.
        args: CrisprHawkSearchInputArgs object containing haplotype reconstruction
            parameters.

    Returns:
        Tuple containing:
            - Dictionary mapping Region objects to lists of Haplotype objects.
            - Boolean indicating if variants are present.
            - Boolean indicating if the VCFs are phased.
    """
    # read input vcf files and fetch variants in each region
    print_verbosity("Reconstructing haplotypes", args.verbosity, VERBOSITYLVL[1])
    start = time()  # track haplotypes reconstruction time
    # initialize haplotypes list with reference sequence haplotype
    haplotypes = initialize_haplotypes(regions, args.debug)
    phased, variants_present = False, False  # default values
    if args.vcfs:  # add variants to regions and solve haplotypes
        variants_present = True  # variants added
        haplotypes, phased = add_variants(
            args.vcfs, regions, haplotypes, args.verbosity, args.debug
        )
    # generate random haplotype IDs
    haplotypes = generate_haplotype_ids(haplotypes)
    print_verbosity(
        f"Haplotypes reconstructed in {time() - start:.2f}s",
        args.verbosity,
        VERBOSITYLVL[2],
    )
    if args.haplotype_table:  # write haplotypes table
        haplotypes_table(haplotypes, args.outdir, args.verbosity, args.debug)
    return haplotypes, variants_present, phased
