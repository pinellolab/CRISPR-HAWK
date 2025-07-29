""" """

from .crisprhawk_error import CrisprHawkHaplotypeError
from .exception_handlers import exception_handler
from .utils import print_verbosity, flatten_list, VERBOSITYLVL
from .region import Region, RegionList
from .variant import VCF, VariantRecord, VTYPES
from .coordinate import Coordinate
from .sequence import Sequence
from .haplotype import Haplotype

from itertools import product
from typing import List, Dict, Tuple
from collections import defaultdict
from time import time

import random
import string
import os

HAPTABCNAMES = ["id", "haplotype", "variants", "samples"]


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
            Exception, "Failed parsing VCF files", os.EX_DATAERR, debug, e  # type: ignore
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


def ishomozygous(haplotypes: List[Haplotype]) -> bool:
    return len({h.sequence.sequence for h in haplotypes}) == 1


def _collapse_haplotypes(
    sequence: str, haplotypes: List[Haplotype], debug: bool
) -> Haplotype:
    hap = Haplotype(
        Sequence(sequence, debug, allow_lower_case=True),
        haplotypes[0].coordinates,
        haplotypes[0].phased,
        0,
        debug,
    )
    samples = "REF" if sequence.isupper() else ",".join({h.samples for h in haplotypes})
    hap.set_samples(samples)
    variants = (
        "NA"
        if sequence.isupper()
        else ",".join(sorted({h.variants for h in haplotypes}))
    )
    hap.set_variants(variants)
    hap.set_posmap(haplotypes[0].posmap)  # same posmap for all collapsed haplotypes
    return hap


def collapse_haplotypes(haplotypes: List[Haplotype], debug: bool) -> List[Haplotype]:
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
    # solve haplotypes for each sample (assumes diploid samples)
    for sample, variants in sample_variants.items():
        hapseqs += _solve_haplotypes_phased(
            refseq, coordinates, phased, variants, sample, debug
        )
    return collapse_haplotypes(hapseqs, debug)  # collapse haplotypes with same sequence


def _split_id(vid: str) -> str:
    chrom, pos, ref_alt = vid.split("-")  # retrieve variant ids fields
    ref, alt = ref_alt.split("/")
    return f"{chrom}-{pos}-{ref}"  # to recover origin of multiallelic sites


def generate_variants_combinations(
    variants: List[VariantRecord], sample: str
) -> List[List[VariantRecord]]:
    variant_groups = {}  # groups of alternative variants
    for variant in variants:
        # retrieve variant id to match multiallelic sites
        vid = _split_id(variant.id[0])
        is_snv = variant.vtype[0] == VTYPES[0]
        is_homozygous = all(sample in call for call in variant.samples[0])
        if is_snv or is_homozygous:  # only one option available
            variant_groups[vid] = [variant]
        elif vid in variant_groups:
            # remove reference allele for alternatives in multiallelic site
            variant_groups[vid] = [v for v in variant_groups[vid] if v is not None]
            variant_groups[vid].append(variant)  # add alt to multiallelic
        else:
            variant_groups[vid] = [variant, None]  # add reference (None)
    return [list(comb) for comb in product(*variant_groups.values())]


def _solve_haplotypes_unphased(
    sequence: str,
    coordinates: Coordinate,
    phased: bool,
    variants: List[VariantRecord],
    sample: str,
    debug: bool,
) -> List[Haplotype]:
    variants_combinations = generate_variants_combinations(variants, sample)
    haps = []  # haplotype list
    for variant_combination in variants_combinations:
        variant_combination = [v for v in variant_combination if v is not None]
        h = Haplotype(Sequence(sequence, debug), coordinates, phased, 0, debug)
        h.add_variants_unphased(variant_combination, sample)  # add variants to sample haplotypes
        haps.append(h)  # insert haplotype to haplotypes list
    return [haps[0]] if ishomozygous(haps) else haps


def solve_haplotypes_unphased(
    sample_variants: Dict[str, List[VariantRecord]],
    hapseqs: List[Haplotype],
    refseq: str,
    coordinates: Coordinate,
    phased: bool,
    debug: bool,
):
    # solve haplotypes for each sample
    for sample, variants in sample_variants.items():
        hapseqs += _solve_haplotypes_unphased(
            refseq, coordinates, phased, variants, sample, debug
        )
    return collapse_haplotypes(hapseqs, debug)


def add_variants(
    vcflist: List[str],
    regions: RegionList,
    haplotypes: Dict[Region, List[Haplotype]],
    verbosity: int,
    debug: bool,
) -> Tuple[Dict[Region, List[Haplotype]], bool]:
    # read VCF files and extract variants located within the region
    vcfs = read_vcf(vcflist, verbosity, debug)
    variants = fetch_variants(vcfs, regions, verbosity, debug)
    
    # for r, vars in variants.items():
    #     print(r)
    #     for v in vars:
    #         print(v, v.samples, v.vtype)
    # print()

    phased = vcfs[regions[0].contig].phased  # assess VCF phasing
    for region in regions:  # reconstruct haplotypes for each region
        if phased:  # phased VCFs
            # recover variants combinations on each chromosome copy
            samples_variants = compute_haplotypes_phased(
                variants[region], vcfs[regions[0].contig].samples
            )
            # solve haplotypes for each sample
            haplotypes[region] = solve_haplotypes_phased(
                samples_variants,
                haplotypes[region],
                region.sequence.sequence,
                region.coordinates,
                phased,
                debug,
            )
        else:  # unphased VCFs
            # Split variants into SNVs and indels
            snvs, indels = [], [] 
            for v in variants[region]:
                if v.vtype[0] == VTYPES[0]:  # assuming VTYPES[0] is "snp"
                    snvs.append(v)
                else:
                    indels.append(v)
            
            # Create haplotype with all SNVs
            if snvs:
                samples_variants_snvs = compute_haplotypes_unphased(
                    snvs, vcfs[regions[0].contig].samples
                )
                snv_haplotypes = solve_haplotypes_unphased(
                    samples_variants_snvs,
                    [],  # start with empty list instead of reference
                    region.sequence.sequence,
                    region.coordinates,
                    phased,
                    debug,
                )
                haplotypes[region].extend(snv_haplotypes)
            
            # Create separate haplotype for each indel
            for indel in indels:
                # Calculate indel window: 23bp upstream and downstream
                indel_pos = indel.position  # assuming this gives the position
                indel_length = len(indel.alt[0]) - len(indel.ref) if indel.alt else 0
                
                # Calculate window boundaries
                window_start = max(region.coordinates.start, indel_pos - 50)
                window_end = min(region.coordinates.stop, indel_pos - (indel_length) + 50)
                
                # Create coordinate object for the indel window
                indel_window_coords = Coordinate(
                    region.contig, 
                    window_start, 
                    window_end,
                    debug
                )

                print(indel)
                print(indel_window_coords)

                # Extract sequence for the indel window
                # Calculate relative positions within the region sequence
                rel_start = window_start - region.coordinates.start - 1
                rel_end = window_end - region.coordinates.start
                indel_window_sequence = region.sequence.sequence[rel_start:rel_end]

                print(region.coordinates, len(region))
                print(rel_start, rel_end)
                print(indel_window_sequence, len(indel_window_sequence))

                # Apply the indel to the sequence before creating haplotype
                indel_rel_pos = indel_pos - window_start  # position relative to window start
                ref_seq = indel.ref[0] if indel.ref else ""
                alt_seq = indel.alt[0] if indel.alt else ""
                
                # Apply indel: replace ref with alt at the relative position
                offset = 0 if indel_length > 0 else abs(indel_length)
                indel_applied_sequence = (
                    indel_window_sequence[:51] + 
                    alt_seq + 
                    indel_window_sequence[51 + 1 + offset:]
                )
                
                # Find SNVs that overlap with this indel window
                overlapping_snvs = [
                    snv for snv in snvs 
                    if window_start <= snv.position <= window_end
                ]
                
                # Combine indel with overlapping SNVs
                indel_variants = [indel] + overlapping_snvs
                
                # Create haplotype for this indel + overlapping SNVs
                samples_variants_indel = compute_haplotypes_unphased(
                    indel_variants, vcfs[regions[0].contig].samples
                )
                
                indel_haplotypes = solve_haplotypes_unphased(
                    samples_variants_indel,
                    [],  # start with empty list
                    indel_applied_sequence,
                    indel_window_coords,
                    phased,
                    debug,
                )
                
                haplotypes[region].extend(indel_haplotypes)
    

            # samples_variants = compute_haplotypes_unphased(
            #     variants[region], vcfs[regions[0].contig].samples
            # )
            # haplotypes[region] = solve_haplotypes_unphased(
            #     samples_variants,
            #     haplotypes[region],
            #     region.sequence.sequence,
            #     region.coordinates,
            #     phased,
            #     debug,
            # )
            # pass
    return haplotypes, phased


def generate_haplotype_ids(
    haplotypes: Dict[Region, List[Haplotype]],
) -> Dict[Region, List[Haplotype]]:
    chars = string.ascii_letters + string.digits  # generate random characters
    for region, haps in haplotypes.items():
        ids = set()
        while len(ids) < len(haps):
            ids.add("hap_" + "".join(random.choices(chars, k=8)))  # generate random ID
        ids = list(ids)
        for i, hap in enumerate(haps):
            hap.set_id(ids[i])  # set haplotype ID
    return haplotypes


def haplotypes_table(
    haplotypes: Dict[Region, List[Haplotype]], outdir: str, verbosity: int, debug: bool
) -> None:
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
            exception_handler(CrisprHawkHaplotypeError, f"Failed writing haplotype table for region {region}", os.EX_IOERR, debug, e)  # type: ignore
    print_verbosity(
        f"Haplotypes table written in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )


def reconstruct_haplotypes(
    vcflist: List[str],
    regions: RegionList,
    store_table: bool,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> Tuple[Dict[Region, List[Haplotype]], bool, bool]:
    # read input vcf files and fetch variants in each region
    print_verbosity("Reconstructing haplotypes", verbosity, VERBOSITYLVL[1])
    start = time()  # track haplotypes reconstruction time
    # initialize haplotypes list with reference sequence haplotype
    haplotypes = initialize_haplotypes(regions, debug)
    phased, variants_present = False, False  # default values
    if vcflist:  # add variants to regions and solve haplotypes
        variants_present = True  # variants added
        haplotypes, phased = add_variants(
            vcflist, regions, haplotypes, verbosity, debug
        )
    # generate random haplotype IDs
    haplotypes = generate_haplotype_ids(haplotypes)
    print_verbosity(
        f"Haplotypes reconstructed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )
    if store_table:  # write haplotypes table
        haplotypes_table(haplotypes, outdir, verbosity, debug)
    exit()
    return haplotypes, variants_present, phased
