""" """

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


def compute_indel_haplotypes_unphased(
    variants: List[VariantRecord], samples: List[str]
) -> Dict[str, List[VariantRecord]]:
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
    hap.set_afs(haplotypes[0].afs)
    hap.set_posmap(
        haplotypes[0].posmap, haplotypes[0].posmap_rev
    )  # same posmap for all collapsed haplotypes
    hap.set_variant_alleles(haplotypes[0].variant_alleles)  # same variant alleles
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
        is_snv = True  # variant.vtype[0] == VTYPES[0]
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
        h.add_variants_unphased(
            variant_combination, sample
        )  # add variants to sample haplotypes
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


def classify_variants(
    variants: List[VariantRecord],
) -> Tuple[List[VariantRecord], List[VariantRecord]]:
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
    # compute and solve snvs-only haplotype (only unphased)
    samples_variants_snvs = compute_haplotypes_unphased(snvs, samples)
    haplotypes_snvs = solve_haplotypes_unphased(
        samples_variants_snvs,
        [],  # start with empty list instead of reference
        refseq,
        coordinates,
        phased,
        debug,
    )
    return haplotypes_snvs


def create_indel_window(
    indel: VariantRecord, region: Region
) -> Tuple[str, Coordinate, int, int, int]:
    # compute indel window; additional 50 bp upstream and downstream
    indel_length = len(indel.alt[0]) - len(indel.ref) if indel.alt else 0
    # compute window boundaries
    window_start = max(region.start, indel.position - 50)
    window_stop = min(region.stop, indel.position - (indel_length) + 50)
    indel_coords = Coordinate(region.contig, window_start, window_stop, 0)
    # extract sequence for the indel window
    startrel = window_start - region.start  # relative window start position
    stoprel = window_stop - region.start + 1  # relative window stop position
    indel_sequence = region.sequence.sequence[startrel:stoprel]
    return indel_sequence, indel_coords, window_start, window_stop, indel_length


def find_overlapping_snvs(
    indel_start: int, indel_stop: int, snvs: List[VariantRecord]
) -> List[VariantRecord]:
    # find snvs overlapping the current indel window
    return [snv for snv in snvs if indel_start <= snv.position <= indel_stop]


def retrieve_indel_samples(indel: VariantRecord) -> Set[str]:
    return set(indel.samples[0][0]) if indel and len(indel.samples) > 0 else set()


def set_haplotypes_samples(
    indel_haplotypes: List[Haplotype], indel_samples: Set[str]
) -> List[Haplotype]:
    # set samples for each haplotype to reflect indel carriers
    for hap in indel_haplotypes:
        if hap.samples != "REF":  # alternative
            # ensure haplotype samples are subset of indel carriers
            hap_samples = set(hap.samples.split(",")) if hap.samples else set()
            final_samples = hap_samples.intersection(indel_samples)
            if final_samples:
                hap.set_samples(",".join(sorted(final_samples)))
    return indel_haplotypes


def create_indels_haplotype_unphased(
    indel: VariantRecord,
    snvs: List[VariantRecord],
    region: Region,
    phased: bool,
    debug: bool,
) -> List[Haplotype]:
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
    # exit()
    return haplotypes, variants_present, phased
