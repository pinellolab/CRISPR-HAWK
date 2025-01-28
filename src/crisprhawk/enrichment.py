"""
"""

from variants import VCF, VariantRecord, VTYPES
from bedfile import RegionList, Region, IndelRegion, Coordinate
from utils import print_verbosity, warning, VERBOSITYLVL, IUPAC_ENCODER
from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkEnrichmentError
from variant_map import VariantMap

from typing import List, Tuple, Union, Dict
from time import time

import os

INDELTYPES = [0, 1]  # indel types -> 0 for insertion, 1 for deletion


def guess_indel_type(ref: str, alt: str, position: int, debug) -> int:
    # guess from ref and alt alleles if the input indel is an insertion
    # (|ref| < |alt|) deletion (|ref| > |alt|)
    if len(ref) == len(alt):  # edge case that should never happen
        exception_handler(
            CrisprHawkEnrichmentError,
            f"Reference and alternative indel alleles at position {position} have the same length ({len(ref)})",
            os.EX_DATAERR,
            debug,
        )
    return INDELTYPES[0] if len(ref) < len(alt) else INDELTYPES[1]


def create_indel_region(
    sequence: str,
    startrel: int,
    stoprel: int,
    region_start: int,
    chrom: str,
    ref: str,
    alt: str,
    guidelen: int,
    pamlen: int,
    indel_type: str,
    debug: bool,
) -> IndelRegion:
    # recover original start and stop coordinates for indel subregion and create
    # coordinate object
    coord = Coordinate(chrom, startrel + region_start, stoprel + region_start)
    # create indel region object for each indel
    return IndelRegion(
        guidelen + pamlen - 1, sequence, coord, ref, alt, indel_type, debug
    )


def adjust_indel_position(
    indel_position: int,
    region_start: int,
    guidelen: int,
    pamlen: int,
    ref: str,
    alt: str,
    debug: bool,
) -> Tuple[int, int, int]:
    # compute indel location in the current region
    posrel = indel_position - 1 - region_start
    indel_length = abs(len(ref) - len(alt))  # compute indel length
    # compute start and stop coordinate for indel subregion
    start, stop = (
        posrel - guidelen - pamlen + 1,
        posrel + guidelen + pamlen - 1,
    )
    # adjust indel stop position according to the indel type
    if guess_indel_type(ref, alt, indel_position, debug) == INDELTYPES[0]:  # insertion
        stop += 1
    else:  # deletion
        stop += indel_length + 1
    return start, stop, posrel


def compute_indel_region(
    region: Region,
    indel: VariantRecord,
    alt: str,
    guidelen: int,
    pamlen: int,
    debug: bool,
) -> IndelRegion:
    # compute start, stop and indel position within the query region
    start, stop, _ = adjust_indel_position(
        indel.position, region.start, guidelen, pamlen, indel.ref, alt, debug
    )
    indelregion = create_indel_region(
        "".join(region[start:stop]),
        start,
        stop,
        region.start,
        indel.contig,
        indel.ref,
        alt,
        guidelen,
        pamlen,
        guess_indel_type(indel.ref, alt, indel.position, debug),
        debug,
    )  # retrieve str from list
    indelregion.enrich()  # insert indel in region
    return indelregion


def annotate_indel_variants(
    variant_map_snps: Dict[int, VariantRecord],
    indel: VariantRecord,
    region: Region,
    alt: str,
    guidelen: int,
    pamlen: int,
    indel_type: int,
    debug: bool,
) -> Dict[int, VariantRecord]:
    # dictionary to map variants to their relative position within the sequence
    variant_map = {}
    # compute start, stop and indel position within the query region
    start, stop, posrel = adjust_indel_position(
        indel.position, region.start, guidelen, pamlen, indel.ref, alt, debug
    )
    # compute region
    indel_length = abs(len(indel.ref) - len(alt))  # compute indel length
    indel_start, indel_stop = posrel, posrel  # compute indel range offset
    if guess_indel_type(indel.ref, alt, indel.position, debug) == INDELTYPES[1]:
        indel_stop += indel_length  # deletion
    for pos, snp in variant_map_snps.items():  # recover snp data within indel region
        if (start <= pos <= stop) and (pos < indel_start or pos > indel_stop):
            vpos = pos - start  # always if vpos < indel_start
            if pos >= indel_start:  # adjust position of variants occurring after indel
                vpos = (
                    vpos + indel_length
                    if indel_type == INDELTYPES[0]
                    else vpos - indel_length
                )
            variant_map[vpos] = snp  # snp flanks indel
    for i in range(len(alt)):
        variant_map[(indel_start - start) + i] = indel  # insert indel data
    return variant_map


def insert_indels(
    indels_map: Dict[Region, List[VariantRecord]],
    variant_maps: Dict[Region, Dict[int, VariantRecord]],
    guidelen: int,
    pamlen: int,
    debug: bool,
) -> Tuple[RegionList, Dict[Region, Dict[int, VariantRecord]]]:
    indelregions = RegionList([], debug)  # list of indel regions
    for r, indels in indels_map.items():
        for indel in indels:
            for altallele in indel.get_altalleles(VTYPES[1]):
                # create a region for each alternative allele in indel
                indelregion = compute_indel_region(
                    r, indel, altallele, guidelen, pamlen, debug
                )
                indelregions.append(indelregion)
                variant_maps[indelregion] = annotate_indel_variants(
                    variant_maps[r],
                    indel,
                    r,
                    altallele,
                    guidelen,
                    pamlen,
                    indelregion.indel_type,
                    debug,
                )
    return indelregions, variant_maps


# def insert_snps(
#     region: Region,
#     snps: List[VariantRecord],
#     no_filter: bool,
#     verbosity: int,
#     debug: bool,
# ) -> Dict[int, VariantRecord]:
#     # dictionary to map variants to their relative position within the sequence
#     variant_map = {}
#     # iterate over variants falling in the input region and add variants
#     for snp in snps:
#         # by default, crispr-hawk discards variants that are not flagged as PASS
#         # on filter, however the user may want to consider tem as well
#         if not no_filter and snp.filter != "PASS":
#             warning(f"Skipping {snp.format()}", verbosity)
#             continue
#         # compute relative variant position (1-based)
#         posrel = snp.position - 1 - region.start
#         # compute iupac char, representing the snp
#         iupac_nt = encode_snp_iupac(region[posrel], snp, debug)
#         if iupac_nt is not None:  # if none, indels -> skip
#             region.enrich(posrel, iupac_nt)  # assign iupac char at position
#             variant_map[posrel] = snp  # map variant to its relative position
#     return variant_map


# def encode_snp_iupac(
#     refnt: str, variant: VariantRecord, debug: bool
# ) -> Union[str, None]:
#     altalleles = "".join(variant.get_altalleles(VTYPES[0]))  # retrieve snps
#     if not altalleles:  # variants are all indels
#         return None
#     if refnt != variant.ref:  # ref alleles must match between vcf and fasta
#         exception_handler(
#             CrisprHawkEnrichmentError,
#             f"Reference allele mismatch between FASTA and VCF file ({refnt} - "
#             f"{variant.ref}) at position {variant.position}",
#             os.EX_DATAERR,
#             debug,
#         )
#     # encode ref and alt as iupac characters
#     return IUPAC_ENCODER["".join([variant.ref, altalleles])]

def load_vcfs(vcflist: List[str], verbosity: int, debug: bool) -> Dict[str, VCF]:
    # load vcf files and map each vcf to its contig (assume on vcf per contig)
    print_verbosity("Loading VCF files", verbosity, VERBOSITYLVL[2])
    start = time()  # track vcf parsing time
    vcfs = {vcf.contig: vcf for vcf in [VCF(f, verbosity, debug) for f in vcflist]}
    print_verbosity(f"Loaded {len(vcfs)} VCFs in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3])
    return vcfs

def adjust_region_coords(region: Region, guidelen: int) -> Tuple[int, int]:
    # remove guide length padding from region's start and stop position to
    # extract variants in that range
    return region.start + guidelen, region.stop - guidelen


def fetch_variants(vcf: VCF, region: Region, guidelen: int, verbosity: int) -> List[VariantRecord]:  
    # recover variants mapped on the query region
    # each region is padded by |guide| nts to avoid missing guides mapped on
    # region's border
    print_verbosity(f"Fetching variants in {region.display()}", verbosity, VERBOSITYLVL[3])
    start = time()  # track variants fecthing time
    variants = vcf.fetch(*adjust_region_coords(region, guidelen))
    # split variants between snps and indels -> they follow different
    # enrichment workflows
    snps = [v for v in variants if VTYPES[0] in v.vtype]  # SNPs
    print_verbosity(f"Fetched {len(snps)} variants in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3])
    return snps

def encode_snp_iupac(
    refnt: str, variant: VariantRecord, debug: bool
) -> Union[str, None]:
    altalleles = "".join(variant.get_altalleles(VTYPES[0]))  # retrieve snps
    if not altalleles:  # variants are all indels
        return None
    if refnt != variant.ref:  # ref alleles must match between vcf and fasta
        exception_handler(
            CrisprHawkEnrichmentError,
            f"Reference allele mismatch between FASTA and VCF file ({refnt} - "
            f"{variant.ref}) at position {variant.position}",
            os.EX_DATAERR,
            debug,
        )
    # encode ref and alt as iupac characters
    return IUPAC_ENCODER["".join([variant.ref, altalleles])]

def insert_snps(
    region: Region, snps: List[VariantRecord], phased: bool, no_filter: bool, verbosity: int, debug: bool,
) -> VariantMap:
    # hashmap variants to their relative position within the sequence
    vmap = VariantMap(phased, debug)  
    # iterate over variants falling in the input region and add variants
    for snp in snps:
        # by default, crispr-hawk discards variants that are not flagged as PASS
        # on filter, however the user may want to consider tem as well
        if not no_filter and snp.filter != "PASS":
            warning(f"Skipping {snp.format()}", verbosity)
            continue
        # compute relative variant position (1-based)
        posrel = snp.position - 1 - region.start
        # compute iupac char, representing the snp
        iupac_nt = encode_snp_iupac(region[posrel], snp, debug)
        if iupac_nt is not None:  # if none, indels -> skip
            region.enrich(posrel, iupac_nt)  # assign iupac char at position
            vmap.insert_variant(posrel, snp)  # map variant to its relative position
    return vmap

def enricher(regions: RegionList, vcfs: Dict[str, VCF], guidelen: int, pamlen: int, no_filter: bool, verbosity: int, debug: bool) -> Dict[Region, VariantMap]:
    variant_maps = {}  # dictionary to map variants in each region
    for region in regions:  # retrieve variants in each region
        print_verbosity(f"Adding variants in {region.display()}", verbosity, VERBOSITYLVL[2])
        snps = fetch_variants(vcfs[region.contig], region, guidelen, verbosity)
        # add variants to regions sequence
        phased = vcfs[region.contig].phased  # assess variants phasing
        variant_maps[region] = insert_snps(region, snps, phased, no_filter, verbosity, debug)
    return variant_maps


# def enricher(
#     regions: RegionList, vcflist: List[str], guidelen: int, pamlen: int, no_filter: bool, verbosity: int, debug: bool,
# ) -> Tuple[Dict[Region, Dict[int, VariantRecord]], bool]:
#     # load vcf files and map each vcf to its contig (assume on vcf per contig)
#     for r in regions:

#         indels_map[r] = indels  # record indels found in the query region
#         # insert variants within region sequence
#         variant_maps[r] = insert_snps(r, snps, no_filter, verbosity, debug)
#     # insert indels within region sequence
#     indelregions, variant_maps = insert_indels(
#         indels_map, variant_maps, guidelen, pamlen, debug
#     )
#     regions.extend(indelregions)  # extend region list to consider indel regions
#     return (
#         regions,
#         variant_maps,
#         vcfs[r.contig].phased,
#     )  # return enriched regions, variant maps, and phasing
