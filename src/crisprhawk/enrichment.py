"""
"""

from variants import VCF, VariantRecord, VTYPES, INDELTYPES
from bedfile import RegionList, Region, IndelRegion, Coordinate
from utils import print_verbosity, warning, VERBOSITYLVL, IUPAC_ENCODER, IUPAC
from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkEnrichmentError
from variant_map import VariantMap

from typing import List, Tuple, Union, Dict
from time import time

import os



def load_vcfs(vcflist: List[str], verbosity: int, debug: bool) -> Dict[str, VCF]:
    # load vcf files and map each vcf to its contig (assume on vcf per contig)
    print_verbosity("Loading VCF files", verbosity, VERBOSITYLVL[2])
    start = time()  # track vcf parsing time
    vcfs = {vcf.contig: vcf for vcf in [VCF(f, verbosity, debug) for f in vcflist]}
    print_verbosity(
        f"Loaded {len(vcfs)} VCFs in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return vcfs


def adjust_region_coords(region: Region, guidelen: int) -> Tuple[int, int]:
    # remove guide length padding from region's start and stop position to
    # extract variants in that range
    return region.start + guidelen, region.stop - guidelen


def fetch_variants(
    vcf: VCF, region: Region, guidelen: int, verbosity: int
) -> Tuple[List[VariantRecord], List[VariantRecord]]:
    # recover variants mapped on the query region
    # each region is padded by |guide| nts to avoid missing guides mapped on
    # region's border
    print_verbosity(
        f"Fetching variants in {region.format(pad=guidelen, string=True)}", verbosity, VERBOSITYLVL[3]
    )
    start = time()  # track variants fecthing time
    variants = vcf.fetch(*adjust_region_coords(region, guidelen))
    # split variants between snps and indels -> they follow different
    # enrichment workflows
    snps = [v for v in variants if VTYPES[0] in v.vtype]  # SNPs
    indels = [v for v in variants if VTYPES[1] in v.vtype]  # indels
    print_verbosity(
        f"Fetched {len(snps) + len(indels)} variants (SNPs: {len(snps)} - INDELs: {len(indels)}) in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return snps, indels


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
    region: Region,
    snps: List[VariantRecord],
    phased: bool,
    no_filter: bool,
    verbosity: int,
    debug: bool,
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
        iupac_nt = encode_snp_iupac(region.sequence[posrel], snp, debug)
        if iupac_nt is not None:  # if none, indels -> skip
            region.enrich(posrel, iupac_nt)  # assign iupac char at position
            vmap.insert_variant(posrel, snp)  # map variant to its relative position
    return vmap

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

def create_indel_region(region: Region, indel: VariantRecord, guidepamlen: int, debug: bool) -> IndelRegion:
    # compute relative indel position
    posrel = indel.position - 1 - region.start
    indel_length = abs(len(indel.ref) - len(indel.alt[0]))  # compute indel length
    # compute start and stop coordinate for indel subregion
    indeltype = guess_indel_type(indel.ref, indel.alt[0], indel.position, debug)
    start, stop = posrel - guidepamlen + 1, posrel + guidepamlen
    if indeltype == INDELTYPES[1]:  # deletion
        stop += indel_length
    refsequence = region.sequence_ref[start:stop]  # extract region sequence
    sequence = region.sequence[start:stop]  # extract region sequence
    coordinates = Coordinate(region.contig, region.start + start, region.start + stop)  # create coordinate object
    # compute offset for insertions and deletions
    offset_ins = len(indel.ref) - 1
    offset_del = len(indel.alt[0]) - 1 
    return IndelRegion(refsequence, sequence, coordinates, offset_ins, offset_del, guidepamlen - 1, indel_length, indeltype, debug)

def compute_vmap(vmap: VariantMap, region: IndelRegion, indel: VariantRecord, offset: int, debug: bool) -> VariantMap:
    vmap_indel = VariantMap(vmap.phased, debug)
    # positions overlapping indel
    for i in range(len(indel.alt[0])):
        vmap_indel.insert_variant(region.indel_pos + i, indel)  # insert indel in vmap
    for i, nt in enumerate(region.sequence):
        if nt in IUPAC[4:]:  # skip ACGT
            position = (region.start + i) - offset
            # adjust position for variant map query
            if i > region.indel_pos:
                position = position - region.indel_length if region.indel_type == INDELTYPES[0] else position + region.indel_length
            vmap_indel.insert_variant(i, vmap[position])
    return vmap_indel


def insert_indels(regions: RegionList, i: int, indels: List[VariantRecord], vmaps: Dict[Region, VariantMap], guidepamlen: int, debug: bool):
    region = regions[i]  # retrieve current region
    for indel_variant in indels:
        for indel in indel_variant.split(VTYPES[1]):
            iregion = create_indel_region(region, indel, guidepamlen, debug)
            iregion.enrich(indel.alt[0])  # insert indel in region 
            regions.append(iregion)  # append indel region to regions list
            vmaps[iregion] = compute_vmap(vmaps[region], iregion, indel, region.start, debug)  # insert indel in vmap
    return vmaps


def enricher(
    regions: RegionList,
    vcfs: Dict[str, VCF],
    guidelen: int,
    pamlen: int,
    no_filter: bool,
    verbosity: int,
    debug: bool,
) -> Dict[Union[Region, IndelRegion], VariantMap]:
    variant_maps = {}  # dictionary to map variants in each region
    regions_num = len(regions)  # regions number
    for i in range(regions_num):
        region = regions[i]  # retrieve region
        print_verbosity(
            f"Adding variants in {region.format(pad=guidelen, string=True)}", verbosity, VERBOSITYLVL[2]
        )
        # retrieve variants in each region
        snps, indels_ = fetch_variants(vcfs[region.contig], region, guidelen, verbosity)        
        # add variants to regions sequence
        phased = vcfs[region.contig].phased  # assess variants phasing
        variant_maps[region] = insert_snps(
            region, snps, phased, no_filter, verbosity, debug
        )
        variant_maps = insert_indels(regions, i, indels_, variant_maps, guidelen + pamlen, debug)
    return variant_maps


