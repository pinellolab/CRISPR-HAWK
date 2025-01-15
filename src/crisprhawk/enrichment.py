"""
"""

from variants import VCF, VariantRecord, Indel, VTYPES
from bedfile import RegionList, Region, Coordinate
from utils import print_verbosity, warning, VERBOSITYLVL, IUPAC_ENCODER
from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkEnrichmentError

from typing import List, Tuple, Union, Dict

import os

INDELTYPES = [0, 1]  # indel types -> 0 for insertion, 1 for deletion


def adjust_region_coords(region: Region, guidelen: int) -> Tuple[int, int]:
    # remove guide length padding from region's start and stop position to
    # extract variants in that range
    return region.start + guidelen, region.stop - guidelen


def fetch_variants(vcf: VCF, region: Region, guidelen: int) -> Tuple[List[VariantRecord], List[VariantRecord]]:
    # recover variants mapped on the query region
    # each region is padded by |guide| nts to avoid missing guides mapped on  
    # region's border
    variants = vcf.fetch(*adjust_region_coords(region, guidelen))
    fetched_vars = len(variants)  # number of fetched variants 
    # split variants between snps and indels -> they follow different  
    # enrichment workflows
    snps = [v for v in variants if VTYPES[0] in v.vtype]  # SNPs
    indels = [v for v in variants if VTYPES[1] in v.vtype]  # indels
    return snps, indels


def enricher(
    regions: RegionList,
    vcflist: List[str],
    guidelen: int,
    pamlen: int,
    no_filter: bool,
    verbosity: int,
    debug: bool,
) -> Tuple[RegionList, Dict[Region, Dict[int, VariantRecord]], bool]:
    # load vcf files and map each vcf to its contig (assume on vcf per contig)
    print_verbosity("Loading VCF files", verbosity, VERBOSITYLVL[2])
    vcfs = {vcf.contig: vcf for vcf in [VCF(f, verbosity, debug) for f in vcflist]}
    print_verbosity(f"Loaded VCF number: {len(vcfs)}", verbosity, VERBOSITYLVL[3])
    variant_maps = {}  # dictionary to map variants in each region
    variants_indels = {}  # dictionary to store indels in each region 
    # retrieve variants in each region
    for r in regions:
        print_verbosity(
            f"Adding variants to region: {r.format(pad=guidelen)}",
            verbosity,
            VERBOSITYLVL[2],
        )
        print_verbosity(
            f"Fetching variants in {r.format(pad=guidelen)}", verbosity, VERBOSITYLVL[3]
        )
        snps, indels = fetch_variants(vcfs[r.contig], r, guidelen)
        print_verbosity(
            f"Fetched {len(snps) + len(indels)} (snps: {len(snps)} - indels: {len(indels)}) variants in {r.format(pad=guidelen)}",
            verbosity,
            VERBOSITYLVL[3],
        )
        # insert variants within region sequence
        variant_maps[r] = insert_snps(r, snps, no_filter, verbosity, debug)
        # record indels found in the query region
        variants_indels[r] = indels
    process_indels(variants_indels, guidelen, pamlen, debug)    # insert_indels(r, indels, guidelen, pamlen, no_filter, verbosity, debug)
    
    exit()
    return (
        regions,
        variant_maps,
        vcfs[r.contig].phased,
    )  # return enriched regions, variant maps, and phasing


def guess_indel_type(ref: str, alt: str, position: int, debug) -> int:
    # guess from ref and alt alleles if the input indel is an insertion 
    # (|ref| < |alt|) deletion (|ref| > |alt|)
    if len(ref) == len(alt):  # edge case that should never happen
        exception_handler(CrisprHawkEnrichmentError, f"Reference and alternative indel alleles at position {position} have the same length ({len(ref)})", os.EX_DATAERR, debug)
    return INDELTYPES[0] if len(ref) < len(alt) else INDELTYPES[1] 


def create_indel_region(sequence: str, startrel: int, stoprel: int, region_start: int, chrom: str, debug: bool) -> Region:
    # recover original start and stop coordinates for indel subregion
    start, stop = startrel + region_start, stoprel + region_start
    # create region object for each indel
    return Region(sequence, Coordinate(chrom, start, stop), debug)


def retrieve_indel_subregion(region: Region, indel: VariantRecord, alt: str, guidelen: int, pamlen: int, debug: bool) -> Region:
    # compute indel location in the current region
    posrel = indel.position - 1 - region.start
    # compute start and stop coordinate for indel subregion
    start, stop = posrel - guidelen, posrel + guidelen + pamlen - 1
    indel_length = abs(len(indel.ref) - len(alt))  # compute indel length    
    # adjust indel stop position according to the indel type
    if guess_indel_type(indel.ref, indel.alt, indel.position, debug) == INDELTYPES[0]:  # insertion
        stop -= indel_length - 1
    else:  # deletion
        stop += indel_length + 1
    return create_indel_region("".join(region[start:stop]), start, stop, region.start, indel.contig, debug)  # retrieve str from list


def process_indels(variant_indels: Dict[Region, List[VariantRecord]], guidelen: int, pamlen: int, debug: bool):
    indelregions = [
        retrieve_indel_subregion(r, indel, altallele, guidelen, pamlen, debug)
        for r, indels in variant_indels.items()
        for indel in indels
        for altallele in indel.get_altalleles(VTYPES[1])
    ]
    for i in indelregions:
        print(i.format())
        print(i._sequence)
            


def insert_snps(
    region: Region,
    snps: List[VariantRecord],
    no_filter: bool,
    verbosity: int,
    debug: bool,
) -> Dict[int, VariantRecord]:
    # dictionary to map variants to their relative position within the sequence
    variant_map = {}
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
            variant_map[posrel] = snp  # map variant to its relative position
    return variant_map

def insert_indels(region: Region, indels: List[VariantRecord], guidelen: int, pamlen: int, no_filter: bool, verbosity: int, debug: bool):
    for indel in indels:
        posrel = indel.position - 1 - region.start
        print(indel.format())
        for alt in indel.get_altalleles(VTYPES[1]):
            x = Indel(region, indel.ref, alt, posrel, guidelen, pamlen, debug)
            print(x.refsequence)
            print(x.indelsequence)
            print()


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



