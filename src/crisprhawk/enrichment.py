"""
"""

from variants import VCF, VariantRecord, Indel, VTYPES
from bedfile import RegionList, Region
from utils import print_verbosity, warning, VERBOSITYLVL, IUPAC_ENCODER
from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkEnrichmentError

from typing import List, Tuple, Union, Dict

import os


def adjust_region_coords(region: Region, guidelen: int) -> Tuple[int, int]:
    # remove guide length padding from region's start and stop position to
    # extract variants in that range
    return region.start + guidelen, region.stop - guidelen


def fetch_variants(vcf: VCF, region: Region, guidelen: int) -> Tuple[List[VariantRecord], List[VariantRecord], int]:
    # recover variants mapped on the query region
    # each region is padded by |guide| nts to avoid missing guides mapped on  
    # region's border
    variants = vcf.fetch(*adjust_region_coords(region, guidelen))
    fetched_vars = len(variants)  # number of fetched variants 
    # split variants between snps and indels -> they follow different  
    # enrichment workflows
    snps = [v for v in variants if VTYPES[0] in v.vtype]  # SNPs
    indels = [v for v in variants if VTYPES[1] in v.vtype]  # indels
    return snps, indels, fetched_vars


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
        snps, indels, fetched_vars = fetch_variants(vcfs[r.contig], r, guidelen)
        print_verbosity(
            f"Fetched {fetched_vars} variants in {r.format(pad=guidelen)}",
            verbosity,
            VERBOSITYLVL[3],
        )
        # insert variants within region sequence
        variant_maps[r] = insert_snps(r, snps, no_filter, verbosity, debug)
        insert_indels(r, indels, guidelen, pamlen, no_filter, verbosity, debug)
    exit()
    return (
        regions,
        variant_maps,
        vcfs[r.contig].phased,
    )  # return enriched regions, variant maps, and phasing


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



