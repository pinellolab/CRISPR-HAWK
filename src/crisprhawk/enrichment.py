"""
"""

from variants import VCF, VariantRecord, VTYPES
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


def enricher(
    regions: RegionList,
    vcflist: List[str],
    guidelen: int,
    no_filter: bool,
    verbosity: int,
    debug: bool,
) -> Tuple[RegionList, Dict[Region, Dict[int, VariantRecord]], bool]:
    # load vcf files and map each vcf to its contig (assume on vcf per contig)
    print_verbosity("Loading VCF files", verbosity, VERBOSITYLVL[2])
    # TODO: check phasing consistency
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
        variants = vcfs[r.contig].fetch(*adjust_region_coords(r, guidelen))
        print_verbosity(
            f"Fetched {len(variants)} variants in {r.format(pad=guidelen)}",
            verbosity,
            VERBOSITYLVL[3],
        )
        # insert variants within region sequence
        variant_maps[r] = insert_variants(r, variants, no_filter, verbosity, debug)
    return regions, variant_maps, vcfs[r.contig].phased  # return enriched regions, variant maps, and phasing


def insert_variants(
    region: Region,
    variants: List[VariantRecord],
    no_filter: bool,
    verbosity: int,
    debug: bool,
) -> Dict[int, VariantRecord]:
    # dictionary to map variants to their relative position within the sequence
    variant_map = {}  
    # iterate over variants falling in the input region and add variants
    for variant in variants:
        # by default, crispr-hawk discards variants that are not flagged as PASS
        # on filter, however the user may want to consider tem as well
        if not no_filter and variant.filter != "PASS":
            warning(variant.format(), verbosity)
            continue
        # compute relative variant position (1-based)
        posrel = variant.position - 1 - region.start
        variant_map[posrel] = variant  # map variant to its relative position
        # compute iupac char, representing the snp
        iupac_nt = encode_snp_iupac(region[posrel], variant, debug)
        if iupac_nt is not None:  # if none, indels -> skip
            region.enrich(posrel, iupac_nt)  # assign iupac char at position
    return variant_map


def encode_snp_iupac(refnt: str, variant: VariantRecord, debug: bool) -> Union[str, None]:
    if variant.vtype == VTYPES[1]:  # variant is indel
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
    return IUPAC_ENCODER["".join([variant.ref, variant.alt])] 
