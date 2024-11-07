"""
"""

from variants import VCF
from bedfile import RegionList, Region
from utils import print_verbosity, warning, VERBOSITYLVL, IUPAC_ENCODER
from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkEnrichmentError

from typing import List, Tuple, Union

import os


def adjust_region_coords(region: Region, guidelen: int) -> Tuple[int, int]:
    # remove guide length padding from region's start and stop position to
    # extract variants in that range
    return region.start + guidelen, region.stop - guidelen


def enricher(regions: RegionList, vcflist: List[str], guidelen: int, no_filter: bool, verbosity: int, debug: bool) -> RegionList:
    # load vcf files and map each vcf to its contig (assume on vcf per contig)
    print_verbosity("Loading VCF files", verbosity, VERBOSITYLVL[2])
    vcfs = {vcf.contig: vcf for vcf in [VCF(f, verbosity, debug) for f in vcflist]}
    print_verbosity(f"Loaded VCF number: {len(vcfs)}", verbosity, VERBOSITYLVL[3])
    # retrieve variants in each region 
    for r in regions:
        print_verbosity(f"Adding variants to region: {r.format(pad=guidelen)}", verbosity, VERBOSITYLVL[2])
        print_verbosity(f"Fetching variants in {r.format(pad=guidelen)}", verbosity, VERBOSITYLVL[3])
        variants = vcfs[r.contig].fetch(*adjust_region_coords(r, guidelen))
        print_verbosity(f"Fetched {len(variants)} variants in {r.format(pad=guidelen)}", verbosity, VERBOSITYLVL[3])
        # insert variants within region sequence
        insert_variants(r, variants, no_filter, verbosity, debug) 
    return regions  # return enriched regions


def insert_variants(region: Region, variants: List[List[str]], no_filter: bool, verbosity: int, debug: bool) -> Region:
    # iterate over veriants falling in the input region and add variants
    for variant in variants:
        # by default, crispr-hawk discards variants that are not flagged as PASS 
        # on filter, however the user may want to consider tem as well
        if not no_filter and variant[6] != "PASS":
            warning(f"Skipping variant {variant[0]}\t{variant[1]}\t{variant[3]}\t{variant[4]}", verbosity)
            continue
        # compute relative variant position (1-based)
        posrel = int(variant[1]) - 1 - region.start  
        ref, alt = variant[3:5]  # retrieve reference and alternative allele
        # compute iupac char, representing the snp
        iupac_nt = encode_snp_iupac(region[posrel], ref, alt, int(variant[1]), debug)
        if iupac_nt is not None:  # if none, indels -> skip
            region.enrich(posrel, iupac_nt)  # assign iupac char at position
    return region


def encode_snp_iupac(refnt: str, ref: str, alt: str, pos: int, debug: bool) -> Union[str, None]:
    if len(ref) > 1:  # deletion, skip
        return 
    alleles_alt = alt.split(",")  # retrieve multiallelic sites
    if len(alleles_alt) == 1 and len(alleles_alt[0]) > 1:  # insertion skip
        return 
    if refnt != ref:  # ref alleles must match between vcf and fasta
        exception_handler(
            CrisprHawkEnrichmentError,
            f"Reference allele mismatch between FASTA and VCF file ({refnt} - "
            f"{ref}) at position {pos}",
            os.EX_DATAERR,
            debug,
        )
    # create iupac string for iupac char encoding - skip indels
    iupac_string = {allele for allele in alleles_alt + [refnt] if len(allele) == 1}
    return IUPAC_ENCODER["".join(iupac_string)]







