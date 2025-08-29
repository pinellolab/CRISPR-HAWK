""" """

from .crisprhawk_argparse import (
    CrisprHawkSearchInputArgs,
    CrisprHawkConverterInputArgs,
    CrisprHawkPrepareDataInputArgs,
)
from .region_constructor import construct_regions
from .haplotypes import reconstruct_haplotypes
from .haplotype import Haplotype
from .region import Region
from .utils import print_verbosity, VERBOSITYLVL, IUPAC_ENCODER
from .search_guides import search
from .annotate import annotate_guides
from .encoder import encode
from .reports import report_guides
from .converter import convert_gnomad_vcf
from .crisprme_data import prepare_data_crisprme
from .bitset import Bitset
from .guide import Guide
from .pam import PAM

from typing import Dict, List
from time import time


def encode_pam(pamseq: str, verbosity: int, debug: bool) -> PAM:
    """Creates and encodes a PAM object from a given sequence.

    This function constructs a PAM object, encodes its sequence, and logs the
    process.

    Args:
        pamseq: The PAM sequence to encode.
        verbosity: The verbosity level for logging.
        debug: Whether to enable debug mode for error handling.

    Returns:
        PAM: The encoded PAM object.
    """
    # construct pam object
    print_verbosity(f"Creating PAM object for PAM {pamseq}", verbosity, VERBOSITYLVL[1])
    start = time()  # encoding start time
    pam = PAM(pamseq, debug)
    pam.encode(verbosity)  # encode pam sequence
    print_verbosity(
        f"PAM object for PAM {pam} created in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[2],
    )
    return pam


def encode_haplotypes(
    haplotypes: Dict[Region, List[Haplotype]], verbosity: int, debug: bool
) -> Dict[Region, List[List[Bitset]]]:
    """Encodes haplotype sequences into lists of Bitset objects for efficient
    guide search.

    This function processes each haplotype sequence for all regions and encodes
    them into bit representations.

    Args:
        haplotypes: A dictionary mapping regions to lists of Haplotype objects.
        verbosity: The verbosity level for logging.
        debug: Whether to enable debug mode for error handling.

    Returns:
        Dict[Region, List[List[Bitset]]]: A dictionary mapping regions to lists
            of encoded haplotype bitsets.
    """
    # encode haplotypes in bit for efficient guide search
    print_verbosity("Encoding haplotypes in bits", verbosity, VERBOSITYLVL[1])
    start = time()  # encoding start time
    haplotypes_bits = {
        region: [encode(hap.sequence.sequence, verbosity, debug) for hap in haps]
        for region, haps in haplotypes.items()
    }  # encode input haplotypes as sequences of bits
    print_verbosity(
        f"Haplotype encoding completed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[2],
    )
    return haplotypes_bits


def guides_search(
    pam: PAM,
    haplotypes: Dict[Region, List[Haplotype]],
    haplotypes_bits: Dict[Region, List[List[Bitset]]],
    guidelen: int,
    right: bool,
    variants_present: bool,
    phased: bool,
    verbosity: int,
    debug: bool,
) -> Dict[Region, List[Guide]]:
    # search guide candidates on encoded haplotypes
    print_verbosity("Searching guides on haplotypes", verbosity, VERBOSITYLVL[1])
    start = time()  # search start time
    guides = {
        region: search(
            pam,
            region,
            haplotype,
            haplotypes_bits[region],
            guidelen,
            right,
            variants_present,
            phased,
            verbosity,
            debug,
        )
        for region, haplotype in haplotypes.items()
    }
    print_verbosity(
        f"Guides search completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )
    return guides


def crisprhawk_search(args: CrisprHawkSearchInputArgs) -> None:
    # extract genomic regions defined in input bed file
    regions = construct_regions(
        args.fasta, args.bedfile, args.fasta_idx, args.verbosity, args.debug
    )
    # reconstruct haplotypes for each input region
    haplotypes, variants_present, phased = reconstruct_haplotypes(
        args.vcfs,
        regions,
        args.haplotype_table,
        args.outdir,
        args.verbosity,
        args.debug,
    )
    # encode pam and haplotype sequences in bit for efficient guides search
    pam = encode_pam(args.pam, args.verbosity, args.debug)
    haplotypes_bits = encode_haplotypes(haplotypes, args.verbosity, args.debug)
    # search guide candidates within input regions
    guides = guides_search(
        pam,
        haplotypes,
        haplotypes_bits,
        args.guidelen,
        args.right,
        variants_present,
        phased,
        args.verbosity,
        args.debug,
    )
    # annotate guide candidates within each region
    guides = annotate_guides(
        guides,
        args.annotations,
        args.gene_annotations,
        pam,
        args.fasta,
        args.estimate_offtargets,
        args.write_offtargets_report,
        args.outdir,
        args.threads,
        args.verbosity,
        args.debug,
    )
    # construct reports
    report_guides(
        guides,
        args.guidelen,
        pam,
        args.right,
        args.annotations,
        args.annotation_colnames,
        args.gene_annotations,
        args.gene_annotation_colnames,
        args.estimate_offtargets,
        args.outdir,
        args.verbosity,
        args.debug,
    )


def crisprhawk_converter(args: CrisprHawkConverterInputArgs) -> None:
    # convert gnomad vcfs; vcfs keep the fundamental data, but add a genotyping
    # field to the vcf; genotypes are assigned based on population-wise allele
    # counts; populations are treated as samples
    convert_gnomad_vcf(
        args.gnomad_vcfs,
        args.joint,
        args.keep,
        args.suffix,
        args.outdir,
        args.threads,
        args.verbosity,
        args.debug,
    )


def crisprhawk_prepare_data_crisprme(args: CrisprHawkPrepareDataInputArgs) -> None:
    prepare_data_crisprme(args.report, args.create_pam, args.outdir, args.debug)
