"""Main workflow and utility functions for the CRISPR-HAWK tool.

This module implements the core logic for CRISPR-HAWK, including guide search,
haplotype encoding, VCF conversion, data preparation for CRISPRme, and CRISPRitz
configuration management. It provides high-level functions that orchestrate the
various steps of the CRISPR-HAWK analysis pipeline.
"""

from .crisprhawk_argparse import (
    CrisprHawkSearchInputArgs,
    CrisprHawkConverterInputArgs,
    CrisprHawkPrepareDataInputArgs,
    CrisprHawkCrispritzConfigInputArgs,
)
from .config_crispritz import CrispritzConfig, config_crispritz
from .region_constructor import construct_regions
from .haplotypes import reconstruct_haplotypes
from .haplotype import Haplotype
from .region import Region
from .utils import print_verbosity, VERBOSITYLVL
from .search_guides import search
from .annotation import annotate_guides
from .encoder import encode
from .reports import report_guides
from .graphical_reports import compute_graphical_reports
from .converter import convert_gnomad_vcf
from .crisprme_data import prepare_data_crisprme
from .bitset import Bitset
from .guide import Guide
from .pam import PAM

from typing import Dict, List
from time import time
import sys


def encode_pam(pamseq: str, right: bool, verbosity: int, debug: bool) -> PAM:
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
    pam = PAM(pamseq, right, debug)
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
    """Searches for guide candidates on encoded haplotypes for each genomic region.

    This function performs guide search using the provided PAM, haplotypes, and
    bit-encoded haplotypes, returning a dictionary of guides per region.

    Args:
        pam: The PAM object used for guide search.
        haplotypes: A dictionary mapping regions to lists of Haplotype objects.
        haplotypes_bits: A dictionary mapping regions to lists of bit-encoded haplotypes.
        guidelen: The length of the guide sequence.
        right: Whether the guide is located downstream of the PAM.
        variants_present: Whether variants are present in the haplotypes.
        phased: Whether the haplotypes are phased.
        verbosity: The verbosity level for logging.
        debug: Whether to enable debug mode for error handling.

    Returns:
        Dict[Region, List[Guide]]: A dictionary mapping regions to lists of found
            Guide objects.
    """
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
    """Executes the CRISPR-HAWK guide search workflow for the provided input arguments.

    This function orchestrates region extraction, haplotype reconstruction, guide
    search, annotation, reporting, and optional graphical report generation.

    Args:
        args (CrisprHawkSearchInputArgs): The parsed and validated input arguments
            for the search workflow.
    """
    # extract genomic regions defined in input bed file
    regions = construct_regions(
        args.fastas, args.bedfile, args.fasta_idx, args.verbosity, args.debug
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
    pam = encode_pam(args.pam, args.right, args.verbosity, args.debug)
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
        args.compute_elevation,
        args.estimate_offtargets,
        args.crispritz_config,
        args.mm,
        args.bdna,
        args.brna,
        args.offtargets_annotations,
        args.offtargets_annotation_colnames,
        args.crispritz_index,
        args.guidelen,
        args.right,
        args.outdir,
        args.threads,
        args.verbosity,
        args.debug,
    )
    # construct reports
    reports = report_guides(
        guides,
        args.guidelen,
        pam,
        args.right,
        args.annotations,
        args.annotation_colnames,
        args.gene_annotations,
        args.gene_annotation_colnames,
        args.estimate_offtargets,
        args.compute_elevation,
        args.outdir,
        args.verbosity,
        args.debug,
    )
    # draw graphical reports
    if args.graphical_reports:
        compute_graphical_reports(reports, args.outdir, args.verbosity, args.debug)


def crisprhawk_converter(args: CrisprHawkConverterInputArgs) -> None:
    """Converts gnomAD VCF files for CRISPR-HAWK analysis using the provided arguments.

    This function processes and converts VCF files, adding genotyping fields and
    assigning genotypes based on population-wise allele counts.

    Args:
        args (CrisprHawkConverterInputArgs): The parsed and validated input arguments
            for the VCF conversion workflow.
    """
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
    """Prepares input data for CRISPRme using the provided arguments.

    This function processes the CRISPR-HAWK report and prepares the necessary files
    for CRISPRme analysis.

    Args:
        args (CrisprHawkPrepareDataInputArgs): The parsed and validated input arguments
            for CRISPRme data preparation.
    """
    prepare_data_crisprme(args.report, args.create_pam, args.outdir, args.debug)


def crisprhawk_crispritz_config(args: CrisprHawkCrispritzConfigInputArgs) -> None:
    """Configures CRISPRitz settings for CRISPR-HAWK using the provided arguments.

    This function manages CRISPRitz configuration, including updating, displaying,
    resetting, and validating the configuration file.

    Args:
        args (CrisprHawkCrispritzConfigInputArgs): The parsed and validated input
            arguments for CRISPRitz configuration.
    """
    config = CrispritzConfig()  # load current config file
    if args.env_name or args.targets_dir:  # change config data
        config_crispritz(config, args.env_name, args.targets_dir)
    if args.show:  # display config
        sys.stdout.write(f"Current config:\n{config.show_config()}\n")
    if args.reset:  # reset config file to default values
        sys.stdout.write("Reverting CRISPRitz config file to default values\n")
        config.reset_to_defaults()
    if args.validate:  # validate current config file
        sys.stdout.write("Validating CRISPRitz config file\n")
        config.validate_config()
        sys.stdout.write("CRISPRitz config file correctly validated\n")
