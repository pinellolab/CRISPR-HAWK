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
from .scoring import scoring_guides
from .search_offtargets import offtargets_search
from .encoder import encode
from .reports import report_guides
from .graphical_reports import compute_graphical_reports
from .converter import convert_gnomad_vcf
from .crisprme_data import prepare_data_crisprme
from .candidate_guides import candidate_guides_analysis
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
    haplotypes: Dict[Region, List[Haplotype]], args: CrisprHawkSearchInputArgs
) -> Dict[Region, List[List[Bitset]]]:
    """Encodes haplotype sequences as bitsets for efficient guide search.

    This function converts each haplotype sequence into a bit-encoded
    representation for all regions, enabling efficient downstream guide search
    operations.

    Args:
        haplotypes: Dictionary mapping Region objects to lists of Haplotype objects.
        args: CrisprHawkSearchInputArgs object containing encoding parameters.

    Returns:
        Dictionary mapping Region objects to lists of bit-encoded haplotypes.
    """
    # encode haplotypes in bit for efficient guide search
    print_verbosity("Encoding haplotypes in bits", args.verbosity, VERBOSITYLVL[1])
    start = time()  # encoding start time
    haplotypes_bits = {
        region: [
            encode(hap.sequence.sequence, args.verbosity, args.debug) for hap in haps
        ]
        for region, haps in haplotypes.items()
    }  # encode input haplotypes as sequences of bits
    print_verbosity(
        f"Haplotype encoding completed in {time() - start:.2f}s",
        args.verbosity,
        VERBOSITYLVL[2],
    )
    return haplotypes_bits


def guides_search(
    pam: PAM,
    haplotypes: Dict[Region, List[Haplotype]],
    haplotypes_bits: Dict[Region, List[List[Bitset]]],
    variants_present: bool,
    phased: bool,
    args: CrisprHawkSearchInputArgs,
) -> Dict[Region, List[Guide]]:
    """Annotates CRISPR guides with variant, functional, gene, and GC content
    information.

    This function processes each region's guides, adding variant, allele frequency,
    reverse complement, GC content, and optional functional and gene annotations,
    returning the updated guides.

    Args:
        guides: Dictionary mapping Region objects to lists of Guide objects.
        args: CrisprHawkSearchInputArgs object containing annotation parameters.

    Returns:
        Dictionary mapping Region objects to lists of annotated Guide objects.
    """
    # search guide candidates on encoded haplotypes
    print_verbosity("Searching guides on haplotypes", args.verbosity, VERBOSITYLVL[1])
    start = time()  # search start time
    guides = {
        region: search(
            pam,
            region,
            haplotype,
            haplotypes_bits[region],
            args.guidelen,
            args.right,
            variants_present,
            phased,
            args.verbosity,
            args.debug,
        )
        for region, haplotype in haplotypes.items()
    }
    print_verbosity(
        f"Guides search completed in {time() - start:.2f}s",
        args.verbosity,
        VERBOSITYLVL[2],
    )
    return guides


def crisprhawk_search(args: CrisprHawkSearchInputArgs) -> None:
    """Executes the main CRISPR-HAWK search workflow using the provided arguments.

    This function orchestrates the guide search, annotation, scoring, off-target
    estimation, report generation, and graphical report creation for the CRISPR-HAWK
    pipeline.

    Args:
        args (CrisprHawkSearchInputArgs): The parsed and validated input arguments
            for the search workflow.
    """
    # extract genomic regions defined in input bed file
    regions = construct_regions(
        args.fastas, args.bedfile, args.fasta_idx, args.verbosity, args.debug
    )
    # reconstruct haplotypes for each input region
    haplotypes, variants_present, phased = reconstruct_haplotypes(regions, args)
    # encode pam and haplotype sequences in bit for efficient guides search
    pam = encode_pam(args.pam, args.right, args.verbosity, args.debug)
    haplotypes_bits = encode_haplotypes(haplotypes, args)
    guides = guides_search(
        pam, haplotypes, haplotypes_bits, variants_present, phased, args
    )  # search guide candidates within input regions
    guides = annotate_guides(guides, args)  # annotate guide candidates
    guides = scoring_guides(guides, pam, args)  # score guide candidates
    if args.estimate_offtargets:  # search off-targets for each guide candidate
        guides = offtargets_search(guides, pam, args)
    reports = report_guides(guides, pam, args)  # construct reports
    if args.graphical_reports:  # draw graphical reports
        compute_graphical_reports(reports, args)
    if args.candidate_guides:  # perform candidate guides analyses
        candidate_guides_analysis(reports, pam, args)


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
