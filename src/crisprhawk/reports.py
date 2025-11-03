"""This module provides functions for constructing, formatting, and storing
CRISPR guide reports.

It includes utilities for processing guide data, handling annotations, collapsing
duplicate entries, and writing output files for genomic regions. The module is
designed to support flexible reporting for different Cas systems and annotation
options, ensuring that guide reports are comprehensive and well-structured.
"""

from .crisprhawk_argparse import CrisprHawkSearchInputArgs
from .crisprhawk_error import CrisprHawkGuidesReportError
from .exception_handlers import exception_handler
from .pam import PAM, CASX, CPF1, SACAS9, SPCAS9, XCAS9
from .guide import Guide
from .region_constructor import PADDING
from .region import Region
from .utils import (
    print_verbosity,
    VERBOSITYLVL,
    GUIDESREPORTPREFIX,
    IUPACTABLE,
    IUPAC,
    STRAND,
)

from typing import List, Dict, Set, Any
from collections import defaultdict
from time import time

import pandas as pd

import os


# REPORTCOLS = [
#     "chr",  # 0
#     "start",  # 1
#     "stop",  # 2
#     "sgRNA_sequence",  # 3
#     "pam",  # 4
#     "pam_class",  # 5
#     "strand",  # 6
#     "score_azimuth",  # 7
#     "score_rs3",  # 8
#     "score_deepcpf1",  # 9
#     "score_cfdon",  # 10
#     "score_elevationon",  # 11
#     "gc_content",  # 12
#     "out_of_frame",  # 13
#     "origin",  # 14
#     "samples",  # 15
#     "variant_id",  # 16
#     "af",  # 17
#     "target",  # 18
#     "haplotype_id",  # 19
#     "offtargets",  # 20
#     "cfd",  # 21
# ]
REPORTCOLS = [
    "chr",  # 0
    "start",  # 1
    "stop",  # 2
    "sgRNA_sequence",  # 3
    "pam",  # 4
    "pam_class",  # 5
    "strand",  # 6
    "score_azimuth",  # 7
    "score_rs3",  # 8
    "score_deepcpf1",  # 9
    "score_cfdon",  # 10
    "score_elevationon",  # 11
    "gc_content",  # 12
    "origin",  # 13
    "samples",  # 14
    "variant_id",  # 15
    "af",  # 16
    "target",  # 17
    "haplotype_id",  # 18
    "offtargets",  # 19
    "cfd",  # 20
]


def compute_pam_class(pam: PAM) -> str:
    """Returns a string representing the PAM class for the given PAM object.

    The output string uses IUPAC notation, converting ambiguous bases to bracketed sets.

    Args:
        pam (PAM): A PAM object representing the protospacer adjacent motif.

    Returns:
        str: A string representation of the PAM class, e.g. 'NGG' -> '[ACGT]GG'.
    """
    # retrieve a string representing the input pam class
    # e.g. NGG -> [ACGT]GG
    return "".join([nt if nt in IUPAC[:4] else f"[{IUPACTABLE[nt]}]" for nt in pam.pam])


def compute_guide_origin(samples: str) -> str:
    """Determines the origin of the guide based on the provided sample string.

    Returns 'ref' if the guide is from the reference genome, otherwise returns 'alt'.

    Args:
        samples (str): A string indicating the sample source.

    Returns:
        str: 'ref' if the sample is from the reference genome, 'alt' otherwise.
    """
    # compute whether the guide came from reference or alternative genomes
    return "ref" if samples == "REF" else "alt"


def compute_strand_orientation(strand: int) -> str:
    """Determines the strand orientation based on the provided strand value.

    Returns '+' for the forward strand and '-' for the reverse strand.

    Args:
        strand (int): An integer representing the strand orientation.

    Returns:
        str: '+' if the strand is forward, '-' if the strand is reverse.
    """
    # retrieve strand orientation
    return "+" if strand == STRAND[0] else "-"  # retrieve strand orientation


def _update_report_fields_spcas9(
    report: Dict[str, List[Any]],
    region_coordinates: str,
    guide: Guide,
    pamclass: str,
    pamlen: int,
    compute_elevation: bool,
) -> Dict[str, List[str]]:
    """Updates the report dictionary with fields specific to the SpCas9 system.

    Adds guide, PAM, scoring, and annotation information for SpCas9 to the report.

    Args:
        report (Dict[str, List[Any]]): The report dictionary to update.
        region_coordinates (str): String representation of the target region coordinates.
        guide (Guide): The guide object containing sequence and annotation data.
        pamclass (str): The PAM class string for the guide.
        pamlen (int): The length of the PAM sequence.
        compute_elevation (bool): Whether to compute and include elevation scores.

    Returns:
        Dict[str, List[str]]: The updated report dictionary with SpCas9-specific fields.
    """
    # update report fields for spcas9 system pam
    report[REPORTCOLS[1]].append(guide.start)  # start and stop position
    report[REPORTCOLS[2]].append(guide.stop)
    report[REPORTCOLS[3]].append(guide.guide)  # guide sequence
    report[REPORTCOLS[4]].append(guide.pam)  # pam guide
    report[REPORTCOLS[5]].append(pamclass)  # extended pam class
    # strand orientation
    report[REPORTCOLS[6]].append(compute_strand_orientation(guide.strand))
    report[REPORTCOLS[7]].append(guide.azimuth_score)  # azimuth score
    report[REPORTCOLS[8]].append(guide.rs3_score)  # rs3 score
    report[REPORTCOLS[10]].append(guide.cfdon_score)  # cfdon score
    if compute_elevation and (
        guide.guidelen + pamlen == 23 and not guide.right
    ):  # elevationon score
        report[REPORTCOLS[11]].append(guide.elevationon_score)
    report[REPORTCOLS[12]].append(guide.gc)  # gc content
    # report[REPORTCOLS[13]].append(guide.ooframe_score)  # out-of-frame score
    report[REPORTCOLS[13]].append(compute_guide_origin(guide.samples))  # genome
    report[REPORTCOLS[14]].append(guide.samples)  # samples list
    report[REPORTCOLS[15]].append(guide.variants)  # variant ids
    report[REPORTCOLS[16]].append(guide.afs_str)  # variants allele frequencies
    report[REPORTCOLS[17]].append(region_coordinates)  # region
    report[REPORTCOLS[18]].append(guide.hapid)  # haplotype id
    return report


def _update_report_fields_cpf1(
    report: Dict[str, List[Any]],
    region_coordinates: str,
    guide: Guide,
    pamclass: str,
    pamlen: int,
    compute_elevation: bool,
) -> Dict[str, List[str]]:
    """Updates the report dictionary with fields specific to the Cpf1 system.

    Adds guide, PAM, scoring, and annotation information for Cpf1 to the report.

    Args:
        report (Dict[str, List[Any]]): The report dictionary to update.
        region_coordinates (str): String representation of the target region coordinates.
        guide (Guide): The guide object containing sequence and annotation data.
        pamclass (str): The PAM class string for the guide.
        pamlen (int): The length of the PAM sequence.
        compute_elevation (bool): Whether to compute and include elevation scores.

    Returns:
        Dict[str, List[str]]: The updated report dictionary with Cpf1-specific fields.
    """
    # update report fields for cpf1 system pam
    report[REPORTCOLS[1]].append(guide.start)  # start and stop position
    report[REPORTCOLS[2]].append(guide.stop)
    report[REPORTCOLS[3]].append(guide.guide)  # guide sequence
    report[REPORTCOLS[4]].append(guide.pam)  # pam guide
    report[REPORTCOLS[5]].append(pamclass)  # extended pam class
    # strand orientation
    report[REPORTCOLS[6]].append(compute_strand_orientation(guide.strand))
    report[REPORTCOLS[9]].append(guide.deepcpf1_score)  # deepcpf1 score
    if compute_elevation and (
        guide.guidelen + pamlen == 23 and not guide.right
    ):  # elevationon score
        report[REPORTCOLS[11]].append(guide.elevationon_score)
    report[REPORTCOLS[12]].append(guide.gc)  # gc content
    # report[REPORTCOLS[13]].append(guide.ooframe_score)  # out-of-frame score
    report[REPORTCOLS[13]].append(compute_guide_origin(guide.samples))  # genome
    report[REPORTCOLS[14]].append(guide.samples)  # samples list
    report[REPORTCOLS[15]].append(guide.variants)  # variant ids
    report[REPORTCOLS[16]].append(guide.afs_str)  # variants allele frequencies
    report[REPORTCOLS[17]].append(region_coordinates)  # region
    report[REPORTCOLS[18]].append(guide.hapid)  # haplotype id
    return report


def _update_report_fields_other(
    report: Dict[str, List[Any]],
    region_coordinates: str,
    guide: Guide,
    pamclass: str,
    pamlen: int,
    compute_elevation: bool,
) -> Dict[str, List[str]]:
    """Updates the report dictionary with fields specific to the remaining Cas
    systems.

    Adds guide, PAM, scoring, and annotation information for remaining Cas systems
    to the report.

    Args:
        report (Dict[str, List[Any]]): The report dictionary to update.
        region_coordinates (str): String representation of the target region coordinates.
        guide (Guide): The guide object containing sequence and annotation data.
        pamclass (str): The PAM class string for the guide.
        pamlen (int): The length of the PAM sequence.
        compute_elevation (bool): Whether to compute and include elevation scores.

    Returns:
        Dict[str, List[str]]: The updated report dictionary with fields.
    """
    # update report fields for other pam
    report[REPORTCOLS[1]].append(guide.start)  # start and stop position
    report[REPORTCOLS[2]].append(guide.stop)
    report[REPORTCOLS[3]].append(guide.guide)  # guide sequence
    report[REPORTCOLS[4]].append(guide.pam)  # pam guide
    report[REPORTCOLS[5]].append(pamclass)  # extended pam class
    # strand orientation
    report[REPORTCOLS[6]].append(compute_strand_orientation(guide.strand))
    if compute_elevation and (
        guide.guidelen + pamlen == 23 and not guide.right
    ):  # elevationon score
        report[REPORTCOLS[11]].append(guide.elevationon_score)
    report[REPORTCOLS[12]].append(guide.gc)  # gc content
    # report[REPORTCOLS[13]].append(guide.ooframe_score)  # out-of-frame score
    report[REPORTCOLS[13]].append(compute_guide_origin(guide.samples))  # genome
    report[REPORTCOLS[14]].append(guide.samples)  # samples list
    report[REPORTCOLS[15]].append(guide.variants)  # variant ids
    report[REPORTCOLS[16]].append(guide.afs_str)  # variants allele frequencies
    report[REPORTCOLS[17]].append(region_coordinates)  # region
    report[REPORTCOLS[18]].append(guide.hapid)  # haplotype id
    return report


def update_report_fields(
    report: Dict[str, List[Any]],
    region_coordinates: str,
    guide: Guide,
    pam: PAM,
    pamclass: str,
    compute_elevation: bool,
) -> Dict[str, List[str]]:
    """Updates the report dictionary with guide-specific fields based on the Cas
    system.

    Selects the appropriate update function for SpCas9, Cpf1, or other Cas systems
    and updates the report.

    Args:
        report (Dict[str, List[Any]]): The report dictionary to update.
        region_coordinates (str): String representation of the target region coordinates.
        guide (Guide): The guide object containing sequence and annotation data.
        pam (PAM): The PAM object representing the Cas system.
        pamclass (str): The PAM class string for the guide.
        compute_elevation (bool): Whether to compute and include elevation scores.

    Returns:
        Dict[str, List[str]]: The updated report dictionary with guide-specific fields.
    """
    if pam.cas_system in [SPCAS9, XCAS9]:  # spcas9 system pam
        return _update_report_fields_spcas9(
            report, region_coordinates, guide, pamclass, len(pam), compute_elevation
        )
    elif pam.cas_system == CPF1:  # cpf1 system pam
        return _update_report_fields_cpf1(
            report, region_coordinates, guide, pamclass, len(pam), compute_elevation
        )
    return _update_report_fields_other(
        report, region_coordinates, guide, pamclass, len(pam), compute_elevation
    )


def update_optional_report_fields(
    report: Dict[str, List[Any]],
    guide: Guide,
    pam: PAM,
    annotations: List[str],
    gene_annotations: List[str],
    estimate_offtargets: bool,
) -> Dict[str, List[str]]:
    """Updates the report dictionary with optional fields such as annotations
    and off-target scores.

    Adds functional and gene annotations, as well as off-target information, to
    the report if specified.

    Args:
        report (Dict[str, List[Any]]): The report dictionary to update.
        guide (Guide): The guide object containing annotation and off-target data.
        pam (PAM): The PAM object representing the Cas system.
        annotations (List[str]): List of functional annotation strings.
        gene_annotations (List[str]): List of gene annotation strings.
        estimate_offtargets (bool): Whether to include off-target information.

    Returns:
        Dict[str, List[str]]: The updated report dictionary with optional fields.
    """
    reportcols = list(report.keys())
    if annotations:
        idx = reportcols.index(REPORTCOLS[18]) + 1  # haplotype_id is last
        for i, annotation in enumerate(guide.funcann):
            report[reportcols[idx + i]].append(annotation)
    if gene_annotations:
        idx = (
            reportcols.index(REPORTCOLS[18]) + 1 + len(annotations)
        )  # haplotype_id is last
        for i, annotation in enumerate(guide.geneann):
            report[reportcols[idx + i]].append(annotation)
    if estimate_offtargets:
        report[REPORTCOLS[19]].append(guide.offtargets)
        if pam.cas_system in [SPCAS9, XCAS9]:  # spcas9 system pam
            report[REPORTCOLS[20]].append(guide.cfd)
    return report


def insert_elevationon_reportcols(guidepam_len: int, right: bool) -> List[str]:
    """Determines whether to include the elevationon score column in the report.

    Returns the elevationon score column if the guide and PAM length is 23 and
    the guide is not on the right strand.

    Args:
        guidepam_len (int): The combined length of the guide and PAM.
        right (bool): Indicates if the guide is on the right strand.

    Returns:
        List[str]: A list containing the elevationon score column if applicable,
            otherwise an empty list.
    """
    return REPORTCOLS[11:12] if guidepam_len == 23 and not right else []


def insert_annotation_reportcols(
    annotations: List[str],
    gene_annotations: List[str],
    anncolnames: List[str],
    gene_anncolnames: List[str],
) -> List[str]:
    """Generates annotation column names for the report based on provided annotations
    and gene annotations.

    Returns a list of annotation and gene annotation column names, using provided
    names or defaulting to generic names.

    Args:
        annotations (List[str]): List of functional annotation strings.
        gene_annotations (List[str]): List of gene annotation strings.
        anncolnames (List[str]): List of column names for annotations.
        gene_anncolnames (List[str]): List of column names for gene annotations.

    Returns:
        List[str]: A list of annotation and gene annotation column names for the report.
    """
    reportcols = [
        anncolnames[i] if anncolnames else f"annotation_{i + 1}"
        for i, _ in enumerate(annotations)
    ]
    reportcols += [
        gene_anncolnames[i] if gene_anncolnames else f"gene_annotation_{i + 1}"
        for i, _ in enumerate(gene_annotations)
    ]
    return reportcols


def insert_offtargets_reportcols(
    estimate_offtargets: bool, cas_system: int, guidepam_len: int, right: bool
) -> List[str]:
    """Determines which off-target related columns to include in the report.

    Returns a list of off-target and CFD score columns based on the Cas system
    and estimation flag.

    Args:
        estimate_offtargets (bool): Whether to include off-target columns.
        cas_system (int): The Cas system identifier.
        guidepam_len (int): The combined length of the guide and PAM.
        right (bool): Indicates if the guide is on the right strand.

    Returns:
        List[str]: A list of off-target and CFD score column names for the report.
    """
    reportcols = []
    if estimate_offtargets:
        reportcols = REPORTCOLS[20:21]
        if cas_system in [SPCAS9, XCAS9]:  # add CFD score
            reportcols += REPORTCOLS[21:]
    return reportcols


def select_reportcols(
    pam: PAM,
    guidelen: int,
    right: bool,
    annotations: List[str],
    annotation_colnames: List[str],
    gene_annotations: List[str],
    gene_annotation_colnames: List[str],
    estimate_offtargets: bool,
) -> List[str]:
    """Selects and returns the appropriate report column names for the given Cas
    system and options.

    Determines the columns to include in the report based on the Cas system,
    guide length, strand orientation, annotations, gene annotations, and
    off-target estimation settings.

    Args:
        pam (PAM): The PAM object representing the Cas system.
        guidelen (int): The length of the guide sequence.
        right (bool): Indicates if the guide is on the right strand.
        annotations (List[str]): List of functional annotation strings.
        annotation_colnames (List[str]): List of column names for annotations.
        gene_annotations (List[str]): List of gene annotation strings.
        gene_annotation_colnames (List[str]): List of column names for gene annotations.
        estimate_offtargets (bool): Whether to include off-target columns.

    Returns:
        List[str]: A list of column names to be included in the report.
    """
    if pam.cas_system in [SPCAS9, XCAS9]:  # spcas9 system pam report columns
        return (
            REPORTCOLS[:9]
            + REPORTCOLS[10:11]
            + insert_elevationon_reportcols(guidelen + len(pam), right)
            + REPORTCOLS[12:19]
            + insert_annotation_reportcols(
                annotations,
                gene_annotations,
                annotation_colnames,
                gene_annotation_colnames,
            )
            + insert_offtargets_reportcols(
                estimate_offtargets, pam.cas_system, guidelen + len(pam), right
            )
        )
    elif pam.cas_system == CPF1:  # cpf1 system pam report
        return (
            REPORTCOLS[:7]
            + REPORTCOLS[9:10]
            + insert_elevationon_reportcols(guidelen + len(pam), right)
            + REPORTCOLS[12:19]
            + insert_annotation_reportcols(
                annotations,
                gene_annotations,
                annotation_colnames,
                gene_annotation_colnames,
            )
            + insert_offtargets_reportcols(
                estimate_offtargets, pam.cas_system, guidelen + len(pam), right
            )
        )
    return (
        REPORTCOLS[:7]
        + insert_elevationon_reportcols(guidelen + len(pam), right)
        + REPORTCOLS[12:19]
        + insert_annotation_reportcols(
            annotations, gene_annotations, annotation_colnames, gene_annotation_colnames
        )
        + insert_offtargets_reportcols(
            estimate_offtargets, pam.cas_system, guidelen + len(pam), right
        )
    )  # all other pams


def process_data(
    region: Region,
    guides: List[Guide],
    pam: PAM,
    annotations: List[str],
    annotation_colnames: List[str],
    gene_annotations: List[str],
    gene_annotation_colnames: List[str],
    estimate_offtargets: bool,
    compute_elevation: bool,
) -> pd.DataFrame:
    """Processes guide data for a genomic region and constructs a report DataFrame.

    Aggregates guide, annotation, and scoring information for a region into a
    structured DataFrame suitable for reporting.

    Args:
        region (Region): The genomic region for which to process guides.
        guides (List[Guide]): List of Guide objects to include in the report.
        pam (PAM): The PAM object representing the Cas system.
        annotations (List[str]): List of functional annotation strings.
        annotation_colnames (List[str]): List of column names for annotations.
        gene_annotations (List[str]): List of gene annotation strings.
        gene_annotation_colnames (List[str]): List of column names for gene annotations.
        estimate_offtargets (bool): Whether to include off-target columns.
        compute_elevation (bool): Whether to compute and include elevation scores.

    Returns:
        pd.DataFrame: A DataFrame containing the processed report data for the region.
    """
    report = {
        cname: []
        for cname in select_reportcols(
            pam,
            guides[0].guidelen,
            guides[0].right,
            annotations,
            annotation_colnames,
            gene_annotations,
            gene_annotation_colnames,
            estimate_offtargets,
        )
    }  # initialize report dictionary
    pamclass = compute_pam_class(pam)  # compute extended pam class
    region_coordinates = str(region.coordinates)  # target region
    for guide in guides:  # iterate over guides and add to report
        report[REPORTCOLS[0]].append(region.contig)  # region contig (chrom)
        # update report with current guide data
        report = update_report_fields(
            report, region_coordinates, guide, pam, pamclass, compute_elevation
        )
        report = update_optional_report_fields(
            report, guide, pam, annotations, gene_annotations, estimate_offtargets
        )
    report = {c: v for c, v in report.items() if v}  # remove empty columns
    return pd.DataFrame(report)  # build dataframe from report data


def construct_report(
    guides: Dict[Region, List[Guide]],
    pam: PAM,
    annotations: List[str],
    annotation_colnames: List[str],
    gene_annotations: List[str],
    gene_annotation_colnames: List[str],
    estimate_offtargets: bool,
    compute_elevation: bool,
) -> Dict[Region, pd.DataFrame]:
    """Constructs a report DataFrame for each genomic region based on provided
    guides and options.

    Processes all guides for each region and returns a dictionary mapping regions
    to their corresponding report DataFrames.

    Args:
        guides (Dict[Region, List[Guide]]): Dictionary mapping regions to lists
            of Guide objects.
        pam (PAM): The PAM object representing the Cas system.
        annotations (List[str]): List of functional annotation strings.
        annotation_colnames (List[str]): List of column names for annotations.
        gene_annotations (List[str]): List of gene annotation strings.
        gene_annotation_colnames (List[str]): List of column names for gene annotations.
        estimate_offtargets (bool): Whether to include off-target columns.
        compute_elevation (bool): Whether to compute and include elevation scores.

    Returns:
        Dict[Region, pd.DataFrame]: A dictionary mapping each region to its
            processed report DataFrame.
    """
    return {
        region: process_data(
            region,
            guides_list,
            pam,
            annotations,
            annotation_colnames,
            gene_annotations,
            gene_annotation_colnames,
            estimate_offtargets,
            compute_elevation,
        )
        for region, guides_list in guides.items()
    }


def _format_elevationon(reportcols: List[str]) -> List[str]:
    """Returns the elevationon score column if present in the report columns.

    Checks if the elevationon score column should be included based on the current
    report columns.

    Args:
        reportcols (List[str]): The list of current report column names.

    Returns:
        List[str]: A list containing the elevationon score column if present,
            otherwise an empty list.
    """
    return REPORTCOLS[11:12] if REPORTCOLS[11] in reportcols else []


def _format_cfd(reportcols: List[str]) -> List[str]:
    """Returns the CFD score column if present in the report columns.

    Checks if the CFD score column should be included based on the current report
    columns.

    Args:
        reportcols (List[str]): The list of current report column names.

    Returns:
        List[str]: A list containing the CFD score column if present, otherwise
            an empty list.
    """
    return REPORTCOLS[21:22] if REPORTCOLS[21] in reportcols else []


def format_reportcols(
    pam: PAM,
    right: bool,
    annotations: List[str],
    gene_annotations: List[str],
    estimate_offtargets: bool,
    reportcols: List[str],
) -> List[str]:
    """Sorts and arranges report column names for output based on Cas system and
    options.

    Determines the order of columns in the final report according to the Cas system,
    strand orientation, annotations, gene annotations, and off-target estimation
    settings.

    Args:
        pam (PAM): The PAM object representing the Cas system.
        right (bool): Indicates if the guide is on the right strand.
        annotations (List[str]): List of functional annotation strings.
        gene_annotations (List[str]): List of gene annotation strings.
        estimate_offtargets (bool): Whether to include off-target columns.
        reportcols (List[str]): The list of current report column names.

    Returns:
        List[str]: A sorted list of column names for the report.
    """
    # coordinates and score columns
    reportcols_sorted = (
        REPORTCOLS[:3] + REPORTCOLS[4:5] + REPORTCOLS[3:4] + REPORTCOLS[5:7]
        if right
        else REPORTCOLS[:7]
    )
    if pam.cas_system in [SPCAS9, XCAS9]:
        reportcols_sorted += REPORTCOLS[7:9] + REPORTCOLS[10:11]
    elif pam.cas_system == CPF1:
        reportcols_sorted += REPORTCOLS[9:10]
    reportcols_sorted += _format_elevationon(reportcols)
    reportcols_sorted += REPORTCOLS[12:17]  # up to target
    if annotations:
        idx = reportcols.index(REPORTCOLS[18]) + 1
        reportcols_sorted += reportcols[idx : idx + len(annotations)]
    if gene_annotations:
        idx = reportcols.index(REPORTCOLS[18]) + 1 + len(annotations)
        reportcols_sorted += reportcols[idx : idx + len(gene_annotations)]
    if estimate_offtargets:
        reportcols_sorted += REPORTCOLS[19:20]
        reportcols_sorted += _format_cfd(reportcols)
    reportcols_sorted += REPORTCOLS[17:19]
    return reportcols_sorted


def format_report(
    report: pd.DataFrame,
    pam: PAM,
    right: bool,
    annotations: List[str],
    gene_annotations: List[str],
    estimate_offtargets: bool,
) -> pd.DataFrame:
    """Formats the report DataFrame for output by sorting and arranging columns.

    Ensures correct data types, sorts by genomic coordinates, and arranges columns
    according to the Cas system and options.

    Args:
        report (pd.DataFrame): The DataFrame containing guide report data.
        pam (PAM): The PAM object representing the Cas system.
        right (bool): Indicates if the guide is on the right strand.
        annotations (List[str]): List of functional annotation strings.
        gene_annotations (List[str]): List of gene annotation strings.
        estimate_offtargets (bool): Whether to include off-target columns.

    Returns:
        pd.DataFrame: The formatted and sorted report DataFrame.
    """
    # force start and stop to int values - they may be treated as float if
    # concatenated with empty dataframe (e.g. no guide found on + or - strand)
    report[REPORTCOLS[1]] = report[REPORTCOLS[1]].astype(int)
    report[REPORTCOLS[2]] = report[REPORTCOLS[2]].astype(int)
    # reset dataframe index and sort by genomic coordinates
    report = report.reset_index(drop=True)
    report = report.sort_values([REPORTCOLS[1], REPORTCOLS[2]], ascending=True)
    # sort report columns
    reportcols = format_reportcols(
        pam,
        right,
        annotations,
        gene_annotations,
        estimate_offtargets,
        report.columns.tolist(),
    )
    report = report[reportcols]
    return report


def store_report(
    report: pd.DataFrame,
    pam: PAM,
    guidesreport: str,
    right: bool,
    annotations: List[str],
    gene_annotations: List[str],
    estimate_offtargets: bool,
    debug: bool,
) -> str:
    """Stores the report DataFrame as a TSV file at the specified path.

    Formats the report, writes it to disk, and handles file-related exceptions
    during the write process.

    Args:
        report (pd.DataFrame): The DataFrame containing guide report data.
        pam (PAM): The PAM object representing the Cas system.
        guidesreport (str): The file path where the report should be saved.
        right (bool): Indicates if the guide is on the right strand.
        annotations (List[str]): List of functional annotation strings.
        gene_annotations (List[str]): List of gene annotation strings.
        estimate_offtargets (bool): Whether to include off-target columns.
        debug (bool): Whether to enable debug mode for exception handling.

    Returns:
        str: The file path where the report was saved.
    """
    try:
        if not report.empty:
            report = format_report(
                report, pam, right, annotations, gene_annotations, estimate_offtargets
            )  # format report
        report.to_csv(guidesreport, sep="\t", index=False)  # store report
        return guidesreport
    except FileNotFoundError as e:
        exception_handler(
            CrisprHawkGuidesReportError,  # type: ignore
            f"Unable to write to {guidesreport}",
            os.EX_OSERR,
            debug,
            e,
        )
    except PermissionError as e:
        exception_handler(
            CrisprHawkGuidesReportError,  # type: ignore
            f"Permission denied to write {guidesreport}",
            os.EX_OSERR,
            debug,
            e,
        )
    except Exception as e:
        exception_handler(
            CrisprHawkGuidesReportError,  # type: ignore
            f"An unexpected error occurred while writing {guidesreport}",
            os.EX_OSERR,
            debug,
            e,
        )


def _polish_samples_phased(samples: str) -> str:
    """Polishes phased genotype sample strings by merging alleles for each sample.

    For phased genotypes, combines alleles for each sample using the maximum value
    for each allele position.

    Args:
        samples (str): A comma-separated string of sample:genotype pairs, where
            genotypes are phased (e.g., "sample1:1|0,sample2:0|1").

    Returns:
        str: A polished, comma-separated string of sample:genotype pairs with
            merged alleles.
    """
    if "|" not in samples:  # unphased genotype, no need for polishing
        return samples
    samplesmap = defaultdict(lambda: [0, 0])  # initialize samples map
    for e in samples.split(","):  # retrive samples with genotypes
        sample, genotype = e.split(":")
        allele1, allele2 = map(int, genotype.split("|"))  # retrieve allele
        # combine using max (equivalent to OR for binary values)
        samplesmap[sample][0] = max(samplesmap[sample][0], allele1)
        samplesmap[sample][1] = max(samplesmap[sample][1], allele2)
    return ",".join([f"{sample}:{a1}|{a2}" for sample, (a1, a2) in samplesmap.items()])


def collapse_samples(samples: pd.Series) -> str:
    """Collapses a pandas Series of sample genotype strings into a single
    polished string.

    Merges all sample genotype entries, sorts and deduplicates them, and polishes
    phased genotypes if present.

    Args:
        samples (pd.Series): A pandas Series containing sample:genotype strings.

    Returns:
        str: A single, polished string of sample:genotype pairs.
    """
    return (
        ""
        if samples.empty
        else _polish_samples_phased(",".join(sorted(set(",".join(samples).split(",")))))
    )


def parse_variant_ids(variant_ids: str) -> Set[str]:
    """Parses a comma-separated string of variant IDs into a set of unique IDs.

    Converts the input string into a set of variant identifiers, or returns an
    empty set if the input is empty.

    Args:
        variant_ids (str): A comma-separated string of variant IDs.

    Returns:
        Set[str]: A set of unique variant IDs.
    """
    return set(variant_ids.split(",")) if variant_ids else set()


def check_variant_ids(variant_ids_list: List[str]) -> str:
    """Checks and merges variant ID sets from a list of variant ID strings.

    Converts each string to a set of variant IDs, ensures uniqueness, and returns
    a sorted, comma-separated string of IDs.

    Args:
        variant_ids_list (List[str]): A list of comma-separated variant ID strings.

    Returns:
        str: A sorted, comma-separated string of unique variant IDs.
    """
    variant_sets = [parse_variant_ids(vid) for vid in variant_ids_list]
    unique_variant_ids_sets = {tuple(sorted(vs)) for vs in variant_sets}
    return ",".join(sorted(unique_variant_ids_sets.pop()))


def collapse_haplotype_ids(hapids: pd.Series) -> str:
    """Collapses a pandas Series of haplotype ID strings into a single, sorted,
    deduplicated string.

    Merges all haplotype ID entries, sorts and deduplicates them, and returns a
    comma-separated string.

    Args:
        hapids (pd.Series): A pandas Series containing haplotype ID strings.

    Returns:
        str: A single, sorted, comma-separated string of unique haplotype IDs.
    """
    return "" if hapids.empty else ",".join(sorted(set(",".join(hapids).split(","))))


def collapse_annotation(anns: pd.Series) -> str:
    """Collapses a pandas Series of annotation strings into a single, deduplicated
    string.

    Merges all annotation entries, removes duplicates, and returns a comma-separated
    string.

    Args:
        anns (pd.Series): A pandas Series containing annotation strings.

    Returns:
        str: A single, comma-separated string of unique annotations.
    """
    return ",".join(set(",".join(anns).split(",")))


def collapse_offtargets(offtargets: pd.Series) -> int:
    """Collapses a pandas Series of off-target counts into a single integer value.

    Returns the unique off-target count from the series, assuming all entries are
    identical.

    Args:
        offtargets (pd.Series): A pandas Series containing off-target count values.

    Returns:
        int: The unique off-target count value.
    """
    return list(set(offtargets))[0]


def collapse_cfd(cfd: pd.Series) -> str:
    """Collapses a pandas Series of CFD score values into a single string value.

    Returns the unique CFD score from the series as a string, assuming all entries
    are identical.

    Args:
        cfd (pd.Series): A pandas Series containing CFD score values.

    Returns:
        str: The unique CFD score value as a string.
    """
    return str(list(set(cfd))[0])


def collapsed_fields(
    pam: PAM,
    annotations: List[str],
    gene_annotations: List[str],
    estimate_offtargets: bool,
    reportcols: List[str],
) -> Dict[str, str]:
    """Defines aggregation functions for collapsing report DataFrame fields.

    Specifies how each report column should be aggregated when collapsing duplicate
    entries, based on the Cas system and provided options.

    Args:
        pam (PAM): The PAM object representing the Cas system.
        annotations (List[str]): List of functional annotation strings.
        gene_annotations (List[str]): List of gene annotation strings.
        estimate_offtargets (bool): Whether to include off-target columns.
        reportcols (List[str]): The list of current report column names.

    Returns:
        Dict[str, str]: A dictionary mapping column names to aggregation functions
            or methods.
    """
    # mandatory report fields
    fields = {
        "pam_class": "first",  # Assuming pam_class is the same across entries
        "origin": "first",  # Assuming origin does not change
        "samples": collapse_samples,  # Merge sample lists
        "variant_id": check_variant_ids,  # Ensure identical variant_id
        "af": "first",  # Merge variants afs
        "target": "first",  # Keep the first target entry
        "haplotype_id": collapse_haplotype_ids,  # Merge haplotype IDs
    }
    # add optional report fields
    if annotations:  # guides functional annotation
        idx = reportcols.index(REPORTCOLS[18]) + 1  # haplotype_id is last
        for colname in reportcols[idx : idx + len(annotations)]:
            fields[colname] = collapse_annotation
    if gene_annotations:
        idx = (
            reportcols.index(REPORTCOLS[18]) + 1 + len(annotations)
        )  # haplotype_id is last
        for colname in reportcols[idx : idx + len(gene_annotations)]:
            fields[colname] = collapse_annotation
    if estimate_offtargets:
        fields["offtargets"] = collapse_offtargets
        if pam.cas_system in [SPCAS9, XCAS9]:
            fields["cfd"] = collapse_cfd
    return fields


def collapse_report_entries(
    report: pd.DataFrame,
    pam: PAM,
    annotations: List[str],
    gene_annotations: List[str],
    estimate_offtargets: bool,
) -> pd.DataFrame:
    """Collapses duplicate entries in the report DataFrame by grouping and
    aggregating fields.

    Groups the report by relevant columns and applies aggregation functions to
    merge duplicate entries based on the Cas system and provided options.

    Args:
        report (pd.DataFrame): The DataFrame containing guide report data.
        pam (PAM): The PAM object representing the Cas system.
        annotations (List[str]): List of functional annotation strings.
        gene_annotations (List[str]): List of gene annotation strings.
        estimate_offtargets (bool): Whether to include off-target columns.

    Returns:
        pd.DataFrame: The collapsed report DataFrame with duplicates merged.
    """
    # Define the columns to group by
    reportcols = report.columns.tolist()
    group_cols = REPORTCOLS[:5]
    if pam.cas_system in [SPCAS9, XCAS9]:
        group_cols += REPORTCOLS[6:9] + REPORTCOLS[10:11] + REPORTCOLS[12:14]
    elif pam.cas_system == CPF1:
        group_cols += REPORTCOLS[6:7] + REPORTCOLS[9:10] + REPORTCOLS[12:14]
    else:
        group_cols += REPORTCOLS[6:7] + REPORTCOLS[12:14]
    if REPORTCOLS[11] in reportcols:  # elevationon computed
        group_cols += REPORTCOLS[11:12]
    if annotations:
        idx = reportcols.index(REPORTCOLS[18]) + 1  # haplotype_id is last
        group_cols += [reportcols[idx + i] for i, _ in enumerate(annotations)]
    if gene_annotations:
        idx = (
            reportcols.index(REPORTCOLS[18]) + 1 + len(annotations)
        )  # haplotype_id is last
        group_cols += [reportcols[idx + i] for i, _ in enumerate(gene_annotations)]
    if estimate_offtargets:
        group_cols.append(REPORTCOLS[19])
        if REPORTCOLS[21] in reportcols:
            group_cols.append(REPORTCOLS[20])
    return report.groupby(group_cols, as_index=False).agg(
        collapsed_fields(
            pam, annotations, gene_annotations, estimate_offtargets, reportcols
        )
    )


def report_guides(
    guides: Dict[Region, List[Guide]], pam: PAM, args: CrisprHawkSearchInputArgs
) -> Dict[Region, str]:
    """Constructs and stores report files for each genomic region and guide set.

    This function processes guides for each region, constructs report DataFrames,
    collapses duplicate entries, and writes the reports to disk. It returns a
    dictionary mapping each region to its report file path.

    Args:
        guides: Dictionary mapping Region objects to lists of Guide objects.
        pam: PAM object representing the Cas system.
        args: CrisprHawkSearchInputArgs object containing report and output parameters.

    Returns:
        Dictionary mapping Region objects to their corresponding report file paths.
    """
    print_verbosity("Constructing reports", args.verbosity, VERBOSITYLVL[1])
    start = time()  # report construction start time
    reports = construct_report(
        guides,
        pam,
        args.annotations,
        args.annotation_colnames,
        args.gene_annotations,
        args.gene_annotation_colnames,
        args.estimate_offtargets,
        args.compute_elevation,
    )  # construct reports
    reports_fnames = {}  # reports TSVs dictionary
    for region, report in reports.items():  # store reports in output folder
        region_name = (
            f"{region.contig}_{region.start + PADDING}_{region.stop - PADDING}"
        )
        guidesreport = os.path.join(
            args.outdir,
            f"{GUIDESREPORTPREFIX}__{region_name}_{pam}_{args.guidelen}.tsv",
        )
        if not report.empty:
            report = collapse_report_entries(
                report,
                pam,
                args.annotations,
                args.gene_annotations,
                args.estimate_offtargets,
            )
        reports_fnames[region] = store_report(
            report,
            pam,
            guidesreport,
            args.right,
            args.annotations,
            args.gene_annotations,
            args.estimate_offtargets,
            args.debug,
        )  # write report
    print_verbosity(
        f"Reports constructed in {time() - start:.2f}s", args.verbosity, VERBOSITYLVL[2]
    )
    return reports_fnames
