"""Provides functions for estimating and reporting CRISPR guide off-targets using
CRISPRitz.

This module handles the creation of input files, execution of off-target searches,
scoring, and annotation of guides with off-target information.
"""

from .crisprhawk_error import (
    CrisprHawkOffTargetsError,
    CrisprHawkCfdScoreError,
    CrisprHawkElevationScoreError,
)
from .config_crispritz import CrispritzConfig, CRISPRITZ
from .exception_handlers import exception_handler
from .utils import (
    print_verbosity,
    suppress_stdout,
    suppress_stderr,
    VERBOSITYLVL,
)
from .region_constructor import PADDING
from .scores.cfdscore.cfdscore import load_mismatch_pam_scores
from .scores.crisprhawk_scores import elevation
from .offtarget import Offtarget
from .bedfile import BedAnnotation
from .region import Region
from .guide import Guide
from .pam import PAM, SPCAS9, XCAS9

from typing import List, Tuple, Dict, Set
from time import time

import pandas as pd
import numpy as np

import subprocess
import os


# off-targets report column names
OTREPCNAMES = [
    "chrom",  # 0
    "position",  # 1
    "strand",  # 2
    "grna",  # 3
    "spacer",  # 4
    "pam",  # 5
    "mm",  # 6
    "bulge_size",  # 7
    "bulg_type",  # 8
    "cfd",  # 9
    "elevation",  # 10
]


def _filter_guides(guides: List[Guide]) -> Set[str]:
    """Returns a set of unique guide sequences in uppercase.
    
    This function ensures that each guide is only searched once by removing duplicates.

    Args:
        guides (List[Guide]): List of Guide objects.

    Returns:
        Set[str]: Set of unique guide sequences in uppercase.
    """
    # avoid searching multiple times the same guide
    return {g.guide.upper() for g in guides}


def _write_guides_file(
    guides: Set[str], pam: PAM, crispritz_dir: str, right: bool, verbosity: int, debug: bool
) -> str:
    """Creates a guides file for off-target estimation using CRISPRitz.
    Formats and writes each guide sequence to a file for downstream analysis.

    Args:
        guides (Set[str]): Set of unique guide sequences.
        pam (PAM): PAM object specifying the PAM sequence.
        crispritz_dir (str): Directory where the guides file will be created.
        right (bool): Boolean indicating PAM orientation.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode.

    Returns:
        str: The path to the created guides file.
    """
    guides_fname = os.path.join(crispritz_dir, "guides.txt")  # guides file
    print_verbosity(
        "Creating guides file for off-target estimation",
        verbosity,
        VERBOSITYLVL[3],
    )
    pamseq = "N" * len(pam)  # pam sequence in guide
    try:
        with open(guides_fname, mode="w") as outfile:
            for guide in guides:  # format each guide for crispritz input
                guide_f = f"{pamseq}{guide}" if right else f"{guide}{pamseq}"
                outfile.write(f"{guide_f}\n")  # write guide
    except (IOError, Exception) as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed writing crispritz guide file {guides_fname}",
            os.EX_DATAERR,
            debug,
            e,
        )
    assert os.stat(guides_fname).st_size > 0
    return guides_fname


def _write_pam_file(
    pam: PAM,
    guidelen: int,
    right: bool,
    crispritz_dir: str,
    verbosity: int,
    debug: bool,
) -> str:
    """Creates a PAM file for off-target estimation using CRISPRitz.
    The PAM sequence is formatted and written to a file for downstream analysis.

    Args:
        pam (PAM): PAM object specifying the PAM sequence.
        guidelen (int): Length of the guide sequence.
        right (bool): Boolean indicating PAM orientation.
        crispritz_dir (str): Directory where the PAM file will be created.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode.

    Returns:
        str: The path to the created PAM file.
    """
    pam_fname = os.path.join(crispritz_dir, "pam.txt")  # pam file
    print_verbosity(
        "Creating PAM file for off-targets estimation", verbosity, VERBOSITYLVL[3]
    )
    try:
        with open(pam_fname, mode="w") as outfile:
            gseq = "N" * guidelen  # guide sequence in pam
            pam_f = f"{pam}{gseq} {-len(pam)}" if right else f"{gseq}{pam} {len(pam)}"
            outfile.write(f"{pam_f}\n")  # write pam
    except (IOError, Exception) as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed writing crispritz PAM file {pam_fname}",
            os.EX_DATAERR,
            debug,
            e,
        )
    assert os.stat(pam_fname).st_size > 0
    return pam_fname


def _prepare_input_data(
    crispritz_config: CrispritzConfig,
    guides: Set[str],
    pam: PAM,
    outdir: str,
    right: bool,
    verbosity: int,
    debug: bool,
) -> Tuple[str, str]:
    """Prepares input files for CRISPRitz off-target estimation.
    Creates the necessary guides and PAM files in the appropriate directory.

    Args:
        crispritz_config (CrispritzConfig): Configuration for CRISPRitz.
        guides (Set[str]): Set of unique guide sequences.
        pam (PAM): PAM object specifying the PAM sequence.
        outdir (str): Output directory for generated files.
        right (bool): Boolean indicating PAM orientation.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode.

    Returns:
        Tuple[str, str]: Paths to the created guides and PAM files.
    """
    if crispritz_config.outdir == ".crispritz_targets":
        # create hidden folder within output directory
        crispritz_dir = os.path.join(outdir, crispritz_config.outdir)
    else:  # use directory reported in the config file
        crispritz_dir = crispritz_config.outdir
    if not os.path.isdir(crispritz_dir):  # stores crispritz targets
        os.makedirs(crispritz_dir)
    # create guides and pam files
    guides_fname = _write_guides_file(guides, pam, crispritz_dir, right, verbosity, debug)
    pam_fname = _write_pam_file(
        pam, len(list(guides)[0]), right, crispritz_dir, verbosity, debug
    )
    return guides_fname, pam_fname


def _format_targets_prefix(
    region: Region, pam: PAM, guidelen: int, crispritz_dir: str
) -> str:
    """Generates a file prefix for CRISPRitz off-target search results.
    The prefix encodes region, PAM, and guide length information for output files.

    Args:
        region (Region): Genomic region for off-target search.
        pam (PAM): PAM object specifying the PAM sequence.
        guidelen (int): Length of the guide sequence.
        crispritz_dir (str): Directory for CRISPRitz output files.

    Returns:
        str: The formatted file prefix for CRISPRitz results.
    """
    return os.path.join(
        crispritz_dir,
        f"{region.contig}_{region.start + PADDING}_{region.stop - PADDING}_{pam.pam}_{guidelen}",
    )


def search(
    crispritz_config: CrispritzConfig,
    crispritz_index: str,
    region: Region,
    pam: PAM,
    guidelen: int,
    guide_fname: str,
    pam_fname: str,
    mm: int,
    bdna: int,
    brna: int,
    threads: int,
    crispritz_dir: str,
    verbosity: int,
    debug: bool,
) -> str:
    """Runs CRISPRitz to search for off-targets in the specified region.
    Executes the CRISPRitz command and returns the path to the results file.

    Args:
        crispritz_config (CrispritzConfig): Configuration for CRISPRitz.
        crispritz_index (str): Path to the CRISPRitz index.
        region (Region): Genomic region for off-target search.
        pam (PAM): PAM object specifying the PAM sequence.
        guidelen (int): Length of the guide sequence.
        guide_fname (str): Path to the guides file.
        pam_fname (str): Path to the PAM file.
        mm (int): Maximum number of mismatches.
        bdna (int): Maximum DNA bulge size.
        brna (int): Maximum RNA bulge size.
        threads (int): Number of threads to use.
        crispritz_dir (str): Directory for CRISPRitz output files.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode.

    Returns:
        str: The path to the CRISPRitz off-targets results file.
    """
    start = time()
    print_verbosity("Searching off-targets with CRISPRitz", verbosity, VERBOSITYLVL[2])
    # define crispritz command
    targets_prefix = _format_targets_prefix(region, pam, guidelen, crispritz_dir)
    crispritz_run = (
        f"{crispritz_config.conda} run -n {crispritz_config.env_name} {CRISPRITZ} "
        f"search {crispritz_index} {pam_fname} {guide_fname} {targets_prefix} -mm "
        f"{mm} -bDNA {bdna} -bRNA {brna} -th {threads} -r"
    )
    try:  # run crispritz to search for off-targets
        with suppress_stdout(), suppress_stderr():
            # suppress stdout and stderr to avoid cluttering the output
            subprocess.check_call(
                crispritz_run,
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
    except subprocess.CalledProcessError as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed to search offtargets for guides in {guide_fname}",
            os.EX_DATAERR,
            debug,
            e,
        )
    targets_fname = f"{targets_prefix}.targets.txt"
    assert os.stat(targets_fname).st_size > 0
    print_verbosity(
        f"Off-targets search with CRISPRitz completed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[2],
    )
    return targets_fname


def _read_offtargets(
    crispritz_targets_file: str, pam: PAM, right: bool, debug: bool
) -> List[Offtarget]:
    """Reads CRISPRitz off-targets from a results file.
    Parses each line and returns a list of Offtarget objects for downstream analysis.

    Args:
        crispritz_targets_file (str): Path to the CRISPRitz targets file.
        pam (PAM): PAM object specifying the PAM sequence.
        right (bool): Boolean indicating PAM orientation.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Offtarget]: List of parsed Offtarget objects.
    """
    try:
        with open(crispritz_targets_file, mode="r") as infile:
            infile.readline()  # skip header
            offtargets = [
                Offtarget(line, pam.pam, right, debug) for line in infile
            ]  # read CRISPRitz off-targets report
    except (IOError, Exception) as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed retrieving CRISPRitz off-targets in {crispritz_targets_file}",
            os.EX_DATAERR,
            debug,
            e,
        )
    return offtargets


def _compute_cfd_score(
    offtargets: List[Offtarget], verbosity: int, debug: bool
) -> List[Offtarget]:
    """Computes the CFD score for a list of off-targets.
    Updates each Offtarget object with its calculated CFD score.

    Args:
        offtargets (List[Offtarget]): List of Offtarget objects to score.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Offtarget]: The list of Offtarget objects with updated CFD scores.
    """
    print_verbosity(
        f"Computing CFD score for {len(offtargets)} off-targets",
        verbosity,
        VERBOSITYLVL[3],
    )
    start = time()
    mmscores, pamscores = load_mismatch_pam_scores(debug)
    try:
        for ot in offtargets:
            ot.compute_cfd(mmscores, pamscores)
    except Exception as e:
        exception_handler(
            CrisprHawkCfdScoreError,
            "Off-targets CFD scoring failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    print_verbosity(
        f"CFD score computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return offtargets


def _compute_elevation_score(
    offtargets: List[Offtarget], verbosity: int, debug: bool
) -> List[Offtarget]:
    """Computes the Elevation score for a list of off-targets.
    Updates each Offtarget object with its calculated Elevation score, skipping
    those with bulges.

    Args:
        offtargets (List[Offtarget]): List of Offtarget objects to score.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Offtarget]: The list of Offtarget objects with updated Elevation scores.
    """
    print_verbosity(
        f"Computing Elevation score for {len(offtargets)} off-targets",
        verbosity,
        VERBOSITYLVL[3],
    )
    start = time()
    # create wildtype and offtarget lists
    wildtypes_list, offtargets_list = [], []
    offtargets_scored, offtargets_filt = [], []
    for ot in offtargets:  # if bulge present skip and keep NA on elevation
        if "-" not in ot.grna and "-" not in ot.spacer:
            wildtypes_list.append(ot.grna.upper())  # wildtypes
            offtargets_list.append(ot.spacer.upper())  # offtargets
            offtargets_scored.append(ot)  # offtargets scored with elevation
        else:
            offtargets_filt.append(ot)
    try:
        for i, score in enumerate(elevation(wildtypes_list, offtargets_list)):
            offtargets_scored[i].elevation = score
    except Exception as e:
        exception_handler(
            CrisprHawkElevationScoreError,
            "Off-targets Elevation scoring failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    print_verbosity(
        f"Elevation score computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return offtargets_filt + offtargets_scored


def report_offtargets(
    crispritz_targets_file: str,
    region: Region,
    pam: PAM,
    guidelen: int,
    compute_elevation: bool,
    right: bool,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> List[Offtarget]:
    """Generates an off-targets report for a genomic region using CRISPRitz results.
    Computes CFD and Elevation scores as needed, writes a tab-separated report,
    and returns the list of Offtarget objects.

    Args:
        crispritz_targets_file (str): Path to the CRISPRitz targets file.
        region (Region): Genomic region for off-target analysis.
        pam (PAM): PAM object specifying the PAM sequence.
        guidelen (int): Length of the guide sequence.
        compute_elevation (bool): Whether to compute Elevation scores.
        right (bool): Boolean indicating PAM orientation.
        outdir (str): Output directory for the report.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Offtarget]: List of Offtarget objects included in the report.
    """
    # write haplotypes table to file
    start = time()  # track haplotypes table writing time
    offtargets = _read_offtargets(crispritz_targets_file, pam, right, debug)
    if pam.cas_system in [SPCAS9, XCAS9]:  # compute CFD score
        offtargets = _compute_cfd_score(offtargets, verbosity, debug)
    if compute_elevation and (
        guidelen + len(pam) == 23 and not right
    ):  # compute elevation score
        offtargets = _compute_elevation_score(offtargets, verbosity, debug)
    print_verbosity("Writing off-targets report", verbosity, VERBOSITYLVL[1])
    report_fname = os.path.join(
        outdir,
        f"offtargets_{region.contig}_{region.start + PADDING}_{region.stop - PADDING}.tsv",
    )
    try:
        with open(report_fname, mode="w") as outfile:
            outfile.write("\t".join(OTREPCNAMES) + "\n")
            outfile.write("\n".join([ot.report_line() for ot in offtargets]))
        ot_table = pd.read_csv(report_fname, sep="\t")
        ot_table = ot_table.sort_values(OTREPCNAMES[:2])
        ot_table.to_csv(report_fname, sep="\t", index=False, na_rep="NA")
    except OSError as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed writing off-targets report for region {region}",
            os.EX_IOERR,
            debug,
            e,
        )
    print_verbosity(
        f"Off-targets report written in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[2],
    )
    return offtargets


def _calculate_offtargets_map(
    offtargets: List[Offtarget], guides: List[Guide]
) -> Dict[str, List[Offtarget]]:
    """Creates a mapping from guide sequences to their associated off-targets.
    Returns a dictionary where each guide maps to a list of its corresponding
    Offtarget objects.

    Args:
        offtargets (List[Offtarget]): List of Offtarget objects.
        guides (List[Guide]): List of Guide objects.

    Returns:
        Dict[str, List[Offtarget]]: Mapping from guide sequence to list of
            Offtarget objects.
    """
    otmap = {g.guide.upper(): [] for g in guides}  # offtargets map
    for ot in offtargets:
        # add each spacer to the corresponding guide (no pam)
        otmap[ot.grna_.upper().replace("-", "")].append(ot)
    return otmap


def _calculate_global_cfd(offtargets: List[Offtarget]) -> float:
    """Calculates the global CFD score for a set of off-targets.
    Returns a summary score representing the overall off-target potential for a guide.

    Args:
        offtargets (List[Offtarget]): List of Offtarget objects.

    Returns:
        float: The calculated global CFD score.
    """
    cfds = [0 if ot.cfd == "NA" else float(ot.cfd) for ot in offtargets]
    return 100 / (100 + sum(cfds))


def annotate_guides_offtargets(
    offtargets: List[Offtarget], guides: List[Guide], verbosity: int
) -> List[Guide]:
    """Annotates each guide with its number of off-targets and global CFD score.
    Updates Guide objects with off-target counts and summary CFD values.

    Args:
        offtargets (List[Offtarget]): List of Offtarget objects.
        guides (List[Guide]): List of Guide objects to annotate.
        verbosity (int): Verbosity level for logging.

    Returns:
        List[Guide]: The list of annotated Guide objects.
    """
    otmap = _calculate_offtargets_map(offtargets, guides)
    start = time()
    print_verbosity("Computing guides global CFD and annotating guides with estimated off-target numbers", verbosity, VERBOSITYLVL[3])
    for guide in guides:
        # set off-targets number and set global CFD
        guide.offtargets = len(otmap[guide.guide.upper()])
        guide.cfd = _calculate_global_cfd(otmap[guide.guide.upper()])
    print_verbosity(
        f"Guides global CFDs and off-target numbers annotation computed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return guides


def search_offtargets(
    guides: List[Guide],
    pam: PAM,
    crispritz_index: str,
    region: Region,
    crispritz_config: CrispritzConfig,
    mm: int,
    bdna: int,
    brna: int,
    guidelen: int,
    compute_elevation: bool,
    right: bool,
    threads: int,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> List[Guide]:
    """Estimates off-targets for a list of CRISPR guides using CRISPRitz.
    Runs the off-target search, computes scores, annotates guides, and returns 
    the updated list.

    Args:
        guides (List[Guide]): List of Guide objects to analyze.
        pam (PAM): PAM object specifying the PAM sequence.
        crispritz_index (str): Path to the CRISPRitz index.
        region (Region): Genomic region for off-target search.
        crispritz_config (CrispritzConfig): Configuration for CRISPRitz.
        mm (int): Maximum number of mismatches.
        bdna (int): Maximum DNA bulge size.
        brna (int): Maximum RNA bulge size.
        guidelen (int): Length of the guide sequence.
        compute_elevation (bool): Whether to compute Elevation scores.
        right (bool): Boolean indicating PAM orientation.
        threads (int): Number of threads to use.
        outdir (str): Output directory for results.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode.

    Returns:
        List[Guide]: The list of Guide objects annotated with off-target information.
    """
    # remove duplicated guides to avoid multiple searches of the same guide
    guides_seqs = _filter_guides(guides)
    guides_fname, pam_fname = _prepare_input_data(
        crispritz_config, guides_seqs, pam, outdir, right, verbosity, debug
    )  # prepare input data for crispritz
    print_verbosity("Estimating off-targets for found guides", verbosity, VERBOSITYLVL[3])
    start = time()
    # search offtargets with crispritz
    targets_fname = search(
        crispritz_config,
        crispritz_index,
        region,
        pam,
        guides[0].guidelen,
        guides_fname,
        pam_fname,
        mm,
        bdna,
        brna,
        threads,
        os.path.dirname(guides_fname),
        verbosity,
        debug,
    )
    offtargets = report_offtargets(
        targets_fname,
        region,
        pam,
        guidelen,
        compute_elevation,
        right,
        outdir,
        verbosity,
        debug,
    )
    guides = annotate_guides_offtargets(offtargets, guides, verbosity)
    print_verbosity(
        f"Off-targets estimation completed in {(time() - start):.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return guides
