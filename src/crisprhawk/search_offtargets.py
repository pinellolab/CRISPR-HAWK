"""
This module provides functions for performing off-target searches for CRISPR guide
RNAs using CRISPRitz.

It includes utilities to estimate and annotate off-targets for guides across genomic
regions, supporting downstream genome editing analysis.
"""

from .config_crispritz import CrispritzConfig
from .offtargets import estimate_offtargets
from .utils import print_verbosity, VERBOSITYLVL
from .region import Region
from .guide import Guide
from .pam import PAM

from typing import Dict, List, Optional
from time import time


def offtargets_search(
    guides: Dict[Region, List[Guide]],
    pam: PAM,
    crispritz_index: str,
    crispritz_config: CrispritzConfig,
    mm: int,
    bdna: int,
    brna: int,
    annotations: List[str],
    anncolnames: List[str],
    guidelen: int,
    compute_elevation: bool,
    right: bool,
    threads: int,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> Dict[Region, List[Guide]]:
    """Performs off-target search for each guide using CRISPRitz.

    This function estimates off-targets for all guides in each region, updating
    the guide objects with off-target information for downstream analysis.

    Args:
        guides (Dict[Region, List[Guide]]): Dictionary mapping regions to lists
            of Guide objects.
        pam (PAM): PAM object specifying the PAM sequence.
        crispritz_index (str): Path to the CRISPRitz index.
        crispritz_config (CrispritzConfig): Configuration for CRISPRitz.
        mm (int): Maximum number of mismatches.
        bdna (int): Maximum DNA bulge size.
        brna (int): Maximum RNA bulge size.
        annotations (List[str]): List of annotation file paths.
        anncolnames (List[str]): List of annotation column names.
        guidelen (int): Length of the guide sequence.
        compute_elevation (bool): Whether to compute Elevation scores.
        right (bool): Boolean indicating PAM orientation.
        threads (int): Number of threads to use.
        outdir (str): Output directory for results.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        Dict[Region, List[Guide]]: The dictionary of regions with guides annotated
            for off-targets.
    """
    # search off-targets for each retrieved guide
    print_verbosity("Searching off-targets", verbosity, VERBOSITYLVL[1])
    start = time()  # offtargets search start time
    for region, guides_list in guides.items():
        guides[region] = estimate_offtargets(
            guides_list,
            pam,
            crispritz_index,
            region,
            crispritz_config,
            mm,
            bdna,
            brna,
            annotations,
            anncolnames,
            guidelen,
            compute_elevation,
            right,
            threads,
            outdir,
            verbosity,
            debug,
        )
    print_verbosity(
        f"Off-targets search completed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[2],
    )
    return guides
