"""
This module provides functions for performing off-target searches for CRISPR guide
RNAs using CRISPRitz.

It includes utilities to estimate and annotate off-targets for guides across genomic
regions, supporting downstream genome editing analysis.
"""

from .crisprhawk_argparse import CrisprHawkSearchInputArgs
from .config_crispritz import CrispritzConfig
from .offtargets import estimate_offtargets
from .utils import print_verbosity, VERBOSITYLVL
from .region import Region
from .guide import Guide
from .pam import PAM

from typing import Dict, List, Optional
from time import time


def offtargets_search(
    guides: Dict[Region, List[Guide]], pam: PAM, args: CrisprHawkSearchInputArgs
) -> Dict[Region, List[Guide]]:
    """Performs off-target search for CRISPR guides using CRISPRitz.

    This function estimates and annotates off-targets for each guide in the provided
    regions, updating the guides with off-target information and returning the
    updated dictionary.

    Args:
        guides: Dictionary mapping Region objects to lists of Guide objects.
        pam: PAM object specifying the protospacer adjacent motif.
        args: CrisprHawkSearchInputArgs object containing search parameters.

    Returns:
        Dictionary mapping Region objects to updated lists of Guide objects with
            off-target information.
    """
    # search off-targets for each retrieved guide
    assert args.crispritz_config  # if here, must be defined
    print_verbosity("Searching off-targets", args.verbosity, VERBOSITYLVL[1])
    start = time()  # offtargets search start time
    for region, guides_list in guides.items():
        guides[region] = estimate_offtargets(
            guides_list,
            pam,
            args.crispritz_index,
            region,
            args.crispritz_config,
            args.mm,
            args.bdna,
            args.brna,
            args.offtargets_annotations,
            args.offtargets_annotation_colnames,
            args.guidelen,
            args.compute_elevation,
            args.right,
            args.threads,
            args.outdir,
            args.verbosity,
            args.debug,
        )
    print_verbosity(
        f"Off-targets search completed in {time() - start:.2f}s",
        args.verbosity,
        VERBOSITYLVL[2],
    )
    return guides
