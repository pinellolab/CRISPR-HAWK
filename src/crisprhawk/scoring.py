"""
This module provides functions for scoring CRISPR guide RNAs using various efficiency
and specificity metrics.

It includes methods to compute Azimuth, RS3, DeepCpf1, Elevation, CFDon, and out-of-frame
scores for guide RNAs.

The module is designed to annotate and enrich guide RNA data with these scores for
downstream genome editing analysis.
"""

from .crisprhawk_error import (
    CrisprHawkAzimuthScoreError,
    CrisprHawkRs3ScoreError,
    CrisprHawkCfdScoreError,
    CrisprHawkDeepCpf1ScoreError,
    CrisprHawkOOFrameScoreError,
)
from .exception_handlers import exception_handler
from .scores import azimuth, rs3, cfdon, deepcpf1, elevationon, ooframe_score
from .utils import flatten_list, print_verbosity, VERBOSITYLVL
from .region import Region
from .guide import Guide, GUIDESEQPAD
from .pam import PAM, SPCAS9, XCAS9, CPF1

from typing import Dict, List, Tuple, Union
from collections import defaultdict
from time import time

import numpy as np

import os


def azimuth_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    """Computes Azimuth scores for a list of guide RNAs.

    This function calculates the Azimuth efficiency score for each guide and updates
    the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with Azimuth scores assigned.
    """
    # create guides np.ndarray required by azimuth; each guide must have 4 nts
    # upstream the guide sequence, and 3 nts downstream the pam
    if not guides:
        return guides
    print_verbosity("Computing Azimuth score", verbosity, VERBOSITYLVL[3])
    start = time()  # azimuth score start time
    guides_seqs = np.array(
        [
            guide.sequence[(GUIDESEQPAD - 4) : (-GUIDESEQPAD + 3)].upper()
            for guide in guides
        ]
    )
    try:  # compute azimuth scores
        azimuth_scores = azimuth(guides_seqs)
    except Exception as e:
        exception_handler(
            CrisprHawkAzimuthScoreError,
            "Azimuth score calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    assert len(azimuth_scores) == len(guides)  # should match
    for i, score in enumerate(azimuth_scores):
        guides[i].azimuth_score = score  # assign score to each guide
    print_verbosity(
        f"Azimuth scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def rs3_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    """Computes RS3 scores for a list of guide RNAs.

    This function calculates the RS3 efficiency score for each guide and updates
    the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with RS3 scores assigned.
    """
    if not guides:
        return guides
    print_verbosity("Computing RS3 score", verbosity, VERBOSITYLVL[3])
    start = time()  # rs3 score start time
    guides_seqs = [
        guide.sequence[(GUIDESEQPAD - 4) : (-GUIDESEQPAD + 3)].upper()
        for guide in guides
    ]
    try:  # compute rs3 scores
        rs3_scores = rs3(guides_seqs)
    except Exception as e:
        exception_handler(
            CrisprHawkRs3ScoreError,
            "RS3 score calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    assert len(rs3_scores) == len(guides)  # should match
    for i, score in enumerate(rs3_scores):
        guides[i].rs3_score = score  # assign score to each guide
    print_verbosity(
        f"RS3 scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def group_guides_position(
    guides: List[Guide], debug: bool
) -> Dict[str, Tuple[Union[None, Guide], List[Guide]]]:
    """Groups guides by their genomic position and strand.

    This function organizes guides into groups based on their start position and
    strand, identifying the reference guide and associated variant guides for each
    group.

    Args:
        guides (List[Guide]): List of Guide objects to group.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        Dict[str, Tuple[Union[None, Guide], List[Guide]]]: A dictionary mapping
            position keys to tuples containing the reference guide and a list of
            guides at that position.
    """

    class _GuideGroup:
        """Helper class to group guides by position and strand.

        This class stores a reference guide and a list of guides for a specific
        genomic position and strand.
        """

        def __init__(self) -> None:
            self._refguide = None
            self._guides = []

        def to_tuple(self) -> Tuple[Union[Guide, None], List[Guide]]:
            return (self._refguide, self._guides)

    pos_guide = defaultdict(_GuideGroup)
    for guide in guides:
        poskey = f"{guide.start}_{guide.strand}"
        if guide.samples == "REF":  # reference guide
            if pos_guide[poskey]._refguide is not None:
                exception_handler(
                    CrisprHawkCfdScoreError,
                    f"Duplicate REF guide at position {guide.start}? CFDon/Elevation-on calculation failed",
                    os.EX_DATAERR,
                    debug,
                )
            pos_guide[poskey]._refguide = guide  # type: ignore
        pos_guide[poskey]._guides.append(guide)
    return {poskey: g.to_tuple() for poskey, g in pos_guide.items()}


def cfdon_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    """Computes CFDon scores for a list of guide RNAs.

    This function calculates the CFDon specificity score for each guide and updates
    the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with CFDon scores assigned.
    """
    print_verbosity("Computing CFDon score", verbosity, VERBOSITYLVL[3])
    start = time()  # cfdon start time
    guide_groups = group_guides_position(guides, debug)  # group guides by positions
    for _, (guide_ref, guides_g) in guide_groups.items():
        try:
            cfdon_scores = cfdon(guide_ref, guides_g, debug)
        except Exception as e:
            exception_handler(
                CrisprHawkCfdScoreError,
                "CFDon score calculation failed",
                os.EX_DATAERR,
                debug,
                e,
            )
        for i, score in enumerate(cfdon_scores):
            guides_g[i].cfdon_score = score
    # revert grouped guides by position into list
    guides = flatten_list([guides_g for _, (_, guides_g) in guide_groups.items()])
    print_verbosity(
        f"CFDon scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def deepcpf1_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    """Computes DeepCpf1 scores for a list of guide RNAs.

    This function calculates the DeepCpf1 efficiency score for each guide and
    updates the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with DeepCpf1 scores assigned.
    """
    if not guides:
        return guides
    print_verbosity("Computing DeepCpf1 score", verbosity, VERBOSITYLVL[3])
    start = time()  # deepcpf1 score start time
    guides_seqs = [
        guide.sequence[(GUIDESEQPAD - 4) : (-GUIDESEQPAD + 3)].upper()
        for guide in guides
    ]
    try:  # compute deepcpf1 scores
        deepcpf1_scores = deepcpf1(guides_seqs)
    except Exception as e:
        exception_handler(
            CrisprHawkDeepCpf1ScoreError,
            "DeepCpf1 score calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    assert len(deepcpf1_scores) == len(guides)  # should match
    for i, score in enumerate(deepcpf1_scores):
        guides[i].deepcpf1_score = score  # assign score to each guide
    print_verbosity(
        f"DeepCpf1 scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def elevationon_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    """Computes Elevation-on scores for a list of guide RNAs.

    This function calculates the Elevation-on efficiency score for each guide and
    updates the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with Elevation-on scores assigned.
    """
    print_verbosity("Computing Elevation-on score", verbosity, VERBOSITYLVL[3])
    start = time()  # cfdon start time
    guide_groups = group_guides_position(guides, debug)  # group guides by positions
    guides = elevationon(guide_groups)
    print_verbosity(
        f"Elevation-on scores computed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return guides


def outofframe_score(
    guides: List[Guide], guidelen: int, right: bool, verbosity: int, debug: bool
) -> List[Guide]:
    """Computes the out-of-frame score for each guide RNA sequence.

    This function calculates the likelihood that a guide induces an out-of-frame
    mutation and updates the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to process.
        guidelen (int): Length of the guide sequence.
        right (bool): Whether the guide is extracted downstream (right side) of
            the PAM.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with out-of-frame scores assigned.
    """
    print_verbosity("Computing out-of-frame score", verbosity, VERBOSITYLVL[3])
    start = time()  # out-of-frame score calculation start time
    try:  # compute out-of-frame score
        idx = GUIDESEQPAD if right else GUIDESEQPAD + guidelen
        # scores = [0] * len(guides)
        scores = ooframe_score(guides, idx)
    except Exception as e:
        exception_handler(
            CrisprHawkOOFrameScoreError,
            "Out-of-frame score calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    for i, score in enumerate(scores):  # set out-of-frame score for each guide
        guides[i].ooframe_score = score
    print_verbosity(
        f"Out-of-frame score computed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return guides


def scoring_guides(
    guides: Dict[Region, List[Guide]],
    pam: PAM,
    compute_elevation: bool,
    guidelen: int,
    right: bool,
    verbosity: int,
    debug: bool,
) -> Dict[Region, List[Guide]]:
    """Scores CRISPR guides using multiple efficiency and specificity metrics.

    This function computes Azimuth, RS3, DeepCpf1, Elevation, and out-of-frame
    scores for each guide, updating the guide objects with the calculated values
    for downstream analysis.

    Args:
        guides (Dict[Region, List[Guide]]): Dictionary mapping regions to lists
            of Guide objects.
        pam (PAM): PAM object specifying the PAM sequence and Cas system.
        compute_elevation (bool): Whether to compute Elevation scores.
        guidelen (int): Length of the guide sequence.
        right (bool): Whether the guide is extracted downstream (right side) of
            the PAM.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        Dict[Region, List[Guide]]: The dictionary of regions with scored Guide
            objects.
    """
    # score guides using azimuth, rs3, deepcpf1, elevation, and out-of-frame scores
    print_verbosity("Scoring guides", verbosity, VERBOSITYLVL[1])
    start = time()  # scoring start time
    for region, guides_list in guides.items():
        if pam.cas_system in [SPCAS9, XCAS9]:  # cas9 system pam
            # score each guide with azimuth
            guides_list = azimuth_score(guides_list, verbosity, debug)
            # score each guide with rs3
            guides_list = rs3_score(guides_list, verbosity, debug)
            # score each guide with CFDon
            guides_list = cfdon_score(guides_list, verbosity, debug)
        if pam.cas_system == CPF1:  # cpf1 system pam
            guides_list = deepcpf1_score(guides_list, verbosity, debug)
        if compute_elevation and (guidelen + len(pam) == 23 and not right):
            # elevation requires 23 bp long sequences, where last 3 bp are pam
            guides_list = elevationon_score(guides_list, verbosity, debug)
        # compute out-of-frame score
        guides_list = outofframe_score(guides_list, guidelen, right, verbosity, debug)
        guides[region] = guides_list  # store scored guides
    print_verbosity(
        f"Scoring completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )
    return guides
