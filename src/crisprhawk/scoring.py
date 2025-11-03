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
from .crisprhawk_argparse import CrisprHawkSearchInputArgs
from .exception_handlers import exception_handler
from .scores import azimuth, rs3, cfdon, deepcpf1, elevationon, ooframe_score
from .utils import calculate_chunks, flatten_list, print_verbosity, VERBOSITYLVL
from .region import Region
from .guide import Guide, GUIDESEQPAD
from .pam import PAM, SPCAS9, XCAS9, CPF1

from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple, Union
from collections import defaultdict
from time import time

import numpy as np

import os


def _extract_guide_sequences(guides: List[Guide]) -> List[str]:
    """Extracts guide RNA sequences with required flanking nucleotides.

    This function returns a list of guide sequences, each including 4 nucleotides
    upstream and 3 nucleotides downstream of the PAM, formatted in uppercase.

    Args:
        guides (List[Guide]): List of Guide objects to extract sequences from.

    Returns:
        List[str]: List of formatted guide sequences with flanking nucleotides.
    """
    return [
        guide.sequence[(GUIDESEQPAD - 4) : (-GUIDESEQPAD + 3)].upper()
        for guide in guides
    ]


def _azimuth(guides_chunks: Tuple[int, List[str]]) -> Tuple[int, List[float]]:
    """Calculates Azimuth scores for a chunk of guide RNAs.

    This function computes Azimuth efficiency scores for a given chunk of guides
    and returns the starting index along with the list of scores.

    Args:
        guides_chunks (Tuple[int, List[Guide]]): A tuple containing the start
            index and a list of Guide objects.

    Returns:
        Tuple[int, List[float]]: The start index and the list of Azimuth scores
            for the guides.
    """
    start_idx, guides = guides_chunks
    scores = azimuth(np.array(guides))
    return start_idx, scores


def _execute_azimuth(
    guides_chunks: List[Tuple[int, List[str]]], size: int, threads: int, debug: bool
) -> List[float]:
    """Executes Azimuth scoring in parallel for guide RNA chunks.

    This function distributes guide RNA chunks across multiple processes to
    compute Azimuth scores efficiently, collecting and returning the scores in
    the original order.

    Args:
        guides_chunks (List[Tuple[int, List[Guide]]]): List of tuples containing
            the start index and guide chunk.
        size (int): Total number of guides.
        threads (int): Number of threads for parallel execution.

    Returns:
        List[float]: List of Azimuth scores for each guide.

    Raises:
        CrisprHawkAzimuthScoreError: If Azimuth score calculation fails for any chunk.
    """
    azimuth_scores = [np.nan] * size  # azimuth scores
    with ProcessPoolExecutor(max_workers=threads) as executor:
        future_to_chunk = {
            executor.submit(_azimuth, chunk): chunk[0] for chunk in guides_chunks
        }
        for future in as_completed(future_to_chunk):
            start_idx = future_to_chunk[future]
            try:
                chunk_start_idx, chunk_scores = future.result()
                for offset, score in enumerate(chunk_scores):
                    azimuth_scores[chunk_start_idx + offset] = score
            except Exception as e:
                exception_handler(
                    CrisprHawkAzimuthScoreError,
                    f"Azimuth score calculation failed for chunk at index {start_idx}",
                    os.EX_DATAERR,
                    debug,
                    e,
                )
    assert all(not np.isnan(s) for s in azimuth_scores)
    assert len(azimuth_scores) == size  # should match
    return azimuth_scores


def azimuth_score(
    guides: List[Guide], threads: int, verbosity: int, debug: bool
) -> List[Guide]:
    """Computes Azimuth scores for a list of guide RNAs.

    This function calculates the Azimuth efficiency score for each guide in
    parallel and updates the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        threads (int): Number of threads to use for parallel computation.
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
    guides_seqs = _extract_guide_sequences(guides)
    guides_seqs_chunks = calculate_chunks(
        guides_seqs, threads
    )  # split guides in chunks
    try:  # compute azimuth scores in parallel
        azimuth_scores = _execute_azimuth(
            guides_seqs_chunks, len(guides), threads, debug
        )
    except Exception as e:
        exception_handler(
            CrisprHawkAzimuthScoreError,
            "Azimuth score parallel execution failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    for i, score in enumerate(azimuth_scores):
        guides[i].azimuth_score = score  # assign score to each guide
    print_verbosity(
        f"Azimuth scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def _rs3(guides_chunks: Tuple[int, List[str]]) -> Tuple[int, List[float]]:
    """Calculates RS3 scores for a chunk of guide RNAs.

    This function computes RS3 efficiency scores for a given chunk of guide
    sequences and returns the starting index along with the list of scores.

    Args:
        guides_chunks (Tuple[int, List[str]]): A tuple containing the start index
            and a list of guide sequences.

    Returns:
        Tuple[int, List[float]]: The start index and the list of RS3 scores for
            the guides.
    """
    start_idx, guides = guides_chunks
    scores = rs3(guides)
    return start_idx, scores


def _execute_rs3(
    guides_chunks: List[Tuple[int, List[str]]], size: int, threads: int, debug: bool
) -> List[float]:
    """Executes RS3 scoring in parallel for guide RNA chunks.

    This function distributes guide RNA chunks across multiple processes to
    compute RS3 scores efficiently, collecting and returning the scores in the
    original order.

    Args:
        guides_chunks (List[Tuple[int, List[str]]]): List of tuples containing
            the start index and guide chunk.
        size (int): Total number of guides.
        threads (int): Number of threads for parallel execution.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[float]: List of RS3 scores for each guide.

    Raises:
        CrisprHawkAzimuthScoreError: If RS3 score calculation fails for any chunk.
    """
    rs3_scores = [np.nan] * size  # rs3 scores
    with ProcessPoolExecutor(max_workers=threads) as executor:
        future_to_chunk = {
            executor.submit(_rs3, chunk): chunk[0] for chunk in guides_chunks
        }
        for future in as_completed(future_to_chunk):
            start_idx = future_to_chunk[future]
            try:
                chunk_start_idx, chunk_scores = future.result()
                for offset, score in enumerate(chunk_scores):
                    rs3_scores[chunk_start_idx + offset] = score
            except Exception as e:
                exception_handler(
                    CrisprHawkAzimuthScoreError,
                    f"RS3 score calculation failed for chunk at index {start_idx}",
                    os.EX_DATAERR,
                    debug,
                    e,
                )
    assert all(not np.isnan(s) for s in rs3_scores)
    assert len(rs3_scores) == size  # should match
    return rs3_scores


def rs3_score(
    guides: List[Guide], threads: int, verbosity: int, debug: bool
) -> List[Guide]:
    """Computes RS3 scores for a list of guide RNAs.

    This function calculates the RS3 efficiency score for each guide in parallel
    and updates the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        threads (int): Number of threads to use for parallel computation.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with RS3 scores assigned.
    """
    if not guides:
        return guides
    print_verbosity("Computing RS3 score", verbosity, VERBOSITYLVL[3])
    start = time()  # rs3 score start time
    guides_seqs = _extract_guide_sequences(guides)
    guides_seqs_chunks = calculate_chunks(
        guides_seqs, threads
    )  # split guides in chunks
    try:  # compute rs3 scores in parallel
        rs3_scores = _execute_rs3(guides_seqs_chunks, len(guides), threads, debug)
    except Exception as e:
        exception_handler(
            CrisprHawkRs3ScoreError,
            "RS3 score calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
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
    guides: Dict[Region, List[Guide]], pam: PAM, args: CrisprHawkSearchInputArgs
) -> Dict[Region, List[Guide]]:
    """Scores CRISPR guides using efficiency and specificity metrics.

    This function computes Azimuth, RS3, DeepCpf1, Elevation, CFDon, and out-of-frame
    scores for each guide, updating the guide objects with the computed values.

    Args:
        guides: Dictionary mapping Region objects to lists of Guide objects.
        pam: PAM object specifying the protospacer adjacent motif and Cas system.
        args: CrisprHawkSearchInputArgs object containing scoring parameters.

    Returns:
        Dictionary mapping Region objects to lists of scored Guide objects.
    """
    # score guides using azimuth, rs3, deepcpf1, elevation, and out-of-frame scores
    print_verbosity("Scoring guides", args.verbosity, VERBOSITYLVL[1])
    start = time()  # scoring start time
    for region, guides_list in guides.items():
        if pam.cas_system in [SPCAS9, XCAS9]:  # cas9 system pam
            # score each guide with azimuth
            guides_list = azimuth_score(
                guides_list, args.threads, args.verbosity, args.debug
            )
            # score each guide with rs3
            guides_list = rs3_score(
                guides_list, args.threads, args.verbosity, args.debug
            )
            # score each guide with CFDon
            guides_list = cfdon_score(guides_list, args.verbosity, args.debug)
        if pam.cas_system == CPF1:  # cpf1 system pam
            guides_list = deepcpf1_score(guides_list, args.verbosity, args.debug)
        if args.compute_elevation and (
            args.guidelen + len(pam) == 23 and not args.right
        ):
            # elevation requires 23 bp long sequences, where last 3 bp are pam
            guides_list = elevationon_score(guides_list, args.verbosity, args.debug)
        # compute out-of-frame score (skipped)
        # guides_list = outofframe_score(
        #     guides_list, args.guidelen, args.right, args.verbosity, args.debug
        # )
        guides[region] = guides_list  # store scored guides
    print_verbosity(
        f"Scoring completed in {time() - start:.2f}s", args.verbosity, VERBOSITYLVL[2]
    )
    return guides
