"""
This module provides functions for scoring CRISPR guide RNAs using various efficiency
and specificity metrics.

It includes methods to compute Azimuth, RS3, DeepCpf1, Elevation, CFDon, and out-of-frame
scores for guide RNAs.

The module is designed to annotate and enrich guide RNA data with these scores for
downstream genome editing analysis.
"""

from .config_utils import ScoringEnvs, CrisprOnConfig
from .crisprhawk_error import (
    CrisprHawkAzimuthScoreError,
    CrisprHawkRs3ScoreError,
    CrisprHawkCfdScoreError,
    CrisprHawkDeepCpf1ScoreError,
    CrisprHawkPlmCrisprScoreError,
    CrisprHawkCRISPRonScoreError,
    CrisprHawksgdesignerScoreError,
)
from .crisprhawk_argparse import CrisprHawkSearchInputArgs
from .exception_handlers import exception_handler
from .scores import (
    azimuth,
    rs3,
    cfdon,
    deepcpf1,
    elevationon,
    plmcrispr,
    crispron,
    sgdesigner,
)
from .utils import calculate_chunks, flatten_list, print_verbosity, VERBOSITYLVL
from .region import Region
from .guide import Guide, GUIDESEQPAD
from .pam import PAM, SPCAS9, XCAS9, CPF1

from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple, Union
from collections import defaultdict
from time import time

import numpy as np
import subprocess
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
    # each guide must have 4 nts upstream the guide sequence, and 3 nts downstream
    # the pam (Azimuth, RS3, and DeepCpf1)
    return [
        guide.sequence[(GUIDESEQPAD - 4) : (-GUIDESEQPAD + 3)].upper()
        for guide in guides
    ]


def _extract_guide_sequences_sgdesigner(guides: List[Guide]) -> List[str]:
    """Extracts guide RNA sequences with required flanking nucleotides.

    This function returns a list of guide sequences, each including
    3 nucleotides downstream of the PAM, formatted in uppercase.

    Args:
        guides (List[Guide]): List of Guide objects to extract sequences from.

    Returns:
        List[str]: List of formatted guide sequences with flanking nucleotides.
    """
    return [
        guide.sequence[(GUIDESEQPAD) : (-GUIDESEQPAD + 3)].upper() for guide in guides
    ]


def _azimuth(guides_chunk: Tuple[int, List[str]]) -> Tuple[int, List[float]]:
    """Calculates Azimuth scores for a chunk of guide RNAs.

    This function computes Azimuth efficiency scores for a given chunk of guides
    and returns the starting index along with the list of scores.

    Args:
        guides_chunk (Tuple[int, List[Guide]]): A tuple containing the start
            index and a list of Guide objects.

    Returns:
        Tuple[int, List[float]]: The start index and the list of Azimuth scores
            for the guides.
    """
    start_idx, guides = guides_chunk
    # create guides np.ndarray required by azimuth
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
    if not guides:
        return guides  # no guide?
    print_verbosity("Computing Azimuth score", verbosity, VERBOSITYLVL[3])
    start = time()  # azimuth score start time
    guides_seqs = _extract_guide_sequences(guides)
    # split guides in chunks
    guides_seqs_chunks = calculate_chunks(guides_seqs, threads)
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


def _rs3(guides_chunk: Tuple[int, List[str]]) -> Tuple[int, List[float]]:
    """Calculates RS3 scores for a chunk of guide RNAs.

    This function computes RS3 efficiency scores for a given chunk of guide
    sequences and returns the starting index along with the list of scores.

    Args:
        guides_chunk (Tuple[int, List[str]]): A tuple containing the start index
            and a list of guide sequences.

    Returns:
        Tuple[int, List[float]]: The start index and the list of RS3 scores for
            the guides.
    """
    start_idx, guides = guides_chunk
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
    # split guides in chunks
    guides_seqs_chunks = calculate_chunks(guides_seqs, threads)
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


def _deepcpf1(guides_chunk: Tuple[int, List[str]]) -> Tuple[int, List[float]]:
    """Calculates DeepCpf1 scores for a chunk of guide RNAs.

    This function computes DeepCpf1 efficiency scores for a given chunk of guide
    sequences and returns the starting index along with the list of scores.

    Args:
        guides_chunk (Tuple[int, List[str]]): A tuple containing the start index
            and a list of guide sequences.

    Returns:
        Tuple[int, List[float]]: The start index and the list of DeepCpf1 scores
            for the guides.
    """
    start_idx, guides = guides_chunk
    scores = deepcpf1(guides)
    return start_idx, scores


def _execute_deepcpf1(
    guides_chunks: List[Tuple[int, List[str]]], size: int, threads: int, debug: bool
) -> List[float]:
    """Executes DeepCpf1 scoring in parallel for guide RNA chunks.

    This function distributes guide RNA chunks across multiple processes to
    compute DeepCpf1 scores efficiently, collecting and returning the scores in
    the original order.

    Args:
        guides_chunks (List[Tuple[int, List[str]]]): List of tuples containing
            the start index and guide chunk.
        size (int): Total number of guides.
        threads (int): Number of threads for parallel execution.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[float]: List of DeepCpf1 scores for each guide.

    Raises:
        CrisprHawkAzimuthScoreError: If DeepCpf1 score calculation fails for
            any chunk.
    """
    deepcpf1_scores = [np.nan] * size  # deepcpf1 scores
    with ProcessPoolExecutor(max_workers=threads) as executor:
        future_to_chunk = {
            executor.submit(_deepcpf1, chunk): chunk[0] for chunk in guides_chunks
        }
        for future in as_completed(future_to_chunk):
            start_idx = future_to_chunk[future]
            try:
                chunk_start_idx, chunk_scores = future.result()
                for offset, score in enumerate(chunk_scores):
                    deepcpf1_scores[chunk_start_idx + offset] = score
            except Exception as e:
                exception_handler(
                    CrisprHawkAzimuthScoreError,
                    f"DeepCpf1 score calculation failed for chunk at index {start_idx}",
                    os.EX_DATAERR,
                    debug,
                    e,
                )
    assert all(not np.isnan(s) for s in deepcpf1_scores)
    assert len(deepcpf1_scores) == size  # should match
    return deepcpf1_scores


def deepcpf1_score(
    guides: List[Guide], threads: int, verbosity: int, debug: bool
) -> List[Guide]:
    """Computes DeepCpf1 scores for a list of guide RNAs.

    This function calculates the DeepCpf1 efficiency score for each guide in
    parallel and updates the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        threads (int): Number of threads to use for parallel computation.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with DeepCpf1 scores assigned.
    """
    if not guides:
        return guides
    print_verbosity("Computing DeepCpf1 score", verbosity, VERBOSITYLVL[3])
    start = time()  # deepcpf1 score start time
    guides_seqs = _extract_guide_sequences(guides)
    # split guides in chunks
    guides_seqs_chunks = calculate_chunks(guides_seqs, threads)
    try:  # compute deepcpf1 scores
        deepcpf1_scores = _execute_deepcpf1(
            guides_seqs_chunks, len(guides), threads, debug
        )
    except Exception as e:
        exception_handler(
            CrisprHawkDeepCpf1ScoreError,
            "DeepCpf1 score calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
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


def _plmcrispr(
    guides_chunk: Tuple[int, List[str]], cas_system: int
) -> Tuple[int, List[float]]:
    """Calculates PLM-CRISPR scores for a chunk of guide RNAs.

    This function computes PLM-CRISPR efficiency scores for a given chunk of
    guide sequences and returns the starting index along with the list of scores.

    Args:
        guides_chunk (Tuple[int, List[str]]): A tuple containing the start index
            and a list of guide sequences.
        cas_system (int): Identifier of the Cas system used for scoring.

    Returns:
        Tuple[int, List[float]]: The start index and the list of PLM-CRISPR
            scores for the guides.
    """
    start_idx, guides = guides_chunk
    scores = plmcrispr(guides, cas_system)
    return start_idx, scores


def _execute_plmcrispr(
    guide_chunks: List[Tuple[int, List[str]]],
    cas_system: int,
    size: int,
    threads: int,
    debug: bool,
) -> List[float]:
    """Executes PLM-CRISPR scoring in parallel for guide RNA chunks.

    This function distributes guide RNA chunks across multiple processes to
    compute PLM-CRISPR efficiency scores, collecting and returning the scores
    in the original guide order.

    Args:
        guide_chunks (List[Tuple[int, List[str]]]): List of tuples containing
            the start index and guide sequence chunk.
        cas_system (int): Identifier of the Cas system used for scoring.
        size (int): Total number of guides.
        threads (int): Number of threads for parallel execution.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[float]: List of PLM-CRISPR scores for each guide.

    Raises:
        CrisprHawkPlmCrisprScoreError: If PLM-CRISPR score calculation fails
            for any chunk.
    """
    plmcrispr_scores = [np.nan] * size  # plm-crispr scores
    with ProcessPoolExecutor(max_workers=1) as executor:
        future_to_chunk = {
            executor.submit(_plmcrispr, chunk, cas_system): chunk[0]
            for chunk in guide_chunks
        }
        for future in as_completed(future_to_chunk):
            start_idx = future_to_chunk[future]
            try:
                chunk_start_idx, chunk_scores = future.result()
                for offset, score in enumerate(chunk_scores):
                    plmcrispr_scores[chunk_start_idx + offset] = score
            except Exception as e:
                exception_handler(
                    CrisprHawkPlmCrisprScoreError,
                    f"PLM-CRISPR score calculation failed for chunk at index {start_idx}",
                    os.EX_DATAERR,
                    debug,
                    e,
                )
    assert all(not np.isnan(s) for s in plmcrispr_scores)
    assert len(plmcrispr_scores) == size  # should match
    return plmcrispr_scores


def _plmcrispr_score(
    guides: List[Guide], cas_system: int, threads: int, verbosity: int, debug: bool
) -> List[Guide]:
    """Computes PLM-CRISPR scores for a list of guide RNAs.

    This function calculates the PLM-CRISPR efficiency score for each guide in
    parallel and updates the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        cas_system (int): Identifier of the Cas system used for scoring.
        threads (int): Number of threads to use for parallel computation.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with PLM-CRISPR scores assigned.
    """
    if not guides:
        return guides  # no guide?
    print_verbosity("Computing PLM-CRISPR score", verbosity, VERBOSITYLVL[3])
    start = time()  # plm-crispr start time
    # retrieve spacer + pam sequence for each input guide
    guides_seqs = [g.guidepam for g in guides]
    try:  # compute plm-crispr scores
        plmcrispr_scores = plmcrispr(guides_seqs, cas_system)
    except Exception as e:
        exception_handler(
            CrisprHawkPlmCrisprScoreError,
            "PLM-CRISPR score calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    for i, score in enumerate(plmcrispr_scores):
        guides[i].plmcrispr_score = score  # assign score to each guide
    print_verbosity(
        f"PLM-CRISPR scores computed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return guides


def _crispron_score(
    guides: List[Guide], config: CrisprOnConfig, verbosity: int, debug: bool
):
    if not guides:
        return guides
    print_verbosity("Computing CRISPRon score", verbosity, VERBOSITYLVL[3])
    start = time()  # crispron score start time
    # crispron requires 4 nt upstream + 20 nt guide + 3 nt PAM + 3 nt downstream
    guides_seqs = _extract_guide_sequences(guides)
    try:  # single threads-only to avoid out of memory execptions
        crispron_scores = crispron(guides_seqs, config.conda, config.env_name)
    except Exception as e:
        exception_handler(
            CrisprHawkCRISPRonScoreError,
            "CRISPRon score calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    for i, score in enumerate(crispron_scores):
        guides[i].crispron_score = score  # assign score to each guide
    print_verbosity(
        f"CRISPRon scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def _sgdesigner(
    guides_chunk: Tuple[int, List[str]],
    env_name: str,
    tmp_parent: str,
) -> Tuple[int, List[float]]:
    """Compute sgDesigner scores for a chunk of 26mers."""
    start_idx, guides = guides_chunk
    scores = sgdesigner(guides, env_name, tmp_parent)
    return start_idx, scores


def _execute_sgdesigner(
    guide_chunks: List[Tuple[int, List[str]]],
    env_name: str,
    tmp_parent: str,
    size: int,
    threads: int,
    debug: bool,
) -> List[float]:
    sgdesigner_scores = [np.nan] * size
    with ProcessPoolExecutor(max_workers=threads) as executor:
        future_to_chunk = {
            executor.submit(_sgdesigner, chunk, env_name, tmp_parent): chunk[0]
            for chunk in guide_chunks
        }
        for future in as_completed(future_to_chunk):
            start_idx = future_to_chunk[future]
            try:
                chunk_start_idx, chunk_scores = future.result()
                for offset, score in enumerate(chunk_scores):
                    sgdesigner_scores[chunk_start_idx + offset] = score
            except Exception as e:
                exception_handler(
                    CrisprHawksgdesignerScoreError,
                    f"sgDesigner score calculation failed for chunk at index {start_idx}",
                    os.EX_DATAERR,
                    debug,
                    e,
                )
    assert all(not np.isnan(s) for s in sgdesigner_scores)
    assert len(sgdesigner_scores) == size
    return sgdesigner_scores


def sgdesigner_score(
    guides: List[Guide],
    env_name: str,
    tmp_parent: str,
    threads: int,
    verbosity: int,
    debug: bool,
) -> List[Guide]:
    """Compute sgDesigner scores for standard 30mer Cas9 guides."""
    if not guides:
        return guides
    print_verbosity("Computing sgDesigner score", verbosity, VERBOSITYLVL[3])
    start = time()

    guides_seqs = _extract_guide_sequences_sgdesigner(guides)

    # sgDesigner needs 26mers: 20 nt guide + 3 nt PAM + 3 nt downstream
    if not all(len(seq) == 26 for seq in guides_seqs):
        exception_handler(
            CrisprHawksgdesignerScoreError,
            "sgDesigner requires 30mer sequences; current guide/PAM configuration does not produce 26mers",
            os.EX_DATAERR,
            debug,
        )
    guides_seqs_chunks = calculate_chunks(guides_seqs, threads)

    try:
        sgdesigner_scores = _execute_sgdesigner(
            guides_seqs_chunks,
            env_name,
            tmp_parent,
            len(guides),
            threads,
            debug,
        )
    except Exception as e:
        exception_handler(
            CrisprHawksgdesignerScoreError,
            "sgDesigner score parallel execution failed",
            os.EX_DATAERR,
            debug,
            e,
        )

    for i, score in enumerate(sgdesigner_scores):
        guides[i].sgdesigner_score = score

    print_verbosity(
        f"sgDesigner scores computed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return guides


def _scoring_guides_cas9(
    guides_list: List[Guide],
    cas_system: int,
    scoring_envs: ScoringEnvs,
    threads: int,
    verbosity: int,
    debug: bool,
) -> List[Guide]:
    # score each guide with azimuth score
    guides_list = azimuth_score(guides_list, threads, verbosity, debug)
    # score each guide with rs3 score
    guides_list = rs3_score(guides_list, threads, verbosity, debug)
    # score each guide with PLM-CRISPR score
    guides_list = _plmcrispr_score(guides_list, cas_system, threads, verbosity, debug)
    # score each guide with CFDon score
    guides_list = cfdon_score(guides_list, verbosity, debug)
    # score each guide with CRISPRon score
    if scoring_envs.crispron_env:
        guides_list = _crispron_score(
            guides_list, scoring_envs.crispron_env, verbosity, debug
        )
    # score each guide with sgDesigner score
    # if scoring_envs.sgdesigner_env:
    #     pass
    return guides_list


def _scoring_guides_cpf1(
    guides_list: List[Guide], threads: int, verbosity: int, debug: bool
):
    # score each guide with deepCpf1 score
    return deepcpf1_score(guides_list, threads, verbosity, debug)


def scoring_guides(
    guides: Dict[Region, List[Guide]],
    pam: PAM,
    scoring_envs: ScoringEnvs,
    args: CrisprHawkSearchInputArgs,
) -> Dict[Region, List[Guide]]:
    # score guides using azimuth, rs3, deepcpf1, elevation, and out-of-frame scores
    print_verbosity("Scoring guides", args.verbosity, VERBOSITYLVL[1])
    start = time()  # scoring start time
    for region, guides_list in guides.items():
        if pam.cas_system in [SPCAS9, XCAS9]:  # cas9 system pam
            guides_list = _scoring_guides_cas9(
                guides_list,
                pam.cas_system,
                scoring_envs,
                args.threads,
                args.verbosity,
                args.debug,
            )
            # score each guide with sgDesigner if environment exists
            if scoring_envs.sgdesigner_env:
                guides_list = sgdesigner_score(
                    guides_list,
                    scoring_envs.sgdesigner_env.env_name,
                    scoring_envs.sgdesigner_env.outdir,
                    args.threads,
                    args.verbosity,
                    args.debug,
                )
        elif pam.cas_system == CPF1:  # cpf1 system pam
            guides_list = _scoring_guides_cpf1(
                guides_list, args.threads, args.verbosity, args.debug
            )
        if args.compute_elevation and (
            args.guidelen + len(pam) == 23 and not args.right
        ):
            # elevation requires 23 bp long sequences, where last 3 bp are pam
            guides_list = elevationon_score(guides_list, args.verbosity, args.debug)
        guides[region] = guides_list  # store scored guides
    print_verbosity(
        f"Scoring completed in {time() - start:.2f}s", args.verbosity, VERBOSITYLVL[2]
    )
    return guides
