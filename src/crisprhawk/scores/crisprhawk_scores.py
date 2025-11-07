"""
Scoring functions for CRISPR guide evaluation.

This module provides wrappers and utilities for computing various CRISPR guide
scores, including Azimuth, RuleSet3, CFD, DeepCpf1, Elevation, and out-of-frame
scores. It integrates multiple models and scoring algorithms to assess guide
efficiency and specificity.
"""

from .azimuth.model_comparison import predict
from rs3.seq import predict_seq
from .cfdscore.cfdscore import compute_cfd, load_mismatch_pam_scores
from .deepCpf1.seqdeepcpf1 import (
    preprocess,
    load_deepcpf1_weights,
    compute_deepcpf1,
    SeqDeepCpf1,
)
from .mhscore.microhomology import calculate_microhomology_score
from .elevation.cmds.predict import Predict
from ..guide import Guide
from ..utils import suppress_stdout, suppress_stderr

from typing import List, Dict, Union, Tuple

import numpy as np


def azimuth(guides: np.ndarray) -> List[float]:
    """Predict Azimuth scores for a set of guides.

    Returns a list of Azimuth scores for the provided guide sequences using the
    Azimuth model.

    Args:
        guides (np.ndarray): An array of guide sequences.

    Returns:
        List[float]: The predicted Azimuth scores for each guide.
    """
    # wrapper for azimuth predict function
    return list(predict(guides))


def rs3(guides: List[str]) -> List[float]:
    """Predict RuleSet3 scores for a list of guide sequences.

    Returns a list of RuleSet3 scores for the provided guide sequences using the
    RuleSet3 model.

    Args:
        guides (List[str]): A list of guide sequences.

    Returns:
        List[float]: The predicted RuleSet3 scores for each guide.
    """
    # wrapper for ruleset3 predict function
    with suppress_stdout(), suppress_stderr():
        rs3scores = predict_seq(guides, sequence_tracr="Hsu2013")
    return list(rs3scores)


def cfdon(
    guide_ref: Union[None, Guide], guides: List[Guide], debug: bool
) -> List[float]:
    """Compute CFD on-target scores for a list of guides.

    Returns a list of CFD (Cutting Frequency Determination) scores for the provided
    guides using the reference guide and loaded mismatch and PAM scoring models.

    Args:
        guide_ref (Union[None, Guide]): The reference guide or None.
        guides (List[Guide]): A list of Guide objects to score.
        debug (bool): Whether to enable debug mode for error handling.

    Returns:
        List[float]: The computed CFD scores for each guide.
    """
    if not guide_ref:
        return [np.nan] * len(guides)
    mmscores, pamscores = load_mismatch_pam_scores(debug)  # load scoring models
    return [
        compute_cfd(guide_ref.guide, sg.guide, sg.pam[-2:], mmscores, pamscores, debug)
        for sg in guides
    ]  # compute cfd score for on-targets


def deepcpf1(guides: List[str]) -> List[float]:
    """Predict DeepCpf1 scores for a list of guide sequences.

    Returns a list of DeepCpf1 scores for the provided guide sequences using the
    DeepCpf1 model.

    Args:
        guides (List[str]): A list of guide sequences.

    Returns:
        List[float]: The predicted DeepCpf1 scores for each guide.
    """
    emb_matrix = preprocess(guides)  # initialize tensor
    model = SeqDeepCpf1()  # initialize seqdeepcpf1 model
    load_deepcpf1_weights(model)  # load models weights
    model.eval()
    scores = compute_deepcpf1(model, emb_matrix)
    return [scores] if isinstance(scores, float) else scores


# TODO: aggregate elevation (currently using a hack)
def elevation(wildtypes: List[str], offtargets: List[str]) -> List[float]:
    """Compute Elevation scores for pairs of wildtype and offtarget sequences.

    Returns a list of Elevation scores for the provided wildtype and offtarget
    sequence pairs using the Elevation model.

    Args:
        wildtypes (List[str]): A list of wildtype guide sequences.
        offtargets (List[str]): A list of offtarget guide sequences.

    Returns:
        List[float]: The predicted Elevation scores for each wildtype-offtarget pair.
    """
    p = Predict()  # initialize elevation predictor
    preds = p.execute(wildtypes, offtargets)  # retrieve elevation scores
    return list(preds["linear-raw-stacker"])


def elevationon(
    guide_groups: Dict[str, Tuple[Union[None, Guide], List[Guide]]],
) -> List[Guide]:
    """Compute Elevation on-target scores for groups of guides.

    Returns a list of Guide objects with their elevationon_score attribute set,
    using the Elevation model for on-target scoring. Guides without a reference
    alternative are assigned a score of NaN.

    Args:
        guide_groups (Dict[str, Tuple[Union[None, Guide], List[Guide]]]): A dictionary
            mapping group keys to tuples of (reference guide or None, list of
            Guide objects).

    Returns:
        List[Guide]: The list of Guide objects with elevationon_score set.
    """
    # optimize input for elevation-on score calculation
    wildtype, offtarget, guides = [], [], []
    for guide_ref, guides_g in guide_groups.values():
        if guide_ref:
            wildtype.extend([guide_ref] * len(guides_g))
            offtarget.extend(guides_g)
        else:
            guides.extend(guides_g)
    # prepare input data for elevation score
    wildtype_ = [g.guidepam.upper() for g in wildtype]
    offtarget_ = [g.guidepam.upper() for g in offtarget]
    scores = elevation(wildtype_, offtarget_)
    for i, score in enumerate(scores):  # assign scores to guides
        offtarget[i].elevationon_score = score
    for g in guides:  # set NA for guides without reference alternative
        g.elevationon_score = np.nan
    return guides + list(offtarget)


def ooframe_score(guides: List[Guide], idx: int) -> List[int]:
    """Compute out-of-frame scores for a list of guides at a given index.

    Returns a list of out-of-frame (ooframe) scores for each guide, calculated
    using the microhomology score at the specified index.

    Args:
        guides (List[Guide]): A list of Guide objects.
        idx (int): The index around which to extract the sequence for scoring.

    Returns:
        List[int]: The out-of-frame scores for each guide.
    """
    guides_seqs = [g.sequence[idx - 30 : idx + 30] for g in guides]
    mhscores = [
        calculate_microhomology_score(gs.upper(), len(gs) // 2) for gs in guides_seqs
    ]
    return [mhscore.ooframe_score for mhscore in mhscores]
