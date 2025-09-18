"""Provides functions for calculating CFD (Cutting Frequency Determination) scores 
for CRISPR guide-target pairs.

This module includes utilities to load mismatch and PAM score models and to compute 
CFD scores using these models, supporting the evaluation of guide RNA specificity 
in genome editing applications.
"""

from crisprhawk.exception_handlers import exception_handler
from crisprhawk.crisprhawk_error import CrisprHawkCfdScoreError
from crisprhawk.utils import dna2rna, reverse_complement

from typing import Dict, Tuple

import pickle
import os

MMSCORES = "mismatch_score.pkl"
PAMSCORES = "pam_scores.pkl"


def load_mismatch_pam_scores(debug: bool) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Loads mismatch and PAM score models for CFD scoring.

    This function loads the mismatch and PAM score dictionaries required for CFD 
    score calculation from precomputed model files.

    Args:
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        Tuple[Dict[str, float], Dict[str, float]]: A tuple containing the mismatch 
            scores and PAM scores dictionaries.

    Raises:
        CrisprHawkCfdScoreError: If there is an error loading the model files.
    """
    modelspath = os.path.join(os.path.abspath(os.path.dirname(__file__)), "models")
    try:  # load mismatches and PAM scores (Doench et al., 2016)
        mmscores = pickle.load(open(os.path.join(modelspath, MMSCORES), mode="rb"))
        pamscores = pickle.load(open(os.path.join(modelspath, PAMSCORES), mode="rb"))
    except OSError as e:
        exception_handler(
            CrisprHawkCfdScoreError,
            "An error occurred while loading CFD model files",
            os.EX_NOINPUT,
            debug,
            e,
        )
    return mmscores, pamscores


def compute_cfd(
    wildtype: str,
    sg: str,
    pam: str,
    mmscores: Dict[str, float],
    pamscores: Dict[str, float],
    debug: bool,
) -> float:
    """Computes the CFD (Cutting Frequency Determination) score for a guide-target 
    pair.

    This function calculates the CFD score by comparing the wildtype and sgRNA
    sequences, applying mismatch and PAM penalties using provided score dictionaries.

    Args:
        wildtype (str): The wildtype DNA sequence.
        sg (str): The sgRNA sequence.
        pam (str): The PAM sequence.
        mmscores (Dict[str, float]): Dictionary of mismatch scores.
        pamscores (Dict[str, float]): Dictionary of PAM scores.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        float: The computed CFD score for the guide-target pair.
    """
    score = 1.0  # initialize cfd score
    wildtype, sg = dna2rna(wildtype), dna2rna(sg)  # convert to RNA sequences
    for i, ntsg in enumerate(sg):
        if i >= 20:  # handle off-targets bulges
            break
        if wildtype[i].upper() == ntsg.upper():
            score *= 1  # no mismatch, score unchanged
            continue
        elif wildtype[i].upper() == "-" or ntsg.upper() == "-":  # handle bulges
            score *= 1
            continue
        # build mismatch dictionary key
        key = (
            f"r{wildtype[i].upper()}:d{reverse_complement(ntsg.upper(), debug)},{i + 1}"
        )
        score *= mmscores[key]
    score *= pamscores[pam.upper()]  # multiply by PAM score
    return score
