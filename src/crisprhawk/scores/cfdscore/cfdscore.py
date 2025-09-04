""" """

from crisprhawk.exception_handlers import exception_handler
from crisprhawk.crisprhawk_error import CrisprHawkCfdScoreError
from crisprhawk.utils import dna2rna, reverse_complement

from typing import Dict, Tuple

import pickle
import os

MMSCORES = "mismatch_score.pkl"
PAMSCORES = "pam_scores.pkl"


def load_mismatch_pam_scores(debug: bool) -> Tuple[Dict[str, float], Dict[str, float]]:
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
