""" """

from .azimuth.model_comparison import predict
from .cfdscore.cfdscore import compute_cfd, load_mismatch_pam_scores
from ..guide import Guide

from typing import List

import numpy as np

HITSCOREMIT = [
    0,
    0,
    0.014,
    0,
    0,
    0.395,
    0.317,
    0,
    0.389,
    0.079,
    0.445,
    0.508,
    0.613,
    0.851,
    0.732,
    0.828,
    0.615,
    0.804,
    0.685,
    0.583,
]


def azimuth(guides: np.ndarray) -> List[float]:
    # wrapper for azimuth predict function
    return list(predict(guides))

def cfdon(guide_ref: Guide, guides: List[Guide], debug: bool) -> List[float]:
    if not guide_ref:
        return [np.nan] * len(guides)
    mmscores, pamscores = load_mismatch_pam_scores(debug)  # load scoring models
    return [compute_cfd(guide_ref.guide, sg.guide, sg.pam[-2:], mmscores, pamscores, debug) for sg in guides]  # compute cfd score for on-targets

# def mit_score(guide: str, sequence: str) -> float:
#     assert len(guide) == len(sequence) == 20
#     mismatch_distances, mm, last_mm_pos = [], 0, 0  # initialize variables
#     score1 = 1.0  # first mit score
#     for i, nt in enumerate(guide):
#         if nt != sequence[i]:  # mismatching position
#             mm += 1
#             if last_mm_pos is not None:
#                 mismatch_distances.append(i - last_mm_pos)
#             score1 *= 1 - HITSCOREMIT[i]
#             last_mm_pos = i
#     # calculate score2 for distribution of mismatches
#     score2 = (
#         1.0 if mm < 2 else 1.0 / (((19 - np.mean(mismatch_distances)) / 19.0) * 4 + 1)
#     )
#     # calculate score3 for mismatch penalty
#     score3 = 1.0 if mm == 0 else 1.0 / (mm**2)
#     return score1 * score2 * score3


# def mit(guides: np.ndarray, sequences: np.ndarray) -> np.ndarray:
#     return np.array([mit_score(g, s) for g, s in zip(guides, sequences)])
