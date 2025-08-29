""" """

from .azimuth.model_comparison import predict
from rs3.seq import predict_seq
from .cfdscore.cfdscore import compute_cfd, load_mismatch_pam_scores
from .deepCpf1.seqdeepcpf1 import (
    preprocess,
    load_deepcpf1_weights,
    compute_deepcpf1,
    SeqDeepCpf1,
)
from .elevation.elevation.cmds.predict import Predict
from ..guide import Guide
from ..utils import suppress_stdout, suppress_stderr

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


def rs3(guides: List[str]) -> List[float]:
    # wrapper for ruleset3 predict function
    with suppress_stdout(), suppress_stderr():
        rs3scores = predict_seq(guides, sequence_tracr="Hsu2013")
    return list(rs3scores)


def cfdon(guide_ref: Guide, guides: List[Guide], debug: bool) -> List[float]:
    if not guide_ref:
        return [np.nan] * len(guides)
    mmscores, pamscores = load_mismatch_pam_scores(debug)  # load scoring models
    return [
        compute_cfd(guide_ref.guide, sg.guide, sg.pam[-2:], mmscores, pamscores, debug)
        for sg in guides
    ]  # compute cfd score for on-targets


def deepcpf1(guides: List[str]) -> List[float]:
    emb_matrix = preprocess(guides)  # initialize tensor
    model = SeqDeepCpf1()  # initialize seqdeepcpf1 model
    load_deepcpf1_weights(model)  # load models weights
    model.eval()
    return compute_deepcpf1(model, emb_matrix)


def elevationon(guide_ref: Guide, guides: List[Guide]) -> List[float]:
    if not guide_ref:
        return [np.nan] * len(guides)
    p = Predict()  # initialize elevation predictor
    wildtype = [guide_ref.guidepam] * len(guides)
    offtarget = [g.guidepam for g in guides]
    scores = p.execute(wildtype, offtarget)
    return [s["CFD"][0][0] for s in scores]  # compute elevation score for on-targets    
