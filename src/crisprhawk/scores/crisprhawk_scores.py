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
from .elevation.cmds.predict import Predict
from ..guide import Guide
from ..utils import suppress_stdout, suppress_stderr

from typing import List, Dict, Union, Tuple

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


def cfdon(
    guide_ref: Union[None, Guide], guides: List[Guide], debug: bool
) -> List[float]:
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


def elevation(wildtypes: List[str], offtargets: List[str]) -> List[float]:
    p = Predict()  # initialize elevation predictor
    preds = p.execute(wildtypes, offtargets)  # retrieve elevation scores
    return [s for s in preds["linear-raw-stacker"]]


def elevationon(
    guide_groups: Dict[str, Tuple[Union[None, Guide], List[Guide]]],
) -> List[Guide]:
    # optimize input for elevation-on score calculation
    wildtype, offtarget, guides = [], [], []
    for _, (guide_ref, guides_g) in guide_groups.items():
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
        offtarget[i].set_elevationon_score(score)
    for g in guides:  # set NA for guides without reference alternative
        g.set_elevationon_score(np.nan)
    return guides + [g for g in offtarget]


# TODO: aggregate elevation
