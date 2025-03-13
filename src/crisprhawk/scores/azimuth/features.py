""" 
"""

from azimuth_error import CrisprHawkAzimuthError

from typing import Dict, Set, Tuple
from itertools import product

import pandas as pd
import numpy as np

NUCFEATS = ["all", "posind", "posdep"]  # nucleotide features
DNA = ["A", "T", "C", "G"]  # DNA alphabet

def derive_alphabet(order: int): 
    return ["".join(e) for e in product(DNA, repeat=order)]

def nucleotide_features(sequence: str, order: int, maxidx: int, feature: str) -> Tuple[pd.Series, pd.Series]:
    if feature not in NUCFEATS:
        raise CrisprHawkAzimuthError(f"Invalid nucleotide feature: {feature}")
    if maxidx is not None:
        maxidx = max(maxidx, len(sequence))
        sequence = sequence[:maxidx]
        assert len(sequence) == 30  # sequence length must be constant (30bp)
    # retrieve nucleotide alphabet for requested dependency order
    alphabet = derive_alphabet(order)
    # construct probability matrices for positional dependence and independence
    features_pos_dep = np.zeros(len(alphabet) * (len(sequence) - (order - 1)))
    features_pos_ind = np.zeros(np.power(len(DNA), order))
    # compute index lists
    index_dep = [f"{nt}_{p}" for p in range(len(sequence) - order + 1) for nt in alphabet]
    index_ind = [f"{nt}" for nt in alphabet]
    # populate positional dependence and independence matrices
    for pos in range(len(sequence) - order + 1):
        nt = sequence[pos : pos + order]
        idx_dep = alphabet.index(nt) + (pos * len(alphabet))
        assert index_dep[idx_dep] == f"{nt}_{pos}"  # check index consistency
        features_pos_dep[idx_dep] = 1.0
        idx_ind = alphabet.index(nt)
        assert index_ind[idx_ind] == nt  # check index consistency
        features_pos_ind[idx_ind] += 1.0
    # check matrix consistency
    if np.any(np.isnan(features_pos_dep)):
        raise CrisprHawkAzimuthError("Positional dependence matrix contains NaN values")
    if np.any(np.isnan(features_pos_ind)):
        raise CrisprHawkAzimuthError("Positional independence matrix contains NaN values")
    if feature ==  NUCFEATS[1]:  # positional indepenedence
        return None, pd.Series(features_pos_ind, index=index_ind)
    elif feature == NUCFEATS[1]:  # positional dependence
        return pd.Series(features_pos_dep, index=index_dep), None
    # both positional dependence and independence
    return pd.Series(features_pos_dep, index=index_dep), pd.Series(features_pos_ind, index=index_ind)


def apply_nucleotide_features(sequences: pd.Series, order: int, include_pos_ind: bool, maxidx: int) -> Tuple[pd.Series, pd.Series]:
    feat_pos_dep = sequences.apply(lambda x: nucleotide_features(x, order, maxidx, NUCFEATS[2])[0], axis=1)  # compute features accouting for positional dependence
    feat_pos_ind = None  # by default do not account for positional independent features
    if include_pos_ind:  # account for positional independence
        feat_pos_ind = sequences.apply(lambda x: nucleotide_features(x, order, maxidx, NUCFEATS[1])[1], axis=1)  # compute features accouting for positional independence    
    return feat_pos_dep, feat_pos_ind

def check_features_consistency(features: Dict[str, pd.Series]) -> None:
    if not features:
        raise CrisprHawkAzimuthError("Empty sequence features table")
    featnums = set([len(ft) for ft in features.values()])  # lengths of feature tables
    if len(featnums) > 1 or list(featnums)[0] < 1:
        raise CrisprHawkAzimuthError("Inconsistent sequence feature table lengths")
    

def retrieve_nucleotide_features(sequences: pd.Series, features: Dict[str, pd.Series], learn_options: Dict[str, float], maxorder: int, maxidx: int) -> Dict[str, pd.Series]:
    for order in range(1, maxorder + 1):
        # retrieve positional dependence and independence features
        feat_pos_dep, feat_pos_ind = apply_nucleotide_features(sequences, order, True, maxidx)
        features[f"nuc_pd_Order{order}"] = feat_pos_dep
        if learn_options["include_pi_nuc_feat"]:  # include positional independence features
            assert feat_pos_ind is not None
            features[f"nuc_pi_Order{order}"] = feat_pos_ind
        check_features_consistency(features)  # check feature table consistency
    return features    

def countgc(sequence: str) -> int:
    return len(sequence[4:24].replace("A", "").replace("T", ""))

def retrieve_gc_features(sequences: pd.Series, features: Dict[str, pd.Series]) -> Dict[str, pd.Series]:
    # compute GC content for each sequence
    gccounts = sequences.apply(lambda  x: countgc(x), axis=1)
    gccounts.name = "GC count"
    features["gc_count"] = gccounts
    # count sequences with GC content above and below 10 nts
    gccounts_above_10 = (gccounts > 10) * 1
    gccounts_above_10.name = "GC > 10"
    features["gc_above_10"] = gccounts_above_10
    gccounts_below_10 = (gccounts < 10) * 1
    gccounts_below_10.name = "GC < 10"
    features["gc_below_10"] = gccounts_below_10
    return features




def featurize(seqtable: pd.DataFrame, learn_options: Dict[str, float], predtable: pd.DataFrame,  geneposition: pd.DataFrame, check_pam: bool, check_len: bool):
    seqlengths = np.unique(seqtable["30mer"].apply(len).values)  # compute lengths of all sequences
    assert len(seqlengths) == 1  # sequences must not have different lengths
    features = {}  # initialize feature dictionary
    if learn_options["nuc_features"]:
        # use spectrum kernels (position-independent) and weighted degree kernels (position-dependent)
        features = retrieve_nucleotide_features(seqtable["30mer"], features, learn_options, learn_options["order"], 30)
    if learn_options["gc_features"]:
        # compute GC content in each sequence
        features = retrieve_gc_features(seqtable["30mer"], features)

