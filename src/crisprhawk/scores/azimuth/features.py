""" 
"""

from .azimuth_error import CrisprHawkAzimuthError

from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from typing import Dict, Tuple, Optional, List
from Bio.SeqUtils import MeltingTemp
from itertools import product

import pandas as pd
import numpy as np

NUCFEATS = ["all", "posind", "posdep"]  # nucleotide features
DNA = ["A", "T", "C", "G"]  # DNA alphabet

def derive_alphabet(order: int, alphabet: List[str]): 
    return ["".join(e) for e in product(alphabet, repeat=order)]

def nucleotide_features(sequence: str, order: int, maxidx: int, feature: str, alphabet: Optional[List[str]] = DNA, prefix: Optional[str] = "") -> Tuple[pd.Series, pd.Series]:
    if feature not in NUCFEATS:
        raise CrisprHawkAzimuthError(f"Invalid nucleotide feature: {feature}")
    if maxidx < len(sequence):
        maxidx = len(sequence)
    sequence = sequence[:maxidx]
    # retrieve nucleotide alphabet for requested dependency order
    alphabet = derive_alphabet(order, alphabet)
    # construct probability matrices for positional dependence and independence
    features_pos_dep = np.zeros(len(alphabet) * (len(sequence) - (order - 1)))
    features_pos_ind = np.zeros(np.power(len(DNA), order))
    # compute index lists
    index_dep = [f"{prefix}{nt}_{p}" for p in range(len(sequence) - order + 1) for nt in alphabet]
    index_ind = [f"{prefix}{nt}" for nt in alphabet]
    # populate positional dependence and independence matrices
    for pos in range(len(sequence) - order + 1):
        nt = sequence[pos : pos + order]
        idx_dep = alphabet.index(nt) + (pos * len(alphabet))
        assert index_dep[idx_dep] == f"{prefix}{nt}_{pos}"  # check index consistency
        features_pos_dep[idx_dep] = 1.0
        idx_ind = alphabet.index(nt)
        assert index_ind[idx_ind] == f"{prefix}{nt}"  # check index consistency
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
    feat_pos_dep = sequences.apply(lambda x: nucleotide_features(x, order, maxidx, NUCFEATS[2])[0])  # compute features accouting for positional dependence
    feat_pos_ind = None  # by default do not account for positional independent features
    if include_pos_ind:  # account for positional independence
        feat_pos_ind = sequences.apply(lambda x: nucleotide_features(x, order, maxidx, NUCFEATS[1])[1])  # compute features accouting for positional independence    
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
    gccounts = sequences.apply(lambda  x: countgc(x))
    gccounts.name = "GC count"
    features["gc_count"] = pd.DataFrame(gccounts)
    # count sequences with GC content above and below 10 nts
    gccounts_above_10 = (gccounts > 10) * 1
    gccounts_above_10.name = "GC > 10"
    features["gc_above_10"] = pd.DataFrame(gccounts_above_10)
    gccounts_below_10 = (gccounts < 10) * 1
    gccounts_below_10.name = "GC < 10"
    features["gc_below_10"] = pd.DataFrame(gccounts_below_10)
    return features

def retrieve_gene_position_features(geneposition: pd.DataFrame, features: Dict[str, pd.Series]) -> Dict[str, pd.Series]:
    for geneset in geneposition.columns.tolist():  # iterate over gene sets
        features[geneset] = pd.DataFrame(geneposition[geneset])
    # add peptide percentage feature data
    features["Percent Peptide <50%"] = pd.DataFrame(features["Percent Peptide"] < 50)
    features["Percent Peptide <50%"]["Percent Peptide <50%"] = features["Percent Peptide <50%"].pop("Percent Peptide")
    return features 

def retrieve_gene_effect_features(predictions: pd.DataFrame, features: Dict[str, pd.Series]) -> Dict[str, pd.Series]:
    genes = predictions["Target gene"]
    encoding = OneHotEncoder()  # encode target genes 
    encoding_labs = LabelEncoder()
    encoding_labs.fit(genes)  # fit genes to encoder
    genes_encoded = np.array(encoding.fit_transform(encoding_labs.transform(genes)[:, None]).todense())
    features["gene effect"] = pd.DataFrame(genes_encoded, columns=[f"gene_{i}" for i in range(genes_encoded.shape[1])], index=genes.index)
    return features


def retrieve_known_pairs_features(predictions: pd.DataFrame, features: Dict[str, pd.Series]) -> Dict[str, pd.Series]:
    features["known pairs"] = pd.DataFrame(predictions["test"])
    return features

def retrieve_nggx_interaction_features(seqtable: pd.Series, features: Dict[str, pd.Series], check_pam: bool) -> Dict[str, pd.Series]:
    nx_features = pd.DataFrame()
    for sequence in seqtable.values:
        if check_pam and sequence[25:27] != "GG":
            raise CrisprHawkAzimuthError(f"Expected NGG PAM, but found {sequence[24:27]}")
        nx = sequence[24] + sequence[27]
        nx_encoding = nucleotide_features(nx, 2, 2, "posdep", prefix="NGGX")[0]
        nx_features = pd.concat([nx_features, nx_encoding], axis=1)
    features["NGGX"] = nx_features.T
    return features

def retrieve_tm_features(seqtable: pd.Series, features: Dict[str, pd.Series], check_pam: bool) -> Dict[str, pd.Series]:
    segments = [(19, 24), (11, 19), (6, 11)]
    featarray = np.ones((seqtable.shape[0], 4))
    for i, sequence in enumerate(seqtable):  
        if check_pam and sequence[25:27] != "GG":
            raise CrisprHawkAzimuthError(f"Expected NGG PAM, but found {sequence[24:27]}")
        rna = False  # do not account for rna melting temp
        featarray[i, 0] = MeltingTemp.Tm_staluc(sequence, rna=rna)  # 30mer melting temp
        featarray[i, 1] = MeltingTemp.Tm_staluc(sequence[segments[0][0]:segments[0][1]], rna=rna)  # 5nts proximal to NGG PAM
        featarray[i, 2] = MeltingTemp.Tm_staluc(sequence[segments[1][0]:segments[1][1]], rna=rna)  # 8nts in middle
        featarray[i, 3] = MeltingTemp.Tm_staluc(sequence[segments[2][0]:segments[2][1]], rna=rna)  # 5nts at start
    columns = [f"Tm global_{rna}", f"5mer_end_{rna}", f"8mer_middle_{rna}", f"5mer_start_{rna}"]
    features["Tm"] = pd.DataFrame(featarray, index=seqtable.index, columns=columns)
    return features

def retrieve_sgrna_score_features(scoretable: pd.Series, features: Dict[str, pd.Series]) -> Dict[str, pd.Series]:
    features["sgRNA Score"] = pd.DataFrame(scoretable)
    return features

def retrieve_drug_features(predictions: pd.DataFrame, features: Dict[str, pd.Series]) -> Dict[str, pd.Series]:
    drugs = predictions.index.get_level_values("drug").tolist()  # drug list
    encoding = OneHotEncoder()
    encoding_labs = LabelEncoder()
    encoding_labs.fit(drugs)  # encode drug labels
    drugs_encoded = np.array(encoding.fit_transform(encoding_labs.transform(drugs)[:, None]).todense())
    features["drug"] = pd.DataFrame(drugs_encoded, columns=[f"drug_{i}" for i in range(drugs_encoded.shape[1])], index=drugs)
    return features

def retrieve_strand_features(strands: pd.Series, features: Dict[str, pd.Series]) -> Dict[str, pd.Series]:
    features["Strand effect"] = pd.DataFrame(strands == "sense") * 1
    return features


def featurize(seqtable: pd.DataFrame, learn_options: Dict[str, float], geneposition: pd.DataFrame, check_pam: bool, check_len: bool) -> Dict[str, pd.Series]:
    seqlengths = np.unique(seqtable["30mer"].apply(len).values)  # compute lengths of all sequences
    assert len(seqlengths) == 1  # sequences must not have different lengths
    features = {}  # initialize feature dictionary
    predictions = pd.DataFrame()  # initialize predictions table (Y in original azimuth)
    if learn_options["nuc_features"]:
        # use spectrum kernels (position-independent) and weighted degree kernels (position-dependent)
        features = retrieve_nucleotide_features(seqtable["30mer"], features, learn_options, learn_options["order"], 30)
    if learn_options["gc_features"]:  # compute GC content in each sequence
        features = retrieve_gc_features(seqtable["30mer"], features)
    if learn_options["include_gene_position"]:  # add gene position features
        features = retrieve_gene_position_features(geneposition, features)
    if learn_options["include_gene_effect"]:  # add gene effect
        features = retrieve_gene_effect_features(predictions, features)
    if learn_options["include_known_pairs"]:  # add known gene-gene pairs
        features = retrieve_known_pairs_features(predictions, features)
    if learn_options["include_NGGX_interaction"]:  # add NGGX interaction features
        features = retrieve_nggx_interaction_features(seqtable["30mer"], features, check_pam)
    if learn_options["include_Tm"]:  # add melting temperature features
        features = retrieve_tm_features(seqtable["30mer"], features, check_pam)
    if learn_options["include_sgRNAscore"]:  # include sgrna scores
        features = retrieve_sgrna_score_features(seqtable["sgRNA Score"], features)
    if learn_options["include_drug"]:  # include drug interaction features
        features = retrieve_drug_features(predictions, features)
    if learn_options["include_strand"]:  # include strand information
        features = retrieve_strand_features(seqtable["Strand"], features)
    check_features_consistency(features)
    return features

def concatenate_features(features: Dict[str, pd.Series]) -> np.ndarray:
    if not features:
        raise CrisprHawkAzimuthError("Empty features dictionary")
    keys = list(features.keys())  # list of features
    featnum = features[keys[0]].shape[0]
    for k, v in features.items():  # check features size consistency
        if v.shape[0] != featnum:
            raise CrisprHawkAzimuthError(f"Mismatching feature dimensions (expected {featnum}, got {v.shape[0]} at column {k})")
    inputs = np.zeros((featnum, 0))  # model test matrix
    for k in keys:  # construct matrix
        inputs = np.hstack((inputs, features[k].values))
    return inputs
    
    

