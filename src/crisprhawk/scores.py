"""
"""

from scores.azimuth.azimuth import load_azimuth_model, override_learn_options, construct_sequences_table
from scores.azimuth.features import featurize, concatenate_features

from typing import List, Optional, Dict

import numpy as np

import pickle
import os

AZIMUTHDIR = "scores/azimuth"
AZIMUTHMODELS = ["V3_model_full.pickle", "V3_model_nopos.pickle"]

def azimuth(
    sequences: List[str], 
    verbosity: int,
    aa_cut_positions: Optional[np.ndarray] = None, 
    peptide_pct: Optional[np.ndarray] = None, 
    model_file: Optional[str] = None, 
    learn_options_dict: Optional[Dict[str, float]] = None, 
    check_pam: Optional[bool] = False, 
    check_len: Optional[bool] = False
) -> np.ndarray:
    # load azimuth model
    model, learn_options = load_azimuth_model(model_file, peptide_pct, aa_cut_positions, verbosity)
    learn_options["V"] = 2  # set V=2 fro azimuth learn parameters
    if learn_options_dict is not None:  # update learn options if provided
        learn_options = override_learn_options(learn_options_dict, learn_options)
    # construct sequence table
    seqtable, geneposition = construct_sequences_table(sequences, peptide_pct, aa_cut_positions)
    # compute sequence features
    features = featurize(seqtable, learn_options, geneposition, check_pam, check_len)
    inputs = concatenate_features(features)  # concatenate features 
    predictions = model.predict(inputs)  # call to scikit-learn, returns a vector of predicted values    
    if all(p in [0, 1] for p in np.unique(predictions)):
        raise Exception
    return predictions


if __name__ == "__main__":
    sequences = np.array(['ACAGCTGATCTCCAGATATGACCATGGGTT', 'CAGCTGATCTCCAGATATGACCATGGGTTT', 'CCAGAAGTTTGAGCCACAAACCCATGGTCA'])
    amino_acid_cut_positions = np.array([2, 2, 4])
    percent_peptides = np.array([0.18, 0.18, 0.35])
    predictions = azimuth(sequences, 3, amino_acid_cut_positions, percent_peptides)
    for i, prediction in enumerate(predictions):
        print(sequences[i], prediction)


    