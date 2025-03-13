"""
"""

from scores.azimuth.azimuth import load_azimuth_model, override_learn_options

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
    print(type(model))
    print(learn_options)
    exit()
    learn_options_dict["V"] = 2  # set V=2 fro azimuth learn parameters
    if learn_options_dict is not None:  # update learn options if provided
        learn_options = override_learn_options(learn_options_dict, learn_options)

if __name__ == "__main__":
    sequences = np.array(['ACAGCTGATCTCCAGATATGACCATGGGTT', 'CAGCTGATCTCCAGATATGACCATGGGTTT', 'CCAGAAGTTTGAGCCACAAACCCATGGTCA'])
    amino_acid_cut_positions = np.array([2, 2, 4])
    percent_peptides = np.array([0.18, 0.18, 0.35])
    azimuth(sequences, 3, amino_acid_cut_positions, percent_peptides)


    