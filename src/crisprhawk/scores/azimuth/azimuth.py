""" 
"""

from .azimuth_error import CrisprHawkAzimuthError
from .azimuth_utils import warning

from sklearn import GradientBoostingRegressor
from typing import List, Optional, Union, Dict, Tuple, Any

# from sklearn.ensemble. import LeastSquaresError

import pandas as pd
import numpy as np

import pickle
import sys
import os

sys.modules["sklearn.ensemble.gradient_boosting"] = sys.modules["sklearn.ensemble"]
# sys.modules["sklearn.ensemble.gradient_boosting.LeastSquaresError"] = LeastSquaresError

AZIMUTHMODELS = ["V3_model_full.pickle", "V3_model_nopos.pickle"]
SEQTABLECOLS = [u"30mer", u"Strand", u"Percent Peptide", u"Amino Acid Cut Position"]

def load_azimuth_model(model_file: Union[str, None], peptide_pct: Union[np.ndarray, None], aa_cut_positions: Union[np.ndarray, None], verbosity: int) -> Tuple[Any, Dict[str, Any]]:
    if model_file is None:  # no input model file provided (most common case)
        # load default azimuth models
        modelsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "models"))
        assert os.path.isdir(modelsdir)
        if np.any(peptide_pct == -1) or (peptide_pct is None and aa_cut_positions is None):
            warning("No azimuth model file specified, using V3 models (no position)", verbosity)
            model_fname = AZIMUTHMODELS[1]  # use V3 model without position
        else:
            warning("No azimuth model file specified, using V3 model (full)", verbosity)
            model_fname = AZIMUTHMODELS[0]  # use V3 model full
        model_file = os.path.join(modelsdir, model_fname)
        assert os.path.isfile(model_file)
    try:
        with open(model_file, mode="rb") as f:
            model, learn_options = pickle.load(f)  # load model and learn parameters
    except (pickle.UnpicklingError, EOFError) as e:  # always fully trace this error
        raise CrisprHawkAzimuthError(f"Error loading azimuth model file: {model_file}") from e
    return model, learn_options

def override_learn_options(learn_options_new: Dict[str, float], learn_options: Dict[str, Any]) -> Dict[str, float]:
    # update all keys in learn options with new values
    for key, value in learn_options_new.items():
        learn_options[key] = value
    return learn_options


def construct_sequences_table(sequences: List[str], peptide_pct: Union[np.ndarray, None], aa_cut_positions: Union[np.ndarray, None]) -> pd.DataFrame:
    seqtable = pd.DataFrame(columns=SEQTABLECOLS[:2], data=zip(sequences, ["NA" for _ in range(len(sequences))]))
    if np.all(peptide_pct != -1) and (peptide_pct is not None and aa_cut_positions is not None):
        geneposition = pd.DataFrame(columns=SEQTABLECOLS[2:], data=zip(peptide_pct, aa_cut_positions))
    else:
        geneposition = pd.DataFrame(columns=SEQTABLECOLS[2:], data=zip(np.ones(len(sequences)) * -1, np.ones(len(sequences)) * -1))
    


def predict(
    sequences: List[str], 
    aa_cut_positions: Optional[np.ndarray] = None, 
    peptide_pct: Optional[np.ndarray] = None, 
    model_file: Optional[str] = None, 
    learn_options_dict: Optional[Dict[str, float]] = None, 
    check_pam: Optional[bool] = False, 
    check_len: Optional[bool] = False
) -> np.ndarray:
    if model_file is None:
        # load default azimuth models
        azimuth_models_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "models"))
        assert os.path.isdir(azimuth_models_dir)



