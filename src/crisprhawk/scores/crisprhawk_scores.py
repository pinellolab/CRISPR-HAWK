"""
"""

from .azimuth.model_comparison import predict

import numpy as np


def azimuth(guides: np.ndarray) -> np.ndarray:
    # wrapper for azimuth predict function
    return predict(guides)


