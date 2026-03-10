"""
This module provides functions for encoding nucleotide sequences into bitset
representations using IUPAC codes.

It enables efficient sequence matching by converting nucleotides and ambiguous
codes into Bitset objects for downstream analysis.
"""

from .exception_handlers import exception_handler
from .utils import VERBOSITYLVL, print_verbosity
from .crisprhawk_error import CrisprHawkIupacTableError

from typing import List, Dict
from time import time

import os

_IUPAC_BITS: Dict[str, int] = {
    "A": 0b0001,
    "C": 0b0010,
    "G": 0b0100,
    "T": 0b1000,
    "N": 0b1111,
    "R": 0b0101,
    "Y": 0b1010,
    "S": 0b0110,
    "W": 0b1001,
    "K": 0b1100,
    "M": 0b0011,
    "B": 0b1110,
    "D": 0b1101,
    "H": 0b1011,
    "V": 0b0111,
}


def _encoder(nt: str, position: int, debug: bool) -> int:
    if nt not in _IUPAC_BITS:
        exception_handler(
            CrisprHawkIupacTableError, 
            f"The nucleotide {nt} at position {position} is not a IUPAC character",
            os.EX_DATAERR,
            debug,
        )
    return _IUPAC_BITS[nt]

def encode(sequence: str, verbosity: int, debug: bool) -> List[int]:
    # encode sequence in bits for efficient matching
    print_verbosity(f"Encoding sequence {sequence} in bits", verbosity, VERBOSITYLVL[3])
    start = time()  # encoding start time
    bits = [_encoder(nt.upper(), i, debug) for i, nt in enumerate(sequence)]
    assert len(bits) == len(sequence)
    print_verbosity(
        f"Encoding completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return bits
