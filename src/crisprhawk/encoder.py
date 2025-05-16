""" """

from exception_handlers import exception_handler
from utils import IUPAC, VERBOSITYLVL, print_verbosity
from crisprhawk_error import CrisprHawkIupacTableError
from bitset import Bitset, SIZE

from typing import List
from time import time

import os


def _encoder(nt: str, position: int, debug: bool) -> Bitset:
    """Encodes a nucleotide character into a Bitset representation using IUPAC 
    codes.
    
    This function converts a nucleotide at a given position into its corresponding 
    bitset, handling ambiguous IUPAC codes.

    Args:
        nt: The nucleotide character to encode.
        position: The position of the nucleotide in the sequence.
        debug: Whether to enable debug mode for the Bitset.

    Returns:
        Bitset: The bitset representation of the nucleotide.

    Raises:
        CrisprHawkIupacTableError: If the nucleotide is not a valid IUPAC character.
    """
    bitset = Bitset(SIZE, debug)  # 4 - bits encoder
    if nt == IUPAC[0]:  # A - 0001
        bitset.set(0)
    elif nt == IUPAC[1]:  # C - 0010
        bitset.set(1)
    elif nt == IUPAC[2]:  # G - 0100
        bitset.set(2)
    elif nt == IUPAC[3]:  # T - 1000
        bitset.set(3)
    elif nt == IUPAC[4]:  # R - 0101 G or A
        bitset.set_bits("0101")
    elif nt == IUPAC[5]:  # Y - 1010 C or T
        bitset.set_bits("1010")
    elif nt == IUPAC[6]:  # S - 0110 C or G
        bitset.set_bits("0110")
    elif nt == IUPAC[7]:  # W  - 1001 A or T
        bitset.set_bits("1001")
    elif nt == IUPAC[8]:  # K - 1100 G or T
        bitset.set_bits("1100")
    elif nt == IUPAC[9]:  # M - 0011 A or C
        bitset.set_bits("0011")
    elif nt == IUPAC[10]:  # B - 1110 --> not A (T or G or C)
        bitset.set_bits("1110")
    elif nt == IUPAC[11]:  # D - 1101 --> not C (A or G or T)
        bitset.set_bits("1101")
    elif nt == IUPAC[12]:  # H - 1011 --> not G (A or C or T)
        bitset.set_bits("1011")
    elif nt == IUPAC[13]:  # V - 0111 --> not T (A or C or G)
        bitset.set_bits("0111")
    elif nt == IUPAC[14]:  # N - 1111 --> any
        bitset.set_bits("1111")
    else:  # default case
        exception_handler(
            CrisprHawkIupacTableError, # type: ignore
            f"The nucleotide {nt} at {position} is not a IUPAC character",
            os.EX_DATAERR,
            debug,
        )
    return bitset


def encode(sequence: str, verbosity: int, debug: bool) -> List[Bitset]:
    """Encodes a nucleotide sequence into a list of Bitset objects for efficient
    matching.

    This function processes each nucleotide in the input sequence and encodes it
    using IUPAC codes.

    Args:
        sequence: The nucleotide sequence to encode.
        verbosity: The verbosity level for logging.
        debug: Whether to enable debug mode for encoding.

    Returns:
        List[Bitset]: A list of Bitset objects representing the encoded sequence.
    """
    # encode sequence in bits for efficient matching
    print_verbosity(f"Encoding sequence {sequence} in bits", verbosity, VERBOSITYLVL[3])
    start = time()  # encoding start time
    bits = [_encoder(nt.upper(), i, debug) for i, nt in enumerate(sequence)]
    assert len(bits) == len(sequence)
    print_verbosity(
        f"Encoding completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return bits
