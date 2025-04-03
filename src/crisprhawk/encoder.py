"""
"""

from exception_handlers import exception_handler
from utils import IUPAC, VERBOSITYLVL, print_verbosity
from crisprhawk_error import CrisprHawkIupacTableError
from bitset import Bitset, SIZE

from typing import List
from time import time

import os


def _encoder(nt: str, position: int, debug: bool) -> Bitset:
    bitset = Bitset(SIZE, debug)  # 4 - bits encoder
    if nt == IUPAC[0]:  # A - 0001
        bitset.set(0)
    elif nt == IUPAC[1]:  # C - 0010
        bitset.set(1)
    elif nt == IUPAC[2]:  # G - 0100
        bitset.set(2)
    elif nt == IUPAC[3]:  # T - 1000
        bitset.set(3)
    elif nt == IUPAC[4]:  # R - 0101 G>A
        bitset.set_bits("0101")
    elif nt == IUPAC[5]:  # Y - 1010 C>T
        bitset.set_bits("1010")
    elif nt == IUPAC[6]:  # S - 0110 C>G
        bitset.set_bits("0110")
    elif nt == IUPAC[7]:  # W  - 1001
        bitset.set_bits("1001")
    elif nt == IUPAC[8]:  # K - 1100
        bitset.set_bits("1100")
    elif nt == IUPAC[9]:  # M - 0011
        bitset.set_bits("0011")
    elif nt == IUPAC[10]:  # B - 1110 --> not A ( T or G or C)
        bitset.set_bits("1110")
    elif nt == IUPAC[11]:  # D - 1101 --> not C
        bitset.set_bits("1101")
    elif nt == IUPAC[12]:  # H - 1011 --> not G
        bitset.set_bits("1011")
    elif nt == IUPAC[13]:  # V - 0111 --> not T
        bitset.set_bits("0111")
    elif nt == IUPAC[14]:  # N - 1111 --> any
        bitset.set_bits("1111")
    else:  # default case
        exception_handler(
            CrisprHawkIupacTableError,
            f"The nucleotide {nt} at {position} is not a IUPAC character",
            os.EX_DATAERR,
            debug
        )
    return bitset


def encode(sequence: str, verbosity: int, debug: bool) -> List[Bitset]:
    # encode sequence in bits for efficient matching
    print_verbosity(f"Encoding sequence {sequence} in bits", verbosity, VERBOSITYLVL[3])
    start = time()  # encoding start time
    bits = [_encoder(nt.upper(), i, debug) for i, nt in enumerate(sequence)]
    assert len(bits) == len(sequence)
    print_verbosity(f"Encoding completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3])
    return bits

