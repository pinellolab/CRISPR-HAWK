""" 
"""

from exception_handlers import exception_handler

from typing import List
from itertools import permutations
from colorama import Fore

import sys
import os

# define static variables shared across software modules
TOOLNAME = "CRISPR-HAWK"  # tool name
COMMAND = "crisprhawk"  # command line call
# define verbosity levels
VERBOSITYLVL = [0, 1, 2, 3]
# complete iupac alphabet
IUPAC = ["A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]
# reverse complement dictionary
RC = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "R": "Y",
    "Y": "R",
    "M": "K",
    "K": "M",
    "H": "D",
    "D": "H",
    "B": "V",
    "V": "B",
    "N": "N",
    "S": "S",
    "W": "W",
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "r": "y",
    "y": "r",
    "m": "k",
    "k": "m",
    "h": "d",
    "d": "h",
    "b": "v",
    "v": "b",
    "n": "n",
    "s": "s",
    "w": "w",
}
# dictionary to encode nucleotides combinations as iupac characters
IUPACTABLE = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "M": "AC",
    "K": "GT",
    "S": "CG",
    "W": "AT",
    "H": "ACT",
    "B": "CGT",
    "V": "ACG",
    "D": "AGT",
    "N": "ACGT",
}
# dictionary to encode nucleotide strings as iupac characters
IUPAC_ENCODER = {
    perm: k
    for k, v in IUPACTABLE.items()
    for perm in {"".join(p) for p in permutations(v)}
}
STRAND = [0, 1]  # strands directions: 0 -> 5'-3'; 1 -> 3'-5'
# report prefix name
GUIDESREPORTPREFIX = "crisprhawk_guides"


# define utils functions
def reverse_complement(sequence: str, debug: bool) -> str:
    try:
        return "".join([RC[nt] for nt in sequence[::-1]])
    except KeyError as e:
        exception_handler(
            ValueError,
            f"Failed reverse complement on {sequence}",
            os.EX_DATAERR,
            debug,
            e,
        )


def warning(message: str, verbosity: int) -> None:
    if verbosity >= VERBOSITYLVL[1]:
        sys.stderr.write(f"{Fore.YELLOW}WARNING: {message}.{Fore.RESET}\n")
    return


def print_verbosity(message: str, verbosity: int, verbosity_threshold: int) -> None:
    if verbosity >= verbosity_threshold:
        sys.stdout.write(f"{message}\n")
    return


def adjust_guide_position(pos: int, pos_rel: List[int]) -> int:
    offset = -1
    for i, p in enumerate(pos_rel):
        if p == pos:
            offset = i
            break
    return offset

def round_score(score: float) -> float:
    # round score to 4 decimal places
    return round(score, 4)