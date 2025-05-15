""" """

from exception_handlers import exception_handler

from typing import Any, List
from itertools import permutations
from colorama import Fore

import sys
import os

# define static variables shared across software modules
TOOLNAME = "CRISPR-HAWK"  # tool name
COMMAND = "crisprhawk"  # command line call
# define verbosity levels
VERBOSITYLVL = [0, 1, 2, 3]
# dna alphabet
DNA = ["A", "C", "G", "T", "N"]
# complete iupac alphabet
IUPAC = DNA + ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"]
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
    """Return the reverse complement of a nucleotide sequence.

    Computes the reverse complement of the input DNA or RNA sequence using the
    RC dictionary. Handles invalid nucleotides by raising an exception.

    Args:
        sequence: The nucleotide sequence to reverse complement.
        debug: Flag to enable debug mode.

    Returns:
        The reverse complement of the input sequence as a string.

    Raises:
        ValueError: If the sequence contains invalid nucleotides.
    """
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
    """Display a warning message if the verbosity level is sufficient.

    Prints a formatted warning message to standard error if the verbosity
    threshold is met.

    Args:
        message: The warning message to display.
        verbosity: The current verbosity level.
    """
    if verbosity >= VERBOSITYLVL[1]:
        sys.stderr.write(f"{Fore.YELLOW}WARNING: {message}.{Fore.RESET}\n")
    return


def print_verbosity(message: str, verbosity: int, verbosity_threshold: int) -> None:
    """Print a message if the verbosity level meets the threshold.

    Outputs the provided message to standard output if the current verbosity is
    greater than or equal to the specified threshold.

    Args:
        message: The message to print.
        verbosity: The current verbosity level.
        verbosity_threshold: The minimum verbosity level required to print the
            message.
    """
    if verbosity >= verbosity_threshold:
        sys.stdout.write(f"{message}\n")
    return


def adjust_guide_position(pos: int, pos_rel: List[int]) -> int:
    """Find the index of a position in a list of relative positions.

    Returns the index of the first occurrence of pos in pos_rel, or -1 if not
    found.

    Args:
        pos: The position to search for.
        pos_rel: The list of relative positions.

    Returns:
        The index of pos in pos_rel, or -1 if not found.
    """
    return next((i for i, p in enumerate(pos_rel) if p == pos), -1)


def round_score(score: float) -> float:
    """Round a score to four decimal places.

    Returns the input score rounded to four decimal places for consistent
    reporting.

    Args:
        score: The score to round.

    Returns:
        The rounded score as a float.
    """
    # round score to 4 decimal places
    return round(score, 4)


def flatten_list(lst: List[List[Any]]) -> List[Any]:
    """Flattens a list of lists into a single list.

    Args:
        lst: The list of lists to flatten.
    Returns:
        A new list containing all the elements of the sublists in a single flattened list.
    """
    return [e for sublist in lst for e in sublist]