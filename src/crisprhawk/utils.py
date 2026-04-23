"""
Utility functions and constants for the CRISPR-HAWK tool.

This module provides helper functions for file and directory management, sequence
manipulation, IUPAC matching, and model extraction. It also defines shared constants
and static variables used across the CRISPR-HAWK software.
"""

from .exception_handlers import exception_handler

from typing import Any, List, Tuple
from itertools import permutations
from colorama import Fore

import subprocess
import contextlib
import shutil
import sys
import io
import os


# ==============================================================================
#
# Define constant variables
#
# ==============================================================================

# define static variables shared across software modules
TOOLNAME = "CRISPR-HAWK"  # tool name
COMMAND = "crisprhawk"  # command line call

# define OS systems
OSSYSTEMS = ["Linux", "Darwin", "Windows"]

# define verbosity levels
VERBOSITYLVL = [0, 1, 2, 3]

# define dna alphabet
DNA = ["A", "C", "G", "T", "N"]

# define complete iupac alphabet
IUPAC = DNA + ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"]

# define reverse complement dictionary
RC = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
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
    "u": "a",
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

# define dictionary to encode nucleotides combinations as iupac characters
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

# define dictionary to encode nucleotide strings as iupac characters
IUPAC_ENCODER = {
    perm: k
    for k, v in IUPACTABLE.items()
    for perm in {"".join(p) for p in permutations(v)}
}

# define strandness
STRAND = [0, 1]  # strands directions: 0 -> 5'-3'; 1 -> 3'-5'


# define reports prefix name
GUIDESREPORTPREFIX = "crisprhawk_guides"
CANDIDATEGUIDESREPORTPREFIX = "crisprhawk_candidate_guides"


# ==============================================================================
#
# Define utilities functions (config files)
#
# ==============================================================================


def reverse_complement(sequence: str, debug: bool) -> str:
    """Return the reverse complement of a nucleotide sequence.

    Computes the reverse complement of the input sequence using the defined nucleotide
    mapping. Raises an error if an invalid character is encountered.

    Args:
        sequence (str): The nucleotide sequence to reverse complement.
        debug (bool): Boolean indicating whether to provide debug information on error.

    Returns:
        str: The reverse complement of the input sequence as a string.
    """
    try:
        return "".join([RC[nt] for nt in sequence[::-1]])
    except KeyError as e:
        exception_handler(
            ValueError,  # type: ignore
            f"Failed reverse complement on {sequence}",
            os.EX_DATAERR,
            debug,
            e,
        )


def warning(message: str, verbosity: int) -> None:
    """Display a warning message if the verbosity level is sufficient.

    Writes a warning message to standard error if the verbosity is at least level 1.

    Args:
        message (str): The warning message to display.
        verbosity (int): The current verbosity level.

    Returns:
        None
    """
    if verbosity >= VERBOSITYLVL[1]:
        sys.stderr.write(f"{Fore.YELLOW}WARNING: {message}.{Fore.RESET}\n")
    return


def print_verbosity(message: str, verbosity: int, verbosity_threshold: int) -> None:
    """Print a message if the verbosity level meets the threshold.

    Writes the message to standard output if the current verbosity is greater
    than or equal to the specified threshold.

    Args:
        message (str): The message to print.
        verbosity (int): The current verbosity level.
        verbosity_threshold (int): The minimum verbosity level required to print
            the message.

    Returns:
        None
    """
    if verbosity >= verbosity_threshold:
        sys.stdout.write(f"{message}\n")
    return


def adjust_guide_position(pos: int, guidelen: int, pamlen: int, right: bool) -> int:
    """Adjust the guide position based on orientation.

    Returns the adjusted position for a guide depending on whether it is on the
    right or left strand.

    Args:
        pos (int): The original position.
        guidelen (int): The length of the guide.
        pamlen (int): The length of the PAM sequence.
        right (bool): Indicates if the guide is on the right strand.

    Returns:
        int: The adjusted guide position.
    """
    return pos if right else pos - guidelen


def round_score(score: float) -> float:
    """Round a score to four decimal places.

    Returns the input score rounded to four decimal places as a float.

    Args:
        score (float): The score to round.

    Returns:
        float: The rounded score.
    """
    return round(score, 4)  # round score to 4 decimal places


def flatten_list(lst: List[List[Any]]) -> List[Any]:
    """Flatten a list of lists into a single list.

    Combines all elements from nested lists into a single flat list.

    Args:
        lst (List[List[Any]]): The list of lists to flatten.

    Returns:
        List[Any]: The flattened list.
    """
    return [e for sublist in lst for e in sublist]


def match_iupac(seq: str, pattern: str) -> bool:
    """Check if a sequence matches a pattern using IUPAC nucleotide codes.

    Compares two sequences and returns True if the sequence matches the pattern
    according to IUPAC codes.

    Args:
        seq (str): The nucleotide sequence to check.
        pattern (str): The IUPAC pattern to match against.

    Returns:
        bool: True if the sequence matches the pattern, False otherwise.
    """
    if len(seq) != len(pattern):
        return False
    seq = seq.upper()  # ensure upper cases
    pattern = pattern.upper()
    return all(snt in list(IUPACTABLE[pnt]) for snt, pnt in zip(seq, pattern))


def dna2rna(sequence: str) -> str:
    """Convert a DNA sequence to its RNA equivalent.

    Replaces all occurrences of 'T' with 'U' and 't' with 'u' in the input sequence.

    Args:
        sequence (str): The DNA sequence to convert.

    Returns:
        str: The RNA sequence.
    """
    return sequence.replace("T", "U").replace("t", "u")


@contextlib.contextmanager
def suppress_stdout():
    """Context manager to suppress standard output.

    Temporarily redirects sys.stdout to an in-memory buffer.
    """
    stdout_channel = sys.stdout
    sys.stdout = io.StringIO()
    try:
        yield
    finally:
        sys.stdout = stdout_channel


@contextlib.contextmanager
def suppress_stderr():
    """Context manager to suppress standard error output.

    Temporarily redirects sys.stderr to an in-memory buffer.
    """
    stderr_channel = sys.stderr
    sys.stderr = io.StringIO()
    try:
        yield
    finally:
        sys.stderr = stderr_channel


def create_folder(dirname: str, exist_ok: bool = False) -> str:
    try:
        os.makedirs(dirname, exist_ok=exist_ok)
    except OSError as e:
        exception_handler(OSError, f"Failed creating folder {dirname}", os.EX_OSERR, True)
    return dirname


def is_lowercase(sequence: str) -> bool:
    """Check if a sequence contains any lowercase characters.

    Returns True if at least one character in the sequence is lowercase, otherwise
    False.

    Args:
        sequence (str): The sequence to check.

    Returns:
        bool: True if the sequence contains lowercase characters, False otherwise.
    """
    return any(c.islower() for c in sequence)


def calculate_chunks(lst: List[Any], threads: int) -> List[Tuple[int, List[Any]]]:
    """Split a list into indexed chunks for parallel processing.

    Returns a list of (start_index, chunk) tuples that partition the input list
    into approximately equal-sized segments based on the number of threads.

    Args:
        lst (List[Any]): The list to be divided into chunks.
        threads (int): The desired number of chunks, typically matching thread count.

    Returns:
        List[Tuple[int, List[Any]]]: A list of tuples where each tuple contains
            the starting index in the original list and the corresponding chunk.
    """
    size = len(lst)  # compute list size
    chunk_size = max(1, size // threads)  # compute chunk sizes
    chunks = []  # create chunks
    for i in range(0, size, chunk_size):
        end_idx = min(i + chunk_size, size)
        chunks.append((i, lst[i:end_idx]))
    return chunks


def remove_file(filename: str) -> None:
    """Remove a file from the filesystem.

    Attempts to delete the specified file. Raises an error if the operation fails.

    Args:
        filename (str): The path to the file to remove.

    Raises:
        OSError: If the file cannot be removed.
    """
    try:
        os.remove(filename)
    except OSError as e:  # always trace this error
        exception_handler(
            OSError, f"Failed to remove file {filename}", os.EX_OSERR, True, e
        )


def remove_file_silent(fname: str) -> None:
    """Remove a file without raising errors if it does not exist.

    Attempts to delete the specified file and silently ignores any OSError that
    occurs during removal, such as a missing file.

    Args:
        fname (str): The path to the file to remove.
    """
    with contextlib.suppress(OSError):
        os.remove(fname)

def rename_file(orig: str, dest: str) -> None:
    try:
        os.rename(orig, dest)
    except OSError as e:
        exception_handler(OSError, f"Failed to rename file {orig}", os.EX_OSERR, True, e)


def remove_folder(folder: str) -> None:
    try:
        shutil.rmtree(folder)
    except OSError as e:  # always trace this error
        exception_handler(
            OSError, f"Failed to remove directory {folder}", os.EX_OSERR, True, e
        )