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
import zipfile
import shutil
import sys
import io
import os

# define static variables shared across software modules
TOOLNAME = "CRISPR-HAWK"  # tool name
COMMAND = "crisprhawk"  # command line call
# define OS systems
OSSYSTEMS = ["Linux", "Darwin", "Windows"]
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
# reports prefix name
GUIDESREPORTPREFIX = "crisprhawk_guides"
CANDIDATEGUIDESREPORTPREFIX = "crisprhawk_candidate_guides"


# define utils functions
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


def create_folder(dirname: str) -> str:
    """Create a new directory with the specified name.

    Ensures the directory exists by creating it if necessary and returns its path.

    Args:
        dirname (str): The name of the directory to create.

    Returns:
        str: The path to the created directory.
    """
    os.makedirs(dirname)
    assert os.path.isdir(dirname)
    return dirname


def remove_folder(dirname: str) -> None:
    """Remove a directory and its contents.

    Attempts to delete the specified directory and all its contents. Raises an
    error if the operation fails.

    Args:
        dirname (str): The path to the directory to remove.

    Raises:
        OSError: If the directory cannot be removed.
    """
    try:
        subprocess.run(["rm", "-rf", dirname], check=True, capture_output=True)
    except subprocess.CalledProcessError as e:  # always trace this error
        raise OSError(f"Failed to clean up folder {dirname}") from e


def remove_file(filename: str) -> None:
    """Remove a file from the filesystem.

    Attempts to delete the specified file. Raises an error if the operation fails.

    Args:
        filename (str): The path to the file to remove.

    Raises:
        OSError: If the file cannot be removed.
    """
    try:
        subprocess.run(["rm", "-rf", filename], check=True, capture_output=True)
    except subprocess.CalledProcessError as e:  # always trace this error
        raise OSError(f"Failed to remove file {filename}") from e


def _uncompress_azimuth_models(azimuthdir: str) -> None:
    """Uncompress the Azimuth model ZIP archive in the specified directory.

    Extracts the Azimuth saved models ZIP file into the given directory and
    removes the archive after extraction. Raises an error if the file is missing
    or extraction fails.

    Args:
        azimuthdir (str): The directory containing the Azimuth models ZIP archive.

    Raises:
        FileNotFoundError: If the Azimuth models ZIP file does not exist.
        zipfile.BadZipFile: If the ZIP file is invalid or extraction fails.
    """
    models_zip = os.path.join(azimuthdir, "saved_models.zip")  # azimuth models
    if not os.path.isfile(models_zip):  # always trace these errors
        raise FileNotFoundError(f"Cannot find Azimuth models: {models_zip}")
    try:
        with zipfile.ZipFile(models_zip, mode="r") as zipref:
            zipref.extractall(path=azimuthdir)  # extract in azimuth directory
    except (zipfile.BadZipFile, RuntimeError) as e:
        raise zipfile.BadZipFile(
            f"An error occurred while unzipping Azimuth models: {models_zip}"
        ) from e
    remove_file(models_zip)


def _uncompress_cfd_models(cfdscoredir: str) -> None:
    """Uncompress the CFD model ZIP archive in the specified directory.

    Extracts the CFD models ZIP file into the given directory and removes the
    archive after extraction. Raises an error if the file is missing or extraction
    fails.

    Args:
        cfdscoredir (str): The directory containing the CFD models ZIP archive.

    Raises:
        FileNotFoundError: If the CFD models ZIP file does not exist.
        zipfile.BadZipFile: If the ZIP file is invalid or extraction fails.
    """
    cfdmodels_zip = os.path.join(cfdscoredir, "models.zip")  # cfd models
    if not os.path.isfile(cfdmodels_zip):  # always trace these errors
        raise FileNotFoundError(f"Cannot find CFD models: {cfdmodels_zip}")
    try:
        with zipfile.ZipFile(cfdmodels_zip, mode="r") as zipref:
            zipref.extractall(path=cfdscoredir)  # extract in cfd directory
    except (zipfile.BadZipFile, RuntimeError) as e:
        raise zipfile.BadZipFile(
            f"An error occurred while unzipping CFD models: {cfdmodels_zip}"
        ) from e
    remove_file(cfdmodels_zip)


def _uncompress_deepcpf1_models(deepcpf1dir: str) -> None:
    """Uncompress the DeepCpf1 model ZIP archive in the specified directory.

    Extracts the DeepCpf1 weights ZIP file into the given directory and removes
    the archive after extraction. Raises an error if the file is missing or
    extraction fails.

    Args:
        deepcpf1dir (str): The directory containing the DeepCpf1 weights ZIP archive.

    Raises:
        FileNotFoundError: If the DeepCpf1 weights ZIP file does not exist.
        zipfile.BadZipFile: If the ZIP file is invalid or extraction fails.
    """
    deepcpf1_weights_zip = os.path.join(deepcpf1dir, "weights.zip")  # deepCpf1 models
    if not os.path.isfile(deepcpf1_weights_zip):  # always trace these errors
        raise FileNotFoundError(f"Cannot find DeepCpf1 models: {deepcpf1_weights_zip}")
    try:
        with zipfile.ZipFile(deepcpf1_weights_zip, mode="r") as zipref:
            zipref.extractall(path=deepcpf1dir)  # extract in DeepCpf1 directory
    except (zipfile.BadZipFile, RuntimeError) as e:
        raise zipfile.BadZipFile(
            f"An error occurred while unzipping DeepCpf1 models: {deepcpf1_weights_zip}"
        ) from e
    remove_file(deepcpf1_weights_zip)


def _uncompress_elevation_models(elevationdir: str) -> None:
    """Uncompress the Elevation model and data ZIP archives in the specified
    directory.

    Extracts the Elevation models and data ZIP files into the given directory and
    removes the archives after extraction. Raises an error if any file is missing
    or extraction fails.

    Args:
        elevationdir (str): The directory containing the Elevation ZIP archives.

    Raises:
        FileNotFoundError: If any Elevation ZIP file does not exist.
        zipfile.BadZipFile: If any ZIP file is invalid or extraction fails.
    """
    _uncompress_elevation_model(
        elevationdir,
        "models.zip",
        "Cannot find Elevation models:",
        "An error occurred while unzipping Elevation models:",
    )
    _uncompress_elevation_model(
        elevationdir,
        "CRISPR.zip",
        "Cannot find Elevation data:",
        "An error occurred while unzipping Elevation data:",
    )


def _uncompress_elevation_model(
    elevationdir, file_zipped: str, errmsg1: str, errmsg2: str
) -> None:
    """Uncompress a specific Elevation ZIP archive in the given directory.

    Extracts the specified ZIP file into the provided directory and removes the
    archive after extraction. Raises an error if the file is missing or extraction
    fails.

    Args:
        elevationdir (str): The directory containing the ZIP archive.
        file_zipped (str): The name of the ZIP file to extract.
        errmsg1 (str): Error message for missing file.
        errmsg2 (str): Error message for extraction failure.

    Raises:
        FileNotFoundError: If the ZIP file does not exist.
        zipfile.BadZipFile: If the ZIP file is invalid or extraction fails.
    """
    elevation_models_zip = os.path.join(elevationdir, file_zipped)
    if not os.path.isfile(elevation_models_zip):  # always trace these errors
        raise FileNotFoundError(f"{errmsg1} {elevation_models_zip}")
    try:
        with zipfile.ZipFile(elevation_models_zip, mode="r") as zipref:
            zipref.extractall(path=elevationdir)  # extract in cfd directory
    except (zipfile.BadZipFile, RuntimeError) as e:
        raise zipfile.BadZipFile(f"{errmsg2}{elevation_models_zip}") from e
    remove_file(elevation_models_zip)


def prepare_package() -> None:
    """Prepare the package by extracting required model and data ZIP files.

    Checks for the presence of required model and data directories. If any are
    missing, extracts the corresponding ZIP archives to ensure all scoring algorithms
    have access to necessary resources.
    """
    # at first run uncompress ZIP files containing models and data used by
    # the scoring algorithms used by crisprhawk
    scoresdir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "scores")
    azimuthdir = os.path.join(scoresdir, "azimuth")  # azimuth
    if not os.path.isdir(azimuthdir):  # always trace these errors
        raise FileNotFoundError("Cannot find Azimuth score modules")
    if not os.path.isdir(os.path.join(azimuthdir, "saved_models")):
        warning("Extracting Azimuth models. This may take some time", 1)
        _uncompress_azimuth_models(azimuthdir)  # uncompress azimuth models
    cfdscoredir = os.path.join(scoresdir, "cfdscore")  # cfd
    if not os.path.isdir(cfdscoredir):  # always trace these errors
        raise FileNotFoundError("Cannot find CFD score modules")
    if not os.path.isdir(os.path.join(cfdscoredir, "models")):
        warning("Extracting CFD models. This may take some time", 1)
        _uncompress_cfd_models(cfdscoredir)  # uncompress CFD models
    deepcpf1dir = os.path.join(scoresdir, "deepCpf1")  # deepCpf1
    if not os.path.isdir(deepcpf1dir):  # always trace these errors
        raise FileNotFoundError("Cannot find DeepCpf1 score modules")
    if not os.path.isdir(os.path.join(deepcpf1dir, "weights")):
        warning("Extracting DeepCpf1 models. This may take some time", 1)
        _uncompress_deepcpf1_models(deepcpf1dir)  # uncompress deepCpf1 models
    elevationdir = os.path.join(scoresdir, "elevation")  # Elevation
    if not os.path.isdir(elevationdir):  # always trace these errors
        raise FileNotFoundError("Cannot find Elevation score modules")
    if not os.path.isdir(os.path.join(elevationdir, "models")) and not os.path.isdir(
        os.path.join(elevationdir, "CRISPR")
    ):
        warning("Extracting Elevation models and data. This may take some time", 1)
        _uncompress_elevation_models(elevationdir)  # uncompress elevation models


def command_exists(command: str) -> bool:
    """Check if a command exists in the system's PATH.

    Returns True if the specified command is found in the system's executable
    search path, otherwise False.

    Args:
        command (str): The command to check for existence.

    Returns:
        bool: True if the command exists, False otherwise.
    """
    return bool(shutil.which(command))


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
    size = len(lst)  # compute list size
    chunk_size = max(1, size // threads)  # compute chunk sizes
    chunks = []  # create chunks
    for i in range(0, size, chunk_size):
        end_idx = min(i + chunk_size, size)
        chunks.append((i, lst[i:end_idx]))
    return chunks
