""" """

from .exception_handlers import exception_handler

from typing import Any, List
from itertools import permutations
from colorama import Fore

import subprocess
import tempfile
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
            ValueError,  # type: ignore
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


def adjust_guide_position(pos: int, guidelen: int, pamlen: int, right: bool) -> int:
    """Adjust the guide position based on orientation.

    Returns the original position if right is True, otherwise subtracts the
    guide length.

    Args:
        pos: The original position.
        guidelen: The length of the guide.
        pamlen: The length of the PAM sequence (unused).
        right: Boolean indicating orientation.

    Returns:
        The adjusted guide position as an integer.
    """
    return pos if right else pos - guidelen


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


def match_iupac(seq: str, pattern: str) -> bool:
    """Check if a nucleotide sequence matches a given IUPAC pattern.

    Compares each nucleotide in the sequence to the corresponding IUPAC code in
    the pattern returning True if all nucleotides are compatible with the pattern.

    Args:
        seq: The nucleotide sequence to check.
        pattern: The IUPAC pattern to match against.

    Returns:
        True if the sequence matches the IUPAC pattern, False otherwise.
    """
    if len(seq) != len(pattern):
        return False
    seq = seq.upper()  # ensure upper cases
    pattern = pattern.upper()
    return all(snt in list(IUPACTABLE[pnt]) for snt, pnt in zip(seq, pattern))


def dna2rna(sequence: str) -> str:
    """Convert a DNA sequence to its RNA equivalent.

    Replaces all occurrences of thymine ('T' or 't') with uracil ('U' or 'u') in
    the input sequence.

    Args:
        sequence: The DNA sequence to convert.

    Returns:
        The RNA sequence as a string.
    """
    return sequence.replace("T", "U").replace("t", "u")


@contextlib.contextmanager
def suppress_stdout():
    """Context manager to suppress standard output.

    Temporarily redirects sys.stdout to an in-memory buffer.

    Returns:
        None
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

    Returns:
        None
    """
    stderr_channel = sys.stderr
    sys.stderr = io.StringIO()
    try:
        yield
    finally:
        sys.stderr = stderr_channel


def create_temp_folder(dirname: str) -> str:
    return tempfile.mkdtemp(prefix=f"crisprhawk_{dirname}_")


def create_folder(dirname: str) -> str:
    os.makedirs(dirname)
    assert os.path.isdir(dirname)
    return dirname


def remove_folder(dirname: str) -> None:
    try:
        subprocess.run(["rm", "-rf", dirname], check=True, capture_output=True)
    except subprocess.CalledProcessError as e:  # always trace this error
        raise OSError(f"Failed to clean up folder {dirname}") from e
    except Exception as e:
        raise Exception(f"Unexpected error while cleaning up folder {dirname}") from e


def remove_file(filename: str) -> None:
    try:
        subprocess.run(["rm", "-rf", filename], check=True, capture_output=True)
    except subprocess.CalledProcessError as e:  # always trace this error
        raise OSError(f"Failed to remove file {filename}") from e
    except Exception as e:
        raise Exception(f"Unexpected error while removing file {filename}") from e


def _uncompress_azimuth_models(azimuthdir: str) -> None:
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
    elevation_models_zip = os.path.join(elevationdir, "models.zip")  # elvation models
    if not os.path.isfile(elevation_models_zip):  # always trace these errors
        raise FileNotFoundError(f"Cannot find Elevation models: {elevation_models_zip}")
    try:
        with zipfile.ZipFile(elevation_models_zip, mode="r") as zipref:
            zipref.extractall(path=elevationdir)  # extract in cfd directory
    except (zipfile.BadZipFile, RuntimeError) as e:
        raise zipfile.BadZipFile(
            f"An error occurred while unzipping Elevation models: {elevation_models_zip}"
        ) from e
    remove_file(elevation_models_zip)
    elevation_data_zip = os.path.join(elevationdir, "CRISPR.zip")  # elvation data
    if not os.path.isfile(elevation_data_zip):  # always trace these errors
        raise FileNotFoundError(f"Cannot find Elevation data: {elevation_data_zip}")
    try:
        with zipfile.ZipFile(elevation_data_zip, mode="r") as zipref:
            zipref.extractall(path=elevationdir)  # extract in cfd directory
    except (zipfile.BadZipFile, RuntimeError) as e:
        raise zipfile.BadZipFile(
            f"An error occurred while unzipping Elevation data: {elevation_data_zip}"
        ) from e
    remove_file(elevation_data_zip)


def prepare_package() -> None:
    # at first run uncompress ZIP files containing models and data used by
    # the scoring algorithms used by crisprhawk
    scoresdir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "scores")
    azimuthdir = os.path.join(scoresdir, "azimuth")  # azimuth
    if not os.path.isdir(azimuthdir):  # always trace these errors
        raise FileNotFoundError(f"Cannot find Azimuth score modules")
    if not os.path.isdir(os.path.join(azimuthdir, "saved_models")):
        warning("Extracting Azimuth models. This may take some time", 1)
        _uncompress_azimuth_models(azimuthdir)  # uncompress azimuth models
    cfdscoredir = os.path.join(scoresdir, "cfdscore")  # cfd
    if not os.path.isdir(cfdscoredir):  # always trace these errors
        raise FileNotFoundError(f"Cannot find CFD score modules")
    if not os.path.isdir(os.path.join(cfdscoredir, "models")):
        warning("Extracting CFD models. This may take some time", 1)
        _uncompress_cfd_models(cfdscoredir)  # uncompress CFD models
    deepcpf1dir = os.path.join(scoresdir, "deepCpf1")  # deepCpf1
    if not os.path.isdir(deepcpf1dir):  # always trace these errors
        raise FileNotFoundError(f"Cannot find DeepCpf1 score modules")
    if not os.path.isdir(os.path.join(deepcpf1dir, "weights")):
        warning("Extracting DeepCpf1 models. This may take some time", 1)
        _uncompress_deepcpf1_models(deepcpf1dir)  # uncompress deepCpf1 models
    elevationdir = os.path.join(scoresdir, "elevation")  # Elevation
    if not os.path.isdir(elevationdir):  # always trace these errors
        raise FileNotFoundError(f"Cannot find Elevation score modules")
    if not os.path.isdir(os.path.join(elevationdir, "models")) and not os.path.isdir(
        os.path.join(elevationdir, "CRISPR")
    ):
        warning("Extracting Elevation models and data. This may take some time", 1)
        _uncompress_elevation_models(elevationdir)  # uncompress elevation models


def command_exists(command: str) -> bool:
    if shutil.which(command):
        return True
    return False


def is_lowercase(sequence: str) -> bool:
    return any(c.islower() for c in sequence)
