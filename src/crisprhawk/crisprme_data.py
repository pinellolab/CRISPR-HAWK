"""
This module provides functions for preparing CRISPRme input data files from
CRISPR-HAWK's guide reports.

It includes utilities to determine PAM orientation, parse guide reports, encode
PAM sequences, and generate the necessary PAM and guide files for downstream analysis.
"""

from .crisprhawk_error import CrisprHawkPrepareDataError
from .exception_handlers import exception_handler
from .utils import IUPAC_ENCODER

from typing import Tuple, List, Union

import re
import os


def is_pam_right(header: str) -> bool:
    """Determines if the PAM sequence occurs upstream or downstream of the guide
    in the header.

    Returns True if the PAM column occurs before the sgRNA sequence column, indicating
    an upstream PAM. Returns False otherwise.

    Args:
        header: The header string from the report file.

    Returns:
        bool: True if PAM is upstream of the guide, False otherwise.
    """
    # if pam column occurs before, pam occurs upstream guide
    sgrna_idx = header.find("sgRNA_sequence")
    pam_idx = header.find("pam")
    assert sgrna_idx != pam_idx
    return pam_idx < sgrna_idx


def read_guides_report(report: str, debug: bool) -> Tuple[bool, List[List[str]]]:
    """Reads a CRISPR-HAWK guide report file and determines the orientation of
    the PAM sequence.

    Returns a tuple containing a boolean indicating if the PAM is upstream and a
    list of guide data fields.

    Args:
        report: Path to the report file.
        debug: Boolean flag to enable debug mode.

    Returns:
        Tuple[bool, List[List[str]]]: A tuple with a boolean for PAM orientation
            and a list of guide data.
    """
    try:
        with open(report, mode="r") as infile:
            # assess if pam occurs upstream guides
            right = is_pam_right(infile.readline().strip())
            fields = [line.strip().split()[:6] for line in infile]  # first 6 cols
    except IOError as e:
        exception_handler(
            CrisprHawkPrepareDataError,
            f"Reading report {report} failed",
            os.EX_IOERR,
            debug,
            e,
        )
    return right, fields


def _replacer(match: re.Match) -> Union[str, None]:
    """Replaces matched IUPAC codes in a PAM class string with their encoded values.

    Returns the encoded value from IUPAC_ENCODER if the key exists, otherwise
    returns the original match.

    Args:
        match: A regular expression match object containing the IUPAC code.

    Returns:
        Union[str, None]: The encoded value or the original match string.
    """
    key = match.group(1)
    return IUPAC_ENCODER.get(key, match.group(0))


def solve_pam(pamclass: str) -> str:
    """Converts IUPAC codes in a PAM class string to their encoded values.

    Returns a PAM string with all IUPAC codes replaced by their corresponding
    encoded values.

    Args:
        pamclass: The PAM class string containing IUPAC codes.

    Returns:
        str: The PAM string with encoded values.
    """
    return re.sub(r"\[([^\]]+)\]", _replacer, pamclass)  # type: ignore


def create_pam_file(
    pamclass: str, right: bool, guidelen: int, outdir: str, debug: bool
) -> None:
    """Creates a PAM file with the specified orientation and guide length.

    Writes a file containing the PAM and guide sequence in the correct orientation
    for CRISPR-HAWK analysis.

    Args:
        pamclass: The PAM class string containing IUPAC codes.
        right: Boolean indicating if the PAM is upstream (True) or downstream (False).
        guidelen: The length of the guide sequence.
        outdir: The output directory for the PAM file.
        debug: Boolean flag to enable debug mode.

    Returns:
        None
    """
    pam = solve_pam(pamclass)  # retrieve original pam with iupac
    try:  # write pam file
        with open(os.path.join(outdir, f"{pam}.txt"), mode="w") as outfile:
            gseq = "N" * guidelen  # add as many Ns as the guide length
            pamfull = (
                f"{pam}{gseq}\t{-len(pam)}" if right else f"{gseq}{pam}\t{len(pam)}"
            )
            outfile.write(f"{pamfull}\n")  # write PAM
    except IOError as e:
        exception_handler(
            CrisprHawkPrepareDataError,
            f"PAM file creation failed for PAM {pamclass}",
            os.EX_IOERR,
            debug,
            e,
        )
    if (
        not os.path.isfile(os.path.join(outdir, f"{pam}.txt"))
        or os.stat(os.path.join(outdir, f"{pam}.txt")).st_size <= 0
    ):
        exception_handler(
            CrisprHawkPrepareDataError,
            f"PAM file empty or not created for PAM {pamclass}",
            os.EX_IOERR,
            debug,
        )


def create_guide_files(
    guides_data: List[List[str]], right: bool, outdir: str, debug: bool
) -> None:
    """Creates guide files for each guide in the provided data with the correct
    PAM orientation.

    Writes a file for each guide containing the guide and PAM sequence in the correct
    orientation for CRISPR-HAWK analysis.

    Args:
        guides_data: List of guide data, each as a list of strings.
        right: Boolean indicating if the PAM is upstream (True) or downstream (False).
        outdir: The output directory for the guide files.
        debug: Boolean flag to enable debug mode.

    Returns:
        None
    """
    pidx = 3 if right else 4  # assess pam column index
    gidx = 4 if right else 3  # assess guide sequence column index
    for gdata in guides_data:
        gfname = os.path.join(
            outdir, f"{gdata[0]}_{gdata[1]}_{gdata[2]}_{gdata[gidx]}_{gdata[pidx]}.txt"
        )
        try:
            with open(gfname, mode="w") as outfile:
                pamseq = "N" * len(gdata[pidx])  # add as many Ns as the pam length
                guidefull = (
                    f"{pamseq}{gdata[gidx]}" if right else f"{gdata[gidx]}{pamseq}"
                )
                outfile.write(f"{guidefull}\n")
        except IOError as e:
            exception_handler(
                CrisprHawkPrepareDataError,
                f"Guide file creation failed for guide {gfname}",
                os.EX_IOERR,
                debug,
                e,
            )
        if not os.path.isfile(gfname) or os.stat(gfname).st_size <= 0:
            exception_handler(
                CrisprHawkPrepareDataError,
                f"Guide file empty or not created for guide {gfname}",
                os.EX_IOERR,
                debug,
            )


def prepare_data_crisprme(
    report: str, create_pam: bool, outdir: str, debug: bool
) -> None:
    """Prepares CRISPRme input data files from a CRISPR-HAWK's guide report.

    Reads the guide report, optionally creates a PAM file, and generates guide
    files for each guide in the report.

    Args:
        report: Path to the guide report file.
        create_pam: Boolean indicating whether to create a PAM file.
        outdir: The output directory for the generated files.
        debug: Boolean flag to enable debug mode.

    Returns:
        None
    """
    # read input guides report
    isright, guides_data = read_guides_report(report, debug)
    if create_pam:  # pam file creation requested
        gdata = guides_data[0]  # retrieve guide's data
        gidx = 4 if isright else 3  # assess guide sequence column index
        create_pam_file(gdata[5], isright, len(gdata[gidx]), outdir, debug)
    # create guide file for each guide in the report
    create_guide_files(guides_data, isright, outdir, debug)
