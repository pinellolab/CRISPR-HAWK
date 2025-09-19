"""
This module provides classes and functions for representing and scoring CRISPR 
off-target sites.

It includes utilities for parsing off-target report lines, computing CFD and elevation 
scores, and formatting output for downstream analysis.
"""
from .exception_handlers import exception_handler
from .scores.cfdscore.cfdscore import compute_cfd
from .bedfile import BedAnnotation
from .guide import Guide
from .utils import round_score

from typing import Dict, Tuple

import numpy as np

import os


class Offtarget:
    """
    Represents an off-target site for CRISPR guide RNAs, including sequence, 
    position, and scoring information. Provides methods to parse off-target report 
    lines, compute CFD and elevation scores, and generate report output.

    Attributes:
        _chrom (str): Chromosome of the off-target site.
        _pos (int): Genomic position of the off-target site.
        _strand (str): Strand of the off-target site.
        _grna_ (str): Raw gRNA sequence (with PAM).
        _grna (str): Formatted gRNA sequence.
        _spacer_ (str): Raw spacer sequence (with PAM).
        _spacer (str): Formatted spacer sequence.
        _pam (str): PAM sequence for the off-target.
        _mm (int): Number of mismatches.
        _bulge_type (str): Type of bulge (if any).
        _bulge_size (int): Size of the bulge.
        _cfd_score (str): CFD score for the off-target (default "NA").
        _elevation_score (str): Elevation score for the off-target (default "NA").
        _debug (bool): Debug mode flag.
    """

    def __init__(self, reportline: str, pam: str, right: bool, debug: bool) -> None:
        """
        Initializes an Offtarget object by parsing a report line and setting relevant 
        attributes. Sets sequence, position, strand, PAM, mismatch, bulge, and scoring 
        information for the off-target site.

        Args:
            reportline (str): The line from the off-target report to parse.
            pam (str): The PAM sequence for the off-target.
            right (bool): Whether the PAM is on the right side of the sequence.
            debug (bool): Flag to enable debug mode.
        """
        self._debug = debug  # store debug mode
        self._parse_reportline(reportline, pam, right)  # read report line content
        self._pam = pam  # set off-target pam
        self._cfd_score = "NA"  # default cfd score
        self._elevation_score = "NA"  # default elevation score

    def __repr__(self) -> str:
        """
        Returns a string representation of the Offtarget object, including its 
        position, spacer, and strand. This helps with debugging and logging by 
        providing a concise summary of the object's key attributes.

        Returns:
            str: String representation of the Offtarget object.
        """
        return (
            f"<{self.__class__.__name__} object; position={self._pos} "
            f"spacer={self._spacer} strand={self._strand}>"
        )

    def _parse_reportline(self, line: str, pam: str, right: bool) -> None:
        """
        Parses a report line to extract and set off-target attributes such as 
        chromosome, position, strand, sequences, mismatches, and bulge information.
        Updates the Offtarget object's internal state with the parsed values for 
        downstream scoring and reporting.

        Args:
            line (str): The report line containing off-target information.
            pam (str): The PAM sequence for the off-target.
            right (bool): Whether the PAM is on the right side of the sequence.
        """
        fields = line.strip().split()
        self._chrom = fields[3]  # set chromosome
        self._pos = int(fields[4])  # set off-target position
        self._strand = fields[6]  # set off-target strand
        self._grna_, self._grna = _format_sequence(
            fields[1], pam, right
        )  # set grna sequence
        self._spacer_, self._spacer = _format_sequence(
            fields[2], _retrieve_pam(fields[2], len(pam), right), right
        )  # set spacer sequence
        self._mm = int(fields[7])  # set mismatches number
        self._bulge_type = fields[0]  # set bulge type
        self._bulge_size = int(fields[8])  # set bulge size

    def compute_cfd(
        self, mmscores: Dict[str, float], pamscores: Dict[str, float]
    ) -> None:
        """
        Computes the CFD (Cutting Frequency Determination) score for the off-target 
        site using provided mismatch and PAM scores. Updates the object's CFD score 
        attribute with the computed value for downstream reporting and analysis.

        Args:
            mmscores (Dict[str, float]): Dictionary of mismatch scores.
            pamscores (Dict[str, float]): Dictionary of PAM scores.

        Returns:
            None
        """
        self._cfd_score = str(
            round_score(
                compute_cfd(
                    self._grna_.upper(),
                    self._spacer_.upper(),
                    self._spacer[-2:],
                    mmscores,
                    pamscores,
                    self._debug,
                )
            )
        )

    

    def report_line(self) -> str:
        """
        Generates a tab-separated string containing all key attributes of the 
        Offtarget object for reporting purposes. Returns a line with chromosome, 
        position, strand, gRNA, spacer, PAM, mismatch count, bulge size, bulge type, 
        CFD score, and elevation score.

        Returns:
            str: Tab-separated string of Offtarget attributes for report output.
        """
        return "\t".join(
            list(
                map(
                    str,
                    [
                        self._chrom,
                        self._pos,
                        self._strand,
                        self._grna,
                        self._spacer,
                        self._pam,
                        self._mm,
                        self._bulge_size,
                        self._bulge_type,
                        self._cfd_score,
                        self._elevation_score,
                    ],
                )
            )
        )

    @property
    def grna(self) -> str:
        return self._grna

    @property
    def grna_(self) -> str:
        return self._grna_

    @property
    def spacer(self) -> str:
        return self._spacer

    @property
    def cfd(self) -> str:
        return self._cfd_score

    @property
    def elevation(self) -> str:
        return self._elevation_score
    
    @elevation.setter
    def elevation(self, value: float) -> None:
        """
        Sets the elevation score for the off-target site, rounding the value or
        setting to "NA" if not a number.
        
        Raises a TypeError if the provided value is not a float.

        Args:
            value (float): The elevation score to set.

        Raises:
            TypeError: If the provided value is not a float.
        """
        if not isinstance(value, float):
            exception_handler(
                TypeError,
                f"Elevation must be a float, got {type(value).__name__} instead",
                os.EX_DATAERR,
                self._debug,
            )
        self._elevation_score = "NA" if np.isnan(value) else str(round_score(value))
    
def _retrieve_pam(sequence: str, length: int, right: bool) -> str:
    """
    Retrieves the PAM sequence from a given DNA sequence based on the specified 
    length and orientation. Returns the PAM sequence from the start or end of the 
    input sequence depending on the 'right' flag.

    Args:
        sequence (str): The DNA sequence containing the PAM.
        length (int): The length of the PAM sequence to retrieve.
        right (bool): If True, retrieves PAM from the start; if False, from the end.

    Returns:
        str: The extracted PAM sequence.
    """
    return sequence[:length] if right else sequence[-length:]


def _format_sequence(sequence: str, pam: str, right: bool) -> Tuple[str, str]:
    """
    Formats a DNA sequence by extracting the spacer and combining it with the PAM 
    sequence according to orientation. Returns a tuple containing the spacer sequence 
    and the full formatted sequence with PAM.

    Args:
        sequence (str): The DNA sequence containing the PAM and spacer.
        pam (str): The PAM sequence to use for formatting.
        right (bool): If True, PAM is placed at the start; if False, at the end.

    Returns:
        Tuple[str, str]: Spacer sequence and formatted sequence with PAM.
    """
    s_ = sequence[len(pam) :] if right else sequence[: -len(pam)]
    s = f"{pam}{s_}" if right else f"{s_}{pam}"
    return s_, s
