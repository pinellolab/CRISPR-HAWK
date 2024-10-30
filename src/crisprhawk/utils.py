""" 
"""

from exception_handlers import exception_handler

import sys
import os

# define static variables shared across software modules
TOOLNAME = "CRISPR-HAWK"  # tool name
COMMAND = "crisprhawk"  # command line call
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
