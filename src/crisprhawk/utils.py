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
# dictionary to encode nucleotide strings as iupac characters
IUPAC_ENCODER = {
    "A": "A",
    "T": "T",
    "C": "C",
    "G": "G",
    "AG": "R",
    "GA": "R",
    "CT": "Y", 
    "TC": "Y",
    "GC": "S", 
    "CG": "S",
    "AT": "W", 
    "TA": "W",
    "GT": "K", 
    "TG": "K",
    "AC": "M", 
    "CA": "M",
    "CGT": "B", 
    "GCT": "B", 
    "TGC": "B", 
    "GTC": "B",
    "CTG": "B", 
    "TCG": "B",
    "AGT": "D", 
    "GAT": "D", 
    "TAG": "D", 
    "ATG": "D", 
    "GTA": "D", 
    "TGA": "D",
    "ACT": "H", 
    "CAT": "H", 
    "TCA": "H", 
    "ATC": "H", 
    "CTA": "H", 
    "TAC": "H",
    "ACG": "V", 
    "CAG": "V", 
    "GAC": "V", 
    "AGC": "V", 
    "CGA": "V", 
    "GCA": "V",
    "ACGT": "N",
    "CAGT": "N",
    "GACT": "N",
    "AGCT": "N",
    "CGAT": "N",
    "GCAT": "N",
    "GCTA": "N",
    "CGTA": "N",
    "TGCA": "N",
    "GTCA": "N",
    "CTGA": "N",
    "TCGA": "N",
    "TAGC": "N",
    "ATGC": "N",
    "GTAC": "N",
    "TGAC": "N",
    "AGTC": "N",
    "GATC": "N",
    "CATG": "N",
    "ACTG": "N",
    "TCAG": "N",
    "CTAG": "N",
    "ATCG": "N",
    "TACG": "N",
}


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
