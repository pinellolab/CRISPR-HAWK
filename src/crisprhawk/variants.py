"""
"""

from crisprhawk_error import CrisprHawkVCFError
from exception_handlers import exception_handler

from typing import Optional, List
from pysam import TabixFile, tabix_index

import os

TBI = "tbi"

class VCF:
    def __init__(self, fname: str, debug: bool, vcfidx: Optional[str] = "") -> None:
        self._fname = fname  # store input filename
        self._debug = debug  # store debug mode flag
        self._vcfidx = self._search_index(vcfidx)  # initialize vcf index
        # initialize TabixFile object with the previously computed index
        self._vcf = TabixFile(self._fname, index=self._vcfidx)


    def _search_index(self, vcfidx: Optional[str] = "") -> str:
        # look for index for the current vcf, if not found compute it
        if not vcfidx:
            if _find_tbi(self._fname):  # index found, store it
                return f"{self._fname}.{TBI}"
            # index not found -> compute it de novo and store it in the same folder 
            # as the input vcf 
            try:
                tabix_index(self._fname, preset="vcf", force=True) 
            except OSError as e:
                exception_handler(CrisprHawkVCFError, f"An error occurred while indexing {self._fname}", os.EX_SOFTWARE, self._debug, e)
            assert _find_tbi(self._fname)
            return f"{self._fname}.{TBI}"
        # precomputed vcf index index must be a non empty file
        if not (os.path.isfile(vcfidx) and os.stat(vcfidx).st_size > 0):
            exception_handler(CrisprHawkVCFError, f"Not existing or empty VCF index {vcfidx}", os.EX_DATAERR, self._debug)

    def fetch(self, contig: str, start: int, stop: int) -> List[List[str]]:
        try:  # extract variants in the input range from vcf file
            # return a list of variants data as list of fields
            return [v.strip().split() for v in self._vcf.fetch(contig, start, stop)]
        except ValueError as e:
            exception_handler(CrisprHawkVCFError, f"Variants extraction failed for coordinates {contig}:{start}-{stop}", os.EX_DATAERR, self._debug, e)



def _find_tbi(vcf: str) -> bool:
    # avoid unexpected crashes due to file location
    vcfindex = f"{os.path.abspath(vcf)}.{TBI}" 
    if os.path.exists(vcfindex):  # index must be a non empty file
        return os.path.isfile(vcfindex) and os.stat(vcfindex).st_size > 0
    return False