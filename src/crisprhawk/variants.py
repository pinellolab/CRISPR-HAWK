"""
"""

from crisprhawk_error import CrisprHawkVCFError
from exception_handlers import exception_handler
from utils import print_verbosity, VERBOSITYLVL

from typing import Optional, List
from pysam import TabixFile, tabix_index

import os

TBI = "tbi"


class VCF:
    def __init__(
        self, fname: str, verbosity: int, debug: bool, vcfidx: Optional[str] = ""
    ) -> None:
        self._fname = fname  # store input filename
        self._debug = debug  # store debug mode flag
        self._verbosity = verbosity  # store verbosity level
        self._vcfidx = self._search_index(vcfidx)  # initialize vcf index
        # initialize TabixFile object with the previously computed index
        self._vcf = TabixFile(self._fname, index=self._vcfidx)
        self._contig = self._vcf.contigs[0]  # assume vcf data about one contig
        self._samples = self._vcf.header[-1].strip().split()[9:]  # recover samples

    def _search_index(self, vcfidx: Optional[str] = "") -> str:
        # look for index for the current vcf, if not found compute it
        if not vcfidx:
            if _find_tbi(self._fname):  # index found, store it
                print_verbosity(
                    f"Tabix index found for {self._fname}",
                    self._verbosity,
                    VERBOSITYLVL[3],
                )
                return f"{self._fname}.{TBI}"
            # index not found -> compute it de novo and store it in the same folder
            # as the input vcf
            print_verbosity(
                f"Tabix index not found. Computing index for {self._fname}",
                self._verbosity,
                VERBOSITYLVL[2],
            )
            try:
                tabix_index(self._fname, preset="vcf", force=True)
            except OSError as e:
                exception_handler(
                    CrisprHawkVCFError,
                    f"An error occurred while indexing {self._fname}",
                    os.EX_SOFTWARE,
                    self._debug,
                    e,
                )
            assert _find_tbi(self._fname)
            print_verbosity(
                f"Tabix index computed for {self._fname}",
                self._verbosity,
                VERBOSITYLVL[3],
            )
            return f"{self._fname}.{TBI}"
        # precomputed vcf index index must be a non empty file
        if not (os.path.isfile(vcfidx) and os.stat(vcfidx).st_size > 0):
            exception_handler(
                CrisprHawkVCFError,
                f"Not existing or empty VCF index {vcfidx}",
                os.EX_DATAERR,
                self._debug,
            )

    def fetch(self, start: int, stop: int) -> List[List[str]]:
        try:  # extract variants in the input range from vcf file
            # return a list of variants data as list of fields
            return [
                v.strip().split() for v in self._vcf.fetch(self._contig, start, stop)
            ]
        except ValueError as e:
            exception_handler(
                CrisprHawkVCFError,
                f"Variants extraction failed for coordinates {self.contig}:{start}-{stop}",
                os.EX_DATAERR,
                self._debug,
                e,
            )

    @property
    def contig(self) -> str:
        if self._contig.startswith("chr"):
            return self._contig
        return f"chr{self._contig}"


def _find_tbi(vcf: str) -> bool:
    # avoid unexpected crashes due to file location
    vcfindex = f"{os.path.abspath(vcf)}.{TBI}"
    if os.path.exists(vcfindex):  # index must be a non empty file
        return os.path.isfile(vcfindex) and os.stat(vcfindex).st_size > 0
    return False
