"""
"""

from crisprhawk_error import CrisprHawkVCFError
from exception_handlers import exception_handler
from utils import print_verbosity, VERBOSITYLVL

from typing import Optional, List, Tuple, Set
from pysam import TabixFile, tabix_index

import os

TBI = "tbi"

class VariantRecord:
    def __init__(self, variant: List[str], samples: List[str], phased: bool, debug: bool) -> None:
        self._debug = debug  # set debug mode flag
        self._chrom = variant[0]  # store chromosome
        self._position = int(variant[1])  # store variant position
        self._ref = variant[3]  # store reference allele
        self._alt = variant[4]  # store alternative allele
        self._vid = self._assign_id(variant[2], self._chrom, self._ref, self._alt, self._position)  # assign variant id
        self._samples = _genotypes_to_samples(variant[9:], samples, phased)

    def _assign_id(self, vid: str, chrom: str, ref: str, alt: str, pos: int) -> str:
        if vid != ".":  # variant id available, return it
            return vid
        # variant id not available, construct the id using chrom, position, ref,
        # and alt (e.g. chrx_100_A_G)
        return f"{chrom}_{pos}_{ref}_{alt}"
    
    
def _genotypes_to_samples(genotypes: List[str], samples: List[str], phased: bool) -> Tuple[Set[str], Set[str]]:
    # define two sets storing samples with variant occurrence on left and right
    # copy respectively
    # if unphased vcf is used only the left set
    sampleshap = (set(), set())  # set avoid repeated values
    for i, gt in enumerate(genotypes):
        if gt[0] != "0":  # left copy
            sampleshap[0].add(samples[i])
        if gt[2] != "0" and phased:
            sampleshap[1].add(samples[i])
        elif gt[2] != "0":
            sampleshap[0].add(samples[i])
    return sampleshap


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
        self._phased = False  # by default treat vcf as unphased

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

    def _is_phased(self) -> None:
        assert hasattr(self, "_vcf")  # otherwise we couldn't establish phasing
        for variant in self._vcf.fetch():  # fecth only the first variant
            gt = variant.strip().split()[9]
            break  # no further iterations required
        # establish from genotype whther the vcf is phased or not
        if "|" in gt:
            self._phased = True


    def fetch(self, start: int, stop: int) -> List[VariantRecord]:
        try:  # extract variants in the input range from vcf file
            self._is_phased()  # assess whether the vcf is phased
            variants = [
                VariantRecord(v.strip().split(), self._samples, self._phased, self._debug)
                for v in self._vcf.fetch(self._contig, start, stop)
            ]
            # variants = [
            #     v.strip().split() for v in self._vcf.fetch(self._contig, start, stop)
            # ]  # return a list of variants data as list of fields
            return variants
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
    
    @property
    def phased(self) -> bool: return self._phased


def _find_tbi(vcf: str) -> bool:
    # avoid unexpected crashes due to file location
    vcfindex = f"{os.path.abspath(vcf)}.{TBI}"
    if os.path.exists(vcfindex):  # index must be a non empty file
        return os.path.isfile(vcfindex) and os.stat(vcfindex).st_size > 0
    return False


