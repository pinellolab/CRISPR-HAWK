"""
"""

from crisprhawk_error import CrisprHawkVCFError
from exception_handlers import exception_handler
from utils import print_verbosity, VERBOSITYLVL

from typing import Optional, List, Tuple, Set
from pysam import TabixFile, tabix_index

import os

TBI = "tbi"  # tabix index file extnsion
VTYPES = ["snp", "indel"]  # variant types


class VariantRecord:
    def __init__(
        self, variant: List[str], samples: List[str], phased: bool, debug: bool
    ) -> None:
        self._debug = debug  # set debug mode flag
        self._chrom = variant[0]  # store chromosome
        self._position = int(variant[1])  # store variant position
        self._ref = variant[3]  # store ref allele
        self._alt = self._retrieve_alt_alleles(variant[4])  # store alt alleles
        self._allelesnum = len(self._alt)  # number of alt alleles
        self._vtype = self._assess_vtype()  # establish whether is a snp or indel
        self._filter = variant[6]  # store filter value
        self._vid = self._assign_id(variant[2])  # assign variant id
        self._samples = _genotypes_to_samples(
            variant[9:], samples, self._allelesnum, phased, self._debug
        )

    def __repr__(self) -> str:
        altalleles = ",".join(self._alt)
        return f"<{self.__class__.__name__} object; variant=\"{self._chrom} {self._position} {self._ref} {altalleles}\">"

    def _retrieve_alt_alleles(self, altalleles: str) -> List[str]:
        # alternative alleles in multiallelic sites are separated by a comma
        try:
            return altalleles.split(",")
        except AttributeError as e:
            exception_handler(
                AttributeError,
                "Alternative alleles must be encoded in a string",
                os.EX_DATAERR,
                self._debug,
                e,
            )
        except Exception as e:
            exception_handler(
                Exception,
                "Unexpected error while parsing alternative alleles",
                os.EX_DATAERR,
                self._debug,
                e,
            )

    def _assess_vtype(self) -> List[str]:
        assert hasattr(self, "_ref")
        assert hasattr(self, "_alt")
        return [_assign_vtype(self._ref, altallele) for altallele in self._alt]

    def _assign_id(self, vid: str) -> List[str]:
        if self._allelesnum == 1:
            # for indels always compute the identifier to avoid confusion
            if (
                vid != "." and self.vtype[0] != VTYPES[1]
            ):  # variant id available, return it
                return [vid]
            # variant id not available, construct the id using chrom, position, ref,
            # and alt (e.g. chrx_100_A_G)
            return [_compute_id(self._chrom, self._position, self._ref, self._alt[0])]
        # if multiallelic site compute the id for each alternative allele
        # avoid potential confusion due to alternative alleles at same position
        # labeled with the same id
        return [
            _compute_id(self._chrom, self._position, self._ref, altallele)
            for altallele in self._alt
        ]

    def format(self) -> str:
        altalleles = ",".join(self.alt)
        return f"{self._chrom}\t{self._position}\t{self._ref}\t{altalleles}"

    def get_altalleles(self, vtype: str) -> List[str]:
        assert vtype in VTYPES
        # return the alternative alleles representing snps or indels
        return [
            altallele
            for i, altallele in enumerate(self._alt)
            if self._vtype[i] == vtype
        ]

    @property
    def filter(self) -> str:
        return self._filter

    @property
    def contig(self) -> str:
        return self._chrom

    @property
    def position(self) -> int:
        return self._position

    @property
    def ref(self) -> str:
        return self._ref

    @property
    def alt(self) -> List[str]:
        return self._alt

    @property
    def vtype(self) -> List[str]:
        return self._vtype

    @property
    def samples(self) -> Tuple[Set[str], Set[str]]:
        return self._samples

    @property
    def id(self) -> List[str]:
        return self._vid


def _assign_vtype(ref: str, alt: str) -> bool:
    if len(ref) != 1 or len(alt) != 1:
        return VTYPES[1]  # indel
    return VTYPES[0]  # snp


def _compute_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    # compute variant id for variants without id, or multiallelic sites
    return f"{chrom}_{pos}_{ref}_{alt}"


def _genotypes_to_samples(
    genotypes: List[str], samples: List[str], allelesnum: int, phased: bool, debug: bool
) -> Tuple[Set[str], Set[str]]:
    # define two sets storing samples with variant occurrence on left and right
    # copy respectively
    # if unphased vcf is used only the left set
    sampleshap = [
        (set(), set()) for _ in range(allelesnum)
    ]  # pairs for each alt allele
    gtsep = "|" if phased else "/"  # define genotype separator char
    for i, gt in enumerate(genotypes):
        try:
            gt1, gt2 = gt.split(":")[0].split(gtsep)
        except ValueError as e:
            exception_handler(
                ValueError, f"Non diploid genotype detected", os.EX_DATAERR, debug, e
            )
        if gt1 != "0" and gt1 != ".":  # left copy
            sampleshap[int(gt1) - 1][0].add(samples[i])
        if phased and (gt2 != "0" and gt2 != "."):  # right copy
            sampleshap[int(gt2) - 1][1].add(samples[i])
        elif gt2 != "0" and gt2 != ".":
            sampleshap[int(gt2) - 1][0].add(samples[i])
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
                VariantRecord(
                    v.strip().split(), self._samples, self._phased, self._debug
                )
                for v in self._vcf.fetch(self._contig, start, stop)
            ]
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
    def phased(self) -> bool:
        return self._phased


def _find_tbi(vcf: str) -> bool:
    # avoid unexpected crashes due to file location
    vcfindex = f"{os.path.abspath(vcf)}.{TBI}"
    if os.path.exists(vcfindex):  # index must be a non empty file
        return os.path.isfile(vcfindex) and os.stat(vcfindex).st_size > 0
    return False
