"""Module for handling variant records and VCF file operations.

This module provides classes and functions for parsing, representing, and fetching
genomic variant data from VCF files. It includes utilities for variant type
assignment, ID computation, genotype parsing, and Tabix indexing.
"""

from .exception_handlers import exception_handler
from .coordinate import Coordinate
from .utils import warning

from typing import Optional, List, Tuple, Set
from pysam import TabixFile, tabix_index

import numpy as np

import sys
import os


TBI = "tbi"  # tabix index suffix
VTYPES = ["snp", "indel"]  # variant types


class VariantRecord:
    """Represents a single variant record from a VCF file.

    This class provides methods for parsing, representing, and manipulating variant
    records, including support for multi-allelic variants and variant type assignment.

    Attributes:
        _debug (bool): Debug mode flag for exception handling.
        _chrom (str): Chromosome name.
        _position (int): Variant position.
        _ref (str): Reference allele.
        _alt (List[str]): List of alternative alleles.
        _allelesnum (int): Number of alternative alleles.
        _vtype (List[str]): List of variant types for each alternative allele.
        _filter (str): Filter status from the VCF.
        _vid (List[str]): List of variant IDs.
        _samples (List[Tuple[Set[str], Set[str]]]): Sample genotype information.
    """

    def __init__(self, debug: bool) -> None:
        """Initialize a VariantRecord object.

        Initializes an empty VariantRecord object. Attributes for the variant
        record will be populated later using the read_vcf_line() method.
        """
        self._debug = debug  # store debug flag

    def __repr__(self) -> str:
        """Return a string representation of the variant record.

        Returns a string representation of the VariantRecord object, including the
        class name, chromosome, position, reference allele, and alternative alleles.

        Returns:
            The string representation of the variant record.
        """
        altalleles = ",".join(self._alt)
        return f'<{self.__class__.__name__} object; variant="{self._chrom} {self._position} {self._ref} {altalleles}">'

    def __str__(self) -> str:
        """Return a string representation of the variant record.

        Returns a string representation of the VariantRecord object, suitable for
        writing to a VCF file. The string includes the chromosome, position,
        reference allele, and alternative alleles, separated by tabs.

        Returns:
            The string representation of the variant record.
        """
        altalleles = ",".join(self._alt)
        return f"{self._chrom}\t{self._position}\t{self._ref}\t{altalleles}"

    def __eq__(self, vrecord: object) -> bool:
        """Check if two VariantRecord objects are equal.

        Two VariantRecord objects are considered equal if they have the same
        chromosome, position, reference allele, and alternative alleles.

        Args:
            vrecord: The other VariantRecord object to compare to.

        Returns:
            True if the two VariantRecord objects are equal, False otherwise.

        Raises:
            AttributeError: If the comparison fails due to missing attributes.
        """
        if not isinstance(vrecord, VariantRecord):
            return NotImplemented
        if not hasattr(vrecord, "_chrom"):  # always trace this error
            raise AttributeError(
                f"Comparison between {self.__class__.__name__} object failed"
            )
        if not hasattr(vrecord, "_position"):  # always trace this error
            raise AttributeError(
                f"Comparison between {self.__class__.__name__} object failed"
            )
        if not hasattr(vrecord, "_ref"):  # always trace this error
            raise AttributeError(
                f"Comparison between {self.__class__.__name__} object failed"
            )
        if not hasattr(vrecord, "_alt"):  # always trace this error
            raise AttributeError(
                f"Comparison between {self.__class__.__name__} object failed"
            )
        return (
            self._chrom == vrecord.contig
            and self._position == vrecord.position
            and self._ref == vrecord.ref
            and self._alt == vrecord.alt
        )

    def __lt__(self, vrecord: "VariantRecord") -> bool:
        """Compare two VariantRecord objects based on position.

        Compares the current VariantRecord object with another VariantRecord object
        based on their positions.

        Args:
            vrecord: The other VariantRecord object to compare to.

        Returns:
            True if the current object's position is less than the other object's
            position, False otherwise.
        """
        return self._position < vrecord.position

    def __gt__(self, vrecord: "VariantRecord") -> bool:
        """Compare two VariantRecord objects based on position.

        Compares the current VariantRecord object with another VariantRecord object
        based on their positions.

        Args:
            vrecord: The other VariantRecord object to compare to.

        Returns:
            True if the current object's position is greater than the other object's
            position, False otherwise.
        """
        return self._position > vrecord.position

    def __hash__(self) -> int:
        """Return a hash value for the variant record.

        Returns a hash value for the VariantRecord object, based on its chromosome,
        position, reference allele, and alternative alleles. This allows
        VariantRecord objects to be used as keys in dictionaries or sets.

        Returns:
            The hash value of the variant record.
        """
        return hash((self._chrom, self._position, self._ref, tuple(self._alt)))

    def _retrieve_alt_alleles(self, altalleles: str) -> List[str]:
        """Retrieve and parse alternative alleles from a string.

        Retrieves alternative alleles from a comma-separated string.

        Args:
            altalleles: A string containing comma-separated alternative alleles.

        Returns:
            A list of alternative alleles.

        Raises:
            AttributeError: If the input is not a string.
            Exception: For any other unexpected error during parsing.
        """
        if not altalleles:
            return []
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
        """Determine the variant type for all alternative alleles.

        Determines if each alternative allele represents a SNP or an indel based on the
        lengths of the reference and alternative alleles.

        Returns:
            A list of variant types (either "snp" or "indel") corresponding to each
            alternative allele.
        """
        assert hasattr(self, "_ref")
        assert hasattr(self, "_alt")
        return [_assign_vtype(self._ref, altallele) for altallele in self._alt]

    def _retrieve_af(self, info: str) -> List[float]:
        i = info.find("AF=")  # find the AF field start index
        if i == -1:  # no AF data in the input VCF
            return [np.nan] * self._allelesnum
        i += 3  # skip 'AF=' in info
        j = info.find(";", i)  # find the next semicolon delimiter
        j = len(info) if j == -1 else j  # patch for no multiple info fields
        afs = list(map(float, info[i:j].split(",")))
        if len(afs) != self._allelesnum:  # one af per af value per allele
            exception_handler(
                ValueError,
                f"AF number does not match the alleles number ({len(afs)} - {self._allelesnum})",
                os.EX_DATAERR,
                self._debug,
            )
        return afs

    def _assign_id(self) -> List[str]:
        """Assign or compute variant IDs.

        Assigns or computes variant IDs based on whether the variant is monoallelic or
        multiallelic. For monoallelic variants, it uses the provided ID or computes one
        if not available. For multiallelic variants, it computes an ID for each
        alternative allele.

        Args:
            vid: The variant ID (may be empty for monoallelic variants).

        Returns:
            A list of variant IDs.
        """
        if self._allelesnum == 1:
            # variant id not available, construct the id using chrom, position, ref,
            # and alt (e.g. chrx-100-A/G)
            return [_compute_id(self._chrom, self._position, self._ref, self._alt[0])]
        # if multiallelic site compute the id for each alternative allele
        # avoid potential confusion due to alternative alleles at same position
        # labeled with the same id
        return [
            _compute_id(self._chrom, self._position, self._ref, altallele)
            for altallele in self._alt
        ]

    def _copy(self, i: int) -> "VariantRecord":
        """Create a copy of the variant record for a specific allele.

        Creates a copy of the current VariantRecord object, representing a single
        alternative allele specified by the index `i`.  This is useful for
        handling multiallelic sites where each alternative allele needs to be
        treated as a separate variant.

        Args:
            i: The index of the alternative allele to copy.

        Returns:
            A new VariantRecord object representing the specified alternative allele.
        """
        # sourcery skip: class-extract-method
        # copy current variant record instance
        vrecord = VariantRecord(self._debug)  # create new instance
        # adjust ref/alt alleles and positions for multiallelic sites
        ref, alt, position = adjust_multiallelic(
            self._ref, self._alt[i], self._position
        )
        vrecord._chrom = self._chrom
        vrecord._position = position
        vrecord._ref = ref
        vrecord._alt = [alt]
        vrecord._allelesnum = 1
        vrecord._vtype = [self._vtype[i]]
        vrecord._filter = self._filter
        vrecord._afs = [self._afs[i]]
        vrecord._vid = [self._vid[i]]
        vrecord._samples = [self._samples[i]]
        return vrecord

    def read_vcf_line(
        self, variant: List[str], samples: List[str], phased: bool
    ) -> None:
        """Read and parse a VCF line.

        Parses a line from a VCF file and populates the attributes of the
        VariantRecord object with the extracted information.

        Args:
            variant: A list of strings representing the fields of the VCF line.
            samples: A list of sample names.
            phased: True if the genotypes are phased, False otherwise.
        """

        self._chrom = variant[0]  # store chromosome
        self._position = int(variant[1])  # store variant position
        self._ref = variant[3]  # store ref allele
        self._alt = self._retrieve_alt_alleles(variant[4])  # store alt alleles
        self._allelesnum = len(self._alt)  # number of alt alleles
        self._vtype = self._assess_vtype()  # establish whether is a snp or indel
        self._filter = variant[6]  # store filter value
        self._afs = self._retrieve_af(variant[7])  # retrieve allele frequencies
        self._vid = self._assign_id()  # assign variant id
        self._samples = _genotypes_to_samples(
            variant[9:], samples, self._allelesnum, phased, self._debug
        )  # recover samples with their genotypes

    def split(self, vtype: Optional[str] = None) -> List["VariantRecord"]:
        """Split a multiallelic variant record by variant type.

        Splits a multiallelic VariantRecord object into a list of VariantRecord objects,
        each representing a single alternative allele of the specified variant type.

        Args:
            vtype: The variant type to select ("snp" or "indel"). If None, all
                variant types are included.

        Returns:
            A list of VariantRecord objects, one for each alternative allele matching
            the specified variant type.
        """
        vtypes_filter = VTYPES if vtype is None else [vtype]
        return [
            self._copy(i)
            for i, _ in enumerate(self._vtype)
            if self._vtype[i] in vtypes_filter
        ]

    def get_altalleles(self, vtype: str) -> List[str]:
        """Retrieve alternative alleles of a specific variant type.

        Returns a list of alternative alleles that match the specified variant type
        (either "snp" or "indel").

        Args:
            vtype: The variant type to select ("snp" or "indel").

        Returns:
            A list of alternative alleles matching the specified type.
        """
        assert vtype in VTYPES
        # return the alternative alleles representing snps or indels
        return [
            altallele
            for i, altallele in enumerate(self._alt)
            if self._vtype[i] == vtype
        ]

    def pytest_initialize(
        self, position: int, ref: str, alt: str, vtype: str, vid: str, afs: List[float]
    ) -> None:
        """Initialize the VariantRecord for pytest with provided values.

        Sets the chromosome, position, reference allele, alternative allele,
        variant type, variant ID, and allele frequencies for testing purposes.

        Args:
            position (int): The variant position.
            ref (str): The reference allele.
            alt (str): The alternative allele.
            vtype (str): The variant type.
            vid (str): The variant ID.
            afs (List[float]): The allele frequencies.
        """
        self._chrom = "chrx"
        self._position = position
        self._ref = ref
        self._alt = [alt]
        self._vtype = [vtype]
        self._vid = [vid]
        self._afs = afs

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
    def afs(self) -> List[float]:
        return self._afs

    @property
    def samples(self) -> List[Tuple[Set[str], Set[str]]]:
        return self._samples

    @property
    def id(self) -> List[str]:
        return self._vid

    @property
    def allelesnum(self) -> int:
        return self._allelesnum


def _assign_vtype(ref: str, alt: str) -> str:
    """Determine the variant type.

    Determines if a variant is an indel or a SNP based on the lengths of the
    reference and alternative alleles.

    Args:
        ref: The reference allele.
        alt: The alternative allele.

    Returns:
        "indel" if the lengths of the reference and alternative alleles are not equal,
        "snp" otherwise.
    """
    return VTYPES[1] if len(ref) != len(alt) else VTYPES[0]


def _compute_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Compute a variant ID.

    Computes a variant ID using the chromosome, position, reference allele, and
    alternative allele, following the IGVF consortium notation.

    Args:
        chrom: The chromosome.
        pos: The position.
        ref: The reference allele.
        alt: The alternative allele.

    Returns:
        The computed variant ID.
    """
    # compute variant id for variants without id, or multiallelic sites
    # use IGVF consortium notation
    return f"{chrom}-{pos}-{ref}/{alt}"


def adjust_multiallelic(ref: str, alt: str, pos: int) -> Tuple[str, str, int]:
    """Adjust reference/alternative alleles and position for multiallelic sites.

    Adjusts the reference and alternative alleles, and the variant position for
    multiallelic sites based on the lengths of the original reference and
    alternative alleles.  This function helps normalize variant representation
    for easier comparison and processing. The function assumes multiallelic
    variants are left-aligned.

    Args:
        ref: The original reference allele.
        alt: The original alternative allele.
        pos: The original variant position.

    Returns:
        A tuple containing the adjusted reference allele, alternative allele, and
        variant position.
    """

    if len(ref) == len(alt):  # likely snp
        ref_new, alt_new = ref[0], alt[0]  # adjust ref/alt alleles
        pos_new = pos  # ref/alt have same length
    elif len(ref) > len(alt):  # deletion
        ref_new = ref[len(alt) - 1 :]  # adjust ref allele
        alt_new = alt[-1]  # adjust alt allele
        pos_new = pos + (len(alt)) - 1  # adjust variant position
    else:  # insertion
        ref_new = ref[-1]  # adjust ref allele
        alt_new = alt[len(ref) - 1 :]  # adjust alt allele
        pos_new = pos + len(ref) - 1  # adjust variant position
    return ref_new, alt_new, pos_new


def _parse_genotype_phased(
    gt_alleles: List[str],
    sample: str,
    sampleshap: List[Tuple[Set[str], Set[str]]],
    debug: bool,
) -> List[Tuple[Set[str], Set[str]]]:
    """Parse phased genotype alleles and update sample sets.

    Updates the sample sets for each alternative allele based on the phased
    genotype information. Assigns the sample to the left or right copy depending
    on the genotype.

    Args:
        gt_alleles: A list containing the alleles for the two haplotypes.
        sample: The sample name.
        sampleshap: A list of tuples of sets, tracking samples for each allele
            and haplotype.
        debug: Whether to enable debug mode for exception handling.

    Returns:
        The updated list of tuples of sets with sample assignments.

    Raises:
        ValueError: If the genotype does not contain exactly two alleles.
    """
    if len(gt_alleles) != 2:
        exception_handler(
            ValueError,
            "Phased genotypes cannot have more than one allele on each copy",
            os.EX_DATAERR,
            debug,
        )
    gt1, gt2 = gt_alleles  # retrieve allele occurring on first and second copy
    if gt1 not in ["0", "."]:  # left copy
        sampleshap[int(gt1) - 1][0].add(sample)
    if gt2 not in ["0", "."]:  # right copy
        sampleshap[int(gt2) - 1][1].add(sample)
    return sampleshap


def _parse_genotype_unphased(
    gt_alleles: List[str],
    sample: str,
    sampleshap: List[Tuple[Set[str], Set[str]]],
) -> List[Tuple[Set[str], Set[str]]]:
    """Parse unphased genotype alleles and update sample sets.

    Updates the sample sets for each alternative allele based on the unphased
    genotype information. Assigns the sample to the appropriate set for each
    allele present in the genotype.

    Args:
        gt_alleles: A list containing the alleles for the two haplotypes or more.
        sample: The sample name.
        sampleshap: A list of tuples of sets, tracking samples for each allele
            and haplotype.

    Returns:
        The updated list of tuples of sets with sample assignments.
    """
    if len(gt_alleles) != 2:  # handle genotypes like 0/1/2
        for gt in gt_alleles:
            if gt not in ["0", "."]:
                sampleshap[int(gt) - 1][0].add(sample)
    else:  # handle genotypes like 0/1
        gt1, gt2 = gt_alleles  # retrieve allele occurring on first and second copy
        if gt1 not in ["0", "."] and gt1 == gt2:  # special case 1/1
            sampleshap[int(gt1) - 1][0].add(sample)
            sampleshap[int(gt2) - 1][1].add(sample)
        else:
            if gt1 not in ["0", "."]:  # 1/0
                sampleshap[int(gt1) - 1][0].add(sample)
            if gt2 not in ["0", "."]:  # 0/1
                sampleshap[int(gt2) - 1][0].add(sample)
    return sampleshap


def _genotypes_to_samples(
    genotypes: List[str], samples: List[str], allelesnum: int, phased: bool, debug: bool
) -> List[Tuple[Set[str], Set[str]]]:
    """Extract sample information from genotypes.

    Parses genotype strings to determine which samples carry each alternative allele.
    Handles both phased and unphased genotypes.

    Args:
        genotypes: A list of genotype strings.
        samples: A list of sample names.
        allelesnum: The number of alternative alleles.
        phased: True if the genotypes are phased, False otherwise.

    Returns:
        A list of tuples containing two lists of sets. The first list contains sets of
        samples with the variant on the left copy (or the only copy if unphased), and
        the second list contains sets of samples with the variant on the right copy
        (only relevant for phased data).

    Raises:
        TypeError: If the genotype string cannot be split.
        Exception: If an unexpected error occurs during genotype parsing.
    """
    # define two sets storing samples with variant occurrence on left and right
    # copy respectively
    # if unphased vcf, is used only the left set
    sampleshap = [(set(), set()) for _ in range(allelesnum)]
    gtsep = "|" if phased else "/"  # define genotype separator char
    for i, gt in enumerate(genotypes):
        try:
            gt_alleles = gt.split(":")[0].split(gtsep)
        except TypeError as e:
            exception_handler(
                TypeError,
                f"Split object is not of {str.__name__} type",
                os.EX_DATAERR,
                debug,
                e,
            )
        except Exception as e:
            exception_handler(
                Exception,
                f"An unexpected error occurred while parsing genotype {gt}",
                os.EX_DATAERR,
                debug,
                e,
            )
        sampleshap = (
            _parse_genotype_phased(gt_alleles, samples[i], sampleshap, debug)
            if phased
            else _parse_genotype_unphased(gt_alleles, samples[i], sampleshap)
        )
    return sampleshap


class VCF:
    """Represents a VCF file for variant fetching and indexing.

    Provides methods to open, index, and fetch variants from a VCF file,
    including support for phased and unphased data. Handles Tabix indexing and
    ensures compatibility with single-contig VCFs.

    Attributes:
        _debug (bool): Debug mode flag for exception handling.
        _verbosity (int): Verbosity level for warnings and errors.
        _fname (str): Path to the VCF file.
        _vcfidx (str): Path to the Tabix index file.
        _vcf (TabixFile): TabixFile object for VCF access.
        _contig (str): Contig name in the VCF file.
        _samples (List[str]): List of sample names in the VCF.
        _phased (bool): Whether the VCF genotypes are phased.
    """

    def __init__(
        self, fname: str, verbosity: int, debug: bool, vcfidx: Optional[str] = ""
    ) -> None:
        """Initialize a VCF object for reading and indexing VCF files.

        Opens a VCF file, checks for or creates a Tabix index, and prepares the
        file for variant fetching. Also determines if the VCF is phased and
        extracts sample and contig information.

        Args:
            fname: The path to the VCF file.
            verbosity: The verbosity level for warnings and errors.
            debug: Whether to enable debug mode for exception handling.
            vcfidx: Optional path to a precomputed Tabix index file.
        """
        self._debug = debug  # store debug flag
        self._verbosity = verbosity  # store verbosity level
        if not os.path.isfile(fname):
            exception_handler(
                FileNotFoundError,
                f"Cannot find input VCF {fname}",
                os.EX_DATAERR,
                self._debug,
            )
        self._fname = fname  # store input filename
        self._vcfidx = self._search_index(vcfidx)  # initialize vcf index
        if not self._vcfidx:  # index not found compute it
            self.index_vcf()
        # initialize TabixFile object with the previously computed index
        self._vcf = TabixFile(self._fname, index=self._vcfidx)
        if len(set(self._vcf.contigs)) != 1:  # assume vcf data about one contig
            exception_handler(
                ValueError,
                f"Input VCF {fname} store variants belonging to multiple contigs",
                os.EX_DATAERR,
                self._debug,
            )
        self._contig = self._vcf.contigs[0]
        self._samples = self._vcf.header[-1].strip().split()[9:]  # recover samples
        self._phased = False  # by default treat vcf as unphased
        self._is_phased()  # check if the input vcf is phased

    def _search_index(self, vcfidx: Optional[str] = "") -> str:
        """Search for or validate a Tabix index for the VCF file.

        Searches for a Tabix index (.tbi) for the associated VCF file if one is not
        provided. If a path to an index is provided, it validates that the index
        exists and is not empty.

        Args:
            vcfidx: An optional path to a Tabix index file.

        Returns:
            The path to the Tabix index file, or an empty string if not found.

        Raises:
            FileNotFoundError: If the provided index file does not exist or is empty.
        """

        # look for index for the current vcf, if not found compute it
        if not vcfidx:
            if find_tbi(self._fname):  # index found, store it
                return f"{self._fname}.{TBI}"
            # index not found -> compute it de novo and store it in the same folder
            # as the input vcf
            sys.stdout.write(f"Tabix index not found for {self._fname}\n")
            return ""
        # precomputed vcf index index must be a non empty file
        if not (os.path.isfile(vcfidx) and os.stat(vcfidx).st_size > 0):
            raise FileNotFoundError(f"Not existing or empty VCF index {vcfidx}")
        return vcfidx

    def index_vcf(self, pytest: Optional[bool] = False) -> None:
        """Create or update the Tabix index for the VCF file.

        Creates or updates the Tabix index (.tbi) for the associated VCF file.
        If an index already exists, it will be overwritten.

        Raises:
            OSError: If an error occurs during indexing.
            RuntimeWarning: If an index already exists.
        """
        if self._vcfidx and not pytest:  # launch warning
            warning("Tabix index already present, forcing update", self._verbosity)
        # compute tabix index if not provided during object initialization
        try:  # create index in the same folder as the input vcf
            tabix_index(self._fname, preset="vcf", force=True)
        except OSError as e:
            exception_handler(
                OSError,
                f"An error occurred while indexing {self._fname}",
                os.EX_DATAERR,
                self._debug,
                e,
            )
        assert find_tbi(self._fname)
        self._vcfidx = f"{self._fname}.{TBI}"

    def _is_phased(self) -> None:
        """Check if the VCF file is phased.

        Checks if the VCF file is phased by examining the genotype of the first variant.
        If the genotype contains a pipe character '|', the VCF is considered phased.
        The result is stored in the _phased attribute.
        """
        assert hasattr(self, "_vcf")  # otherwise we couldn't establish phasing
        for variant in self._vcf.fetch():  # fecth only the first variant
            gt = variant.strip().split()[9]
            # establish from genotype whther the vcf is phased or not
            if "|" in gt:
                self._phased = True
            break  # no further iterations required

    def fetch(self, coordinate: Coordinate) -> List[VariantRecord]:
        """Fetch variants within a specified genomic interval.

        Fetches variants from the VCF file that overlap with the specified genomic
        interval. Handles potential mismatches between chromosome prefixes and checks for
        phasing.

        Args:
            coordinate: A Coordinate object specifying the genomic interval.

        Returns:
            A list of VariantRecord objects representing the variants within the
            specified interval.

        Raises:
            ValueError: If the VCF and coordinate contigs mismatch or if an invalid
                reference or position is provided.
            IndexError: If an attempt is made to fetch data outside the bounds of the
                indexed file.
            Exception: For any other unexpected error during fetching.
        """

        if self._contig != coordinate.contig:
            # may be just caused by a prefix in contig
            vcfcontig = self._contig.replace("chr", "")
            coordcontig = coordinate.contig.replace("chr", "")
            if vcfcontig != coordcontig:
                exception_handler(
                    ValueError,
                    f"Mismatching VCF and coordinate contigs ({self._contig} - {coordinate.contig})",
                    os.EX_DATAERR,
                    self._debug,
                )
        try:  # extract variants in the input range from vcf file
            return [
                _create_variant_record(
                    v.strip().split(), self._samples, self._phased, self._debug
                )
                for v in self._vcf.fetch(
                    self._contig, coordinate.start, coordinate.stop
                )
            ]
        except ValueError as e:
            exception_handler(
                ValueError,
                f"Invalid reference or position provided ({self._contig}\t{coordinate.start}\t{coordinate.stop})",
                os.EX_DATAERR,
                self._debug,
                e,
            )
        except IndexError as e:
            exception_handler(
                IndexError,
                f"Tried to fetch data outside the bounds of the indexed file ({self._contig}\t{coordinate.start}\t{coordinate.stop})",
                os.EX_DATAERR,
                self._debug,
                e,
            )
        except Exception as e:
            exception_handler(
                Exception,
                f"An unexpected error occurred ({self._contig}\t{coordinate.start}\t{coordinate.stop})",
                os.EX_DATAERR,
                self._debug,
                e,
            )

    @property
    def contig(self) -> str:
        return self._contig if self._contig.startswith("chr") else f"chr{self._contig}"

    @property
    def phased(self) -> bool:
        return self._phased

    @property
    def samples(self) -> List[str]:
        return self._samples


def find_tbi(vcf: str) -> bool:
    """Check if a Tabix index exists for a VCF file.

    Checks if a Tabix index (.tbi) exists for the given VCF file and is a non-empty
    file.

    Args:
        vcf: The path to the VCF file.

    Returns:
        True if the index exists and is a non-empty file, False otherwise.
    """
    # avoid unexpected crashes due to file location
    vcfindex = f"{os.path.abspath(vcf)}.{TBI}"
    if os.path.exists(vcfindex):  # index must be a non empty file
        return os.path.isfile(vcfindex) and os.stat(vcfindex).st_size > 0
    return False


def _create_variant_record(
    variant: List[str], samples: List[str], phased: bool, debug: bool
) -> VariantRecord:
    """Create and populate a VariantRecord from VCF line data.

    Instantiates a VariantRecord object and populates it with information parsed from a VCF line. Returns the fully populated VariantRecord.

    Args:
        variant: A list of strings representing the fields of the VCF line.
        samples: A list of sample names.
        phased: True if the genotypes are phased, False otherwise.
        debug: Whether to enable debug mode for exception handling.

    Returns:
        A populated VariantRecord object.
    """
    vrecord = VariantRecord(debug)  # create variant record instance
    vrecord.read_vcf_line(variant, samples, phased)  # read vcf line
    return vrecord
