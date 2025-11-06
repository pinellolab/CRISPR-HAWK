"""
This module provides classes and utility functions for representing and manipulating
haplotypes.

It includes the Haplotype and HaplotypeIndel classes for managing genomic regions
with sequence variants, as well as helper functions for variant sorting, chain
computation, and IUPAC encoding.
"""

from .crisprhawk_error import CrisprHawkHaplotypeError, CrisprHawkIupacTableError
from .exception_handlers import exception_handler
from .region import Region
from .sequence import Sequence
from .coordinate import Coordinate
from .variant import VariantRecord, VTYPES
from .utils import match_iupac, IUPAC_ENCODER, IUPACTABLE

from typing import List, Dict, Tuple

import os


class Haplotype(Region):
    """Represents a haplotype, which is a specific combination of alleles or sequence
    variants in a genomic region.

    Provides methods to add, modify, and manage variants, phasing, and sample
    information for a haplotype.

    Attributes:
        _debug (bool): Debug flag for verbose output.
        _size (int): Size of the haplotype sequence.
        _variant_alleles (dict): Maps variant positions to their alleles.
        _variants (str): String representation of haplotype variants.
        _afs (dict): Maps variant IDs to their allele frequencies.
        _samples (str): Sample information for the haplotype.
        _phased (bool): Indicates if the haplotype is phased.
        _chromcopy (int): Chromosome copy number.
        _posmap (dict): Maps relative positions to genomic positions.
        _posmap_reverse (dict): Maps genomic positions to relative positions.
        _id (str): Haplotype identifier.
        _sequence (Sequence): The sequence object representing the haplotype.
        _coordinates (Coordinate): The coordinate object for the haplotype region.
    """

    def __init__(
        self,
        sequence: Sequence,
        coord: Coordinate,
        phased: bool,
        chromcopy: int,
        debug: bool,
    ) -> None:
        """Initializes a Haplotype object with sequence, coordinates, phasing,
        chromosome copy, and debug flag.

        Sets up the haplotype's sequence, region, and internal state for variant
        and sample management.

        Args:
            sequence (Sequence): The sequence object representing the haplotype.
            coord (Coordinate): The coordinate object for the haplotype region.
            phased (bool): Indicates if the haplotype is phased.
            chromcopy (int): Chromosome copy number.
            debug (bool): Debug flag for verbose output.
        """
        self._debug = debug  # store debug flag
        super().__init__(sequence, coord)  # genomic sequence and region coordinates
        self._size = len(sequence)  # haplotype size
        self._variant_alleles = {}  # map alleles to each variant
        self._variants = "NA"  # haplotype variants
        self._afs = {}  # map allele frequency to each variant
        self._samples = "REF"  # haplotype samples
        self._phased = phased  # haplotype phasing
        self._chromcopy = chromcopy  # chromosome copy
        self._initialize_posmap([], self._coordinates.start)  # initialize position map

    def __str__(self) -> str:
        """Returns a string representation of the haplotype showing sample and
        sequence.

        Provides a concise summary of the haplotype's sample information and its
        nucleotide sequence.

        Returns:
            str: A string in the format 'samples: sequence'.
        """
        return f"{self._samples}: {self._sequence.sequence}"

    def _initialize_posmap(self, chains: List[int], start: int) -> None:
        """Initializes the position mapping for the haplotype based on indel chains
        and start position.

        Updates the internal position maps to reflect the current structure of the
        haplotype sequence.

        Args:
            chains (List[int]): List of indel chain lengths.
            start (int): The starting genomic position for the haplotype.
        """
        self._posmap = {
            p: start + i for i, p in enumerate(range(self._size + sum(chains)))
        }
        self._posmap_reverse = {pos: posrel for posrel, pos in self._posmap.items()}

    def _update_sequence(self, start: int, stop: int, alt: str) -> None:
        """Updates the haplotype sequence by replacing a region with an alternate
        allele.

        Modifies the internal sequence to reflect the specified variant.

        Args:
            start (int): Start index of the region to replace.
            stop (int): Stop index (exclusive) of the region to replace.
            alt (str): Alternate allele sequence to insert.
        """
        self._sequence._sequence_raw = (
            self._sequence._sequence_raw[:start]
            + list(alt.lower())
            + self._sequence._sequence_raw[stop:]
        )

    def substring(self, start: int, stop: int) -> str:
        """Returns a substring of the haplotype sequence between the specified
        start and stop positions.

        Extracts and concatenates the nucleotide sequence from the given range.

        Args:
            start (int): Start index of the substring.
            stop (int): Stop index (exclusive) of the substring.

        Returns:
            str: The substring of the haplotype sequence.
        """
        return "".join(self._sequence._sequence_raw[start:stop])

    def _update_posmap(self, posrel: int, chain: int) -> None:
        """Updates the position mapping after a variant is inserted or deleted.

        Adjusts the internal position maps to reflect changes in the haplotype
        sequence structure.

        Args:
            posrel (int): The relative position in the haplotype sequence.
            chain (int): The length change caused by the variant (positive for
                insertions, negative for deletions).
        """
        if chain < 0:
            for pos in range(posrel + 1, max(self._posmap.keys()) + 1):
                self._posmap[pos] -= chain  # self._posmap[pos] + abs(chain)
        if chain > 0:
            for pos in range(posrel + 1, max(self._posmap.keys()) + 1):
                self._posmap[pos] = (
                    self._posmap[posrel]
                    if pos < posrel + chain + 1
                    else self._posmap[pos] - chain
                )
        self._posmap_reverse = {pos: posrel for posrel, pos in self._posmap.items()}

    def _update_variant_alleles(self, pos: int, stop: int, offset: int) -> None:
        """Updates the variant alleles mapping after an indel event in the haplotype.

        This method removes variant alleles within the affected region and shifts
        positions of alleles occurring after the indel by the specified offset.

        Args:
            pos (int): The position where the indel starts.
            stop (int): The position where the indel ends (exclusive).
            offset (int): The offset to apply to positions after the indel.

        Returns:
            None
        """
        self._variant_alleles = {
            p: alleles
            for p, alleles in self._variant_alleles.items()
            if p <= pos or p >= stop
        }
        self._variant_alleles = {
            (p + offset if p > pos else p): alleles
            for p, alleles in self._variant_alleles.items()
        }

    def _insert_variant_phased(
        self, position: int, ref: str, alt: str, chain: int, offset: int
    ) -> None:
        """Inserts a variant into the haplotype sequence for phased data.

        Updates the haplotype by applying the given variant and adjusting the
        sequence and position mapping.

        Args:
            position (int): Genomic position of the variant.
            ref (str): Reference allele sequence.
            alt (str): Alternate allele sequence.
            chain (int): Length change caused by the variant (positive for insertions,
                negative for deletions).
            offset (int): Cumulative offset from previous variants.
        """
        posrel = self._posmap_reverse[position]
        posrel_stop = posrel + abs(chain) + 1 if chain < 0 else posrel + 1
        if posrel_stop > self._size:
            posrel_stop = (self._size + offset) - 1
        refnt = self.substring(posrel, posrel_stop)
        if refnt != ref and refnt.isupper():
            raise ValueError(
                f"Mismatching reference alleles in VCF and reference sequence "
                f"at position {position} ({refnt} - {ref})"
            )
        self._update_sequence(posrel, posrel_stop, alt)
        self._update_posmap(posrel, chain)

    def add_variants_phased(self, variants: List[VariantRecord], sample: str) -> None:
        """Adds a list of phased variants to the haplotype and updates its sequence
        and sample information.

        Applies each variant to the haplotype, reconstructs the sequence, and updates
        sample, variant, and allele frequency data.

        Args:
            variants (List[VariantRecord]): List of VariantRecord objects to add
                to the haplotype.
            sample (str): Sample identifier to associate with the haplotype.
        """
        if not self._phased:
            exception_handler(
                ValueError,
                "Unphased haplotype, unable to add phased variants",
                os.EX_DATAERR,
                True,
            )  # always trace this error
        variants = _sort_variants(variants)
        chains = _compute_chains(variants)  # retrieve original sequence positions
        # position map to handle indels' presence
        self._initialize_posmap(chains, self._coordinates.start)
        for i, variant in enumerate(variants):  # add variants to haplotype
            self._insert_variant_phased(
                variant.position,
                variant.ref,
                variant.alt[0],
                chains[i],
                sum(chains[:i]),
            )
        # add chromcopy as suffix to haplotype's samples
        suffix = "1|0" if self._chromcopy == 0 else "0|1"
        self._sequence = Sequence(
            "".join(self._sequence._sequence_raw), self._debug, allow_lower_case=True
        )  # reconstruct haplotype sequence
        self._samples = f"{sample}:{suffix}" if self._phased else sample
        self._variants = ",".join([v.id[0] for v in variants])
        self._afs = {v.id[0]: v.afs[0] for v in variants}

    def _insert_variant_unphased(
        self, position: int, ref: str, alt: str, vtype: str, chain: int, offset: int
    ) -> None:
        """Inserts a variant into the haplotype sequence for unphased data.

        Applies the given variant to the haplotype, updating the sequence, variant
        alleles, and position mapping.

        Args:
            position (int): Genomic position of the variant.
            ref (str): Reference allele sequence.
            alt (str): Alternate allele sequence.
            vtype (str): Variant type (e.g., SNV or indel).
            chain (int): Length change caused by the variant (positive for insertions,
                negative for deletions).
            offset (int): Cumulative offset from previous variants.

        Raises:
            ValueError: If the reference allele does not match the haplotype sequence
                at the specified position.
        """
        try:
            posrel = self._posmap_reverse[position]
        except KeyError:  # position may be deleted by previous deletion
            return
        posrel_stop = posrel + abs(chain) + 1 if chain < 0 else posrel + 1
        if posrel_stop > self._size:
            posrel_stop = (self._size + offset) - 1
        refnt = self.substring(posrel, posrel_stop)  # retrieve ref nt in haplotype
        if not match_iupac(ref, refnt):  # check ref alleles coherence
            raise ValueError(
                f"Mismatching reference alleles in VCF and reference sequence at position {position} ({refnt} - {ref})"
            )
        if posrel in self._variant_alleles:  # to solve haplotypes
            if position == self._variant_alleles[posrel][0][2]:
                self._variant_alleles[posrel].append((ref, alt, position))
        else:
            self._variant_alleles[posrel] = [(ref, alt, position)]
        if vtype == VTYPES[0]:  # if snv encode as iupac
            alt = _encode_iupac(refnt, alt, position, self._debug)
        if vtype == VTYPES[1]:  # if indel update variant alleles map positions
            self._update_variant_alleles(posrel, posrel_stop, len(alt) - len(ref))
        # update haplotype sequence and positions map
        self._update_sequence(posrel, posrel_stop, alt)
        self._update_posmap(posrel, chain)

    def add_variants_unphased(self, variants: List[VariantRecord], sample: str) -> None:
        """Adds a list of unphased variants to the haplotype and updates its
        sequence and sample information.

        Applies each unphased variant to the haplotype, reconstructs the sequence,
        and updates sample, variant, and allele frequency data.

        Args:
            variants (List[VariantRecord]): List of VariantRecord objects to add
                to the haplotype.
            sample (str): Sample identifier to associate with the haplotype.
        """
        variants = _sort_variants(variants)
        chains = _compute_chains(variants)
        self._initialize_posmap(chains, self._coordinates.start)
        for i, variant in enumerate(variants):
            self._insert_variant_unphased(
                variant.position,
                variant.ref,
                variant.alt[0],
                variant.vtype[0],
                chains[i],
                sum(chains[:i]),
            )
        self._sequence = Sequence(
            "".join(self._sequence._sequence_raw), self._debug, allow_lower_case=True
        )
        self._samples = sample
        self._variants = ",".join([v.id[0] for v in variants])
        self._afs = {v.id[0]: v.afs[0] for v in variants}

    def homozygous_samples(self) -> None:
        """Sets the haplotype samples to homozygous by updating their phasing value.

        Changes the phasing value for all samples to indicate homozygosity,
        supporting diploid cases.

        Raises:
            CrisprHawkHaplotypeError: If the haplotype is a reference haplotype.
        """
        # if samples are homozygous, change their phasing value (support diploid)
        if self._samples == "REF":
            exception_handler(
                CrisprHawkHaplotypeError,
                "REF haplotype cannot be homozygous",
                os.EX_DATAERR,
                self._debug,
            )
        self._samples = ",".join(
            [
                f"{s[0]}:1|1"
                for sample in self._samples.split(",")
                for s in [sample.split(":")]
            ]
        )

    def set_afs(self, afs: Dict[str, float]) -> None:
        """Sets the allele frequencies for the haplotype.

        Assigns the provided dictionary of allele frequencies to the haplotype.

        Args:
            afs (Dict[str, float]): Dictionary mapping variant IDs to their
                allele frequencies.
        """
        self._afs = afs  # set allele frequencies to haplotype

    def set_posmap(self, posmap: Dict[int, int], posmap_rev: Dict[int, int]) -> None:
        """Sets the position mapping dictionaries for the haplotype.

        Assigns the provided position map and its reverse to the haplotype for
        coordinate translation.

        Args:
            posmap (Dict[int, int]): Dictionary mapping relative positions to
                genomic positions.
            posmap_rev (Dict[int, int]): Dictionary mapping genomic positions to
                relative positions.
        """
        self._posmap = posmap  # set position map to haplotype
        self._posmap_reverse = posmap_rev

    def set_variant_alleles(
        self, variant_alleles: Dict[int, List[Tuple[str, str, int]]]
    ) -> None:
        """Sets the variant alleles mapping for the haplotype.

        Assigns the provided dictionary of variant alleles to the haplotype.

        Args:
            variant_alleles (Dict[int, List[Tuple[str, str, int]]]): Dictionary mapping
                positions to (ref, alt) allele tuples.
        """
        self._variant_alleles = variant_alleles  # set variant alleles

    @property
    def samples(self) -> str:
        return self._samples

    @samples.setter
    def samples(self, value: str) -> None:
        """Sets the samples attribute for the haplotype.

        Validates that the provided value is a string and assigns it to the
        haplotype's samples.

        Args:
            value (str): The sample string to assign.

        Raises:
            CrisprHawkHaplotypeError: If the provided value is not a string.
        """
        if not isinstance(value, str):
            exception_handler(
                CrisprHawkHaplotypeError,
                f"Samples must be a string, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._samples = value  # set samples to haplotype

    @property
    def variants(self) -> str:
        return self._variants

    @variants.setter
    def variants(self, value: str) -> None:
        """Sets the variants attribute for the haplotype.

        Validates that the provided value is a string and assigns it to the
        haplotype's variants.

        Args:
            value (str): The variants string to assign.

        Raises:
            CrisprHawkHaplotypeError: If the provided value is not a string.
        """
        if not isinstance(value, str):
            exception_handler(
                CrisprHawkHaplotypeError,
                f"Variants must be a string, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._variants = value  # set variants to haplotype

    @property
    def afs(self) -> Dict[str, float]:
        return self._afs

    @property
    def phased(self) -> bool:
        return self._phased

    @property
    def posmap(self) -> Dict[int, int]:
        return self._posmap

    @property
    def posmap_rev(self) -> Dict[int, int]:
        return self._posmap_reverse

    @property
    def id(self) -> str:
        return self._id

    @id.setter
    def id(self, value: str) -> None:
        """Sets the haplotype ID attribute.

        Validates that the provided value is a string and assigns it as the
        haplotype's ID.

        Args:
            value (str): The haplotype ID to assign.

        Raises:
            CrisprHawkHaplotypeError: If the provided value is not a string.
        """
        if not isinstance(value, str):
            exception_handler(
                CrisprHawkHaplotypeError,
                f"Haplotype id must be a string, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._id = value  # set haplotype ID

    @property
    def variant_alleles(self) -> Dict[int, List[Tuple[str, str, int]]]:
        return self._variant_alleles


def _sort_variants(variants: List[VariantRecord]) -> List[VariantRecord]:
    """Sorts a list of variants so that SNPs appear before indels.

    Returns a new list with SNPs sorted before indels, each group sorted individually.

    Args:
        variants (List[VariantRecord]): List of VariantRecord objects to sort.

    Returns:
        List[VariantRecord]: Sorted list of variants with SNPs first, then indels.
    """
    # sort variants set to have snps before indels
    snps, indels = [], []
    for variant in variants:
        if variant.vtype[0] == VTYPES[0]:  # snp
            snps.append(variant)
        else:  # indel
            indels.append(variant)
    return sorted(snps) + sorted(indels)


def _compute_chains(variants: List[VariantRecord]) -> List[int]:
    """Computes the length difference between alternate and reference alleles for
    each variant.

    Returns a list of length changes (insertions or deletions) for the provided
    variants.

    Args:
        variants (List[VariantRecord]): List of VariantRecord objects to process.

    Returns:
        List[int]: List of length differences for each variant (alt length - ref
            length).
    """
    return [len(v.alt[0]) - len(v.ref) for v in variants]


def _encode_iupac(ref: str, alt: str, position: int, debug: bool) -> str:
    """Encodes a reference and alternate allele pair as an IUPAC character.

    Returns the IUPAC code for the given alleles or raises an error if the
    combination is invalid.

    Args:
        ref (str): Reference allele.
        alt (str): Alternate allele.
        position (int): Genomic position of the variant.
        debug (bool): Debug flag for verbose output.

    Returns:
        str: IUPAC character representing the allele combination.

    Raises:
        CrisprHawkIupacTableError: If the allele combination cannot be encoded
            as IUPAC.
    """
    try:
        return IUPAC_ENCODER["".join({IUPACTABLE[ref.upper()], alt})]
    except KeyError as e:
        exception_handler(
            CrisprHawkIupacTableError,
            f"An error occurred while encoding {ref}>{alt} at position {position} as IUPAC character",
            os.EX_DATAERR,
            debug,
            e,
        )


class HaplotypeIndel(Haplotype):
    """Represents a haplotype with additional support for indel offset and position
    tracking.

    Extends the Haplotype class to manage indel-specific attributes such as offset
    and indel position.

    Attributes:
        _offset (int): The offset value for the haplotype indel.
        _indel_position (int): The position of the indel in the haplotype.
    """

    def __init__(
        self,
        sequence: Sequence,
        coord: Coordinate,
        phased: bool,
        chromcopy: int,
        debug: bool,
    ) -> None:
        super().__init__(sequence, coord, phased, chromcopy, debug)
        self._offset = 0
        self._indel_position = -1

    @property
    def offset(self) -> int:
        return self._offset

    @offset.setter
    def offset(self, value: int) -> None:
        """Sets the offset value for the haplotype indel.

        Validates that the provided value is an integer and assigns it as the offset.

        Args:
            value (int): The offset value to assign.

        Raises:
            CrisprHawkHaplotypeError: If the provided value is not an integer.
        """
        if not isinstance(value, int):
            exception_handler(
                CrisprHawkHaplotypeError,
                f"Offset must be a int, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._offset = value

    @property
    def indel_position(self) -> int:
        return self._indel_position

    @indel_position.setter
    def indel_position(self, value: int) -> None:
        """Sets the indel position for the haplotype.

        Validates that the provided value is an integer and assigns it as the
        indel position.

        Args:
            value (int): The indel position to assign.

        Raises:
            CrisprHawkHaplotypeError: If the provided value is not an integer.
        """
        if not isinstance(value, int):
            exception_handler(
                CrisprHawkHaplotypeError,
                f"Indel position must be a int, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._indel_position = value
