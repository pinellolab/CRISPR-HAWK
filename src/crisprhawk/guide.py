"""Provides classes and utilities for representing and manipulating CRISPR guide
sequences.

This module defines the Guide class for encapsulating CRISPR guide properties,
annotations, and scoring, as well as the GuideIterator for iterating over guide
sequences. It supports guide sequence management, annotation, and scoring for
CRISPR analysis workflows.
"""

from .exception_handlers import exception_handler
from .crisprhawk_error import CrisprHawkGuideError
from .utils import round_score, RC

from typing import Dict, Union, List

import numpy as np

import os

GUIDESEQPAD = 50  # upstream and downstream sequence padding for guides scoring


class Guide:
    """Represents a CRISPR guide and provides methods for manipulating and
    annotating guide properties.

    This class encapsulates guide sequence, position, direction, and associated
    scores and annotations for CRISPR analysis.

    Attributes:
        _debug (bool): Debug mode flag.
        _guidelen (int): Guide length.
        _pamlen (int): PAM length.
        _start (int): Guide start position.
        _stop (int): Guide stop position.
        _sequence (str): Guide sequence as a string.
        _right (bool): Guide position on the right side of PAM.
        _pamseq (str): PAM sequence.
        _guideseq (str): Guide sequence (without PAM).
        _direction (int): Guide direction.
        _samples (str): Samples carrying guide variants.
        _variants (str): Variants overlapping guide.
        _afs (Dict[str, float]): Allele frequencies of variants overlapping guide.
        _hapid (str): Haplotype ID.
        _guide_id (str): Unique guide identifier.
        _azimuth_score (str): Azimuth score for guide.
        _rs3_score (str): RS3 score for guide.
        _cfdon_score (str): CFDOn score for guide.
        _elevationon_score (str): ElevationOn score for guide.
        _deepcpf1_score (str): DeepCpf1 score for guide.
        _ooframe_score (str): Out-of-frame score for guide.
        _cfd (str): CFD score for guide.
        _gc (str): GC content of guide.
        _funcann (List[str]): Functional annotations for guide.
        _geneann (List[str]): Gene annotations for guide.
        _offtargets_num (int): Number of off-targets for guide.
        _afs_ (str): Allele frequencies as a string.
    """

    def __init__(
        self,
        position_start: int,
        position_stop: int,
        sequence: str,
        guidelen: int,
        pamlen: int,
        direction: int,
        samples: str,
        variants: str,
        afs: Dict[str, float],
        debug: bool,
        right: bool,
        hapid: str,
    ) -> None:
        """Initializes a Guide object with sequence, position, direction, and
        annotation data.

        This constructor sets up the guide's sequence, genomic coordinates, direction,
        and initializes all relevant scores and annotations.

        Args:
            position_start: The start position of the guide.
            position_stop: The stop position of the guide.
            sequence: The nucleotide sequence of the guide.
            guidelen: The length of the guide.
            pamlen: The length of the PAM sequence.
            direction: The direction of the guide (strand).
            samples: Samples carrying guide variants.
            variants: Variants overlapping the guide.
            afs: Allele frequencies of variants overlapping the guide.
            debug: Flag to enable debug mode.
            right: Indicates if the guide is on the right side of the PAM.
            hapid: Haplotype ID for the guide.

        Returns:
            None
        """
        self._debug = debug  # store debug mode
        self._guidelen = guidelen  # guide length
        self._pamlen = pamlen  # pam lenght
        self._start = position_start  # guide start position
        self._stop = position_stop  # guide stop position
        self._sequence = sequence  # sequence as string
        self._right = right  # guide position on the right side of pam
        self._compute_pamguide_sequences()  # compute pam and guide sequences
        self._direction = direction  # guide direction
        self._samples = samples  # samples carrying guide variants
        self._variants = variants  # variants overlapping guide
        self._afs = afs  # allele frequencies of variants overlapping guide
        self._hapid = hapid  # haplotype ID
        self._compute_guide_id()  # compute unique guide ID
        self._initialize_scores()  # initialize scores to NAs
        self._initialize_annotations()  # initialize annotations to NAs

    def __repr__(self) -> str:
        """Returns a string representation of the Guide object.

        This method provides a concise summary of the guide's start, stop, sequence,
        and direction for debugging and display purposes.

        Returns:
            str: A string describing the Guide object.
        """
        return (
            f"<{self.__class__.__name__} object; start={self._start} "
            f"stop={self._stop} sequence={self._sequence} direction={self._direction}>"
        )

    def __len__(self) -> int:
        """Returns the length of the guide sequence.

        This method provides the total number of nucleotides in the guide's sequence.

        Returns:
            int: The length of the guide sequence.
        """
        return len(self._sequence)

    def __getitem__(self, idx: Union[int, slice]) -> str:
        """Returns a substring of the guide sequence at the specified index or slice.

        This method allows indexed or sliced access to the guide's sequence, raising
        an error if the index is out of range.

        Args:
            idx: An integer index or slice specifying the desired portion of the
                sequence.

        Returns:
            str: The substring of the guide sequence at the given index or slice.

        Raises:
            CrisprHawkGuideError: If the index is out of range.
        """
        assert hasattr(self, "_sequence")
        try:
            return "".join(self._sequence[idx])
        except IndexError as e:
            exception_handler(
                CrisprHawkGuideError,  # type: ignore
                f"Index {idx} out of range",
                os.EX_DATAERR,
                self._debug,
                e,
            )

    def __iter__(self) -> "GuideIterator":
        """Returns an iterator over the guide sequence.

        This method enables iteration over each nucleotide in the guide's sequence.

        Returns:
            GuideIterator: An iterator for the guide sequence.
        """
        return GuideIterator(self)

    def _compute_pamguide_sequences(self) -> None:
        """Computes and sets the PAM and guide sequences for the Guide object.

        This method extracts the PAM and guide sequences from the full guide
        sequence, accounting for padding and orientation.
        """
        assert hasattr(self, "_sequence")  # check if sequence is available
        sequence = self._sequence[GUIDESEQPAD:-GUIDESEQPAD]  # remove padding
        if self._right:  # guide on the right side of pam
            self._pamseq = sequence[: self._pamlen]
            self._guideseq = sequence[self._pamlen :]
        else:  # guide on the left side of pam
            self._pamseq = sequence[-self._pamlen :]
            self._guideseq = sequence[: -self._pamlen]

    def _compute_guide_id(self) -> None:
        """Generates and sets a unique identifier for the Guide object.

        This method creates a guide ID string using the guide's start, stop,
        direction, haplotype ID, and guide sequence.
        """
        assert (
            hasattr(self, "_start")
            and hasattr(self, "_stop")
            and hasattr(self, "_direction")
            and hasattr(self, "_hapid")
            and hasattr(self, "_guideseq")
        )
        self._guide_id = f"{self._start}_{self._stop}_{self._direction}_{self._hapid}_{self._guideseq}"  # unique identifier for the guide

    def _initialize_scores(self) -> None:
        """Initializes all scoring attributes for the Guide object to 'NA'.

        This method sets default values for all guide scoring metrics, ensuring
        they are marked as not available until computed.
        """
        # initialize scores for guide to NA
        self._azimuth_score = "NA"
        self._rs3_score = "NA"
        self._cfdon_score = "NA"
        self._elevationon_score = "NA"
        self._deepcpf1_score = "NA"
        self._ooframe_score = "NA"
        self._cfd = "NA"

    def _initialize_annotations(self) -> None:
        """Initializes all annotation attributes for the Guide object to default
        values.

        This method sets the guide's GC content to 'NA' and initializes the functional
        and gene annotation lists as empty.
        """
        # initialize annotations for guide to NA
        self._gc = "NA"
        self._offtargets_num = "NA"
        self._funcann = []
        self._geneann = []

    def reverse_complement(self) -> None:
        """Computes the reverse complement of the guide sequence and updates orientation.

        This method reverses the guide sequence, computes its complement, and toggles
        the guide's orientation relative to the PAM.
        """
        assert hasattr(self, "_sequence")
        # compute reverse complement sequence
        self._sequence = "".join([RC[nt] for nt in self._sequence[::-1]])
        self._right = not self._right  # update guide direction
        self._compute_pamguide_sequences()  # adjust pam and guide sequences

    @property
    def start(self) -> int:
        return self._start

    @property
    def stop(self) -> int:
        return self._stop

    @property
    def strand(self) -> int:
        return self._direction

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def samples(self) -> str:
        return self._samples

    @property
    def variants(self) -> str:
        return self._variants

    @variants.setter
    def variants(self, value: str) -> None:
        """Sets the variants property for the Guide object.

        This method updates the guide's variants information, validating that the
        input is a non-empty string.

        Args:
            value: A string representing the variants overlapping the guide.

        Raises:
            CrisprHawkGuideError: If the value is not a non-empty string.
        """
        print(f"value: {value}")
        if not isinstance(value, str) or not value.strip():
            exception_handler(
                CrisprHawkGuideError,
                f"Variants must be a non-empty string, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._variants = value  # set variants

    @property
    def afs(self) -> Dict[str, float]:
        return self._afs

    @property
    def afs_str(self) -> str:
        return self._afs_

    @afs_str.setter
    def afs_str(self, value: List[str]) -> None:
        """Sets the allele frequencies string for the Guide object.

        This method updates the guide's allele frequencies string, converting a
        list of values to a comma-separated string or 'NA' if not available.

        Args:
            value: A list of allele frequency strings.

        Returns:
            None
        """
        self._afs_ = (
            "NA"
            if not value or (len(set(value)) == 1 and value[0] == "NA")
            else ",".join(value)
        )

    @property
    def pam(self) -> str:
        return self._pamseq

    @property
    def pamlen(self) -> int:
        return self._pamlen

    @property
    def guide(self) -> str:
        return self._guideseq

    @property
    def guidepam(self) -> str:
        if self._right:
            return self._pamseq + self._guideseq
        return self._guideseq + self._pamseq

    @property
    def guidelen(self) -> int:
        return self._guidelen

    @property
    def right(self) -> bool:
        return self._right

    @property
    def guide_id(self) -> str:
        return self._guide_id

    @property
    def azimuth_score(self) -> str:
        return self._azimuth_score

    @azimuth_score.setter
    def azimuth_score(self, value: float) -> None:
        """Sets the azimuth score for the Guide object.

        This method validates and updates the guide's azimuth score, converting
        it to a string or 'NA' if not available.

        Args:
            value: The azimuth score as a float.

        Raises:
            CrisprHawkGuideError: If the value is not a float or is not valid.
        """
        if not isinstance(value, float):
            exception_handler(
                CrisprHawkGuideError,
                f"Azimuth score must be a float, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._azimuth_score = "NA" if np.isnan(value) else str(round_score(value))

    @property
    def rs3_score(self) -> str:
        return self._rs3_score

    @rs3_score.setter
    def rs3_score(self, value: float) -> None:
        """Sets the RS3 score for the Guide object.

        This method validates and updates the guide's RS3 score, converting it to
        a string or 'NA' if not available.

        Args:
            value: The RS3 score as a float.

        Raises:
            CrisprHawkGuideError: If the value is not a float or is not valid.
        """
        if not isinstance(value, float):
            exception_handler(
                CrisprHawkGuideError,
                f"RS3 score must be a float, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._rs3_score = "NA" if np.isnan(value) else str(round_score(value))

    @property
    def deepcpf1_score(self) -> str:
        return self._deepcpf1_score

    @deepcpf1_score.setter
    def deepcpf1_score(self, value: float) -> None:
        """Sets the DeepCpf1 score for the Guide object.

        This method validates and updates the guide's DeepCpf1 score, converting
        it to a string or 'NA' if not available.

        Args:
            value: The DeepCpf1 score as a float.

        Raises:
            CrisprHawkGuideError: If the value is not a float or is not valid.
        """
        if not isinstance(value, float):
            exception_handler(
                CrisprHawkGuideError,
                f"DeepCpf1 score must be a float, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._deepcpf1_score = "NA" if np.isnan(value) else str(round_score(value))

    @property
    def cfdon_score(self) -> str:
        return self._cfdon_score

    @cfdon_score.setter
    def cfdon_score(self, value: float) -> None:
        """Sets the CFD-on score for the Guide object.

        This method validates and updates the guide's CFD-on score, converting it
        to a string or 'NA' if not available.

        Args:
            value: The CFD-on score as a float.

        Raises:
            CrisprHawkGuideError: If the value is not a float or is not valid.
        """
        if not isinstance(value, float):
            exception_handler(
                CrisprHawkGuideError,
                f"CFD-on score must be a float, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._cfdon_score = "NA" if np.isnan(value) else str(round_score(value))

    @property
    def elevationon_score(self) -> str:
        return self._elevationon_score

    @elevationon_score.setter
    def elevationon_score(self, value: float) -> None:
        """Sets the Elevation-on score for the Guide object.

        This method validates and updates the guide's Elevation-on score, converting
        it to a string or 'NA' if not available.

        Args:
            value: The Elevation-on score as a float.

        Raises:
            CrisprHawkGuideError: If the value is not a float or is not valid.
        """
        if not isinstance(value, float):
            exception_handler(
                CrisprHawkGuideError,
                f"Elevation-on score must be a float, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._elevationon_score = "NA" if np.isnan(value) else str(round_score(value))

    @property
    def ooframe_score(self) -> str:
        return self._ooframe_score

    @ooframe_score.setter
    def ooframe_score(self, value: int) -> None:
        """Sets the out-of-frame score for the Guide object.

        This method validates and updates the guide's out-of-frame score, converting
        it to a string.

        Args:
            value: The out-of-frame score as an integer.

        Raises:
            CrisprHawkGuideError: If the value is not an int or is not valid.
        """
        if not isinstance(value, int):
            exception_handler(
                CrisprHawkGuideError,
                f"Out-of-frame score must be an int, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._ooframe_score = str(value)

    @property
    def gc(self) -> str:
        return self._gc

    @gc.setter
    def gc(self, value: float) -> None:
        """Sets the GC content for the Guide object.

        This method validates and updates the guide's GC content, converting it
        to a string.

        Args:
            value: The GC content as a float.

        Raises:
            CrisprHawkGuideError: If the value is not a float.
        """
        if not isinstance(value, float):
            exception_handler(
                CrisprHawkGuideError,
                f"Out-of-frame score must be a float, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._gc = str(value)

    @property
    def hapid(self) -> str:
        return self._hapid

    @property
    def funcann(self) -> List[str]:
        if not hasattr(self, "_funcann"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _funcann attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return self._funcann

    @funcann.setter
    def funcann(self, value: str) -> None:
        if not isinstance(value, str):
            exception_handler(
                CrisprHawkGuideError,
                f"Annotation score must be a string, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._funcann.append(value)

    @property
    def geneann(self) -> List[str]:
        if not hasattr(self, "_geneann"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _geneann attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return self._geneann

    @geneann.setter
    def geneann(self, value) -> None:
        """Sets the gene annotation for the Guide object.

        This method validates and appends a gene annotation string to the guide's
        gene annotation list.

        Args:
            value: The gene annotation as a string.

        Raises:
            CrisprHawkGuideError: If ```python the value is not a string.
        """
        if not isinstance(value, str):
            exception_handler(
                CrisprHawkGuideError,
                f"Gene annotation score must be a string, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._geneann.append(value)

    @property
    def offtargets(self) -> str:
        if not hasattr(self, "_offtargets_num"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _offtargets_num attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return self._offtargets_num

    @offtargets.setter
    def offtargets(self, value: int) -> None:
        """Sets the number of off-targets for the Guide object.

        This method validates and updates the number of off-targets associated
        with the guide.

        Args:
            value: The number of off-targets as an integer.

        Raises:
            CrisprHawkGuideError: If the value is not an int or is not valid.
        """
        if not isinstance(value, int) or not value:
            exception_handler(
                CrisprHawkGuideError,
                f"Off-targets number must be an int, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._offtargets_num = str(value)

    @property
    def cfd(self) -> str:
        if not hasattr(self, "_cfd"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _cfd attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return self._cfd

    @cfd.setter
    def cfd(self, value: float) -> None:
        """Sets the CFD score for the Guide object.

        This method validates and updates the guide's CFD score, converting it to
        a string or 'NA' if not available.

        Args:
            value: The CFD score as a float.

        Raises:
            CrisprHawkGuideError: If the value is not a float or is not valid.
        """
        if not isinstance(value, float):
            exception_handler(
                CrisprHawkGuideError,
                f"CFD score must be a float, got {type(value).__name__} instead",
                os.EX_DATAERR,
                True,
            )
        self._cfd = "NA" if np.isnan(value) else str(round_score(value))


class GuideIterator:
    """Iterator for traversing the sequence of a Guide object.

    This class enables sequential access to each nucleotide in a Guide's sequence.

    Attributes:
        _guide (Guide): The Guide object to iterate over.
        _index (int): The current index in the guide sequence.
    """

    def __init__(self, guide: Guide) -> None:
        """Initializes the GuideIterator with a Guide object.

        This constructor sets up the iterator to traverse the sequence of the provided
        Guide.

        Args:
            guide: The Guide object to iterate over.

        Raises:
            AttributeError: If the Guide object does not have a '_guide' attribute.
        """
        if not hasattr(guide, "_sequence"):  # always trace this error
            exception_handler(
                AttributeError,  # type: ignore
                f"Missing _sequence attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        self._sequence = guide  # guide object to iterate over
        self._index = 0  # iterator index used over the guide sequence

    def __next__(self) -> str:
        """Returns the next nucleotide in the Guide sequence during iteration.

        This method retrieves the next nucleotide from the guide's sequence or
        raises StopIteration when finished.

        Returns:
            str: The next nucleotide in the guide sequence.

        Raises:
            StopIteration: If the end of the guide sequence is reached.
        """
        if self._index < len(self._sequence):
            result = self._sequence[self._index]
            self._index += 1  # go to next position in the guide sequence
            return result
        raise StopIteration  # stop iteration over guide object
