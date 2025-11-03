"""
This module provides classes and functions for identifying and scoring microhomology
patterns in nucleotide sequences.

This module adapts the method described in Bae et al., Nat Methods 2014.

It includes utilities for matching sequences with IUPAC ambiguity codes, finding
and scoring microhomology patterns, and calculating microhomology scores for guide
sequences. The module is designed to support genome editing applications by evaluating
the likelihood of microhomology-mediated end joining events.
"""

from itertools import product
from typing import Union, List, Tuple

import math

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


class MicrohomologyPattern:
    """
    Represents a microhomology pattern within a nucleotide sequence.

    This class stores the sequence and positional information for a microhomology
    pattern, including the left and right start/stop positions and the pattern
    length. It provides equality, hashing, and string representation methods for
    use in microhomology scoring.

    Args:
        sequence (str): The microhomology sequence.
        left_start (int): Start index of the pattern on the left side.
        left_stop (int): Stop index of the pattern on the left side.
        right_start (int): Start index of the pattern on the right side.
        right_stop (int): Stop index of the pattern on the right side.
        length (int): Length of the microhomology pattern.

    Attributes:
        _sequence (str): The microhomology sequence.
        _left_start (int): Start index of the pattern on the left side.
        _left_stop (int): Stop index of the pattern on the left side.
        _right_start (int): Start index of the pattern on the right side.
        _right_stop (int): Stop index of the pattern on the right side.
        _length (int): Length of the microhomology pattern.
        _hash (int): Precomputed hash value for the pattern.
    """

    def __init__(
        self,
        sequence: str,
        left_start: int,
        left_stop: int,
        right_start: int,
        right_stop: int,
        length: int,
    ) -> None:
        self._sequence = sequence
        self._left_start = left_start
        self._left_stop = left_stop
        self._right_start = right_start
        self._right_stop = right_stop
        self._length = length
        self._hash = (
            self._calculate_hash()
        )  # precompute hash for performance (immutable object)

        def __repr__(self) -> str:
            return (
                f"<{self.__class__.__name__} object; sequence={self._sequence} "
                f"left_start={self._left_start} left_stop={self._left_stop} "
                f"right_start={self._right_start} right_stop={self._right_stop} "
                f"length={self._length}>"
            )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, MicrohomologyPattern):
            return NotImplemented
        return (
            self._sequence == other.sequence
            and self._left_start == other.left_start
            and self._left_stop == other.left_stop
            and self._right_start == other.right_start
            and self._right_stop == other.right_stop
            and self._length == other.length
        )

    def __hash__(self) -> int:
        return self._hash

    def _calculate_hash(self) -> int:
        """
        Compute a hash value for the microhomology pattern.

        Returns a hash based on the sequence and positional attributes of the
        pattern, enabling use in sets and as dictionary keys.

        Returns:
            int: The hash value for the microhomology pattern.
        """
        return hash(
            (
                self._sequence,
                self._left_start,
                self._left_stop,
                self._right_start,
                self._right_stop,
                self._length,
            )
        )

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def left_start(self) -> int:
        return self._left_start

    @property
    def left_stop(self) -> int:
        return self._left_stop

    @property
    def right_start(self) -> int:
        return self._right_start

    @property
    def right_stop(self) -> int:
        return self._right_stop

    @property
    def length(self) -> int:
        return self._length


class MicrohomologyResult:
    """
    Stores the results of microhomology scoring for a guide sequence.

    This class holds the total microhomology score, the out-of-frame score, and
    the list of deletion patterns with their associated scores. It provides
    property accessors for each result component.

    Args:
        mh_score (int): The total microhomology score.
        ooframe_score (int): The out-of-frame score as a percentage of the total
            score.
        deletion_patterns (List[Tuple[float, str]]): List of tuples containing
            the score and the corresponding deletion sequence.

    Attributes:
        _mh_score (int): The total microhomology score.
        _ooframe_score (int): The out-of-frame score.
        _deletion_patterns (List[Tuple[float, str]]): The deletion patterns and
            their scores.
    """

    def __init__(
        self,
        mh_score: int,
        ooframe_score: int,
        deletion_patterns: List[Tuple[float, str]],
    ) -> None:
        self._mh_score = mh_score
        self._ooframe_score = ooframe_score
        self._deletion_patterns = deletion_patterns

    @property
    def mh_score(self) -> int:
        return self._mh_score

    @property
    def ooframe_score(self) -> int:
        return self._ooframe_score

    @property
    def deletion_patterns(self) -> List[Tuple[float, str]]:
        return self._deletion_patterns


def _match(seq1: str, seq2: str) -> bool:
    """
    Determines if two nucleotide sequences match, considering IUPAC ambiguity codes.

    Compares each position in the two sequences and returns True if all positions
    have at least one nucleotide in common according to IUPAC codes, otherwise False.

    Args:
        seq1 (str): The first nucleotide sequence.
        seq2 (str): The second nucleotide sequence.

    Returns:
        bool: True if the sequences match at all positions, False otherwise.

    Raises:
        ValueError: If an invalid IUPAC character is encountered.
    """
    seq1, seq2 = seq1.upper(), seq2.upper()
    for i, (nt1, nt2) in enumerate(zip(seq1, seq2)):
        if nt1 not in IUPACTABLE or nt2 not in IUPACTABLE:
            raise ValueError(
                f"Invalid IUPAC characters ({nt1} - {nt2}) at position {i}"
            )
        ntiupac1, ntiupac2 = set(IUPACTABLE[nt1]), set(IUPACTABLE[nt2])
        if not ntiupac1.intersection(ntiupac2):  # check possible overlap
            return False
    return True


def _find_microhomology_patterns(
    sequence: str, start: int, stop: int
) -> List[MicrohomologyPattern]:
    """
    Identify all microhomology patterns in a nucleotide sequence.

    Searches for and returns all microhomology patterns of length 2 up to start-1
    within the specified region of the sequence. Each pattern is represented as
    a MicrohomologyPattern object.

    Args:
        sequence (str): The nucleotide sequence to search.
        start (int): The starting index for pattern search.
        stop (int): The number of positions to search after the start index.

    Returns:
        List[MicrohomologyPattern]: A list of found microhomology patterns.
    """
    patterns = []  # list of microhomology patterns founf for current guide
    for k in reversed(list(range(2, start))):
        # search for patterns of length k (from 2 to start - 1)
        for j, i in product(range(start, start + stop - k + 1), range(start - k + 1)):
            leftseq = sequence[i : i + k]
            rightseq = sequence[j : j + k]
            if _match(leftseq, rightseq):
                deletion_length = j - i
                patterns.append(
                    MicrohomologyPattern(leftseq, i, i + k, j, j + k, deletion_length)
                )
    return patterns


def _compute_pattern_scores(
    patterns: List[MicrohomologyPattern],
    sequence: str,
    length_weight: Union[None, float],
) -> Tuple[List[Tuple[float, str]], float, float]:
    """
    Compute scores for microhomology patterns in a sequence.

    Calculates the score for each microhomology pattern based on length, GC content,
    and position, and classifies them as in-frame or out-of-frame. Returns the
    deletion patterns with their scores, and the total in-frame and out-of-frame scores.

    Args:
        patterns (List[MicrohomologyPattern]): The list of microhomology patterns.
        sequence (str): The nucleotide sequence containing the patterns.
        length_weight (Union[None, float]): The length weighting factor for scoring.

    Returns:
        Tuple[List[Tuple[float, str]], float, float]: A tuple containing the list
            of (score, deletion sequence) pairs, the total in-frame score, and the
            total out-of-frame score.
    """
    deletion_patterns = []  # initialize deletion patterns
    iframe_score, ooframe_score = 0.0, 0.0  # initialize scores
    for pattern in patterns:
        # compute length factor (exponential decay)
        assert length_weight  # should never be none
        length_factor = round(1 / math.exp((pattern.length) / length_weight), 3)
        # compute GC content bonus
        gc_count = pattern.sequence.count("G") + pattern.sequence.count("C")
        at_count = len(pattern.sequence) - gc_count
        # compute score: AT bases = 1x, GC bases = 2x
        score = 100 * length_factor * (at_count + (gc_count * 2))
        if pattern.length % 3 == 0:  # categorize as in-frame or out-of-frame
            iframe_score += score
        else:
            ooframe_score += score
        # create deletion sequence (replace deleted region with dashes)
        deletion_seq = (
            sequence[: pattern.left_stop]
            + ("-" * pattern.length)
            + sequence[pattern.right_stop :]
        )
        deletion_patterns.append((float(score), deletion_seq))
    return deletion_patterns, iframe_score, ooframe_score


def _remove_duplicate_patterns(
    patterns: List[MicrohomologyPattern],
) -> List[MicrohomologyPattern]:
    """
    Remove duplicate microhomology patterns from a list.

    Filters out duplicate microhomology patterns based on their positional
    attributes, ensuring only unique patterns are retained in the result.

    Args:
        patterns (List[MicrohomologyPattern]): The list of microhomology patterns
            to filter.

    Returns:
        List[MicrohomologyPattern]: The list of unique microhomology patterns.
    """
    patterns_unique = []
    for i, pattern in enumerate(patterns):
        n = sum(
            pattern.left_start >= patterns[j].left_start
            and pattern.left_stop <= patterns[j].left_stop
            and pattern.right_start >= patterns[j].right_start
            and pattern.right_stop <= patterns[j].right_stop
            and (pattern.left_start - patterns[j].left_start)
            == (pattern.right_start - patterns[j].right_start)
            and (pattern.left_stop - patterns[j].left_stop)
            == (pattern.right_stop - patterns[j].right_stop)
            for j in range(i)
        )
        if n == 0:
            patterns_unique.append(pattern)
    return patterns_unique


def calculate_microhomology_score(
    guideseq: str, start: int, lweight: float = 20.0
) -> MicrohomologyResult:
    """
    Calculate the microhomology score for a guide sequence.

    Identifies all unique microhomology patterns in the guide sequence, computes
    their scores, and returns a MicrohomologyResult containing the total score,
    out-of-frame score, and deletion patterns.

    Args:
        guideseq (str): The guide sequence (must be uppercase).
        start (int): The starting index for microhomology pattern search.
        lweight (Optional[float]): The length weighting factor for scoring (default: 20.0).

    Returns:
        MicrohomologyResult: The result object containing scores and deletion patterns.
    """
    assert guideseq.isupper()  # must be upper case (handle alternative guides)
    # find all microhomology patterns for current guide
    patterns = _find_microhomology_patterns(guideseq, start, len(guideseq) - start)
    if not patterns:  # nothing to score here
        return MicrohomologyResult(0, 0, [])
    # remove duplicate patterns and compute scores
    patterns_unique = _remove_duplicate_patterns(patterns)
    deletion_patterns, iframe_score, ooframe_score = _compute_pattern_scores(
        patterns_unique, guideseq, lweight
    )
    # compute final scores
    mh_score = int(iframe_score + ooframe_score)
    ooframe_score = int((ooframe_score * 100) / mh_score)
    return MicrohomologyResult(mh_score, ooframe_score, deletion_patterns)
