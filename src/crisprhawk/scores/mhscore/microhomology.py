""" """

from typing import Optional, Union, List, Tuple

import math


class MicrohomologyPattern:
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


def _find_microhomology_patterns(
    sequence: str, start: int, stop: int
) -> List[MicrohomologyPattern]:
    patterns = []  # list of microhomology patterns founf for current guide
    for k in reversed(
        list(range(2, start))
    ):  # search for patterns of length k (from 2 to start - 1)
        for j in range(start, start + stop - k + 1):  # check pos after breakpoints
            for i in range(start - k + 1):  # check positions before breakpoints
                leftseq = sequence[i : i + k]
                rightseq = sequence[j : j + k]
                if leftseq == rightseq:
                    deletion_length = j - i
                    patterns.append(
                        MicrohomologyPattern(
                            leftseq, i, i + k, j, j + k, deletion_length
                        )
                    )
    return patterns


def _compute_pattern_scores(
    patterns: List[MicrohomologyPattern],
    sequence: str,
    length_weight: Union[None, float],
) -> Tuple[List[Tuple[float, str]], float, float]:
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
    patterns_unique = []
    for i, pattern in enumerate(patterns):
        n = 0
        for j in range(i):
            if (
                (pattern.left_start >= patterns[j].left_start)
                and (pattern.left_stop <= patterns[j].left_stop)
                and (pattern.right_start >= patterns[j].right_start)
                and (pattern.right_stop <= patterns[j].right_stop)
                and (pattern.left_start - patterns[j].left_start)
                == (pattern.right_start - patterns[j].right_start)
                and (pattern.left_stop - patterns[j].left_stop)
                == (pattern.right_stop - patterns[j].right_stop)
            ):
                n += 1
        if n == 0:
            patterns_unique.append(pattern)
    return patterns_unique


def calculate_microhomology_score(
    guideseq: str, start: int, lweight: Optional[float] = 20.0
):
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
