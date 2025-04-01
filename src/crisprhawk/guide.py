"""
"""

from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkGuideError
from utils import RC

from typing import List, Union, Set
import os

GUIDESEQPAD = 10  # upstream and downstream sequence padding for guides scoring

class Guide:
    def __init__(
        self,
        position: int,
        sequence: List[str],
        guidelen: int,
        pamlen: int,
        direction: int,
        debug: bool,
        right: bool,
    ) -> None:
        self._guidelen = guidelen  # guide length
        self._pamlen = pamlen  # pam lenght
        self._position = position if right else position - guidelen  # guide position
        self._sequence_raw = list(sequence)  # sequence as list of chars
        self._sequence = sequence  # sequence as string
        self._direction = direction  # guide direction
        self._debug = debug  # store debug mode
        self._samples = "NA"  # samples carrying guide variants
        self._variants = "NA"  # variants overlapping guide
        self._right = right  # guide position on the right side of pam
        self._compute_pamguide_sequences()  # compute pam and guide sequences

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} object; position={self._position} sequence={self._sequence} direction={self._direction}>"

    def __len__(self) -> int:
        return len(self._sequence_raw)

    def __getitem__(self, idx: Union[int, slice]) -> str:
        assert hasattr(self, "_sequence")
        try:
            return "".join(self._sequence_raw[idx])
        except IndexError as e:
            exception_handler(
                CrisprHawkGuideError,
                f"Index {idx} out of range",
                os.EX_DATAERR,
                self._debug,
                e,
            )

    def __iter__(self) -> "GuideIterator":
        return GuideIterator(self)

    def _compute_pamguide_sequences(self) -> None:
        assert hasattr(self, "_sequence_raw")  # check if sequence is available
        if self._right:  # guide on the right side of pam
            self._pamseq = self._sequence[: self._pamlen]
            self._guideseq = self._sequence[self._pamlen :]
        else:  # guide on the left side of pam
            self._pamseq = self._sequence[-self._pamlen :]
            self._guideseq = self._sequence[: -self._pamlen]

    def reverse_complement(self) -> None:
        assert hasattr(self, "_sequence")
        self._sequence_raw = [RC[nt] for nt in self._sequence[::-1]]
        self._sequence = "".join(self._sequence_raw)  # guide reverse complement
        self._right = not self._right  # update guide direction
        self._compute_pamguide_sequences()  # adjust pam and guide sequences

    def set_position(self, position: int) -> None:
        self._position = position  # set guide position

    def set_samples(self, samples: Set[str]) -> None:
        self._samples = ",".join(samples)  # set samples carrying guide variants

    def set_variants(self, variants: Set[str]) -> None:
        self._variants = ",".join(variants)  # set variants overlapping guide

    @property
    def position(self) -> int:
        return self._position

    @property
    def strand(self) -> int:
        return self._direction

    @property
    def sequence(self) -> str:
        return self._sequence[GUIDESEQPAD:-GUIDESEQPAD]

    @property
    def samples(self) -> str:
        return self._samples

    @property
    def variants(self) -> str:
        return self._variants

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
    def guidelen(self) -> int:
        return self._guidelen

    @property
    def right(self) -> bool:
        return self._right

    @property
    def debug(self) -> bool:
        return self._debug


class GuideIterator:
    def __init__(self, guide: Guide) -> None:
        if not hasattr(guide, "_sequence_raw"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _sequence_raw attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        self._guide = guide  # guide object to iterate over
        self._index = 0  # iterator index used over the guide sequence

    def __next__(self) -> str:
        if self._index < len(self._guide):
            result = self._guide[self._index]
            self._index += 1  # go to next position in the guide sequence
            return result
        raise StopIteration  # stop iteration over guide object
