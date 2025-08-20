""" """

from .exception_handlers import exception_handler
from .crisprhawk_error import CrisprHawkGuideError
from .utils import round_score, RC

from typing import List, Union, Set
import os

GUIDESEQPAD = 10  # upstream and downstream sequence padding for guides scoring


class Guide:
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
        debug: bool,
        right: bool,
        hapid: str,
    ) -> None:
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
        self._hapid = hapid  # haplotype ID
        self._compute_guide_id()  # compute unique guide ID
        self._initialize_scores()  # initialize scores to NAs

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} object; start={self._start} stop={self._stop} sequence={self._sequence} direction={self._direction}>"

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, idx: Union[int, slice]) -> str:
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
        return GuideIterator(self)

    def _compute_pamguide_sequences(self) -> None:
        assert hasattr(self, "_sequence")  # check if sequence is available
        sequence = self._sequence[GUIDESEQPAD:-GUIDESEQPAD]  # remove padding
        if self._right:  # guide on the right side of pam
            self._pamseq = sequence[: self._pamlen]
            self._guideseq = sequence[self._pamlen :]
        else:  # guide on the left side of pam
            self._pamseq = sequence[-self._pamlen :]
            self._guideseq = sequence[: -self._pamlen]

    def _compute_guide_id(self) -> None:
        assert (
            hasattr(self, "_start")
            and hasattr(self, "_stop")
            and hasattr(self, "_direction")
            and hasattr(self, "_hapid")
            and hasattr(self, "_guideseq")
        )
        self._guide_id = f"{self._start}_{self._stop}_{self._direction}_{self._hapid}_{self._guideseq}"  # unique identifier for the guide

    def _initialize_scores(self) -> None:
        # initialize scores for guide to NA
        self._azimuth_score = "NA"

    def reverse_complement(self) -> None:
        assert hasattr(self, "_sequence")
        # compute reverse complement sequence
        self._sequence = "".join([RC[nt] for nt in self._sequence[::-1]])
        self._right = not self._right  # update guide direction
        self._compute_pamguide_sequences()  # adjust pam and guide sequences

    def set_variants(self, variants: str) -> None:
        self._variants = variants  # set variants

    def set_azimuth_score(self, score: float) -> None:
        if not isinstance(score, float):
            exception_handler(
                TypeError,
                f"Expected azimuth score of type {float.__name__}, got {type(score).__name__}",
                os.EX_DATAERR,
                self._debug,
            )
        self._azimuth_score = str(round_score(score))

    def set_rs3_score(self, score: float) -> None:
        if not isinstance(score, float):
            exception_handler(
                TypeError,
                f"Expected rs3 score of type {float.__name__}, got {type(score).__name__}",
                os.EX_DATAERR,
                self._debug,
            )
        self._rs3_score = str(round_score(score))

    def set_deepcpf1_score(self, score: float) -> None:
        if not isinstance(score, float):
            exception_handler(
                TypeError,
                f"Expected deepCpf1 score of type {float.__name__}, got {type(score).__name__}",
                os.EX_DATAERR,
                self._debug,
            )
        self._deepcpf1_score = str(round_score(score))

    def set_cfdon_score(self, score: float) -> None:
        if not isinstance(score, float):
            exception_handler(
                TypeError,
                f"Expected cfdon score of type {float.__name__}, got {type(score).__name__}",
                os.EX_DATAERR,
                self._debug,
            )
        self._cfdon_score = str(round_score(score))

    def set_func_ann(self, annotation: str) -> None:
        self._funcann = annotation

    def set_gene_ann(self, annotation: str) -> None:
        self._geneann = annotation

    def set_offtargets(self, offtargets_num: int) -> None:
        self._offtargets_num = offtargets_num

    def set_cfd(self, cfd: float) -> None:
        self._cfd = cfd

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

    @property
    def rs3_score(self) -> str:
        return self._rs3_score

    @property
    def deepcpf1_score(self) -> str:
        return self._deepcpf1_score

    @property
    def cfdon_score(self) -> str:
        return self._cfdon_score

    @property
    def hapid(self) -> str:
        return self._hapid

    @property
    def funcann(self) -> str:
        if not hasattr(self, "_funcann"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _funcann attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return self._funcann

    @property
    def geneann(self) -> str:
        if not hasattr(self, "_geneann"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _geneann attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return self._geneann

    @property
    def offtargets(self) -> int:
        if not hasattr(self, "_offtargets_num"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _offtargets_num attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return self._offtargets_num

    @property
    def cfd(self) -> float:
        if not hasattr(self, "_cfd"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _cfd attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return self._cfd


class GuideIterator:
    def __init__(self, guide: Guide) -> None:
        if not hasattr(guide, "_guide"):  # always trace this error
            exception_handler(
                AttributeError,  # type: ignore
                f"Missing _guide attribute on {self.__class__.__name__}",
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
