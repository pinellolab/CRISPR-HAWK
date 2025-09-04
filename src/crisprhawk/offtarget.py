""" """

from .exception_handlers import exception_handler
from .scores.cfdscore.cfdscore import compute_cfd
from .bedfile import BedAnnotation
from .guide import Guide
from .utils import round_score

from typing import Dict, Tuple

import numpy as np

import os


class Offtarget:

    def __init__(self, reportline: str, pam: str, right: bool, debug: bool) -> None:
        self._debug = debug  # store debug mode
        self._parse_reportline(reportline, pam, right)  # read report line content
        self._pam = pam  # set off-target pam
        self._cfd_score = "NA"  # default cfd score
        self._elevation_score = "NA"  # default elevation score

    def __repr__(self) -> str:
        return (
            f"<{self.__class__.__name__} object; position={self._pos} "
            f"spacer={self._spacer} strand={self._strand}>"
        )

    def _parse_reportline(self, line: str, pam: str, right: bool) -> None:
        fields = line.strip().split()
        self._chrom = fields[3]  # set chromosome
        self._pos = int(fields[4])  # set off-target position
        self._strand = fields[6]  # set off-target strand
        self._grna_, self._grna = _format_sequence(
            fields[1], pam, right
        )  # set grna sequence
        self._spacer_, self._spacer = _format_sequence(
            fields[2], _retrieve_pam(fields[2], len(pam), right), right
        )  # set spacer sequence
        self._mm = int(fields[7])  # set mismatches number
        self._bulge_type = fields[0]  # set bulge type
        self._bulge_size = int(fields[8])  # set bulge size

    def compute_cfd(
        self, mmscores: Dict[str, float], pamscores: Dict[str, float]
    ) -> None:
        self._cfd_score = str(
            round_score(
                compute_cfd(
                    self._grna_.upper(),
                    self._spacer_.upper(),
                    self._spacer[-2:],
                    mmscores,
                    pamscores,
                    self._debug,
                )
            )
        )

    def set_elevation(self, score: float) -> None:
        if not isinstance(score, float):
            exception_handler(
                TypeError,
                f"Expected elevation-on score of type {float.__name__}, got {type(score).__name__}",
                os.EX_DATAERR,
                self._debug,
            )
        if not np.isnan(score):
            self._elevation_score = str(round_score(score))

    # def annotate_functional(self, funcann: BedAnnotation) -> None:
    #     annotation = funcann.fetch_features(self._chrom, self._start, self._stop, 9)
    #     self._funcann = annotation

    # def annotate_gene(self, geneann: BedAnnotation) -> None:
    #     annotation = geneann.fetch_features(self._chrom, self._start, self._stop, 22)
    #     self._geneann = annotation

    def report_line(self) -> str:
        return "\t".join(
            list(
                map(
                    str,
                    [
                        self._chrom,
                        self._pos,
                        self._strand,
                        self._grna,
                        self._spacer,
                        self._pam,
                        self._mm,
                        self._bulge_size,
                        self._bulge_type,
                        self._cfd_score,
                        self._elevation_score,
                    ],
                )
            )
        )

    @property
    def grna(self) -> str:
        return self._grna

    @property
    def grna_(self) -> str:
        return self._grna_

    @property
    def spacer(self) -> str:
        return self._spacer

    @property
    def cfd(self) -> str:
        return self._cfd_score

    @property
    def elevation(self) -> str:
        return self._elevation_score


def _retrieve_pam(sequence: str, length: int, right: bool) -> str:
    return sequence[:length] if right else sequence[-length:]


def _format_sequence(sequence: str, pam: str, right: bool) -> Tuple[str, str]:
    s_ = sequence[len(pam) :] if right else sequence[: -len(pam)]
    s = f"{pam}{s_}" if right else f"{s_}{pam}"
    return s_, s
