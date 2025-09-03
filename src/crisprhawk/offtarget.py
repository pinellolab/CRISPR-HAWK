""" """

from .exception_handlers import exception_handler
from .scores.cfdscore.cfdscore import compute_cfd
from .bedfile import BedAnnotation
from .guide import Guide
from .utils import round_score

from typing import Dict

import numpy as np

import os


class Offtarget:

    def __init__(self, reportline: str, pam: str, debug: bool) -> None:
        self._debug = debug  # store debug mode
        self._parse_reportline(reportline)  # read report line content
        self._pam = pam  # set off-target pam
        self._cfd_score = "NA"  # default cfd score
        self._elevation_score = "NA"  # default elevation score        

    def __repr__(self) -> str:
        return (
            f"<{self.__class__.__name__} object; position={self._pos} "
            f"spacer={self._spacer} strand={self._strand}>"
        )


    def _parse_reportline(self, line: str) -> None:
        fields = line.strip().split()
        self._chrom = fields[3]  # set chromosome
        self._pos = int(fields[4])  # set off-target position
        self._strand = fields[6]  # set off-target strand
        self._grna = fields[1]  # set grna sequence
        self._spacer = fields[2]  # set spacer sequence
        self._mm = int(fields[7])  # set mismatches number
        self._bulge_size = int(fields[8])  # set bulge size


    def compute_cfd(
        self, wildtype: Guide, mmscores: Dict[str, float], pamscores: Dict[str, float]
    ) -> None:
        self._cfd = compute_cfd(
            wildtype.guide,
            self._sequence,
            self._pam[-2:],
            mmscores,
            pamscores,
            self._debug,
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
            self._elevationon_score = str(round_score(score))


    def annotate_functional(self, funcann: BedAnnotation) -> None:
        annotation = funcann.fetch_features(self._chrom, self._start, self._stop, 9)
        self._funcann = annotation

    def annotate_gene(self, geneann: BedAnnotation) -> None:
        annotation = geneann.fetch_features(self._chrom, self._start, self._stop, 22)
        self._geneann = annotation

    def report_line(self, guide: Guide) -> str:
        return "\t".join(
            list(
                map(
                    str,
                    [
                        guide.guidepam,
                        self._chrom,
                        self._start,
                        self._stop,
                        self._strand,
                        # spacer,
                        self._pam,
                        self._mm,
                        self._cfd,
                        self._elevation,
                    ],
                )
            )
        )

    @property
    def cfd(self) -> float:
        return self._cfd
