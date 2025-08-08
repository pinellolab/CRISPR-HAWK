""" """

from .scores.cfdscore.cfdscore import compute_cfd
from .bedfile import BedAnnotation
from .guide import Guide

from typing import Dict

import numpy as np


class Offtarget:

    def __init__(
        self,
        chrom: str,
        start: int,
        stop: int,
        strand: str,
        sequence: str,
        pam: str,
        mm: int,
        debug: bool,
    ) -> None:
        self._debug = debug  # store debug mode
        self._chrom = chrom  # set chromosome
        self._start = start  # set start position
        self._stop = stop  # set stop position
        self._strand = strand  # set off-target strand
        self._sequence = sequence  # set off-target sequence
        self._pam = pam  # set off-target pam
        self._mm = mm  # set mismatches number
        self._cfd = np.nan  # default cfd score
        self._funcann = "NA"  # default functional annotation
        self._geneann = "NA"  # default gene annotation

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} object; start={self._start} stop={self._stop} sequence={self._sequence + self._pam} strand={self._strand}>"

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

    def annotate_functional(self, funcann: BedAnnotation) -> None:
        annotation = funcann.fetch_features(self._chrom, self._start, self._stop, 9)
        self._funcann = annotation

    def annotate_gene(self, geneann: BedAnnotation) -> None:
        annotation = geneann.fetch_features(self._chrom, self._start, self._stop, 22)
        self._geneann = annotation

    def report_line(self, guide: Guide) -> str:
        # align guide and spacer, reporting mismatching nts in lower case
        assert len(self._sequence) == len(guide.guide)
        spacer = "".join(
            [
                nt if nt == guide.guide[i] else nt.lower()
                for i, nt in enumerate(self._sequence)
            ]
        )
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
                        spacer,
                        self._pam,
                        self._mm,
                        self._cfd,
                        self._funcann,
                        self._geneann,
                    ],
                )
            )
        )

    @property
    def cfd(self) -> float:
        return self._cfd
