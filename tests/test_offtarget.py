from crisprhawk.offtarget import Offtarget, _retrieve_pam, _format_sequence

import pytest


class DummyOfftarget(Offtarget):
    def __init__(self, reportline, pam, right, debug):
        # Avoid calling parent __init__ to skip actual parsing
        self._debug = debug
        self._pam = pam
        self._chrom = "chr1"
        self._pos = 12345
        self._strand = "+"
        self._grna_ = "ACGTACGTACGT"
        self._grna = "ACGTACGTACGT"
        self._spacer_ = "ACGTACGT"
        self._spacer = "ACGTACGT"
        self._mm = 2
        self._bulge_type = "DNA"
        self._bulge_size = 1
        self._cfd_score = "0.85"
        self._elevation_score = "0.90"


def test_retrieve_pam_right():
    seq = "AGGTCGATCGG"
    pam = _retrieve_pam(seq, 3, True)
    assert pam == seq[:3]


def test_retrieve_pam_left():
    seq = "AGGTCGATCGG"
    pam = _retrieve_pam(seq, 3, False)
    assert pam == seq[-3:]


def test_format_sequence_right():
    seq = "AGGTCGATCGG"
    pam = "AGG"
    spacer, formatted = _format_sequence(seq, pam, True)
    assert formatted.startswith(pam)
    assert spacer == seq[len(pam) :]


def test_format_sequence_left():
    seq = "AGGTCGATCGG"
    pam = "CGG"
    spacer, formatted = _format_sequence(seq, pam, False)
    assert formatted.endswith(pam)
    assert spacer == seq[: -len(pam)]


def test_offtarget_repr():
    ot = DummyOfftarget("", "AGG", True, False)
    rep = repr(ot)
    assert "position=12345" in rep
    assert "spacer=ACGTACGT" in rep
    assert "strand=+" in rep


def test_offtarget_report_line():
    ot = DummyOfftarget("", "AGG", True, False)
    line = ot.report_line()
    assert isinstance(line, str)
    assert "chr1" in line
    assert "12345" in line
    assert "ACGTACGT" in line
    assert "0.85" in line
    assert "0.90" in line


def test_offtarget_elevation_setter():
    ot = DummyOfftarget("", "AGG", True, False)
    ot.elevation = 0.75
    assert ot.elevation == "0.75"


def test_offtarget_elevation_setter_nan():
    import numpy as np

    ot = DummyOfftarget("", "AGG", True, False)
    ot.elevation = float("nan")
    assert ot.elevation == "NA"


def test_offtarget_elevation_setter_typeerror():
    ot = DummyOfftarget("", "AGG", True, False)
    with pytest.raises(SystemExit):
        ot.elevation = "not_a_float"
