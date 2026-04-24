import pytest
import os
from crisprhawk.pam import (
    PAM,
    CASX,
    CPF1,
    SACAS9,
    SPCAS9,
    XCAS9,
    CASXPAM,
    CPF1PAM,
    SACAS9PAM,
    SPCAS9PAM,
    XCAS9PAM,
)
from crisprhawk.crisprhawk_error import CrisprHawkPamError


def test_valid_pam_initialization():
    """Test standard SpCas9 PAM initialization."""
    pam = PAM("NGG", right=False, debug=False)
    assert pam.pam == "NGG"
    assert pam.pamrc == "CCN"
    assert pam.cas_system == SPCAS9


def test_cas_system_detection():
    """Verify that different PAMs map to the correct Cas systems."""
    # Note: Logic in pam.py uses 'right' parameter for orientation
    assert PAM("TTCN", right=False, debug=False).cas_system == CASX
    assert PAM("TTN", right=True, debug=False).cas_system == CPF1
    assert PAM("NNGRRT", right=False, debug=False).cas_system == SACAS9
    assert PAM("NGG", right=False, debug=False).cas_system == SPCAS9
    assert PAM("NGK", right=False, debug=False).cas_system == XCAS9


def test_invalid_pam_sequence(monkeypatch):
    """Verify that invalid nucleotides trigger an exception."""
    # Since exception_handler might exit, we mock it to raise for testing
    from crisprhawk import pam

    monkeypatch.setattr(
        pam, "exception_handler", lambda *args: pytest.fail("Exception handler called")
    )

    # Alternatively, if exception_handler is expected to raise the error:
    def mock_handler(exc_type, msg, *args):
        raise exc_type(msg)

    monkeypatch.setattr(pam, "exception_handler", mock_handler)

    with pytest.raises(ValueError):
        PAM("XYZ", right=False, debug=True)


def test_pam_encoding_and_integers():
    """Verify that encode() produces integers and not Bitset objects."""
    pam = PAM("NGG", right=False, debug=False)
    pam.encode(verbosity=0)

    assert isinstance(pam.bits, int)
    assert isinstance(pam.bitsrc, int)
    assert isinstance(pam.bits_list, list)
    assert all(isinstance(b, int) for b in pam.bits_list)


def test_pam_equality_and_string():
    """Test utility methods of the PAM class."""
    pam1 = PAM("NGG", right=False, debug=False)
    pam2 = PAM("ngg", right=False, debug=False)
    assert pam1 == pam2
    assert str(pam1) == "NGG"
    assert "NGG" in repr(pam1)
    assert len(pam1) == 3
