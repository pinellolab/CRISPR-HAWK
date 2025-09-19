from crisprhawk.pam import PAM, CASXPAM, CPF1PAM, SACAS9PAM, SPCAS9PAM, XCAS9PAM, CASX, CPF1, SACAS9, SPCAS9, XCAS9
from crisprhawk.crisprhawk_error import CrisprHawkPamError

import pytest


def test_valid_pam_initialization():
    pam = PAM("NGG", right=False, debug=False)
    assert pam.pam == "NGG"
    assert isinstance(pam.pamrc, str)
    assert pam.cas_system == SPCAS9

def test_invalid_pam_raises():
    with pytest.raises(ValueError):
        PAM("XYZ", right=False, debug=True)

def test_pam_equality():
    pam1 = PAM("NGG", right=False, debug=False)
    pam2 = PAM("NGG", right=False, debug=False)
    pam3 = PAM("NGA", right=False, debug=False)
    assert pam1 == pam2
    assert pam1 != pam3

def test_pam_repr_and_str():
    pam = PAM("NGG", right=False, debug=False)
    assert repr(pam).startswith("<PAM object; sequence=NGG>")
    assert str(pam) == "NGG"

def test_pam_encoding_and_bits():
    pam = PAM("NGG", right=False, debug=False)
    pam.encode(verbosity=0)
    bits = pam.bits
    bitsrc = pam.bitsrc
    assert isinstance(bits, list)
    assert isinstance(bitsrc, list)
    assert len(bits) == len(pam.pam)
    assert len(bitsrc) == len(pam.pamrc)

def test_cas_system_detection():
    assert PAM("TTCN", right=False, debug=False).cas_system == CASX
    assert PAM("TTN", right=True, debug=False).cas_system == CPF1
    assert PAM("NNGRRT", right=False, debug=False).cas_system == SACAS9
    assert PAM("NGG", right=False, debug=False).cas_system == SPCAS9
    assert PAM("NGK", right=False, debug=False).cas_system == XCAS9

def test_bits_property_error():
    pam = PAM("NGG", right=False, debug=False)
    with pytest.raises(AttributeError):
        _ = pam.bits

def test_bitsrc_property_error():
    pam = PAM("NGG", right=False, debug=False)
    with pytest.raises(AttributeError):
        _ = pam.bitsrc
