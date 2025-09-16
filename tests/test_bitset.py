from crisprhawk.bitset import Bitset
from crisprhawk.crisprhawk_error import CrisprHawkBitsetError

import pytest

class DummyException(Exception):
    pass

def dummy_exception_handler(exc_type, message, code, debug):
    raise exc_type(message)

def test_bitset_init_valid(monkeypatch):
    # Patch exception_handler to raise immediately
    monkeypatch.setattr("crisprhawk.bitset.exception_handler", dummy_exception_handler)
    b = Bitset(4, debug=False)
    assert b.size == 4
    assert b.bits == 0
    assert str(b) == "0000"

def test_bitset_init_invalid(monkeypatch):
    monkeypatch.setattr("crisprhawk.bitset.exception_handler", dummy_exception_handler)
    with pytest.raises(CrisprHawkBitsetError):
        Bitset(0, debug=True)

def test_set_and_test(monkeypatch):  # sourcery skip: extract-duplicate-method
    monkeypatch.setattr("crisprhawk.bitset.exception_handler", dummy_exception_handler)
    b = Bitset(4, debug=False)
    b.set(2)
    assert b.test(2) is True
    assert str(b) == "0100"
    b.set(0)
    assert b.test(0) is True
    assert str(b) == "0101"

def test_set_out_of_bounds(monkeypatch):
    monkeypatch.setattr("crisprhawk.bitset.exception_handler", dummy_exception_handler)
    b = Bitset(4, debug=True)
    with pytest.raises(CrisprHawkBitsetError):
        b.set(4)

def test_reset_and_test(monkeypatch):
    monkeypatch.setattr("crisprhawk.bitset.exception_handler", dummy_exception_handler)
    b = Bitset(4, debug=False)
    b.set(1)
    b.set(2)
    b.reset(1)
    assert b.test(1) is False
    assert b.test(2) is True
    assert str(b) == "0100"

def test_reset_out_of_bounds(monkeypatch):
    monkeypatch.setattr("crisprhawk.bitset.exception_handler", dummy_exception_handler)
    b = Bitset(4, debug=True)
    with pytest.raises(CrisprHawkBitsetError):
        b.reset(4)

def test_set_bits(monkeypatch):
    monkeypatch.setattr("crisprhawk.bitset.exception_handler", dummy_exception_handler)
    b = Bitset(4, debug=False)
    b.set_bits("1010")
    assert str(b) == "1010"
    b.set_bits("0000")
    assert str(b) == "0000"

def test_set_bits_invalid(monkeypatch):
    monkeypatch.setattr("crisprhawk.bitset.exception_handler", dummy_exception_handler)
    b = Bitset(4, debug=True)
    with pytest.raises(CrisprHawkBitsetError):
        b.set_bits("10a0")

def test_and_operator(monkeypatch):
    monkeypatch.setattr("crisprhawk.bitset.exception_handler", dummy_exception_handler)
    b1 = Bitset(4, debug=False)
    b2 = Bitset(4, debug=False)
    b1.set_bits("1100")
    b2.set_bits("1010")
    b3 = b1 & b2
    assert isinstance(b3, Bitset)
    assert str(b3) == "1000"

def test_and_operator_size_mismatch(monkeypatch):
    monkeypatch.setattr("crisprhawk.bitset.exception_handler", dummy_exception_handler)
    b1 = Bitset(4, debug=True)
    b2 = Bitset(3, debug=True)
    with pytest.raises(CrisprHawkBitsetError):
        _ = b1 & b2

def test_to_bool(monkeypatch):
    monkeypatch.setattr("crisprhawk.bitset.exception_handler", dummy_exception_handler)
    b = Bitset(4, debug=False)
    assert b.to_bool() is False
    b.set(1)
    assert b.to_bool() is True

def test_str_and_repr(monkeypatch):
    monkeypatch.setattr("crisprhawk.bitset.exception_handler", dummy_exception_handler)
    b = Bitset(4, debug=False)
    b.set_bits("1010")
    assert str(b) == "1010"
    assert "Bitset" in repr(b)
    assert "value=1010" in repr(b)
    assert "size=4" in repr(b)