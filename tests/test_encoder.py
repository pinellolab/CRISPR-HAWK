from crisprhawk import encoder

import pytest


class DummyBitset:
    def __init__(self, size, debug):
        self.size = size
        self.debug = debug
        self.bits = None
        self.set_calls = []
        self.set_bits_calls = []

    def set(self, idx):
        self.set_calls.append(idx)
        self.bits = idx

    def set_bits(self, bits):
        self.set_bits_calls.append(bits)
        self.bits = bits


def test_encoder_valid_nucleotides(monkeypatch):
    # Patch Bitset to DummyBitset for test
    monkeypatch.setattr(encoder, "Bitset", DummyBitset)
    monkeypatch.setattr(encoder, "SIZE", 4)
    monkeypatch.setattr(encoder, "IUPAC", "ACGTNRYSWKMBDHV")
    # Test all valid IUPAC codes
    for idx, nt in enumerate(encoder.IUPAC):
        b = encoder._encoder(nt, 0, False)
        assert isinstance(b, DummyBitset)
        # Check that set or set_bits was called
        if idx < 4:
            assert b.set_calls or b.set_bits_calls
        else:
            assert b.set_bits_calls


def test_encoder_invalid_nucleotide(monkeypatch):
    monkeypatch.setattr(encoder, "Bitset", DummyBitset)
    monkeypatch.setattr(encoder, "SIZE", 4)
    monkeypatch.setattr(encoder, "IUPAC", "ACGTNRYSWKMBDHV")

    # Patch exception_handler to raise
    def raise_exc(*args, **kwargs):
        raise Exception("Invalid nucleotide")

    monkeypatch.setattr(encoder, "exception_handler", raise_exc)
    with pytest.raises(Exception):
        encoder._encoder("Z", 0, False)


def test_encode(monkeypatch):
    # Patch Bitset to DummyBitset for test
    monkeypatch.setattr(encoder, "Bitset", DummyBitset)
    monkeypatch.setattr(encoder, "SIZE", 4)
    monkeypatch.setattr(encoder, "IUPAC", "ACGTNRYSWKMBDHV")
    monkeypatch.setattr(encoder, "print_verbosity", lambda *a, **k: None)
    monkeypatch.setattr(encoder, "VERBOSITYLVL", [0, 1, 2, 3])
    seq = "ACGT"
    bits = encoder.encode(seq, verbosity=0, debug=False)
    assert len(bits) == 4
    for b in bits:
        assert isinstance(b, DummyBitset)
