from crisprhawk import encoder
import pytest
from unittest.mock import MagicMock


def test_encoder_valid_nucleotides():
    """Verify that _encoder returns correct bitmasks for standard nucleotides."""
    # A: 0b0001, C: 0b0010, G: 0b0100, T: 0b1000, N: 0b1111
    assert encoder._encoder("A", 0, False) == 1
    assert encoder._encoder("C", 1, False) == 2
    assert encoder._encoder("G", 2, False) == 4
    assert encoder._encoder("T", 3, False) == 8
    assert encoder._encoder("N", 4, False) == 15


def test_encoder_all_iupac_characters():
    """Verify all IUPAC characters defined in the module return expected integer values."""
    from crisprhawk.encoder import _IUPAC_BITS

    for nt, bitmask in _IUPAC_BITS.items():
        assert encoder._encoder(nt, 0, False) == bitmask


def test_encoder_invalid_nucleotide(monkeypatch):
    """Verify that an invalid character triggers the exception_handler."""
    # Mock exception_handler to raise an exception for testing purposes
    mock_handler = MagicMock(side_effect=ValueError("Mocked Error"))
    monkeypatch.setattr(encoder, "exception_handler", mock_handler)

    with pytest.raises(ValueError, match="Mocked Error"):
        encoder._encoder("Z", 0, False)

    mock_handler.assert_called_once()


def test_encode_full_sequence():
    """Verify the high-level encode function processes sequences correctly."""
    sequence = "ACGTN"
    # Expected: [1, 2, 4, 8, 15]
    result = encoder.encode(sequence, verbosity=0, debug=False)

    assert result == [1, 2, 4, 8, 15]
    assert len(result) == len(sequence)


def test_encode_case_insensitivity():
    """Verify that the encode function handles lowercase input by converting it to uppercase."""
    sequence = "acgtn"
    result = encoder.encode(sequence, verbosity=0, debug=False)
    assert result == [1, 2, 4, 8, 15]


def test_encode_empty_sequence():
    """Verify that encoding an empty string returns an empty list."""
    assert encoder.encode("", verbosity=0, debug=False) == []
