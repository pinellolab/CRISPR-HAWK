"""
"""

from crisprhawk_error import CrisprHawkBitsetError
from exception_handlers import exception_handler

import os

SIZE = 4  # bit size


class Bitset(object):

    def __init__(self, size: int, debug: bool) -> None:
        if size < 1:
            exception_handler(
                CrisprHawkBitsetError,
                f"Forbidden Bitset size ({size})",
                os.EX_DATAERR,
                debug,
            )
        self._size = size  # number of bits to encode
        self._bits = 0  # bit initialized to 0
        self._debug = debug

    def __str__(self) -> str:
        # bin(self._bits) converts an integer to a binary string
        # [2:] remove the first two characters from the output of bin(self._bits)
        # .zfill(self._size) add as many zeros as required to reach self._size
        # Example: bits = 5, size =  8 bits  --> 101 --> 00000101 final
        return bin(self._bits)[2:].zfill(self._size)

    def __repr__(self) -> str:
        # value is the binary string, size is the len of the binary string
        return f"<{self.__class__.__name__} object; value={self}, size={self._size}>"

    def __and__(self, bitset: "Bitset") -> "Bitset":
        if self._size != bitset.size:
            exception_handler(
                CrisprHawkBitsetError,
                f"{self.__class__.__name__} objects must have the same size for AND operator",
                os.EX_DATAERR,
                self._debug,
            )
        result = Bitset(self._size, self._debug)  # allocate bits for AND result
        result._bits = self._bits & bitset.bits  # perform AND between bitsets
        return result

    @property
    def size(self) -> int:
        return self._size

    @property
    def bits(self) -> int:
        return self._bits

    def set(self, index: int) -> None:
        if index >= self._size:
            exception_handler(
                CrisprHawkBitsetError,
                f"Index {index} out of bounds, unable to set bit",
                os.EX_DATAERR,
                self._debug,
            )
        # bitwise OR operation -> sets 1 at the input position and shifts the
        # 1s to the left (<<)
        # EXAMPLE: if index = 3, 0000 -> 1000
        self._bits |= 1 << index

    def reset(self, index: int) -> None:
        if index >= self._size:
            exception_handler(
                CrisprHawkBitsetError,
                f"Index {index} out of bounds, unable to reset bit",
                os.EX_DATAERR,
                self._debug,
            )
        # reset bit at position index
        self._bits &= ~(1 << index)  # ~ is not operator

    def set_bits(self, bits: str) -> None:
        if any(bit not in "01" for bit in bits):
            exception_handler(
                CrisprHawkBitsetError,
                f"{bits} is not a bit string",
                os.EX_DATAERR,
                self._debug,
            )
        bitstring_size = len(bits)
        for i, bit in enumerate(bits):
            if bit == "0":  # force bit reset
                self.reset(bitstring_size - 1 - i)
            else:  # force set bit (bit == 1)
                self.set(bitstring_size - 1 - i)
        assert str(self) == bits

    def to_bool(self) -> bool:
        return bool(self._bits)  # cast bitset to bool

    def test(self, index: int) -> bool:
        if index >= self._size:
            exception_handler(
                CrisprHawkBitsetError,
                f"Index {index} out of bounds, unable to test bit",
                os.EX_DATAERR,
                self._debug,
            )
        return bool(self._bits & (1 << index))  # test if bit at position index
