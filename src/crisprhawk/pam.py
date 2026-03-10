"""Provides the PAM class for representing and encoding Protospacer Adjacent
Motif sequences.

This module defines the PAM class, which validates, stores, and encodes PAM
sequences and their reverse complements for efficient sequence matching.
"""

from .exception_handlers import exception_handler
from .crisprhawk_error import CrisprHawkPamError
from .utils import reverse_complement, IUPAC
from .encoder import encode

from typing import List

import os

# list PAMs for each cas system
CASXPAM = ["TTCN"]
CPF1PAM = [
    "TTN",
    "TTTN",
    "TYCV",
    "TATV",
    "TTTV",
    "TTTR",
    "ATTN",
    "TTTA",
    "TCTA",
    "TCCA",
    "CCCA",
    "YTTV",
    "TTYN",
]
SACAS9PAM = ["NNGRRT", "NNNRRT"]
SPCAS9PAM = ["NGG", "NGA", "NRG", "NGC"]
XCAS9PAM = ["NGK", "NGN", "NNG"]

# list cas systems
CASX = 0
CPF1 = 1
SACAS9 = 2
SPCAS9 = 3
XCAS9 = 4


class PAM:
    def __init__(self, pamseq: str, right: bool, debug: bool):
        """Initializes a PAM object with a given sequence and debug mode.

        This method validates the PAM sequence, stores it in uppercase, and
        computes its reverse complement.

        Args:
            pamseq: The PAM sequence to be represented.
            debug: Whether to enable debug mode for error handling.
        """
        self._debug = debug  # store debig mode flag
        if any(nt.upper() not in IUPAC for nt in pamseq):
            exception_handler(
                ValueError, f"Invalid PAM sequence {pamseq}", os.EX_DATAERR, self._debug  # type: ignore
            )
        self._sequence = pamseq.upper()  # store pam sequence
        self._sequence_rc = reverse_complement(pamseq, debug)  # reverse complement
        self._assess_cas_system(right)  # assess pam's cas system

    def __len__(self) -> int:
        """Returns the length of the PAM sequence.

        This method allows the PAM object to be used with the built-in len()
        function.

        Returns:
            int: The length of the PAM sequence.
        """
        return len(self._sequence)

    def __eq__(self, pam: object) -> bool:
        """Checks equality between this PAM object and another.

        Compares the stored PAM sequence with another PAM object's sequence to
        determine equality.

        Args:
            pam: The object to compare with this PAM instance.

        Returns:
            bool: True if the sequences are equal and the object is a PAM instance,
                False otherwise.
        """
        return self._sequence == pam.pam if isinstance(pam, PAM) else NotImplemented

    def __repr__(self) -> str:
        """Returns a string representation of the PAM object for debugging.

        This method provides a detailed string useful for developers to inspect
        the PAM object.

        Returns:
            str: A string representation of the PAM object.
        """
        return f"<{self.__class__.__name__} object; sequence={self._sequence}>"

    def __str__(self) -> str:
        """Returns the PAM sequence as a string.

        This method allows the PAM object to be converted to its sequence string
        representation.

        Returns:
            str: The PAM sequence.
        """
        return f"{self._sequence}"

    def _assess_cas_system(self, right: bool) -> None:
        self._cas_system = -1  # unknown cas system pam
        if self._sequence in CASXPAM:  # casx system pam
            self._cas_system = CASX
        elif self._sequence in CPF1PAM and right:  # cpf1 cas system pam
            self._cas_system = CPF1
        elif self._sequence in SACAS9PAM:  # sacas9 system pam
            self._cas_system = SACAS9
        elif self._sequence in SPCAS9PAM and not right:  # spcas9 system pam
            self._cas_system = SPCAS9
        elif self._sequence in XCAS9PAM and not right:  # xcas9 pam
            self._cas_system = XCAS9

    def encode(self, verbosity: int) -> None:
        try:  # encode in bit fwd and rev pam sequence
            self._sequence_bits = encode(self._sequence, verbosity, self._debug)
            self._sequence_rc_bits = encode(self._sequence_rc, verbosity, self._debug)
            # packpam bit-based representation
            self._packed_bits = _pack_bits(self._sequence_bits)
            self._packed_bitsrc = _pack_bits(self._sequence_rc_bits)
        except ValueError as e:
            exception_handler(
                CrisprHawkPamError,  # type: ignore
                "PAM bit encoding failed",
                os.EX_DATAERR,
                self._debug,
                e,
            )

    @property
    def pam(self) -> str:

        return self._sequence

    @property
    def pamrc(self) -> str:
        return self._sequence_rc

    @property
    def bits(self) -> int:
        return self._packed_bits

    @property
    def bitsrc(self) -> int:
        return self._packed_bitsrc
    
    @property
    def bits_list(self) -> List[int]:
        return self._sequence_bits

    @property
    def cas_system(self) -> int:
        return self._cas_system


def _pack_bits(bits: List[int]) -> int:
    packed = 0
    for b in bits:
        packed = (packed << 4) | b
    return packed