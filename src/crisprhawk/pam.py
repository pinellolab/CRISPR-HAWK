"""Provides the PAM class for representing and encoding Protospacer Adjacent 
Motif sequences.

This module defines the PAM class, which validates, stores, and encodes PAM 
sequences and their reverse complements for efficient sequence matching.
"""

from .exception_handlers import exception_handler
from .crisprhawk_error import CrisprHawkPamError
from .utils import reverse_complement, IUPAC
from .encoder import encode
from .bitset import Bitset

from typing import Optional, List

import os


class PAM:
    def __init__(self, pamseq: str, debug: bool):
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
                ValueError, f"Invalid PAM sequence {pamseq}", os.EX_DATAERR, self._debug # type: ignore
            )
        self._sequence = pamseq.upper()  # store pam sequence
        self._sequence_rc = reverse_complement(pamseq, debug)  # reverse complement

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

    def encode(self, verbosity: int) -> None:
        """Encodes the PAM sequence and its reverse complement into bit 
        representations.

        This method prepares the PAM object for efficient sequence matching by 
        encoding both the forward and reverse complement sequences.

        Args:
            verbosity: The verbosity level for logging.

        Raises:
            CrisprHawkPamError: If encoding the PAM sequence fails.
        """
        try:  # encode in bit fwd and rev pam sequence
            self._sequence_bits = encode(self._sequence, verbosity, self._debug)
            self._sequence_rc_bits = encode(self._sequence_rc, verbosity, self._debug)
        except ValueError as e:
            exception_handler(
                CrisprHawkPamError, # type: ignore
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
    def bits(self) -> List[Bitset]:
        if not hasattr(self, "_sequence_bits"):  # always trace these errors
            exception_handler(
                AttributeError, # type: ignore
                f"Missing _sequence_bits attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return self._sequence_bits

    @property
    def bitsrc(self) -> List[Bitset]:
        if not hasattr(self, "_sequence_rc_bits"):  # always trace these errors
            exception_handler(
                AttributeError, # type: ignore
                f"Missing _sequence_rc_bits attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return self._sequence_rc_bits
