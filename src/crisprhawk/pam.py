"""
"""

from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkPamError
from utils import reverse_complement, IUPAC
from encoder import encode
from bitset import Bitset

from typing import Optional, List

import os

class PAM:
    def __init__(self, pamseq: str, debug: bool):
        self._debug = debug  # store debig mode flag
        if any(nt.upper() not in IUPAC for nt in pamseq): 
            exception_handler(ValueError, f"Invalid PAM sequence {pamseq}", os.EX_DATAERR, self._debug)
        self._sequence = pamseq.upper()  # store pam sequence
        self._sequence_rc = reverse_complement(pamseq, debug)  # reverse complement

    def __len__(self) -> int:
        return len(self._sequence)
    
    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} object; sequence={self._sequence}>"

    def __str__(self) -> str:
        return f"{self._sequence}"

    def encode(self, verbosity: Optional[int] = 0) -> None:
        try:  # encode in bit fwd and rev pam sequence
            self._sequence_bits = encode(self._sequence, verbosity, self._debug)  
            self._sequence_rc_bits = encode(self._sequence_rc, verbosity, self._debug)
        except ValueError as e:
            exception_handler(
                CrisprHawkPamError,
                f"PAM bit encoding failed",
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
            exception_handler(AttributeError, f"Missing _sequence_bits attribute on {self.__class__.__name__}", os.EX_DATAERR, True)
        return self._sequence_bits
    
    @property
    def bitsrc(self) -> List[Bitset]:
        if not hasattr(self, "_sequence_rc_bits"):  # always trace these errors
            exception_handler(AttributeError, f"Missing _sequence_rc_bits attribute on {self.__class__.__name__}", os.EX_DATAERR, True)
        return self._sequence_rc_bits
