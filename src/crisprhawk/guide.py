"""
"""

from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkGuideError
from utils import RC

from typing import List, Union
import os

class Guide:
    def __init__(self, position: int, sequence: List[str], direction: int, debug: bool) -> None:
        self._position = position  # guide position
        self._sequence_raw = sequence  # sequence as list of chars
        self._sequence = "".join(sequence)  # sequence as string
        self._direction = direction  # guide direction 
        self._debug = debug  # store debug mode
        self._samples = "NA"  # samples carrying guide variants
        self._variants = "NA"  # variants overlapping guide

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} object; position={self._position} sequence={self._sequence} direction={self._direction}>"
        
    def __len__(self) -> int:
        return len(self._sequence_raw)
    
    def __getitem__(self, idx: Union[int, slice]) -> str:
        assert hasattr(self, "_sequence")
        try:
            return "".join(self._sequence_raw[idx])
        except IndexError as e:
            exception_handler(
                CrisprHawkGuideError,
                f"Index {idx} out of range",
                os.EX_DATAERR,
                self._debug,
                e
            )
    
    def reverse_complement(self) -> None:
        assert hasattr(self, "_sequence")
        self._sequence_raw = [RC[nt] for nt in self._sequence[::-1]]
        self._sequence = "".join(self._sequence_raw)  # guide reverse complement

    
    @property
    def position(self) -> int:
        return self._position
    
    @property
    def strand(self) -> int:
        return self._direction
    
    @property
    def samples(self) -> str:
        return self._samples

    @property
    def variants(self) -> str:
        return self._variants