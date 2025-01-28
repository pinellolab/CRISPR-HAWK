""" 
"""

from crisprhawk_error import CrisprHawkVariantMapError
from exception_handlers import exception_handler
from variants import VariantRecord

import os

class VariantMap:
    def __init__(self, phased: bool, debug: bool) -> None:
        self._map = {}  # variants hashmap (pos->var)
        self._positions = []  # list of variant positions
        self._variants = []  # list of hashed variants
        self._phased = phased  # variants phasing
        self._debug = debug  # debug mode

    def __repr__(self) -> str:
        values = ", ".join([f"{pos}: {var}" for pos, var in self._map.items()])
        return f"<{self.__class__.__name__} object; phased={self._phased}; [{values}]>"
    
    def __getitem__(self, position: int) -> VariantRecord:
        try:
            return self._map[position]
        except KeyError as e:
            exception_handler(CrisprHawkVariantMapError, f"Variants not found at position {position} in current map", os.EX_DATAERR, self._debug, e)

    def __contains__(self, position: int) -> bool:
        return position in self._positions

    def _update_map(self) -> None:
        self._positions = list(self._map.keys())
        self._variants = list(self._map.values())
        assert len(self._positions) == len(self._variants)

    def insert_variant(self, position: int, variant: VariantRecord) -> None:
        if position in self._positions:
            exception_handler(KeyError, f"Position {position} already present in variants hashmap", os.EX_DATAERR, self._debug)
        self._map[position] = variant  # insert variants
        self._update_map()  # update variants map info

    @property
    def phased(self) -> bool:
        return self._phased