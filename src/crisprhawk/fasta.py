"""
"""

from crisprhawk_error import CrisprHawkFastaError
from exception_handlers import exception_handler
from bedfile import Coordinate

from typing import Optional

import pysam
import os

FAI = "fai"  # fasta index extension format

class Sequence:
    def __init__(self, sequence: str) -> None:
        self._sequence = sequence.upper()  # force sequence nucleotides to upper cases

    @property
    def sequence(self) -> str: return self._sequence

class Region(Sequence):
    def __init__(self, sequence, coord: Coordinate):
        super().__init__(sequence)
        self._contig = coord.contig  # store contig name from which the sequence was extracted
        self._start = coord.start  # store start position of sequence
        self._stop = coord.stop  # store stop position of sequence

    @property
    def contig(self) -> str: return self._contig
    @property
    def start(self) -> int: return self._start
    @property
    def stop(self) -> int: return self._stop

class Fasta:
    def __init__(self, fname: str, debug: bool, faidx: Optional[str] = "") -> None:
        self._fname = fname  # store input file name
        self._debug = debug  # store debug mode flag
        self._faidx = self._search_index(faidx)  # initialize fasta index
        self._fasta = pysam.FastaFile(self._fname, filepath_index=self._faidx)  # initialize FastaFile object with the previously computed index
        self._contigs = self._fasta.references

    def _search_index(self, faidx: Optional[str] = "") -> str:
        # look for index file for the current fasta file, if not found compute it
        if not faidx:  # index not provided from input arguments -> search it
            if _find_fai(self._fname):  # index found, return it
                return f"{self._fname}.{FAI}"
            # not found the index -> compute it de nove and store it in the same 
            # folder of the indexed fasta
            try:
                pysam.faidx(self._fname)  # index fasta using samtools
            except OSError as e:
                exception_handler(CrisprHawkFastaError, f"An error occurred while indexing {self._fname}", os.EX_SOFTWARE, self._debug, e)
        fastaindex = f"{self._fname}.{FAI}"
        assert _find_fai(fastaindex)
        return fastaindex
    
    def extract(self, coord: Coordinate) -> Sequence:
        if coord.contig not in self._contigs:  # conting not available in fasta -> raise error
            exception_handler(CrisprHawkFastaError, f"Input contig ({coord.contig}) not available in {self._fname}", os.EX_DATAERR, self._debug)
        try:  # extract sequence in the input range from fasta file
            return Sequence(self._fasta.fetch(coord.contig, coord.start, coord.stop).strip())
        except CrisprHawkFastaError as e:  # failed extraction
            exception_handler(CrisprHawkFastaError, f"Sequence extraction failed for coordinates ({coord})", os.EX_DATAERR, self._debug, e)
    

def _find_fai(fastafile: str) -> bool:
    # search index for the input fasta file, assumes that the index is located
    # in the same folder as the indexed fasta
    fastaindex = f"{os.path.abspath(fastafile)}.{FAI}"  # avoid unexpected crashes
    if os.path.exists(fastaindex):
        return os.path.isfile(fastaindex) and os.stat(fastaindex).st_size > 0
    return False

