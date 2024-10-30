"""
"""

from crisprhawk_error import (
    CrisprHawkFastaError,
    CrisprHawkPamError,
    CrisprHawkIupacTableError,
)
from exception_handlers import exception_handler
from utils import IUPAC, IUPACTABLE, reverse_complement
from bitset import Bitset, SIZE

from typing import Optional, Union, List, Tuple
from itertools import product

import pysam
import os

FAI = "fai"  # fasta index extension format


class Sequence:
    def __init__(self, sequence: str, debug: bool) -> None:
        self._sequence = sequence.upper()  # force sequence nucleotides to upper cases
        # cast str to list for fast index access
        self._sequence_raw = list(self._sequence)
        self._debug = debug  # store debug mode flag

    def __len__(self) -> int:
        return len(self._sequence)

    def __str__(self) -> str:
        return self._sequence

    def __getitem__(self, idx: Union[int, slice]) -> str:
        if not hasattr(self, "_sequence_raw"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _sequence_raw attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        try:
            return self._sequence_raw[idx]
        except IndexError as e:
            exception_handler(
                CrisprHawkFastaError,
                f"Index {idx} out of range",
                os.EX_DATAERR,
                self._debug,
                e,
            )

    def encode(self) -> None:
        # encode sequence in bits for efficient matching
        self._sequence_bits = [
            _encoder(nt, i, self._debug) for i, nt in enumerate(self._sequence_raw)
        ]
        assert len(self._sequence_bits) == len(self)

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def bits(self) -> List[Bitset]:
        return self._sequence_bits


class SequenceIterator:
    def __init__(self, sequence: Sequence) -> None:
        if not hasattr(self, "_sequence_raw"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _sequence_raw attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        self._sequence = sequence  # sequence object to iterate over
        self._index = 0  # iterator index used over the sequence

    def __next__(self) -> str:
        if self._index < len(self._sequence):
            self._index += 1  # go to next position in the sequence
            return self._sequence[self._index - 1]
        raise StopIteration  # stop iteration over sequence object


class PAM:
    def __init__(self, pamseq: str, debug: bool):
        self._sequence = Sequence(pamseq, debug)  # initialize forward pam
        self._sequence_rc = Sequence(
            reverse_complement(pamseq, debug), debug
        )  # initialize pam rc
        self._debug = debug  # store debug mode flag

    def __len__(self) -> int:
        return len(self._sequence)

    def __str__(self) -> str:
        return f"{self._sequence}"

    def encode(self) -> None:
        try:  # encode in bit fwd and rev pam sequence
            self._sequence_rc.encode()  # reverse strand
            self._sequence.encode()  # forward strand
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
    def bits(self) -> Tuple[str]:
        if not hasattr(self.pam, "_sequence_bits"):  # always trace this erro
            exception_handler(
                AttributeError,
                f"Missing _sequence_bits attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return (self._sequence.bits, self._sequence_rc.bits)  # return arrays of bits


class Fasta:
    def __init__(self, fname: str, debug: bool, faidx: Optional[str] = "") -> None:
        self._fname = fname  # store input file name
        self._debug = debug  # store debug mode flag
        self._faidx = self._search_index(faidx)  # initialize fasta index
        self._fasta = pysam.FastaFile(
            self._fname, filepath_index=self._faidx
        )  # initialize FastaFile object with the previously computed index
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
                exception_handler(
                    CrisprHawkFastaError,
                    f"An error occurred while indexing {self._fname}",
                    os.EX_SOFTWARE,
                    self._debug,
                    e,
                )
        assert _find_fai(self._fname)
        return f"{self._fname}.{FAI}"

    def fetch(self, contig: str, start: int, stop: int) -> str:
        if contig not in self._contigs:  # conting not available in fasta
            exception_handler(
                CrisprHawkFastaError,
                f"Input contig ({contig}) not available in {self._fname}",
                os.EX_DATAERR,
                self._debug,
            )
        try:  # extract sequence in the input range from fasta file
            return self._fasta.fetch(contig, start, stop).strip()
        except ValueError as e:  # failed extraction
            exception_handler(
                CrisprHawkFastaError,
                f"Sequence extraction failed for coordinates ({contig}:{start}-{stop})",
                os.EX_DATAERR,
                self._debug,
                e,
            )


def _find_fai(fastafile: str) -> bool:
    # search index for the input fasta file, assumes that the index is located
    # in the same folder as the indexed fasta
    fastaindex = f"{os.path.abspath(fastafile)}.{FAI}"  # avoid unexpected crashes
    if os.path.exists(fastaindex):
        return os.path.isfile(fastaindex) and os.stat(fastaindex).st_size > 0
    return False


def _encoder(nt: str, position: int, debug: bool) -> Bitset:
    bitset = Bitset(SIZE, debug)  # 4 - bits encoder
    if nt == IUPAC[0]:  # A - 0001
        bitset.set(0)
    elif nt == IUPAC[1]:  # C - 0010
        bitset.set(1)
    elif nt == IUPAC[2]:  # G - 0100
        bitset.set(2)
    elif nt == IUPAC[3]:  # T - 1000
        bitset.set(3)
    elif nt == IUPAC[4]:  # R - 0101 G>A
        bitset.set_bits("0101")
    elif nt == IUPAC[5]:  # Y - 1010 C>T
        bitset.set_bits("1010")
    elif nt == IUPAC[6]:  # S - 0110 C>G
        bitset.set_bits("0110")
    elif nt == IUPAC[7]:  # W  - 1001
        bitset.set_bits("1001")
    elif nt == IUPAC[8]:  # K - 1100
        bitset.set_bits("1100")
    elif nt == IUPAC[9]:  # M - 0011
        bitset.set_bits("0011")
    elif nt == IUPAC[10]:  # B - 1110 --> not A ( T or G or C)
        bitset.set_bits("1110")
    elif nt == IUPAC[11]:  # D - 1101 --> not C
        bitset.set_bits("1101")
    elif nt == IUPAC[12]:  # H - 1011 --> not G
        bitset.set_bits("1011")
    elif nt == IUPAC[13]:  # V - 0111 --> not T
        bitset.set_bits("0111")
    elif nt == IUPAC[14]:  # N - 1111 --> any
        bitset.set_bits("1111")
    else:  # default case
        exception_handler(
            CrisprHawkFastaError,
            f"The nucleotide {nt} at {position} is not a IUPAC character",
        )
    return bitset


def explode_iupac_sequence(iupac_sequence: List[str], debug) -> List[str]:
    try:  # decode iupac nucleotides using iupac table
        sequence_combinations = [IUPACTABLE[nt] for nt in iupac_sequence]
        # explode the input iupac sequence in all their combinations
        return ["".join(s) for s in product(*sequence_combinations)]
    except ValueError as e:
        exception_handler(
            CrisprHawkIupacTableError,
            "IUPAC sequence decoding failed",
            os.EX_DATAERR,
            debug,
            e,
        )
