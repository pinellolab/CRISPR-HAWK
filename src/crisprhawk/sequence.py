""" """

from exception_handlers import exception_handler
from coordinate import Coordinate
from utils import IUPAC, warning

from typing import Optional, Union

import pysam
import sys
import os

FAI = "fai"  # fasta index extension format


class Sequence:
    """Represents a DNA sequence.

    Stores a DNA sequence and provides methods for accessing and manipulating it.
    The sequence is stored as a string and as a list of characters for efficient indexing.

    Attributes:
        _sequence: The DNA sequence as a string.
        _sequence_raw: The DNA sequence as a list of characters.
        _sequence_bits: The encoded DNA sequence as a list of bits.
    """

    def __init__(self, sequence: str, debug: bool, allow_lower_case: Optional[bool] = False) -> None:
        """Initialize a Sequence object with a DNA string.

        Validates and stores the DNA sequence, converting it to uppercase and preparing it for efficient access.

        Args:
            sequence: The DNA sequence as a string.
            debug: Flag to enable debug mode.

        Raises:
            ValueError: If the input string contains non-DNA characters.
        """
        self._debug = debug  # store debug flag
        sequence = sequence if allow_lower_case else sequence.upper()  # force sequence nucleotides to upper case
        if any(nt.upper() not in IUPAC for nt in set(sequence)):
            exception_handler(
                ValueError,
                "The input string is not a DNA string",
                os.EX_DATAERR,
                self._debug,
            )
        self._sequence = sequence  # store sequence
        # cast str to list for fast index access
        self._sequence_raw = list(self._sequence)

    def __eq__(self, sequence: "Sequence") -> bool:
        """Compare two Sequence objects for equality.

        Compares the DNA sequences of two Sequence objects.

        Args:
            sequence: The Sequence object to compare to.

        Returns:
            True if the sequences are identical, False otherwise.
        """
        return self._sequence == sequence.sequence

    def __len__(self) -> int:
        """Return the length of the sequence.

        Returns the length of the DNA sequence.

        Returns:
            The length of the sequence.
        """
        return len(self._sequence)

    def __str__(self) -> str:
        """Return the DNA sequence as a string.

        Returns the DNA sequence stored in the object.

        Returns:
            The DNA sequence.
        """
        return self._sequence

    def __iter__(self) -> "SequenceIterator":
        """Return an iterator over the nucleotides in the sequence.

        Returns an iterator that allows iterating over the nucleotides in the Sequence
        object.

        Returns:
            An iterator over the nucleotides.
        """
        return SequenceIterator(self)

    def __getitem__(self, idx: Union[int, slice]) -> str:
        """Return the nucleotide at the given index or slice.

        Returns the nucleotide(s) at the specified index or slice in the DNA sequence.

        Args:
            idx: The index or slice to retrieve.

        Returns:
            The nucleotide or a slice of the sequence.

        Raises:
            AttributeError: If the _sequence_raw attribute is missing.
            IndexError: If the index is out of range.
        """
        if not hasattr(self, "_sequence_raw"):  # always trace this error
            raise AttributeError(
                f"Missing _sequence_raw attribute on {self.__class__.__name__}"
            )
        try:
            return self._sequence_raw[idx]
        except IndexError as e:
            raise IndexError(f"Index {idx} out of range") from e

    @property
    def sequence(self) -> str:
        return self._sequence


class SequenceIterator:
    """Iterator for Sequence objects.

    Iterates over the nucleotides in a Sequence object.

    Attributes:
        sequence: The Sequence object to iterate over.
        _index: The current index of the iterator.
    """

    def __init__(self, sequence: Sequence) -> None:
        """Initialize the SequenceIterator.

        Args:
            sequence: The Sequence object to iterate over.

        Raises:
            AttributeError: If the sequence object is missing the _sequence_raw
            attribute.
        """
        if not hasattr(sequence, "_sequence_raw"):  # always trace this error
            raise AttributeError(
                f"Missing _sequence_raw attribute on {self.__class__.__name__}"
            )
        self._sequence = sequence  # sequence object to iterate over
        self._index = 0  # iterator index used over the sequence

    def __next__(self) -> str:
        """Return the next nucleotide in the sequence.

        Returns the next nucleotide in the Sequence object and advances the iterator.


        Returns:
            The next nucleotide in the sequence.

        Raises:
            StopIteration: If there are no more nucleotides to iterate over.
        """
        if self._index < len(self._sequence):
            result = self._sequence[self._index]
            self._index += 1  # go to next position in the sequence
            return result
        raise StopIteration  # stop iteration over sequence object


class Fasta:
    """Handles FASTA file operations.

    Provides methods for interacting with FASTA files, including indexing, fetching
    sequences, and managing associated data.

    Attributes:
        _fname (str): The path to the FASTA file.
        _faidx (str): The path to the FASTA index file.
        _fasta (pysam.FastaFile): The pysam FastaFile object.
        _contigs (list): A list of contig names in the FASTA file.
    """

    def __init__(
        self, fname: str, verbosity: int, debug: bool, faidx: Optional[str] = ""
    ) -> None:
        """Initialize a Fasta object for handling FASTA file operations.

        Sets up the Fasta object by validating the input FASTA file, searching for or creating an index, and initializing internal data structures.

        Args:
            fname: The path to the FASTA file.
            verbosity: The verbosity level for warnings and messages.
            debug: Flag to enable debug mode.
            faidx: Optional path to a FASTA index file.

        Raises:
            FileNotFoundError: If the FASTA file or provided index file does not exist.
        """
        self._debug = debug  # store debug flag
        self._verbosity = verbosity  # store verbosity level
        if not os.path.isfile(fname):
            exception_handler(
                FileNotFoundError,
                f"Cannot find input FASTA {fname}",
                os.EX_DATAERR,
                self._debug,
            )
        self._fname = fname  # store input file name
        self._faidx = self._search_index(faidx)  # initialize fasta index
        if not self._faidx:  # index not found compute it
            self.index_fasta()
        # initialize FastaFile object with the previously computed index
        # (class is a wrapper for pysam FastaFile class)
        self._fasta = pysam.FastaFile(self._fname, filepath_index=self._faidx)
        self._contigs = self._fasta.references  # add contig names

    def __repr__(self):
        """Return a string representation of the Fasta object.

        Returns a string representation of the Fasta object, including the class name.

        Returns:
            A string representation of the Fasta object.
        """
        return f"<{self.__class__.__name__} object>"

    def _search_index(self, faidx: Optional[str] = "") -> str:
        """Search for or validate a FASTA index.

        Searches for a FASTA index (.fai) for the associated FASTA file if one is not
        provided. If a path to an index is provided, it validates that the index
        exists and is not empty.

        Args:
            faidx: An optional path to a FASTA index file.

        Returns:
            The path to the FASTA index file, or an empty string if not found.

        Raises:
            FileNotFoundError: If the provided index file does not exist or is empty.
        """
        # look for index file for the current fasta file, if not found compute it
        if not faidx:  # index not provided from input arguments -> search it
            if _find_fai(self._fname):  # index found, return it
                return f"{self._fname}.{FAI}"
            # not found the index, print message and return empty string
            sys.stdout.write(f"FASTA index not found for {self._fname}\n")
            return ""
        # input fasta index index must not be empty
        if not (os.path.isfile(faidx) and os.stat(faidx).st_size > 0):
            exception_handler(
                FileNotFoundError,
                f"Not existing or empty FASTA index {faidx}",
                os.EX_DATAERR,
                self._debug,
            )
        return faidx

    def index_fasta(self) -> None:
        """Create or update the FASTA index.

        Creates or updates the FASTA index (.fai) for the associated FASTA file.
        If an index already exists, it will be overwritten.

        Raises:
            RuntimeError: If an error occurs during indexing.
        """
        if self._faidx:  # launch warning
            warning("FASTA index already present, forcing update", self._verbosity)
        # compute fasta index if not provided during object creation or found
        try:  # create index in the same folder of base fasta
            sys.stdout.write(f"Generating FASTA index...\n")
            pysam.faidx(self._fname)  # index fasta using samtools
        except OSError as e:
            exception_handler(
                RuntimeError,
                f"An error occurred while indexing {self._fname}",
                os.EX_SOFTWARE,
                self._debug,
                e,
            )
        assert _find_fai(self._fname)
        self._faidx = f"{self._fname}.{FAI}"

    def fetch(self, coord: Coordinate) -> Sequence:
        """Fetch the DNA sequence for the given coordinate.

        Fetches the DNA sequence from the FASTA file for the specified genomic
        coordinate.

        Args:
            coord: The genomic coordinate to fetch the sequence for.

        Returns:
            A Sequence object containing the fetched sequence.

        Raises:
            ValueError: If the contig is not found in the FASTA file or if sequence
                extraction fails.
        """
        if coord.contig not in self._contigs:  # conting not available in fasta
            exception_handler(
                ValueError,
                f"Input contig ({coord.contig}) not available in {self._fname}",
                os.EX_DATAERR,
                self._debug,
            )
        try:  # extract sequence in the input range from fasta file
            return Sequence(
                self._fasta.fetch(coord.contig, coord.start - 1, coord.stop).strip(),
                self._debug,
            )
        except ValueError as e:  # failed extraction
            exception_handler(
                ValueError,
                f"Sequence extraction failed for coordinates {coord.contig}:{coord.start}-{coord.stop}",
                os.EX_DATAERR,
                self._debug,
                e,
            )

    @property
    def fname(self) -> str:
        return self._fname


def _find_fai(fastafile: str) -> bool:
    """Check if a FASTA index file exists for the given FASTA file.

    Checks if a FASTA index file (.fai) exists for the given FASTA file in the same
    directory.

    Args:
        fastafile: The path to the FASTA file.

    Returns:
        True if the index file exists and is not empty, False otherwise.
    """
    # search index for the input fasta file, assumes that the index is located
    # in the same folder as the indexed fasta
    fastaindex = f"{os.path.abspath(fastafile)}.{FAI}"  # avoid unexpected crashes
    if os.path.exists(fastaindex):
        return os.path.isfile(fastaindex) and os.stat(fastaindex).st_size > 0
    return False
