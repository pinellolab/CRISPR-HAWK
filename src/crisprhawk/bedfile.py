"""Provides classes and functions for working with BED files.

This module defines the Bed class for representing genomic regions from a BED file,
along with associated iterator and parsing functions. It allows for reading,
manipulating, and extracting regions from BED files, and integrating with FASTA
files for sequence retrieval.
"""

from .coordinate import Coordinate
from .sequence import Fasta
from .exception_handlers import exception_handler
from .region import Region, RegionList
from .utils import warning


from typing import List, Union, Optional, Dict
from pysam import TabixFile, tabix_index

import sys
import os

TBI = "tbi"  # tabix index file format


class Bed:
    """Represents genomic regions from a BED file.

    Stores genomic coordinates from a BED file and provides methods for accessing,
    manipulating, and extracting regions.

    Attributes:
        _fname (str): The path to the BED file.
        _coordinates (List[Coordinate]): A list of Coordinate objects representing the genomic regions.
    """

    def __init__(self, bedfile: str, padding: int, debug: bool) -> None:
        """Initialize a Bed object with genomic regions from a BED file.

        Loads genomic coordinates from the specified BED file and prepares them for region extraction and manipulation.

        Args:
            bedfile: The path to the BED file.
            padding: The padding to add to the start and stop coordinates.
            debug: Flag to enable debug mode.

        Raises:
            FileNotFoundError: If the BED file does not exist.
        """
        self._debug = debug  # store debug flag
        if not os.path.isfile(bedfile):
            exception_handler(
                FileNotFoundError,
                f"Cannot find input BED file {bedfile}",
                os.EX_DATAERR,
                self._debug,
            )
        self._fname = bedfile  # store input file name
        # read input bed file content and store a list of coordinates
        self._coordinates = self._read(padding)

    def __repr__(self) -> str:
        """Return a string representation of the Bed object.

        The representation includes the class name and the number of stored regions.

        Returns:
            A string representation of the Bed object.
        """
        return f"<{self.__class__.__name__} object; stored regions={len(self)}>"

    def __len__(self) -> int:
        """Return the number of regions in the Bed object.

        Returns the number of regions stored in the _coordinates attribute.

        Returns:
            The number of regions.

        Raises:
            AttributeError: If the _coordinates attribute is missing.
        """
        if not hasattr(self, "_coordinates"):  # always trace this error
            raise AttributeError(
                f"Missing _coordinates attribute on {self.__class__.__name__}"
            )
        return len(self._coordinates)

    def __iter__(self) -> "BedIterator":
        """Return an iterator over the regions in the Bed object.

        Returns an iterator that allows iterating over the regions in the Bed object.

        Returns:
            An iterator over the regions in the Bed object.
        """
        return BedIterator(self)

    def __getitem__(
        self, idx: Union[int, slice]
    ) -> Union[Coordinate, List[Coordinate]]:
        """Return the regions coordinate at the given index or slice.

        Returns the regions coordinate at the given index or slice from the
        _coordinates list.

        Args:
            idx: The index or slice of the desired regions coordinate.

        Returns:
            The regions coordinate at the given index or slice.

        Raises:
            AttributeError: If the _coordinates attribute is missing.
            IndexError: If the index is out of range.
            TypeError: If the index is not an integer or a slice.
        """
        if not hasattr(self, "_coordinates"):  # always trace this error
            raise AttributeError(
                f"Missing _ccordinates attribute on {self.__class__.__name__}"
            )
        try:  # access _coordinates list to return the corresponding coordinate
            return self._coordinates[idx]
        except IndexError as e:
            raise IndexError(f"Index {idx} out of range") from e
        except TypeError as e:
            raise TypeError(
                f"Invalid index type ({type(idx).__name__}), expected {int.__name__} or {slice.__name__}"
            ) from e

    def _read(self, padding: int) -> List[Coordinate]:
        """Read the BED file and return a list of Coordinate objects.

        Reads the BED file specified by the filename attribute, parses each line,
        and returns a list of Coordinate objects representing the genomic intervals.

        Args:
            padding: The padding to add to the start and stop coordinates.

        Returns:
            A list of Coordinate objects.

        Raises:
            FileNotFoundError: If the BED file is not found.
            PermissionError: If there is a permission error reading the BED file.
            IOError: If there is an I/O error reading the BED file.
            Exception: If any other error occurs while reading the BED file.
        """
        coordinates = []  # list of coordinates from the input bed (Coordinate objs)
        try:
            with open(self._fname, mode="r") as infile:  # begin Bedfile parsing
                # valid bed files must have at least three columns: chromosome, start,
                # and end coordinates; separator is not necessarily tab
                coordinates.extend(
                    _parse_bed_line(line, i + 1, padding, self._debug)
                    for i, line in enumerate(infile)
                    if not line.startswith("#") and line.strip()
                )
        except FileNotFoundError as e:  # bed file not found
            exception_handler(
                FileNotFoundError,
                f"Unable to find {self._fname}",
                os.EX_DATAERR,
                self._debug,
                e,
            )
        except PermissionError as e:  # permission error on reading
            exception_handler(
                PermissionError,
                f"Permission denied when trying reading {self._fname}",
                os.EX_DATAERR,
                self._debug,
                e,
            )
        except IOError as e:  # i/o error on read
            exception_handler(
                IOError,
                f"I/O error while reading {self._fname}",
                os.EX_DATAERR,
                self._debug,
                e,
            )
        except Exception as e:  # generic exception caught
            exception_handler(
                Exception,
                f"An unexpected error occurred while reading {self._fname}",
                os.EX_DATAERR,
                self._debug,
                e,
            )  # sourcery skip: raise-specific-error
        return coordinates

    def extract_regions(self, fastas: Dict[str, Fasta]) -> RegionList:
        if not hasattr(self, "_coordinates") or self._coordinates is None:
            raise AttributeError("Missing coordinates list, cannot extract sequences")
        return RegionList(
            [Region(fastas[c.contig].fetch(c), c) for c in self._coordinates]
        )


class BedIterator:
    """Iterator for Bed objects.

    Iterates over the regions in a Bed object.

    Attributes:
        bed: The Bed object to iterate over.
        _index: The current index of the iterator.
    """

    def __init__(self, bed: Bed) -> None:
        """Initialize the BedIterator.

        Initializes the iterator with the given Bed object.

        Args:
            bed: The bed object to iterate over.

        Raises:
            AttributeError: If the bed object is missing the _coordinates attribute.
        """
        if not hasattr(bed, "_coordinates"):  # always trace this error
            raise AttributeError(
                f"Missing _ccordinates attribute on {bed.__class__.__name__}"
            )
        self._bed = bed  # bed object to be iterated
        self._index = 0  # iterator index used over the coordinates list

    def __next__(self) -> Union[Coordinate, List[Coordinate]]:
        """Return the next region coordinate in the Bed object.

        Returns the next region coordinate in the Bed object, and advances the
        iterator.

        Returns:
            The next region coordinate in the Bed object.

        Raises:
            StopIteration: If there are no more regions to iterate over.
        """
        if self._index < len(self._bed):  # not reached end of bed coordinates
            result = self._bed[self._index]
            self._index += 1  # go to next position in the list
            return result
        raise StopIteration  # stop iteration over bed object


def _parse_bed_line(
    bedline: str, linenum: int, padding: int, debug: bool
) -> Coordinate:
    """Parse a single line from a BED file and return a Coordinate object.

    Extracts chromosome, start, and stop positions from a BED file line, applies padding, and returns a Coordinate object. Handles errors for missing or invalid fields.

    Args:
        bedline: The line from the BED file to parse.
        linenum: The line number in the BED file (for error reporting).
        padding: The padding to add to the start and stop coordinates.
        debug: Flag to enable debug mode.

    Returns:
        A Coordinate object representing the parsed region.

    Raises:
        ValueError: If there are fewer than three columns or stop < start.
        TypeError: If start or stop values are not integers.
    """
    columns = bedline.strip().split()  # recover bed fields for current line
    # minimum fields required: chrom, start, stop
    # (see https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
    if len(columns) < 3:
        exception_handler(
            ValueError,
            f"Less than three columns at line {linenum}",
            os.EX_DATAERR,
            debug,
        )
    try:  # initialize chrom, start and stop fields
        chrom, start, stop = columns[0], int(columns[1]), int(columns[2])
    except ValueError as e:  # raise if start or stop are not valid int
        exception_handler(
            TypeError,
            f"Start/stop values at line {linenum} are not {int.__name__}",
            os.EX_DATAERR,
            debug,
            e,
        )
    if stop < start:  # ensure valid genomic range
        exception_handler(
            ValueError,
            f"Stop < start coordinate ({stop} < {start}) at line {linenum}",
            os.EX_DATAERR,
            debug,
        )
    # if required, pad the input region sequence up and downstream
    return Coordinate(chrom, start, stop, padding)


class BedAnnotation:

    def __init__(self, fname: str, verbosity: int, debug: bool) -> None:
        self._debug = debug  # store debug flag
        self._verbosity = verbosity  # store verbosity level
        self._fname = fname  # store input file name
        self._bedidx = self._search_index()  # initialize bed index
        if not self._bedidx:  # index not found, compute it
            self.index_bed()
        # initialize TabixFile object with the previously computed index
        self._bed = TabixFile(self._fname, index=self._bedidx)

    def _search_index(self, bedidx: Optional[str] = "") -> str:
        """Search for or validate a Tabix index for the BED file.

        Searches for a Tabix index (.tbi) for the associated BED file if one is not
        provided. If a path to an index is provided, it validates that the index
        exists and is not empty.

        Args:
            bedidx: An optional path to a Tabix index file.

        Returns:
            The path to the Tabix index file, or an empty string if not found.

        Raises:
            FileNotFoundError: If the provided index file does not exist or is empty.
        """

        # look for index for the current vcf, if not found compute it
        if not bedidx:
            if _find_tbi(self._fname):  # index found, store it
                return f"{self._fname}.{TBI}"
            # index not found -> compute it de novo and store it in the same folder
            # as the input vcf
            sys.stdout.write(f"Tabix index not found for {self._fname}\n")
            return ""
        # precomputed vcf index index must be a non empty file
        if not (os.path.isfile(bedidx) and os.stat(bedidx).st_size > 0):
            raise FileNotFoundError(f"Not existing or empty VCF index {bedidx}")
        return bedidx

    def index_bed(self, pytest: Optional[bool] = False) -> None:
        """Create or update the Tabix index for the BED file.

        Creates or updates the Tabix index (.tbi) for the associated BED file.
        If an index already exists, it will be overwritten.

        Raises:
            OSError: If an error occurs during indexing.
            RuntimeWarning: If an index already exists.
        """
        if self._bedidx and not pytest:  # launch warning
            warning("Tabix index already present, forcing update", self._verbosity)
        # compute tabix index if not provided during object initialization
        try:  # create index in the same folder as the input vcf
            tabix_index(self._fname, preset="bed", force=True)
        except OSError as e:
            exception_handler(
                OSError,
                f"An error occurred while indexing {self._fname}",
                os.EX_DATAERR,
                self._debug,
                e,
            )
        assert _find_tbi(self._fname)
        self._bedidx = f"{self._fname}.{TBI}"

    def fetch_features(
        self, contig: str, start: int, stop: int
    ) -> Union[List[str], None]:
        if contig not in self._bed.contigs:
            return None
        return [e.strip() for e in self._bed.fetch(contig, start, stop)]

    @property
    def contigs(self) -> List[str]:
        return self._bed.contigs


def _find_tbi(bedfile: str) -> bool:
    """Check if a Tabix index exists for a BED file.

    Checks if a Tabix index (.tbi) exists for the given BED file and is a non-empty
    file.

    Args:
        bedfile: The path to the BED file.

    Returns:
        True if the index exists and is a non-empty file, False otherwise.
    """
    # avoid unexpected crashes due to file location
    bedindex = f"{os.path.abspath(bedfile)}.{TBI}"
    if os.path.exists(bedindex):  # index must be a non empty file
        return os.path.isfile(bedindex) and os.stat(bedindex).st_size > 0
    return False
