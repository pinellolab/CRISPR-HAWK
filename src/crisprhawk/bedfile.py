"""
"""

from crisprhawk_error import CrisprHawkBedError
from exception_handlers import exception_handler
from sequences import Sequence, Fasta

from typing import List, Union

import os


class Coordinate:
    def __init__(self, contig: str, start: int, stop: int) -> None:
        # initialize a genomic coordinate interval
        # used when a bed file is provided in the input to select the
        # genomic regions to consider when searching gRNAs
        self._contig = contig  # set contig name
        self._start = start  # set start coordinate
        self._stop = stop  # set stop coordinate

    def __str__(self):
        return f"{self._contig}:{self._start}-{self._stop}"

    @property
    def contig(self) -> str:
        return self._contig

    @property
    def start(self) -> int:
        return self._start

    @property
    def stop(self) -> int:
        return self._stop


def _parse_bed_line(bedline: str, linenum: int, debug: bool) -> Coordinate:
    columns = bedline.strip().split()  # recover bed fields for current line
    # minimum fields required: chrom, start, stop
    # (see https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
    if len(columns) < 3:
        exception_handler(
            CrisprHawkBedError,
            f"Less than 3 columns at line {linenum}",
            os.EX_DATAERR,
            debug,
        )
    try:  # initialize chrom, start and stop fields
        chrom, start, stop = columns[0], int(columns[1]), int(columns[2])
    except ValueError as e:  # raise if start or stop are not valid int
        exception_handler(
            e,
            f"Start/stop values at line {linenum} are not {int.__name__}",
            os.EX_DATAERR,
            debug,
        )
    if stop < start:  # ensure valid genomic range
        exception_handler(
            CrisprHawkBedError,
            f"Stop < start coordinate ({stop} < {start}) at line {linenum}",
            os.EX_DATAERR,
            debug,
        )
    return Coordinate(chrom, start, stop)  # initialize coordinate object

class Region(Sequence):
    def __init__(self, sequence: str, coord: Coordinate, debug: bool):
        super().__init__(sequence, debug)  # init parent sequence object
        # store contig name, start and stop coordinates of the extracted region
        self._coordinates = coord

    def __str__(self):
        return f">{self._coordinates}\n{super().__str__()}"

    @property
    def contig(self) -> str:
        return self._coordinates.contig

    @property
    def start(self) -> int:
        return self._coordinates.start

    @property
    def stop(self) -> int:
        return self._coordinates.stop


class Bed:
    def __init__(self, bedfile: str, debug: bool) -> None:
        self._fname = bedfile  # store input file name
        self._debug = debug  # store debug mode flag
        # read input bed file content and store a list of coordinates
        self._coordinates = self._read()
    
    def __len__(self) -> int:
        if not hasattr(self, "_coordinates"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _coordinates attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return len(self._coordinates)

    def __getitem__(self, idx: Union[int, slice]) -> Coordinate:
        if not hasattr(self, "_coordinates"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _ccordinates attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        try:  # access _coordinates list to return the corresponding coordinate
            return self._coordinates[idx]
        except IndexError as e:
            exception_handler(
                CrisprHawkBedError,
                f"Index {idx} out of range",
                os.EX_DATAERR,
                self._debug,
                e,
            )
        except TypeError as e:
            exception_handler(
                CrisprHawkBedError,
                f"Invalid index type ({type(idx).__name__}), expected {int.__name__} or {slice.__name__}",
            )

    def _read(self) -> List[Coordinate]:
        coordinates = []  # list of coordinates from the input bed (Coordinate objs)
        try:
            with open(self._fname, mode="r") as infile:  # begin Bedfile parsing
                # valid bed files must have at least three columns: chromosome, start,
                # and end coordinates; separator is not necessarily tab
                for i, line in enumerate(infile):
                    if (
                        line.startswith("#") or not line.strip()
                    ):  # skip comments and empy lines
                        continue
                    # add coordinate to coordinates list
                    coordinates.append(_parse_bed_line(line, i + 1, self._debug))
        except FileNotFoundError as e:  # bed file not found
            exception_handler(
                CrisprHawkBedError,
                f"Unable to find {self._fname}",
                os.EX_OSFILE,
                self._debug,
                e,
            )
        except PermissionError as e:  # permission error on reading
            exception_handler(
                CrisprHawkBedError,
                f"Permission denied when trying reading {self._fname}",
                os.EX_OSERR,
                self._debug,
                e,
            )
        except IOError as e:  # i/o error on read
            exception_handler(
                CrisprHawkBedError,
                f"I/O error while reading {self._fname}",
                os.EX_IOERR,
                self._debug,
                e,
            )
        except Exception as e:  # generic exception caught
            exception_handler(
                CrisprHawkBedError,
                f"An unexpected error occurred while reading {self._fname}",
                os.EX_OSERR,
                self._debug,
                e,
            )
        return coordinates
    
    def extract(self, fasta: Fasta) -> List[Region]:
        if not hasattr(self, "_coordinates"):
            exception_handler(CrisprHawkBedError, "Missing coordinates list, cannot extract sequences", os.EX_DATAERR, self._debug)
        return [
            Region(fasta.fetch(c.contig, c.start, c.stop), c, self._debug)
            for c in self._coordinates
        ]


class BedIterator:
    def __init__(self, bed: Bed) -> None:
        if not hasattr(bed, "_coordinates"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _ccordinates attribute on {bed.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        self._bed = bed  # bed object to be iterated
        self._index = 0  # iterator index used over the coordinates list

    def __next__(self) -> Coordinate:
        if self._index < len(self._bed):  # not reached end of bed coordinates
            self._index += 1  # go to next position in the list
            return self._bed[self._index - 1]
        raise StopIteration  # stop iteration over bed object
    
