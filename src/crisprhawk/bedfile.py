"""
"""

from crisprhawk_error import CrisprHawkBedError, CrisprHawkEnrichmentError
from exception_handlers import exception_handler
from sequences import Sequence, Fasta

from typing import List, Union, Optional, final

import os


class Coordinate:
    def __init__(self, contig: str, start: int, stop: int) -> None:
        # initialize a genomic coordinate interval
        # used when a bed file is provided in the input to select the
        # genomic regions to consider when searching gRNAs
        self._contig = contig  # set contig name
        self._start = start  # set start coordinate
        self._stop = stop  # set stop coordinate

    def __str__(self) -> str:
        return f"{self._contig}:{self._start}-{self._stop}"

    def format(self, pad: Optional[int] = 0) -> str:
        return "_".join(
            list(map(str, [self._contig, self._start + pad, self._stop - pad]))
        )

    @property
    def contig(self) -> str:
        return self._contig

    @property
    def start(self) -> int:
        return self._start

    @property
    def stop(self) -> int:
        return self._stop


def _parse_bed_line(
    bedline: str, linenum: int, guidelen: int, debug: bool
) -> Coordinate:
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
    # handle pam occurrences on the edges of the region by padding up and downstream
    # the input region by |g| nts, where g is the guide
    return Coordinate(chrom, start - guidelen, stop + guidelen)


class Region(Sequence):
    def __init__(self, sequence: str, coord: Coordinate, debug: bool):
        super().__init__(sequence, debug)  # init parent sequence object
        # store contig name, start and stop coordinates of the extracted region
        self._coordinates = coord

    def __str__(self):
        return f">{self._coordinates}\n{super().__str__()}"

    def format(self, pad: Optional[int] = 0) -> str:
        return self._coordinates.format(pad)

    def enrich(self, pos: int, iupac_nt: str) -> None:
        try:
            self._sequence_raw[pos] = iupac_nt
        except IndexError as e:
            exception_handler(
                f"Enrichment failed for region {self.format()}",
                CrisprHawkEnrichmentError,
                os.EX_DATAERR,
                self._debug,
                e,
            )

    @final
    def contains(self, indelregion: "IndelRegion") -> bool:
        if not isinstance(indelregion, IndelRegion):
            exception_handler(
                TypeError,
                f"Full overlap check on input regiin can only be done on {IndelRegion.__name__}",
                os.EX_DATAERR,
                self._debug,
            )
        # check whether the test region is contained in the current region
        return (
            self.contig == indelregion.contig
            and self.start <= indelregion.start
            and self.stop >= indelregion.stop
        )

    @property
    def contig(self) -> str:
        return self._coordinates.contig

    @property
    def start(self) -> int:
        return self._coordinates.start

    @property
    def stop(self) -> int:
        return self._coordinates.stop


class IndelRegion(Region):
    def __init__(
        self,
        indelpos: int,
        sequence: str,
        coord: Coordinate,
        ref: str,
        alt: str,
        indel_type: str,
        debug: bool,
    ) -> None:
        super().__init__(sequence, coord, debug)  # initialize indel region
        self._refseq = sequence  # keep ref seq to simplify haplotype tracking
        self._indelpos = indelpos  # store indel position
        self._ref = ref  # indel reference allele
        self._alt = alt  # indel alt allele
        self._indel_type = indel_type  # insertion or deletion
        self._indel_len = abs(len(alt) - len(ref))  # compute indel length

    def _update_sequence(self) -> None:
        assert hasattr(self, "_sequence_raw")
        # update pointed sequence after changes on raw sequence data
        self._sequence = "".join(self._sequence_raw)

    def enrich(self) -> None:
        try:
            # enrich sequence replacing reference allele with indel allele
            self._sequence_raw = (
                self._sequence_raw[: self._indelpos]
                + list(self._alt)
                + self._sequence_raw[self._indelpos + len(self._ref) :]
            )
            self._sequence = "".join(self._sequence_raw)  # adjust pointed sequence
        except IndexError as e:
            exception_handler(
                f"Indel enrichement failed for region {self.format()}",
                CrisprHawkEnrichmentError,
                os.EX_DATAERR,
                self._debug,
                e,
            )

    @property
    def refseq(self) -> str:
        return self._refseq

    @property
    def indelpos(self) -> int:
        return self._indelpos

    @property
    def indel_len(self) -> int:
        return self._indel_len

    @property
    def indel_type(self) -> int:
        return self._indel_type


class RegionList:
    def __init__(self, regions: List[Region], debug: bool) -> None:
        self._regions = regions
        self._debug = debug

    def __len__(self) -> int:
        return len(self._regions)

    def __iter__(self) -> "RegionListIterator":
        return RegionListIterator(self)

    def __getitem__(self, idx: Union[int, slice]) -> Region:
        if not hasattr(self, "_regions"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _regions attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        try:
            return self._regions[idx]
        except IndexError as e:
            exception_handler(
                CrisprHawkBedError,
                f"Index {idx} out of range",
                os.EX_DATAERR,
                self._debug,
                e,
            )

    def extend(self, regions: "RegionList") -> None:
        if not isinstance(regions, self.__class__):
            exception_handler(
                TypeError,
                f"Cannot extend {self.__class__.__name__} with objects of type {type(regions).__name__}",
                os.EX_DATAERR,
                self._debug,
            )
        for r in regions:  # extend regions list
            self._regions.append(r)

    def append(self, region: Union[Region, IndelRegion]) -> None:
        if not isinstance(region, Region) and not isinstance(region, IndelRegion):
            exception_handler(
                TypeError,
                f"Cannot append to {self.__class__.__name__} objects of type {type(region).__name__}",
                os.EX_DATAERR,
                self._debug,
            )
        self._regions.append(region)

    def format(self, sep: Optional[str] = "\n", pad: Optional[int] = 0) -> str:
        return sep.join([r.format(pad) for r in self._regions])


class RegionListIterator:
    def __init__(self, regions: RegionList) -> None:
        self._regions = regions  # region list object to iterate over
        self._index = 0  # iterator index used over the list

    def __next__(self) -> Region:
        if self._index < len(self._regions):
            result = self._regions[self._index]
            self._index += 1  # go to next position in list
            return result
        raise StopIteration  # stop iteration over regions list


class Bed:
    def __init__(self, bedfile: str, guidelen: int, debug: bool) -> None:
        self._fname = bedfile  # store input file name
        self._debug = debug  # store debug mode flag
        # read input bed file content and store a list of coordinates
        self._coordinates = self._read(guidelen)

    def __len__(self) -> int:
        if not hasattr(self, "_coordinates"):  # always trace this error
            exception_handler(
                AttributeError,
                f"Missing _coordinates attribute on {self.__class__.__name__}",
                os.EX_DATAERR,
                True,
            )
        return len(self._coordinates)

    def __iter__(self) -> "BedIterator":
        return BedIterator(self)

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

    def _read(self, guidelen: int) -> List[Coordinate]:
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
                    coordinates.append(
                        _parse_bed_line(line, i + 1, guidelen, self._debug)
                    )
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

    def extract(self, fasta: Fasta) -> RegionList:
        if not hasattr(self, "_coordinates"):
            exception_handler(
                CrisprHawkBedError,
                "Missing coordinates list, cannot extract sequences",
                os.EX_DATAERR,
                self._debug,
            )
        return RegionList(
            [
                Region(fasta.fetch(c.contig, c.start, c.stop), c, self._debug)
                for c in self._coordinates
            ],
            self._debug,
        )


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
            result = self._bed[self._index]
            self._index += 1  # go to next position in the list
            return result
        raise StopIteration  # stop iteration over bed object
