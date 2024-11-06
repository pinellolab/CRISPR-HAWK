"""
"""

from crisprhawk_error import CrisprHawkBedError, CrisprHawkEnrichmentError
from exception_handlers import exception_handler
from sequences import Sequence, Fasta
from utils import IUPAC_ENCODER

from typing import List, Union, Optional

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

    def enrich(self, variant: List[str], no_filter: bool) -> None:
        # this check ids made by default, however the user may decide to skip it
        # by using --no-filter option
        if not no_filter and variant[6] != "PASS":  # skip variant if not pass
            return
        vpos = int(variant[1])  # variant position (1-based)
        vposrel = (
            vpos - 1 - self.start
        )  # compute relative variant position (bed coords 0-based)
        ref, alt = variant[3], variant[4]  # retrieve ref and alt alleles
        # compute iupac nt to represent current variant
        iupac_nt = _encode_snp_iupac(
            self._sequence_raw[vposrel], ref, alt, vpos, self._debug
        )
        if iupac_nt is None:  # indel
            return
        self._sequence_raw[vposrel] = iupac_nt  # replace ref allele with iupac

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

    def extract(self, fasta: Fasta) -> List[Region]:
        if not hasattr(self, "_coordinates"):
            exception_handler(
                CrisprHawkBedError,
                "Missing coordinates list, cannot extract sequences",
                os.EX_DATAERR,
                self._debug,
            )
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
            result = self._bed[self._index]
            self._index += 1  # go to next position in the list
            return result
        raise StopIteration  # stop iteration over bed object


def _encode_snp_iupac(
    refnt: str, ref: str, alt: str, pos: int, debug: bool
) -> Union[str, None]:
    if len(ref) != 1:  # deletion, skip
        return
    alleles_alt = alt.split(",")  # retrieve multiallelic alt alleles
    if len(alleles_alt) == 1 and len(alleles_alt[0]) > 1:  # insertion, skip
        return
    if refnt != ref:  # ref alleles must match between vcf and fasta
        exception_handler(
            CrisprHawkEnrichmentError,
            f"Reference allele mismatch between FASTA and VCF file ({refnt} - {ref}) at position {pos}",
            os.EX_DATAERR,
            debug,
        )
    # create iupac string to use for the iupac char encoding
    iupac_string = {
        allele for allele in alleles_alt + [refnt] if len(allele) == 1
    }  # skip indels and focus on just snps
    return IUPAC_ENCODER[iupac_string]
