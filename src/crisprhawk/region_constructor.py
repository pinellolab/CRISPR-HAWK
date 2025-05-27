""" """

from .exception_handlers import exception_handler
from .utils import print_verbosity, VERBOSITYLVL
from .region import RegionList
from .sequence import Fasta
from .bedfile import Bed


from time import time

import os

PADDING = 100  # padding length for regions


def read_fasta(fastafile: str, fasta_idx: str, verbosity: int, debug: bool) -> Fasta:
    """Read a FASTA file and return a Fasta object.

    Loads the specified FASTA file and returns a Fasta object for sequence
    operations. Handles errors and prints progress messages.

    Args:
        fastafile: The path to the FASTA file.
        fasta_idx: The path to the FASTA index file.
        verbosity: The verbosity level for logging.
        debug: Flag to enable debug mode.

    Returns:
        A Fasta object representing the loaded FASTA file.

    Raises:
        Exception: If the FASTA file cannot be parsed.
    """
    # read input fasta and bed and construct Region object for each genomic region
    print_verbosity(f"Reading input Fasta file {fastafile}", verbosity, VERBOSITYLVL[3])
    start = time()  # track processing time
    try:  # create fasta object
        fasta = Fasta(fastafile, verbosity, debug, faidx=fasta_idx)
    except Exception as e:
        exception_handler(
            Exception, # type: ignore
            f"Failed parsing Fasta file ({fastafile})",
            os.EX_DATAERR,
            debug,
            e,
        )
    print_verbosity(
        f"Fasta file parsed in {(time() - start):.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return fasta


def read_bed(bedfile: str, verbosity: int, debug: bool) -> Bed:
    """Read a BED file and return a Bed object.

    Loads the specified BED file and returns a Bed object for region extraction
    and manipulation. Handles errors and prints progress messages.

    Args:
        bedfile: The path to the BED file.
        verbosity: The verbosity level for logging.
        debug: Flag to enable debug mode.

    Returns:
        A Bed object representing the loaded BED file.

    Raises:
        Exception: If the BED file cannot be parsed.
    """
    # read input bed and construct Region object for each genomic region
    print_verbosity(f"Parsing input BED file {bedfile}", verbosity, VERBOSITYLVL[3])
    start = time()  # track processing time
    try:  # create bed object
        bed = Bed(bedfile, PADDING, debug)  # guidelen used to pad regions
    except Exception as e:
        exception_handler(
            Exception, f"Failed parsing BED file ({bedfile})", os.EX_DATAERR, debug, e # type: ignore
        )
    print_verbosity(
        f"Parsed {len(bed)} regions in {(time() - start):.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return bed


def extract_regions(bed: Bed, fasta: Fasta, verbosity: int, debug: bool) -> RegionList:
    """Extract genomic regions from a FASTA file using BED coordinates.

    Retrieves and returns a RegionList object containing sequences from the
    FASTA file for each region specified in the BED file. Handles errors and
    prints progress messages.

    Args:
        bed: The Bed object containing genomic coordinates.
        fasta: The Fasta object representing the FASTA file.
        verbosity: The verbosity level for logging.
        debug: Flag to enable debug mode.

    Returns:
        A RegionList object containing the extracted regions.

    Raises:
        Exception: If region extraction fails due to missing coordinates.
    """
    # extract genomic sequences foe each input region
    print_verbosity(
        f"Extracting regions from {fasta.fname}", verbosity, VERBOSITYLVL[3]
    )
    start = time()  # track processing time
    try:  # extract regions from input fasta file
        regions = bed.extract_regions(fasta)
    except AttributeError as e:
        exception_handler(
            Exception, # type: ignore
            f"Failed region extraction from {fasta.fname}",
            os.EX_DATAERR,
            debug,
            e,
        )
    print_verbosity(f"Extracted regions:\n{str(regions)}", verbosity, VERBOSITYLVL[3])
    print_verbosity(
        f"Regions extracted in {(time() - start):.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return regions


def construct_regions(
    fastafile: str,
    bedfile: str,
    fasta_idx: str,
    verbosity: int,
    debug: bool,
) -> RegionList:
    """Construct and return a list of genomic regions from input files.

    Reads the input FASTA and BED files, extracts genomic regions, and returns a
    RegionList object containing the regions. Handles progress reporting and
    error management.

    Args:
        fastafile: The path to the FASTA file.
        bedfile: The path to the BED file.
        fasta_idx: The path to the FASTA index file.
        verbosity: The verbosity level for logging.
        debug: Flag to enable debug mode.

    Returns:
        A RegionList object containing the constructed genomic regions.
    """
    # read input fasta and bed and construct Region object for each genomic region
    print_verbosity("Retrieving input genomic regions", verbosity, VERBOSITYLVL[1])
    start = time()  # track processing time
    fasta = read_fasta(fastafile, fasta_idx, verbosity, debug)
    bed = read_bed(bedfile, verbosity, debug)
    regions = extract_regions(bed, fasta, verbosity, debug)
    print_verbosity(
        f"Genomic regions retrieved in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[2],
    )
    return regions
