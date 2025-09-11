""" """

from .exception_handlers import exception_handler
from .utils import print_verbosity, VERBOSITYLVL
from .region import RegionList
from .sequence import Fasta
from .bedfile import Bed

from typing import List, Dict, Tuple

from time import time

import os

PADDING = 100  # padding length for regions


def read_fasta(
    fastafiles: List[str], fasta_idx: str, verbosity: int, debug: bool
) -> Dict[str, Fasta]:
    # read input fasta and bed and construct Region object for each genomic region
    print_verbosity("Reading input Fasta files", verbosity, VERBOSITYLVL[3])
    start = time()  # track processing time
    try:
        fastas = [
            Fasta(fastafile, verbosity, debug, faidx=fasta_idx)
            for fastafile in fastafiles
        ]
        fastas_dict = {fasta.contig: fasta for fasta in fastas}
    except Exception as e:
        exception_handler(
            Exception,
            f"Failed parsing Fasta files",
            os.EX_DATAERR,
            debug,
            e,
        )
    print_verbosity(
        f"Fasta files parsed in {(time() - start):.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return fastas_dict


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
            Exception, f"Failed parsing BED file ({bedfile})", os.EX_DATAERR, debug, e  # type: ignore
        )
    print_verbosity(
        f"Parsed {len(bed)} regions in {(time() - start):.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return bed


def extract_regions(
    bed: Bed, fastas: Dict[str, Fasta], verbosity: int, debug: bool
) -> RegionList:
    # extract genomic sequences foe each input region
    print_verbosity(f"Extracting regions from FASTA files", verbosity, VERBOSITYLVL[3])
    start = time()  # track processing time
    try:  # extract regions from input fasta file
        regions = bed.extract_regions(fastas)
    except (AttributeError, KeyError, Exception) as e:
        exception_handler(
            Exception,
            f"Failed region extraction from FASTA files",
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
    fastafiles: List[str],
    bedfile: str,
    fasta_idx: str,
    verbosity: int,
    debug: bool,
) -> RegionList:
    # read input fasta and bed and construct Region object for each genomic region
    print_verbosity("Retrieving input genomic regions", verbosity, VERBOSITYLVL[1])
    start = time()  # track processing time
    fastas = read_fasta(fastafiles, fasta_idx, verbosity, debug)
    bed = read_bed(bedfile, verbosity, debug)
    regions = extract_regions(bed, fastas, verbosity, debug)
    print_verbosity(
        f"Genomic regions retrieved in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[2],
    )
    return regions
