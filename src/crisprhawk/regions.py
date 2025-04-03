"""
"""

from exception_handlers import exception_handler
from utils import print_verbosity, VERBOSITYLVL

from hapsolver import Fasta, Bed, RegionList
from time import time 

import os

PADDING = 100  # padding length for regions

def read_fasta(fastafile: str, fasta_idx: str, verbosity: int, debug: bool) -> Fasta:
    # read input fasta and bed and construct Region object for each genomic region
    print_verbosity(f"Reading input Fasta file {fastafile}", verbosity, VERBOSITYLVL[3])
    start = time()  # track processing time
    try:  # create fasta object
        fasta = Fasta(fastafile, faidx=fasta_idx)  
    except (FileNotFoundError, OSError) as e:
        exception_handler(Exception, f"Failed parsing Fasta file ({fastafile})", os.EX_DATAERR, debug, e)
    print_verbosity(f"Fasta file parsed in {(time() - start):.2f}s", verbosity, VERBOSITYLVL[3])
    return fasta

def read_bed(bedfile: str, verbosity: int, debug: bool) -> Bed:
    # read input bed and construct Region object for each genomic region
    print_verbosity(f"Parsing input BED file {bedfile}", verbosity, VERBOSITYLVL[3])
    start = time()  # track processing time
    try:  # create bed object
        bed = Bed(bedfile, PADDING)  # guidelen used to pad regions
    except (FileNotFoundError, PermissionError, IOError, Exception) as e:
        exception_handler(Exception, f"Failed parsing BED file ({bedfile})", os.EX_DATAERR, debug, e)
    print_verbosity(f"Parsed {len(bed)} regions in {(time() - start):.2f}s", verbosity, VERBOSITYLVL[3])
    return bed

def extract_regions(bed: Bed, fasta: Fasta, verbosity: int, debug: bool) -> RegionList:
    # extract genomic sequences foe each input region
    print_verbosity(f"Extracting regions from {fasta.fname}", verbosity, VERBOSITYLVL[3])
    start = time()  # track processing time
    try:  # extract regions from input fasta file
        regions = bed.extract_regions(fasta)
    except AttributeError as e:
        exception_handler(Exception, f"Failed region extraction from {fasta.fname}", os.EX_DATAERR, debug, e)
    print_verbosity(f"Extracted regions:\n{str(regions)}", verbosity, VERBOSITYLVL[3])
    print_verbosity(f"Regions extracted in {(time() - start):.2f}s", verbosity, VERBOSITYLVL[3])
    return regions

def construct_regions(fastafile: str, bedfile: str, fasta_idx: str, guidelen: int, verbosity: int, debug: bool) -> RegionList:
    # read input fasta and bed and construct Region object for each genomic region
    print_verbosity("Retrieving input genomic regions", verbosity, VERBOSITYLVL[1])
    start = time()  # track processing time
    fasta = read_fasta(fastafile, fasta_idx, verbosity, debug)
    bed = read_bed(bedfile, verbosity, debug)
    regions = extract_regions(bed, fasta, verbosity, debug)
    print_verbosity(f"Genomic regions retrieved in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2])
    return regions



