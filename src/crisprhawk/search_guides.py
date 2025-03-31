"""
"""

from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkBitsetError
from utils import print_verbosity, VERBOSITYLVL
from guide import Guide
from bitset import Bitset
from pam import PAM
from regions import PADDING

from typing import List, Tuple
from hapsolver import Region, Haplotype

import os

def match(bitset1: Bitset, bitset2: Bitset, position: int, debug: bool) -> bool:
    # bitwise matching operation for input bitsets
    # assumes the two bitsets have the same length
    try:
        return all((ntbit & bitset2[i]).to_bool() for i, ntbit in enumerate(bitset1))
    except ValueError as e:
        exception_handler(CrisprHawkBitsetError, f"PAM bitwise matching failed at position {position}", os.EX_DATAERR, debug, e,)

def scan_haplotype(pam: PAM, haplotype: List[Bitset], start: int, stop: int, debug: bool) -> Tuple[List[int], List[int]]:
    # lists storing hits for input pam on forward and reverse strands
    matches_fwd, matches_rev = [], []
    for pos in range(start, stop + 1):  # start scanning haplotype for guides
        if match(pam.bits, haplotype[pos:(pos + len(pam))], pos, debug):
            matches_fwd.append(pos)  # hit on forward strand
        if match(pam.bitsrc, haplotype[pos:(pos + len(pam))], pos, debug):
            matches_rev.append(pos)  # hit on negative strand
    return matches_fwd, matches_rev


def pam_search(pam: PAM, region: Region, haplotypes: List[Haplotype], haplotypes_bits: List[List[Bitset]], guidelen: int, verbosity: int, debug: bool):
    
    for i, hap in enumerate(haplotypes):
        print_verbosity(f"Searching PAM occurrences in haplotype {hap.samples}", verbosity, VERBOSITYLVL[3])
        # define scan stop position for each haplotype
        scan_stop = hap.position_map[(len(region) - PADDING)][-1] - len(pam) + 1
        scan_start = PADDING  # start position for guides search
        # scan haplotype for pam occurrences 
        print(hap.sequence[scan_start:scan_stop + 3], hap.samples)
        print(scan_start, scan_stop)
        matches_fwd, matches_rev = scan_haplotype(pam, haplotypes_bits[i], scan_start, scan_stop, debug)
        print(matches_fwd, matches_rev)
        print_verbosity(f"Found {len(matches_fwd) + len(matches_rev)} PAM occurrences", verbosity, VERBOSITYLVL[3])


    
    

def search(pam: PAM, region: Region, haplotypes: List[Haplotype], haplotypes_bits: List[List[Bitset]], guidelen: int, right: bool, verbosity: int, debug: bool) -> List[Guide]:
    # search pam occurrences on forward and reverse strand of the input sequence
    print_verbosity(f"Searching guide candidates in {region.coordinates}", verbosity, VERBOSITYLVL[2])
    pam_search(pam, region, haplotypes, haplotypes_bits, guidelen, verbosity, debug)






