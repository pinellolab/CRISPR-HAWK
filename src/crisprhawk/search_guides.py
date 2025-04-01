"""
"""

from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkBitsetError
from utils import print_verbosity, adjust_guide_position, VERBOSITYLVL
from guide import Guide, GUIDESEQPAD
from bitset import Bitset
from pam import PAM
from regions import PADDING

from typing import List, Tuple
from hapsolver import Region, Haplotype, flatten_list
from time import time

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
    for pos in range(start, stop):  # start scanning haplotype for guides
        if match(pam.bits, haplotype[pos:(pos + len(pam))], pos, debug):
            matches_fwd.append(pos)  # hit on forward strand
        if match(pam.bitsrc, haplotype[pos:(pos + len(pam))], pos, debug):
            matches_rev.append(pos)  # hit on negative strand
    return matches_fwd, matches_rev


def pam_search(pam: PAM, region: Region, haplotypes: List[Haplotype], haplotypes_bits: List[List[Bitset]], verbosity: int, debug: bool) -> List[Tuple[List[int], List[int]]]:
    pam_hits = []  # list of pam hits
    for i, hap in enumerate(haplotypes):
        print_verbosity(f"Searching PAM occurrences in haplotype {hap.samples}", verbosity, VERBOSITYLVL[3])
        # define scan stop position for each haplotype
        scan_stop = hap.position_map[(len(region) - PADDING)][-1] - len(pam) + 1
        scan_start = hap.position_map[PADDING][-1]  # start position for guides search
        # scan haplotype for pam occurrences 
        matches_fwd, matches_rev = scan_haplotype(pam, haplotypes_bits[i], scan_start, scan_stop, debug)
        print_verbosity(f"Found {len(matches_fwd) + len(matches_rev)} PAM occurrences ({len(matches_fwd)} on 5'-3'; {len(matches_rev)} on 3'-5')", verbosity, VERBOSITYLVL[3])
        pam_hits.append((matches_fwd, matches_rev))  # store pam hits for each haplotype
    return pam_hits

def extract_guide_sequence(haplotype: Haplotype, position: int, pamlen: int, guidelen: int, right: bool) -> str:
    if right:
        return "".join(haplotype[position - GUIDESEQPAD:position + guidelen + pamlen + GUIDESEQPAD])
    return "".join(haplotype[position - guidelen - GUIDESEQPAD: position + pamlen + GUIDESEQPAD])

def retrieve_guides(pam_hits: List[int], haplotype: Haplotype, guidelen: int, pamlen: int, direction: int, right: bool, verbosity: int, debug: bool) -> List[Guide]:
    guides = []  # list of haplotype-specific guides
    print_verbosity(f"Retrieving guides from haplotype {haplotype.samples}", verbosity, VERBOSITYLVL[3])
    start = time()
    for pos in pam_hits:
        # retrieve guide sequence from haplotype
        guideseq = extract_guide_sequence(haplotype, pos, pamlen, guidelen, right)
        if haplotype.samples != "REF" and guideseq[GUIDESEQPAD:-GUIDESEQPAD].isupper():  # reference guide
            continue
        position = haplotype.position_map_rev[pos][0] 
        position += adjust_guide_position(pos, haplotype.position_map[position])
        guides.append(Guide(position, guideseq, guidelen, pamlen, direction, debug, right))
    print_verbosity(f"Guides retrieved in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3])
    return guides
        
    

def search(pam: PAM, region: Region, haplotypes: List[Haplotype], haplotypes_bits: List[List[Bitset]], guidelen: int, right: bool, verbosity: int, debug: bool) -> List[Guide]:
    print_verbosity(f"Searching guide candidates in {region.coordinates}", verbosity, VERBOSITYLVL[2])
    # search pam occurrences on forward and reverse strand of the input sequence
    pam_hits = pam_search(pam, region, haplotypes, haplotypes_bits, verbosity, debug)
    # retrieve guide candidates sequence for each haplotype
    guides = flatten_list(
        [
            retrieve_guides(pam_hits[i][strand], haplotype, guidelen, len(pam), strand, (not right if strand == 1 else right), verbosity, debug)
            for i, haplotype in enumerate(haplotypes)
            for strand in [0, 1]
        ]
    )
    return guides





