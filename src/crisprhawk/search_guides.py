"""
"""

from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkBitsetError
from sequences import PAM
from bedfile import Region
from bitset import Bitset

from typing import Tuple, List

import os


def match(bitset1: Bitset, bitset2: Bitset) -> bool:
    # bitwise matching operation for input bitsets
    # assumes the two bitsets have the same length
    return all(
        (ntbit & bitset2[i]).to_bool() for i, ntbit in enumerate(bitset1) 
    )

def pam_search(pam: PAM, region: Region, debug: bool) -> Tuple[List[int], List[int]]:
    # assert that the input region has already bin encoded in bits
    assert hasattr(region, "_sequence_bits")
    pam.encode()  # encode pam linear time search
    stop_position = len(region) - len(pam) + 1  # matching limit on the right
    # lists storing hits for input pam on forward and reverse strands
    matches_fwd, matches_rev = [], [] 
    for pos in range(stop_position):  # start pam search 
        try:
            if match(region.bits[pos:(pos + len(pam))], pam.bits[0]):
                matches_fwd.append(pos)
        except ValueError as e:
            exception_handler(CrisprHawkBitsetError, "PAM bitwise matching failed on forward strand", os.EX_DATAERR, debug, e)
        try:
            if match(region.bits[pos:(pos + len(pam))], pam.bits[1]):
                matches_rev.append(pos)
        except ValueError as e:
            exception_handler(CrisprHawkBitsetError, "PAM bitwise matching failed on reverse strand", os.EX_DATAERR, debug, e)
    return matches_fwd, matches_rev

def retrieve_guides(matches_fwd: List[int], matches_rev: List[int], pam: PAM, region: Region, guidelen: int, right: bool) -> Tuple[List[str], List[str]]:
    if right:
        guides_fwd = [
            "".join(region[pos:pos + guidelen + len(pam)]) 
            for pos in matches_fwd 
            if guidelen + len(pam) <= len(region) - pos
        ]
        guides_rev = [
            "".join(region[pos - guidelen - len(pam):pos]) 
            for pos in matches_rev 
            if guidelen + len(pam) <= pos
        ]
    else: 
        guides_fwd = [
            "".join(region[pos - guidelen:pos])  
            for pos in matches_fwd 
            if guidelen <= pos
        ]
        guides_rev = [
            "".join(region[pos + len(pam):pos + guidelen + len(pam)]) 
            for pos in matches_rev
            if guidelen <= len(region) - pos
        ]
    return guides_fwd, guides_rev
        




def search(pam: PAM, region: Region, guidelen: int, right: bool, debug: bool):
    # search pam occurrences on forward and reverse strand of the input sequence
    matches_fwd, matches_rev = pam_search(pam, region, debug)
    # recover guide candidates found on forward and reverse strands
    guides_fwd, guides_rev = retrieve_guides(matches_fwd, matches_rev, pam, region, guidelen, right)
    print(len(guides_fwd), len(guides_rev))


            