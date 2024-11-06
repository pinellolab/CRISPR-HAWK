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
    return all((ntbit & bitset2[i]).to_bool() for i, ntbit in enumerate(bitset1))


def pam_search(
    pam: PAM, region: Region, guidelen: int, debug: bool
) -> Tuple[List[int], List[int]]:
    # assert that the input region has already bin encoded in bits
    assert hasattr(region, "_sequence_bits")
    pam.encode()  # encode pam linear time search
    # matching limit on the right - ignore pad nts
    stop_position = len(region) - guidelen - len(pam) + 1
    # lists storing hits for input pam on forward and reverse strands
    matches_fwd, matches_rev = [], []
    # start pam search from real region start - ignore pad nts
    for pos in range(guidelen, stop_position):
        try:
            if match(region.bits[pos : (pos + len(pam))], pam.bits[0]):
                matches_fwd.append(pos)
        except ValueError as e:
            exception_handler(
                CrisprHawkBitsetError,
                "PAM bitwise matching failed on forward strand",
                os.EX_DATAERR,
                debug,
                e,
            )
        try:
            if match(region.bits[pos : (pos + len(pam))], pam.bits[1]):
                matches_rev.append(pos)
        except ValueError as e:
            exception_handler(
                CrisprHawkBitsetError,
                "PAM bitwise matching failed on reverse strand",
                os.EX_DATAERR,
                debug,
                e,
            )
    return matches_fwd, matches_rev


def extract_guide(
    region: Region, pos: int, guidelen: int, pamlen: int, right: bool
) -> str:
    if right:
        return region[pos : pos + guidelen + pamlen]
    return region[pos - guidelen : pos + pamlen]


def valid_position(pos: int, guidelen: int, regionlen: int, right: bool) -> bool:
    if right:
        return guidelen <= regionlen - pos
    return guidelen <= pos


def retrieve_guides(
    matches_fwd: List[int],
    matches_rev: List[int],
    pamlen: int,
    region: Region,
    guidelen: int,
    right: bool,
) -> Tuple[List[List[str]], List[List[str]]]:

    regionlen = len(region)  # length of the input region
    # retrieve guides from the forward strand
    guides_fwd = [
        extract_guide(region, pos, guidelen, pamlen, right)
        for pos in matches_fwd
        if valid_position(pos, guidelen, regionlen, right)
    ]
    # retrieve guides from the reverse strand
    guides_rev = [
        extract_guide(region, pos, guidelen, pamlen, not right)
        for pos in matches_rev
        if valid_position(pos, guidelen, regionlen, not right)
    ]
    return guides_fwd, guides_rev


def search(
    pam: PAM, region: Region, guidelen: int, right: bool, debug: bool
) -> Tuple[Tuple[List[int], List[int]], Tuple[List[List[str]], List[List[str]]]]:
    # search pam occurrences on forward and reverse strand of the input sequence
    matches_fwd, matches_rev = pam_search(pam, region, guidelen, debug)
    # recover guide candidates found on forward and reverse strands
    guides_fwd, guides_rev = retrieve_guides(
        matches_fwd, matches_rev, len(pam), region, guidelen, right
    )
    return (matches_fwd, matches_rev), (guides_fwd, guides_rev)
